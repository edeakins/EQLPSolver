/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2021 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/*    Authors: Julian Hall, Ivet Galabova, Qi Huangfu, Leona Gottwald    */
/*    and Michael Feldmeier                                              */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "mip/HighsMipSolver.h"

#include "lp_data/HighsLpUtils.h"
#include "lp_data/HighsModelUtils.h"
#include "mip/HighsCliqueTable.h"
#include "mip/HighsCutPool.h"
#include "mip/HighsDomain.h"
#include "mip/HighsImplications.h"
#include "mip/HighsLpRelaxation.h"
#include "mip/HighsMipSolverData.h"
#include "mip/HighsPseudocost.h"
#include "mip/HighsSearch.h"
#include "mip/HighsSeparation.h"
#include "presolve/HPresolve.h"
#include "presolve/HighsPostsolveStack.h"
#include "presolve/PresolveComponent.h"
#include "util/HighsCDouble.h"

HighsMipSolver::HighsMipSolver(const HighsOptions& options, const HighsLp& lp,
                               const HighsSolution& solution, bool submip)
    : options_mip_(&options),
      model_(&lp),
      solution_objective_(kHighsInf),
      submip(submip),
      rootbasis(nullptr),
      pscostinit(nullptr),
      clqtableinit(nullptr),
      implicinit(nullptr) {}

HighsMipSolver::~HighsMipSolver() = default;

void HighsMipSolver::run() {
  modelstatus_ = HighsModelStatus::kNotset;
  // std::cout << options_mip_->presolve << std::endl;
  timer_.start(timer_.solve_clock);

  mipdata_ = decltype(mipdata_)(new HighsMipSolverData(*this));
  mipdata_->init();
  mipdata_->runPresolve();
  if (modelstatus_ != HighsModelStatus::kNotset) {
    highsLogUser(options_mip_->log_options, HighsLogType::kInfo,
                 "Presolve: %s\n",
                 utilModelStatusToString(modelstatus_).c_str());
    if (modelstatus_ == HighsModelStatus::kOptimal) {
      mipdata_->lower_bound = 0;
      mipdata_->upper_bound = 0;
      mipdata_->transformNewIncumbent(std::vector<double>());
    }
    cleanupSolve();
    return;
  }

  mipdata_->runSetup();
restart:
  if (modelstatus_ == HighsModelStatus::kNotset) {
    mipdata_->evaluateRootNode();
  }
  if (mipdata_->nodequeue.empty()) {
    cleanupSolve();
    return;
  }

  std::shared_ptr<const HighsBasis> basis;
  HighsSearch search{*this, mipdata_->pseudocost};
  mipdata_->debugSolution.registerDomain(search.getLocalDomain());
  HighsSeparation sepa(*this);

  search.setLpRelaxation(&mipdata_->lp);
  sepa.setLpRelaxation(&mipdata_->lp);

  mipdata_->lower_bound = mipdata_->nodequeue.getBestLowerBound();

  highsLogUser(options_mip_->log_options, HighsLogType::kInfo,
               "\nstarting tree search\n");
  mipdata_->printDisplayLine();
  search.installNode(mipdata_->nodequeue.popBestBoundNode());
  int64_t numStallNodes = 0;
  int64_t lastLbLeave = 0;
  int64_t numQueueLeaves = 0;
  HighsInt numHugeTreeEstim = 0;
  int64_t numNodesLastCheck = mipdata_->num_nodes;
  double treeweightLastCheck = 0.0;
  while (search.hasNode()) {
    // set iteration limit for each lp solve during the dive to 10 times the
    // average nodes

    HighsInt iterlimit = 10 * mipdata_->lp.getAvgSolveIters();
    iterlimit = std::max(HighsInt{10000}, iterlimit);

    mipdata_->lp.setIterationLimit(iterlimit);

    // perform the dive and put the open nodes to the queue
    size_t plungestart = mipdata_->num_nodes;
    bool limit_reached = false;
    bool heuristicsCalled = false;
    while (true) {
      if (!heuristicsCalled && mipdata_->moreHeuristicsAllowed()) {
        search.evaluateNode();
        if (search.currentNodePruned()) {
          ++mipdata_->num_leaves;
          search.flushStatistics();
          break;
        }

        heuristicsCalled = true;

        if (mipdata_->incumbent.empty())
          mipdata_->heuristics.randomizedRounding(
              mipdata_->lp.getLpSolver().getSolution().col_value);

        if (mipdata_->incumbent.empty())
          mipdata_->heuristics.RENS(
              mipdata_->lp.getLpSolver().getSolution().col_value);
        else
          mipdata_->heuristics.RINS(
              mipdata_->lp.getLpSolver().getSolution().col_value);

        mipdata_->heuristics.flushStatistics();
      }

      if (mipdata_->domain.infeasible()) break;

      search.dive();
      ++mipdata_->num_leaves;

      search.flushStatistics();
      if (mipdata_->checkLimits()) {
        limit_reached = true;
        break;
      }

      if (!search.backtrack()) break;

      if (search.getCurrentEstimate() >= mipdata_->upper_limit) break;

      if (mipdata_->num_nodes - plungestart >=
          std::min(100., mipdata_->num_nodes * 0.1))
        break;

      if (mipdata_->dispfreq != 0) {
        if (mipdata_->num_leaves - mipdata_->last_displeave >=
            std::min(mipdata_->dispfreq,
                     1 + int64_t(0.01 * mipdata_->num_leaves)))
          mipdata_->printDisplayLine();
      }

      // printf("continue plunging due to good esitmate\n");
    }
    search.openNodesToQueue(mipdata_->nodequeue);
    mipdata_->lower_bound = std::min(mipdata_->upper_bound,
                                     mipdata_->nodequeue.getBestLowerBound());

    if (limit_reached) break;

    if (mipdata_->dispfreq != 0) {
      if (mipdata_->num_leaves - mipdata_->last_displeave >=
          std::min(mipdata_->dispfreq,
                   1 + int64_t(0.01 * mipdata_->num_leaves)))
        mipdata_->printDisplayLine();
    }

    // the search datastructure should have no installed node now
    assert(!search.hasNode());

    // propagate the global domain
    mipdata_->domain.propagate();
    mipdata_->pruned_treeweight += mipdata_->nodequeue.pruneInfeasibleNodes(
        mipdata_->domain, mipdata_->feastol);

    // if global propagation detected infeasibility, stop here
    if (mipdata_->domain.infeasible() || mipdata_->nodequeue.empty()) {
      mipdata_->nodequeue.clear();
      mipdata_->pruned_treeweight = 1.0;
      mipdata_->lower_bound = std::min(kHighsInf, mipdata_->upper_bound);
      break;
    }

    // if global propagation found bound changes, we update the local domain
    if (!mipdata_->domain.getChangedCols().empty()) {
      highsLogDev(options_mip_->log_options, HighsLogType::kInfo,
                  "added %" HIGHSINT_FORMAT " global bound changes\n",
                  (HighsInt)mipdata_->domain.getChangedCols().size());
      mipdata_->cliquetable.cleanupFixed(mipdata_->domain);
      for (HighsInt col : mipdata_->domain.getChangedCols())
        mipdata_->implications.cleanupVarbounds(col);

      mipdata_->domain.setDomainChangeStack(std::vector<HighsDomainChange>());
      search.resetLocalDomain();

      mipdata_->domain.clearChangedCols();
      mipdata_->removeFixedIndices();
    }

    if (!submip) {
      double currNodeEstim =
          numNodesLastCheck + (mipdata_->num_nodes - numNodesLastCheck) *
                                  double(1.0 - mipdata_->pruned_treeweight) /
                                  std::max(double(mipdata_->pruned_treeweight -
                                                  treeweightLastCheck),
                                           mipdata_->feastol);

      if (currNodeEstim >=
          1000 * (mipdata_->numRestarts + 1) * mipdata_->num_nodes) {
        ++numHugeTreeEstim;
        // if (!submip)
        //   printf("%" HIGHSINT_FORMAT " (nodeestim: %.1f)\n",
        //   numHugeTreeEstim, currNodeEstim);
      } else {
        numHugeTreeEstim = 0;
        treeweightLastCheck = double(mipdata_->pruned_treeweight);
        numNodesLastCheck = mipdata_->num_nodes;
      }

      if (numHugeTreeEstim >= 50) {
        mipdata_->performRestart();
        goto restart;
      }
    }

    // remove the iteration limit when installing a new node
    // mipdata_->lp.setIterationLimit();

    // loop to install the next node for the search
    while (!mipdata_->nodequeue.empty()) {
      // printf("popping node from nodequeue (length = %" HIGHSINT_FORMAT ")\n",
      // (HighsInt)nodequeue.size());
      assert(!search.hasNode());

      if (numQueueLeaves - lastLbLeave >= 10) {
        search.installNode(mipdata_->nodequeue.popBestBoundNode());
        lastLbLeave = numQueueLeaves;
      } else {
        search.installNode(mipdata_->nodequeue.popBestNode());
        if (search.getCurrentLowerBound() == mipdata_->lower_bound)
          lastLbLeave = numQueueLeaves;
      }

      ++numQueueLeaves;

      if (search.getCurrentEstimate() >= mipdata_->upper_limit) {
        ++numStallNodes;
        if (options_mip_->mip_max_stall_nodes != kHighsIInf &&
            numStallNodes >= options_mip_->mip_max_stall_nodes) {
          limit_reached = true;
          modelstatus_ = HighsModelStatus::kIterationLimit;
          break;
        }
      } else
        numStallNodes = 0;

      assert(search.hasNode());

      // set the current basis if available
      if (basis) {
        mipdata_->lp.setStoredBasis(basis);
        mipdata_->lp.recoverBasis();
      }

      // we evaluate the node directly here instead of performing a dive
      // because we first want to check if the node is not fathomed due to
      // new global information before we perform separation rounds for the node
      search.evaluateNode();

      // if the node was pruned we remove it from the search and install the
      // next node from the queue
      if (search.currentNodePruned()) {
        search.backtrack();
        ++mipdata_->num_leaves;
        ++mipdata_->num_nodes;
        search.flushStatistics();

        if (mipdata_->domain.infeasible()) {
          mipdata_->nodequeue.clear();
          mipdata_->pruned_treeweight = 1.0;
          mipdata_->lower_bound = std::min(kHighsInf, mipdata_->upper_bound);
          break;
        }

        if (mipdata_->checkLimits()) {
          limit_reached = true;
          break;
        }

        mipdata_->lower_bound = std::min(
            mipdata_->upper_bound, mipdata_->nodequeue.getBestLowerBound());

        if (mipdata_->dispfreq != 0) {
          if (mipdata_->num_leaves - mipdata_->last_displeave >=
              std::min(mipdata_->dispfreq,
                       1 + int64_t(0.01 * mipdata_->num_leaves)))
            mipdata_->printDisplayLine();
        }
        continue;
      }

      // the node is still not fathomed, so perform separation
      sepa.separate(search.getLocalDomain());

      // after separation we store the new basis and proceed with the outer loop
      // to perform a dive from this node
      if (mipdata_->lp.getStatus() != HighsLpRelaxation::Status::kError &&
          mipdata_->lp.getStatus() != HighsLpRelaxation::Status::kNotSet)
        mipdata_->lp.storeBasis();

      basis = mipdata_->lp.getStoredBasis();
      if (!basis || !isBasisConsistent(mipdata_->lp.getLp(), *basis)) {
        HighsBasis b = mipdata_->firstrootbasis;
        b.row_status.resize(mipdata_->lp.numRows(), HighsBasisStatus::kBasic);
        basis = std::make_shared<const HighsBasis>(std::move(b));
        mipdata_->lp.setStoredBasis(basis);
      }

      break;
    }

    if (limit_reached) break;
  }

  cleanupSolve();
}

void HighsMipSolver::cleanupSolve() {
  timer_.start(timer_.postsolve_clock);
  bool havesolution = solution_objective_ != kHighsInf;
  dual_bound_ = mipdata_->lower_bound + model_->offset_;
  primal_bound_ = mipdata_->upper_bound + model_->offset_;
  node_count_ = mipdata_->num_nodes;

  if (modelstatus_ == HighsModelStatus::kNotset) {
    if (havesolution)
      modelstatus_ = HighsModelStatus::kOptimal;
    else
      modelstatus_ = HighsModelStatus::kInfeasible;
  }

  model_ = orig_model_;
  timer_.stop(timer_.postsolve_clock);
  timer_.stop(timer_.solve_clock);
  std::string solutionstatus = "-";

  if (havesolution) {
    bool feasible =
        bound_violation_ <= options_mip_->mip_feasibility_tolerance &&
        integrality_violation_ <= options_mip_->mip_feasibility_tolerance &&
        row_violation_ <= options_mip_->mip_feasibility_tolerance;
    solutionstatus = feasible ? "feasible" : "infeasible";
  }
  highsLogUser(options_mip_->log_options, HighsLogType::kInfo,
               "\nSolving report\n"
               "  Status            %s\n"
               "  Primal bound      %.12g\n"
               "  Dual bound        %.12g\n"
               "  Solution status   %s\n",
               utilModelStatusToString(modelstatus_).c_str(), primal_bound_,
               dual_bound_, solutionstatus.c_str());
  if (solutionstatus != "-")
    highsLogUser(options_mip_->log_options, HighsLogType::kInfo,
                 "                    %.12g (objective)\n"
                 "                    %.12g (bound viol.)\n"
                 "                    %.12g (int. viol.)\n"
                 "                    %.12g (row viol.)\n",
                 solution_objective_, bound_violation_, integrality_violation_,
                 row_violation_);
  highsLogUser(options_mip_->log_options, HighsLogType::kInfo,
               "  Timing            %.2f (total)\n"
               "                    %.2f (presolve)\n"
               "                    %.2f (postsolve)\n"
               "  Nodes             %llu\n"
               "  LP iterations     %llu (total)\n"
               "                    %llu (strong br.)\n"
               "                    %llu (separation)\n"
               "                    %llu (heuristics)\n",
               timer_.read(timer_.solve_clock),
               timer_.read(timer_.presolve_clock),
               timer_.read(timer_.postsolve_clock),
               (long long unsigned)mipdata_->num_nodes,
               (long long unsigned)mipdata_->total_lp_iterations,
               (long long unsigned)mipdata_->sb_lp_iterations,
               (long long unsigned)mipdata_->sepa_lp_iterations,
               (long long unsigned)mipdata_->heuristic_lp_iterations);

  assert(modelstatus_ != HighsModelStatus::kNotset);
}
