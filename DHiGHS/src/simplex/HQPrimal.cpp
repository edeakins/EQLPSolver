/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file simplex/HQPrimal.cpp
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#include "simplex/HQPrimal.h"
#include "io/HighsIO.h"
#include "lp_data/HConst.h"
#include "simplex/HSimplex.h"
#include "simplex/SimplexTimer.h"
#include "simplex/HighsSimplexInterface.h"
#include "util/HighsRandom.h"
#include "util/HighsUtils.h"

#include <cassert>
#include <cstdio>
#include <iostream>

using std::runtime_error;

double fTranTime = 0;
double Chuzc1Time = 0;
double Chuzc2Time = 0;
double updatePrimalTime = 0;
double collectPrimalInfsTime = 0;
double bTranTime = 0;
double updateDualTime = 0;
double updateFactorTime = 0;
double invertTime = 0;
double computePrimalTime = 0;
double computePrimalObjTime = 0;
double computeDualTime = 0;
double computePrimalInfsTime = 0;
double computeDualInfsTime = 0;
double reportRebTime = 0;
bool liftStart = false;
int ivHint = 0;

HighsStatus HQPrimal::solve() {
  HighsOptions& options = workHMO.options_;
  HighsSolutionParams& scaled_solution_params = workHMO.scaled_solution_params_;
  HighsSimplexInfo& simplex_info = workHMO.simplex_info_;
  HighsSimplexLpStatus& simplex_lp_status = workHMO.simplex_lp_status_;
  workHMO.scaled_model_status_ = HighsModelStatus::NOTSET;
  // Assumes that the LP has a positive number of rows, since
  // unconstrained LPs should be solved in solveLpSimplex
  bool positive_num_row = workHMO.simplex_lp_.numRow_ > 0;
  assert(positive_num_row);
  if (!positive_num_row) {
    HighsLogMessage(options.logfile, HighsMessageType::ERROR,
		    "HPrimal::solve called for LP with non-positive (%d) number of constraints",
		    workHMO.simplex_lp_.numRow_);
    return HighsStatus::Error;
  }
  HighsTimer& timer = workHMO.timer_;
  invertHint = INVERT_HINT_NO;

  // Setup aspects of the model data which are needed for solve() but better
  // left until now for efficiency reasons.
  // ToDo primal simplex version
  // setup_for_solve(workHMO);

  // Set SolveBailout to be true if control is to be returned immediately to
  // calling function
  // ToDo Move to Simplex
  //  SolveBailout = false;

  // Initialise working environment
  // Does LOTS, including initialisation of edge weights. Should only
  // be called if model dimension changes
  // ToDo primal simplex version
  // init();
  // initParallel();

  // ToDo primal simplex version
  // initialise_cost(workHMO, 1); //  model->initCost(1);
  assert(simplex_lp_status.has_fresh_invert);
  if (!simplex_lp_status.has_fresh_invert) {
    printf(
        "ERROR: Should enter with fresh INVERT - unless no_invert_on_optimal "
        "is set\n");
  }
  // Consider initialising edge weights - create Primal variants
  //
#ifdef HiGHSDEV
  //  printf("simplex_lp_status.has_dual_steepest_edge_weights 2 = %d;
  //  dual_edge_weight_mode = %d; DualEdgeWeightMode::STEEPEST_EDGE = %d\n",
  //	 simplex_lp_status.has_dual_steepest_edge_weights,
  //dual_edge_weight_mode, DualEdgeWeightMode::STEEPEST_EDGE);cout<<flush;
  //  printf("Edge weights known? %d\n",
  //  !simplex_lp_status.has_dual_steepest_edge_weights);cout<<flush;
#endif
  /*
  if (!simplex_lp_status.has_dual_steepest_edge_weights) {
    // Edge weights are not known
    // Set up edge weights according to dual_edge_weight_mode and
  initialise_dual_steepest_edge_weights
    // Using dual Devex edge weights
    // Zero the number of Devex frameworks used and set up the first one
    devex_index.assign(solver_num_tot, 0);
    initialiseDevexFramework();
    // Indicate that edge weights are known
    simplex_lp_status.has_dual_steepest_edge_weights = true;
  }
  */

  // ToDo Determine primal simplex phase from initial primal values
  //
  /*
  compute_primal(workHMO);
  compute_primal_infeasible_in_??(workHMO, &dualInfeasCount);
  solvePhase = ??InfeasCount > 0 ? 1 : 2;
  */
  solvePhase = 0;  // Frig to skip while (solvePhase) {*}

  // Check that the model is OK to solve:
  //
  // Level 0 just checks the flags
  //
  // Level 1 also checks that the basis is OK and that the necessary
  // data in work* is populated.
  //
  // Level 2 (will) checks things like the nonbasic duals and basic
  // primal values
  //
  // Level 3 (will) checks expensive things like the INVERT and
  // steepeest edge weights
  //
  // ToDo Write primal simplex equivalent
  /*
  bool ok = ok_to_solve(workHMO, 1, solvePhase);
  if (!ok) {printf("NOT OK TO SOLVE???\n"); cout << flush;}
  assert(ok);
  */
#ifdef HiGHSDEV
  //  reportSimplexLpStatus(simplex_lp_status, "Before HQPrimal major solving
  //  loop");
#endif
  // The major solving loop

  // Initialise the iteration analysis. Necessary for strategy, but
  // much is for development and only switched on with HiGHSDEV
  // ToDo Move to simplex and adapt so it's OK for primal and dual
  //  iterationAnalysisInitialise();

  while (solvePhase) {
    /*
    int it0 = scaled_solution_params.simplex_iteration_count;
    switch (solvePhase) {
      case 1:
        timer.start(simplex_info.clock_[SimplexPrimalPhase1Clock]);
        solvePhase1();
        timer.stop(simplex_info.clock_[SimplexPrimalPhase1Clock]);
        simplex_info.primal_phase1_iteration_count +=
    (scaled_solution_params.simplex_iteration_count - it0); break; case 2:
        timer.start(simplex_info.clock_[SimplexPrimalPhase2Clock]);
        solvePhase2();
        timer.stop(simplex_info.clock_[SimplexPrimalPhase2Clock]);
        simplex_info.primal_phase2_iteration_count +=
    (scaled_solution_params.simplex_iteration_count - it0); break; case 4: break; default:
        solvePhase = 0;
        break;
    }
    // Jump for primal
    if (solvePhase == 4) break;
    // Possibly bail out
    if (SolveBailout) break;
    */
  }
  solvePhase = 2;
  if (workHMO.options_.simplex_strategy == SIMPLEX_STRATEGY_SWAP)
    solvePhase = 4;
  if (workHMO.options_.simplex_strategy == SIMPLEX_STRATEGY_UNFOLD)
    solvePhase = 3;
  if (workHMO.scaled_model_status_ != HighsModelStatus::REACHED_TIME_LIMIT) {
    if (solvePhase == 2) {
      int it0 = scaled_solution_params.simplex_iteration_count;

      timer.start(simplex_info.clock_[SimplexPrimalPhase2Clock]);
      solvePhase2();
      timer.stop(simplex_info.clock_[SimplexPrimalPhase2Clock]);

      simplex_info.primal_phase2_iteration_count +=
          (scaled_solution_params.simplex_iteration_count - it0);
    }
    else if (solvePhase == 3){
      int it0 = scaled_solution_params.simplex_iteration_count;

      timer.start(simplex_info.clock_[SimplexPrimalPhase2Clock]);
      solvePhase3();
      timer.stop(simplex_info.clock_[SimplexPrimalPhase2Clock]);

      simplex_info.primal_phase2_iteration_count +=
          (scaled_solution_params.simplex_iteration_count - it0);
    }
    else if (solvePhase == 4){
      int it0 = scaled_solution_params.simplex_iteration_count;

      timer.start(simplex_info.clock_[SimplexPrimalPhase2Clock]);
      solvePhaseSwap();
      timer.stop(simplex_info.clock_[SimplexPrimalPhase2Clock]);

      simplex_info.primal_phase2_iteration_count +=
          (scaled_solution_params.simplex_iteration_count - it0);
    }
  }
  /*
  // ToDo Adapt ok_to_solve to be used by primal
  bool ok = ok_to_solve(workHMO, 1, solvePhase);// model->OKtoSolve(1,
  solvePhase); if (!ok) {printf("NOT OK After Solve???\n"); cout << flush;}
  assert(ok);
  */
  // buildTableau();
  return HighsStatus::OK;
}

// This function collects the reduced Amatrix after a master iteration of 
// the unfolding procedure
void HQPrimal::buildTableau(){
  int _numRow_ = workHMO.lp_.numRow_;// - workHMO.lp_.numResiduals_;
  int _numCol_ = workHMO.lp_.numCol_;
  HighsTableau& tableau = workHMO.tableau_;
  HVector rowAp;
  int nnz = 0;
  rowAp.setup(row_ap.array.size());
  // tableau.ARtableauStart.push_back(0);
  std::cout << "[" << std::endl;
  for (int i = 0; i < _numRow_; ++i){
    if (!workHMO.lp_.degenSlack_[i]) continue;
    // tableau.ARreducedRHS.push_back(workHMO.simplex_info_.baseValue_[i]);
    row_ep.clear();
    row_ep.count = 1;
    row_ep.index[0] = i;
    row_ep.array[i] = 1;
    row_ep.packFlag = true;
    workHMO.factor_.btran(row_ep, analysis->row_ep_density);
    computeTableauRowFull(workHMO, row_ep, rowAp);
    std::cout << "[ ";
    for (int j = 0; j < _numCol_; ++j){
      if (j == _numCol_ - 1){
        std::cout << rowAp.array[j] << " ],"; 
        break;
      }
      std::cout << rowAp.array[j] << ", ";
      // if (fabs(rowAp.array[j]) > 1e-10){
      //   nnz++;
      //   tableau.ARtableauIndex.push_back(j);
      //   tableau.ARtableauValue.push_back(rowAp.array[j]);
      // }
    }
    std::cout << std::endl;
    std::cin.get();
    // tableau.ARtableauStart.push_back(nnz);
  }
  std::cout << "]" << std::endl;
  std::cin.get();
}

void HQPrimal::solvePhaseSwap() {
  HighsSimplexInfo& simplex_info = workHMO.simplex_info_;
  HighsSimplexLpStatus& simplex_lp_status = workHMO.simplex_lp_status_;
  HighsTimer& timer = workHMO.timer_;
  printf("HQPrimal::solvePhase3\n");
  // When starting a new phase the (updated) primal objective function
  // value isn't known. Indicate this so that when the value
  // computed from scratch in build() isn't checked against the the
  // updated value
  simplex_lp_status.has_primal_objective_value = 0;
  // Set invertHint so that it's assigned when first tested
  invertHint = INVERT_HINT_NO;
  // Set solvePhase=2 so it's set if solvePhase2() is called directly
  solvePhase = 4;
  // Set up local copies of model dimensions
  solver_num_col = workHMO.simplex_lp_.numCol_;
  solver_num_row = workHMO.simplex_lp_.numRow_;
  solver_num_tot = solver_num_col + solver_num_row;

  analysis = &workHMO.simplex_analysis_;

  // Setup update limits
  simplex_info.update_limit =
      min(100 + solver_num_row / 100,
          1000);  // TODO: Consider allowing the dual limit to be used
  simplex_info.update_count = 0;

  // Setup local vectors
  col_aq.setup(solver_num_row);
  row_ep.setup(solver_num_row);
  row_ap.setup(solver_num_col);

  // Setup local dep graph vectors
  depStart.resize(workHMO.lp_.numResiduals_ + 1);

  ph1SorterR.reserve(solver_num_row);
  ph1SorterT.reserve(solver_num_row);

#ifdef HiGHSDEV
  printf("HQPrimal::solvePhase3 - WARNING: Setting analysis->col_aq_density = 0\n");
  printf("HQPrimal::solvePhase3 - WARNING: Setting analysis->row_ep_density = 0\n");
#endif
  analysis->col_aq_density = 0;
  analysis->row_ep_density = 0;

  devexReset();

  no_free_columns = true;
  for (int iCol = 0; iCol < solver_num_tot; iCol++) {
    if (highs_isInfinity(-workHMO.simplex_info_.workLower_[iCol])) {
      if (highs_isInfinity(workHMO.simplex_info_.workUpper_[iCol])) {
        // Free column
        no_free_columns = false;
        break;
      }
    }
  }
#ifdef HiGHSDEV
  if (no_free_columns) {
    printf("Model has no free columns\n");
  } else {
    printf("Model has free columns\n");
  }
#endif

  // // Setup other buffers

  // HighsPrintMessage(workHMO.options_.output, workHMO.options_.message_level, ML_DETAILED,
  //       "primal-phase3-start\n");
  // // Main solving structure
  // pivotArtificial.assign(workHMO.lp_.artificialVariables.size(), false);
  // for (;;) {
  //   timer.start(simplex_info.clock_[IteratePrimalRebuildClock]);
  //   primalRebuild();
  //   timer.stop(simplex_info.clock_[IteratePrimalRebuildClock]);

  //   if (isPrimalPhase1) {
  //     for (;;) {
  //       /* Primal phase 1 choose column */
  //       phase1ChooseColumnArtificial();
  //       std::cout << "columnIn: " << columnIn << std::endl;
  //       if (columnIn == -1) {
  //         printf("==> Primal phase 1 choose column failed.\n");
  //         invertHint = INVERT_HINT_CHOOSE_COLUMN_FAIL;
  //         break;
  //       }

  //       /* Primal phsae 1 choose row */
  //       phase1ChooseRow();
  //       std::cout << "rowOut: " << rowOut << std::endl;
  //       if (rowOut == -1) {
  //         printf("Primal phase 1 choose row failed.\n");
  //         exit(0);
  //       }
  //       std::cin.get();
  //       /* Primal phase 1 update */
  //       phase1Update();
  //       if (invertHint) {
  //         break;
  //       }
  //     }
  //     /* Go to the next rebuild */
  //     if (invertHint) {
  //       /* Stop when the invert is new */
  //       if (simplex_lp_status.has_fresh_rebuild) {
  //         break;
  //       }
  //       continue;
  //     }
  //   }
  // }
  rSwapped.resize(workHMO.lp_.numCol_);
  sSwapped.resize(workHMO.lp_.numRow_);
  swapDegenerate();
  // buildTableau();

  if (workHMO.scaled_model_status_ == HighsModelStatus::REACHED_TIME_LIMIT) {
    return;
  }

  if (columnIn == -1) {
    HighsPrintMessage(workHMO.options_.output, workHMO.options_.message_level, ML_DETAILED,
          "primal-optimal\n");
    HighsPrintMessage(workHMO.options_.output, workHMO.options_.message_level, ML_DETAILED,
          "problem-optimal\n");
    workHMO.scaled_model_status_ = HighsModelStatus::OPTIMAL;
  } else {
    HighsPrintMessage(workHMO.options_.output, workHMO.options_.message_level, ML_MINIMAL,
          "primal-unbounded\n");
    workHMO.scaled_model_status_ = HighsModelStatus::PRIMAL_UNBOUNDED;
  }
  computeDualObjectiveValue(workHMO);
}

void HQPrimal::solvePhase3() {
  HighsSimplexInfo& simplex_info = workHMO.simplex_info_;
  HighsSimplexLpStatus& simplex_lp_status = workHMO.simplex_lp_status_;
  HighsTimer& timer = workHMO.timer_;
  printf("HQPrimal::solvePhase3\n");
  // When starting a new phase the (updated) primal objective function
  // value isn't known. Indicate this so that when the value
  // computed from scratch in build() isn't checked against the the
  // updated value
  simplex_lp_status.has_primal_objective_value = 0;
  // Set invertHint so that it's assigned when first tested
  invertHint = INVERT_HINT_NO;
  // Set solvePhase=2 so it's set if solvePhase2() is called directly
  solvePhase = 3;
  // Set up local copies of model dimensions
  solver_num_col = workHMO.simplex_lp_.numCol_;
  solver_num_row = workHMO.simplex_lp_.numRow_;
  solver_num_tot = solver_num_col + solver_num_row;

  analysis = &workHMO.simplex_analysis_;

  // Setup update limits
  simplex_info.update_limit =
      min(100 + solver_num_row / 100,
          1000);  // TODO: Consider allowing the dual limit to be used
  simplex_info.update_count = 0;

  // Setup local vectors
  col_aq.setup(solver_num_row);
  row_ep.setup(solver_num_row);
  row_ap.setup(solver_num_col);

  // Setup local dep graph vectors
  depStart.resize(workHMO.lp_.numResiduals_ + 1);

  ph1SorterR.reserve(solver_num_row);
  ph1SorterT.reserve(solver_num_row);

#ifdef HiGHSDEV
  printf("HQPrimal::solvePhase3 - WARNING: Setting analysis->col_aq_density = 0\n");
  printf("HQPrimal::solvePhase3 - WARNING: Setting analysis->row_ep_density = 0\n");
#endif
  analysis->col_aq_density = 0;
  analysis->row_ep_density = 0;

  devexReset();

  no_free_columns = true;
  for (int iCol = 0; iCol < solver_num_tot; iCol++) {
    if (highs_isInfinity(-workHMO.simplex_info_.workLower_[iCol])) {
      if (highs_isInfinity(workHMO.simplex_info_.workUpper_[iCol])) {
        // Free column
        no_free_columns = false;
        break;
      }
    }
  }
#ifdef HiGHSDEV
  if (no_free_columns) {
    printf("Model has no free columns\n");
  } else {
    printf("Model has free columns\n");
  }
#endif

  // // Setup other buffers

  // HighsPrintMessage(workHMO.options_.output, workHMO.options_.message_level, ML_DETAILED,
  //       "primal-phase3-start\n");
  // // Main solving structure
  // pivotArtificial.assign(workHMO.lp_.artificialVariables.size(), false);
  // for (;;) {
  //   timer.start(simplex_info.clock_[IteratePrimalRebuildClock]);
  //   primalRebuild();
  //   timer.stop(simplex_info.clock_[IteratePrimalRebuildClock]);

  //   if (isPrimalPhase1) {
  //     for (;;) {
  //       /* Primal phase 1 choose column */
  //       phase1ChooseColumnArtificial();
  //       std::cout << "columnIn: " << columnIn << std::endl;
  //       if (columnIn == -1) {
  //         printf("==> Primal phase 1 choose column failed.\n");
  //         invertHint = INVERT_HINT_CHOOSE_COLUMN_FAIL;
  //         break;
  //       }

  //       /* Primal phsae 1 choose row */
  //       phase1ChooseRow();
  //       std::cout << "rowOut: " << rowOut << std::endl;
  //       if (rowOut == -1) {
  //         printf("Primal phase 1 choose row failed.\n");
  //         exit(0);
  //       }
  //       std::cin.get();
  //       /* Primal phase 1 update */
  //       phase1Update();
  //       if (invertHint) {
  //         break;
  //       }
  //     }
  //     /* Go to the next rebuild */
  //     if (invertHint) {
  //       /* Stop when the invert is new */
  //       if (simplex_lp_status.has_fresh_rebuild) {
  //         break;
  //       }
  //       continue;
  //     }
  //   }
  // }
  unfold();
  // buildTableau();

  if (workHMO.scaled_model_status_ == HighsModelStatus::REACHED_TIME_LIMIT) {
    return;
  }

  if (columnIn == -1) {
    HighsPrintMessage(workHMO.options_.output, workHMO.options_.message_level, ML_DETAILED,
          "primal-optimal\n");
    HighsPrintMessage(workHMO.options_.output, workHMO.options_.message_level, ML_DETAILED,
          "problem-optimal\n");
    workHMO.scaled_model_status_ = HighsModelStatus::OPTIMAL;
  } else {
    HighsPrintMessage(workHMO.options_.output, workHMO.options_.message_level, ML_MINIMAL,
          "primal-unbounded\n");
    workHMO.scaled_model_status_ = HighsModelStatus::PRIMAL_UNBOUNDED;
  }
  computeDualObjectiveValue(workHMO);
}

void HQPrimal::solvePhase2() {
  HighsSimplexInfo& simplex_info = workHMO.simplex_info_;
  HighsSimplexLpStatus& simplex_lp_status = workHMO.simplex_lp_status_;
  HighsTimer& timer = workHMO.timer_;
  printf("HQPrimal::solvePhase2\n");
  // When starting a new phase the (updated) primal objective function
  // value isn't known. Indicate this so that when the value
  // computed from scratch in build() isn't checked against the the
  // updated value
  simplex_lp_status.has_primal_objective_value = 0;
  // Set invertHint so that it's assigned when first tested
  invertHint = INVERT_HINT_NO;
  // Set solvePhase=2 so it's set if solvePhase2() is called directly
  solvePhase = 2;
  // Set up local copies of model dimensions
  solver_num_col = workHMO.simplex_lp_.numCol_;
  solver_num_row = workHMO.simplex_lp_.numRow_;
  solver_num_tot = solver_num_col + solver_num_row;

  analysis = &workHMO.simplex_analysis_;

  // Setup update limits
  simplex_info.update_limit =
      min(100 + solver_num_row / 100,
          1000);  // TODO: Consider allowing the dual limit to be used
  simplex_info.update_count = 0;

  // Setup local vectors
  col_aq.setup(solver_num_row);
  row_ep.setup(solver_num_row);
  row_ap.setup(solver_num_col);

  ph1SorterR.reserve(solver_num_row);
  ph1SorterT.reserve(solver_num_row);

#ifdef HiGHSDEV
  printf("HQPrimal::solvePhase2 - WARNING: Setting analysis->col_aq_density = 0\n");
  printf("HQPrimal::solvePhase2 - WARNING: Setting analysis->row_ep_density = 0\n");
#endif
  analysis->col_aq_density = 0;
  analysis->row_ep_density = 0;

  devexReset();

  no_free_columns = true;
  for (int iCol = 0; iCol < solver_num_tot; iCol++) {
    if (highs_isInfinity(-workHMO.simplex_info_.workLower_[iCol])) {
      if (highs_isInfinity(workHMO.simplex_info_.workUpper_[iCol])) {
        // Free column
        no_free_columns = false;
        break;
      }
    }
  }
#ifdef HiGHSDEV
  if (no_free_columns) {
    printf("Model has no free columns\n");
  } else {
    printf("Model has free columns\n");
  }
#endif

  // Setup other buffers

  HighsPrintMessage(workHMO.options_.output, workHMO.options_.message_level, ML_DETAILED,
		    "primal-phase2-start\n");
  // Main solving structure
  for (;;) {
    timer.start(simplex_info.clock_[IteratePrimalRebuildClock]);
    primalRebuild();
    timer.stop(simplex_info.clock_[IteratePrimalRebuildClock]);

    if (isPrimalPhase1) {
      for (;;) {
        /* Primal phase 1 choose column */
        phase1ChooseColumn();
        if (columnIn == -1) {
          printf("==> Primal phase 1 choose column failed.\n");
          invertHint = INVERT_HINT_CHOOSE_COLUMN_FAIL;
          break;
        }

        /* Primal phsae 1 choose row */
        phase1ChooseRow();
        if (rowOut == -1) {
          printf("Primal phase 1 choose row failed.\n");
          exit(0);
        }

        /* Primal phase 1 update */
        phase1Update();
        if (invertHint) {
          break;
        }
      }
      /* Go to the next rebuild */
      if (invertHint) {
        /* Stop when the invert is new */
        if (simplex_lp_status.has_fresh_rebuild) {
          break;
        }
        continue;
      }
    }

    for (;;) {
      primalChooseColumn();
      if (columnIn == -1) {
        invertHint = INVERT_HINT_POSSIBLY_OPTIMAL;
        break;
      }
      primalChooseRow();
      if (rowOut == -1) {
        invertHint = INVERT_HINT_POSSIBLY_PRIMAL_UNBOUNDED;
        break;
      }
      primalUpdate();
      if (invertHint) {
        break;
      }
    }

    double currentRunHighsTime = timer.readRunHighsClock();
    if (currentRunHighsTime > workHMO.options_.time_limit) {
      workHMO.scaled_model_status_ = HighsModelStatus::REACHED_TIME_LIMIT;
      break;
    }
    // If the data are fresh from rebuild(), break out of
    // the outer loop to see what's ocurred
    // Was:	if (simplex_info.update_count == 0) break;
    if (simplex_lp_status.has_fresh_rebuild) {
#ifdef HiGHSDEV
      if (num_flip_since_rebuild)
        printf("Consider doing a primal rebuild if flips have occurred\n");
#endif
      //      if (num_flip_since_rebuild == 0)
      break;
    }
  }

  if (workHMO.scaled_model_status_ == HighsModelStatus::REACHED_TIME_LIMIT) {
    return;
  }

  if (columnIn == -1) {
    HighsPrintMessage(workHMO.options_.output, workHMO.options_.message_level, ML_DETAILED,
		      "primal-optimal\n");
    HighsPrintMessage(workHMO.options_.output, workHMO.options_.message_level, ML_DETAILED,
		      "problem-optimal\n");
    workHMO.scaled_model_status_ = HighsModelStatus::OPTIMAL;
  } else {
    HighsPrintMessage(workHMO.options_.output, workHMO.options_.message_level, ML_MINIMAL,
		      "primal-unbounded\n");
    workHMO.scaled_model_status_ = HighsModelStatus::PRIMAL_UNBOUNDED;
  }
  computeDualObjectiveValue(workHMO);
}

void HQPrimal::primalRebuild() {
  HighsSimplexInfo& simplex_info = workHMO.simplex_info_;
  HighsSimplexLpStatus& simplex_lp_status = workHMO.simplex_lp_status_;
  HighsTimer& timer = workHMO.timer_;
  // Move this to Simplex class once it's created
  //  simplex_method.record_pivots(-1, -1, 0);  // Indicate REINVERT

  // Rebuild workHMO.factor_ - only if we got updates
  int sv_invertHint = invertHint;
  invertHint = INVERT_HINT_NO;
  // Possibly Rebuild workHMO.factor_
  bool reInvert = simplex_info.update_count > 0;
  // if (liftStart) reInvert = true;
  if (!invert_if_row_out_negative) {
    // Don't reinvert if columnIn is negative [equivalently, if sv_invertHint ==
    // INVERT_HINT_POSSIBLY_OPTIMAL]
    if (sv_invertHint == INVERT_HINT_POSSIBLY_OPTIMAL) {
      assert(columnIn == -1);
      reInvert = false;
    }
  }
  double init;
  // std::cout << "reinvert: " << reInvert << std::endl;
  // std::cin.get();
  // std::cout << "reIvert variable: " << reInvert << std::endl;
  // std::cin.get();
  if (reInvert) {
    // std::cout << "reinvert called" << std::endl;
    // std::cin.get();
    double init = timer.readRunHighsClock();
    timer.start(simplex_info.clock_[InvertClock]);
    int rankDeficiency = compute_factor(workHMO);
    timer.stop(simplex_info.clock_[InvertClock]);
    double invertTime = timer.readRunHighsClock() - init;
    // std::cout << "invert time: " << invertTime << std::endl;
    if (rankDeficiency) {
      throw runtime_error("Primal reInvert: singular-basis-matrix");
    }
    simplex_info.update_count = 0;
  }
  timer.start(simplex_info.clock_[ComputeDualClock]);
  compute_dual(workHMO);
  timer.stop(simplex_info.clock_[ComputeDualClock]);

  timer.start(simplex_info.clock_[ComputePrimalClock]);
  compute_primal(workHMO);
  timer.stop(simplex_info.clock_[ComputePrimalClock]);

  // Primal objective section
  timer.start(simplex_info.clock_[ComputePrObjClock]);
  computePrimalObjectiveValue(workHMO);
  timer.stop(simplex_info.clock_[ComputePrObjClock]);

  double primal_objective_value = simplex_info.primal_objective_value;
#ifdef HiGHSDEV
  // Check the objective value maintained by updating against the
  // value when computed exactly - so long as there is a value to
  // check against
  if (simplex_lp_status.has_primal_objective_value) {
    double absPrimalObjectiveError = fabs(
        simplex_info.updated_primal_objective_value - primal_objective_value);
    double rlvPrimalObjectiveError =
        absPrimalObjectiveError / max(1.0, fabs(primal_objective_value));
    // TODO Investigate these Primal objective value errors
    if (rlvPrimalObjectiveError >= 1e-8) {
      HighsLogMessage(workHMO.options_.logfile, HighsMessageType::WARNING,
                      "Primal objective value error |rel| = %12g (%12g)",
                      absPrimalObjectiveError, rlvPrimalObjectiveError);
    }
  }
#endif
  simplex_info.updated_primal_objective_value = primal_objective_value;

  timer.start(simplex_info.clock_[ComputePrIfsClock]);
  computePrimalInfeasible(workHMO);
  timer.stop(simplex_info.clock_[ComputePrIfsClock]);

  timer.start(simplex_info.clock_[ComputeDuIfsClock]);
  computeDualInfeasible(workHMO);
  timer.stop(simplex_info.clock_[ComputeDuIfsClock]);

  /* Whether to switch to primal phase 1 */
  isPrimalPhase1 = 0;
  HighsSolutionParams& scaled_solution_params = workHMO.scaled_solution_params_;
  if (scaled_solution_params.num_primal_infeasibilities > 0) {
    isPrimalPhase1 = 1;
    phase1ComputeDual();
  }

  timer.start(simplex_info.clock_[ReportRebuildClock]);
  reportRebuild(sv_invertHint);
  timer.stop(simplex_info.clock_[ReportRebuildClock]);
#ifdef HiGHSDEV
  if (simplex_info.analyse_rebuild_time) {
    int iClock = simplex_info.clock_[IteratePrimalRebuildClock];
    int totalRebuilds = timer.clock_num_call[iClock];
    double totalRebuildTime = timer.read(iClock);
    printf(
        "Primal     rebuild %d (%1d) on iteration %9d: Total rebuild time %g\n",
        totalRebuilds, sv_invertHint, scaled_solution_params.simplex_iteration_count,
        totalRebuildTime);
  }
#endif
  num_flip_since_rebuild = 0;
  // Data are fresh from rebuild
  simplex_lp_status.has_fresh_rebuild = true;
}

void HQPrimal::buildSlackRSubMatrix(){
  int iCol, iRow, rowPointer, rowIdx = -1, colIdx = -1;
  HighsLp& workLp = workHMO.lp_;
  int numRow = workLp.numRow_;// workLp.numDegenRow_ + workLp.numResiduals_;
  int numCol = workLp.numCol_; // + workLp.numRow_ - workLp.numResiduals_;
  spSubMatRowPerm.assign(workLp.numRow_, -1);
  spSubMatRowUnperm.assign(workLp.numRow_, -1);
  for (iCol = 0; iCol < workLp.numCol_; ++iCol){
    for (rowPointer = workLp.Astart_[iCol]; rowPointer < workLp.Astart_[iCol + 1]; ++rowPointer){
      iRow = workLp.Aindex_[rowPointer];
      if (workLp.degenSlack_[iRow] || iRow >= workLp.numRow_ - workLp.numResiduals_){
        if (spSubMatRowPerm[iRow] == -1){ 
          spSubMatRowPerm[iRow] = ++rowIdx;
          spSubMatRowUnperm[rowIdx] = iRow;
        }
        spTriplets.push_back(Eigen::Triplet<double>(iRow, iCol, workLp.Avalue_[rowPointer]));
      }
    }
  }
  // for (iRow = 0; iRow < workLp.numRow_ - workLp.numResiduals_; ++iRow){
  //   if (rowIdxTrunc[iRow] > -1)
  //     spTriplets.push_back(Eigen::Triplet<double>(rowIdxTrunc[iRow], iRow + workLp.numCol_, -1));
  // }
  spSubMat.resize(numRow, numCol);
  spSubMat.setFromTriplets(spTriplets.begin(), spTriplets.end());
}

void HQPrimal::buildDependencyMatrix(){
  int col, row, idx;
  double mValue;
  depMatStart.push_back(0);
  HighsLp& workLp = workHMO.lp_;
  HighsBasis& workBasis = workHMO.basis_;
  std::vector<int> freq(workLp.numCol_, 0);
  int rowCnt = 0;
  int pairings = 0;
  for (int row = 0; row < workLp.numRow_ - workLp.numResiduals_; ++row){
    rowCnt = 0;
    if (!workLp.degenSlack_[row]){ 
      depMatStart.push_back(depMatIndex.size());
      continue;
    }
    // if (row >= workLp.numRow_ - workLp.numResiduals_) continue;
    row_ep.clear();
    row_ap.clear();
    row_ep.count = 1;
    row_ep.index[0] = row;
    row_ep.array[row] = 1;
    row_ep.packFlag = true;
    workHMO.factor_.btran(row_ep, analysis->row_ep_density);
    computeTableauRowFull(workHMO, row_ep, row_ap);
    for (int colCnt = 0; colCnt < row_ap.count; ++colCnt){
      if (!(row_ap.index[colCnt] >= workLp.numCol_ - workLp.numResiduals_)) continue;
      if (std::fabs(row_ap.array[row_ap.index[colCnt]]) < HIGHS_OCMATCHING_TINY) continue;
      // depMatValue.push_back(row_ap.array[row_ap.index[colCnt]]);
      depMatIndex.push_back(row_ap.index[colCnt]);
    }
    depMatStart.push_back(depMatIndex.size());
  }
}

void HQPrimal::matchingHeuristic(){
  int iRow, iCol, iR, innerCnt = 0, outerCnt = 0, pairCnt = 1;
  HighsLp& workLp = workHMO.lp_;
  // while(outerCnt != pairCnt){
    pairCnt = outerCnt;
    for (iRow = 0; iRow < workLp.numRow_ - workLp.numResiduals_; ++iRow){
      if (slackPaired[iRow]) continue;
      innerCnt = 0;
      for (iCol = depMatStart[iRow]; iCol < depMatStart[iRow + 1]; ++iCol){
        iR = depMatIndex[iCol];
        if (!rColPaired[iR] && !rInNeighborhood[iR] && !innerCnt){
          rColPaired[iR] = true;
          rColPairing[iR] = iRow;
          slackPaired[iRow] = true;
          slackPairing[iRow] = iR;
          innerCnt++;
        }
        else{
          rInNeighborhood[iR] = true;
        }
      }
      outerCnt += innerCnt;
    }
    depMatStart.clear();
    depMatIndex.clear();
    changeDegenSlack();
    buildNewBasis();
    buildDependencyMatrix();
    for (iCol = workLp.numCol_ - workLp.numResiduals_; iCol < workLp.numCol_; ++iCol)
      if (rInNeighborhood[iCol]) rInNeighborhood[iCol] = false;
  // }
}

void HQPrimal::changeDegenSlack(){
  int iRow;
  HighsLp& workLp = workHMO.lp_;
  for (iRow = 0; iRow < workLp.numRow_ - workLp.numResiduals_; ++iRow){
    if (slackPaired[iRow]) workLp.degenSlack_[iRow] = false;
  }
}

void HQPrimal::buildReduceLUFactor(){
  HighsTimer& timer = workHMO.timer_;
  HighsSimplexInfo& simplex_info = workHMO.simplex_info_;
  spSubMatRowPerm.assign(workHMO.lp_.numRow_, -1);
  spSubMatRowUnperm.assign(workHMO.lp_.numRow_, -1);
  spSubMatColPerm.assign(workHMO.lp_.numCol_, -1);
  spSubMatColUnperm.assign(workHMO.lp_.numCol_, -1);
  degenRowPair.assign(workHMO.lp_.numRow_, -1);
  rColPair.assign(workHMO.lp_.numCol_, -1);
  zeros.resize(workHMO.lp_.numS_);
  std::vector<double> temp(workHMO.lp_.numResiduals_, 0);
  for (int i = 0; i < workHMO.lp_.numResiduals_; ++i){
    printBinvAr.push_back(temp);
  }
  rReduceAstart.push_back(0);
  rZeroStart.push_back(0);
  rReduceNumRow = rReduceNumCol = workHMO.lp_.numResiduals_;
  for (int i = 0; i < workHMO.lp_.numResiduals_; ++i){
    computeReduceResiduals(workHMO.lp_.residuals_[i]);
  }
  // // std::cout << "nnz for r sub matrix: " << rReduceNnz << std::endl;
  // // std::cin.get();
  // timer.start(simplex_info.clock_[InvertPreClock]);
  // rReduceFactor.setup(rReduceNumCol, rReduceNumRow, &rReduceAstart[0],
  //                &rReduceAindex[0], &rReduceAvalue[0],
  //                &rSwapBasis[0], &workHMO.simplex_analysis_);
  // rReduceRankDeficiency = rReduceFactor.build();
  // timer.stop(simplex_info.clock_[InvertPreClock]);
  // // std::cout << "Time to do R sub Matrix Factor: " << timer.clock_time[simplex_info.clock_[InvertPreClock]] << std::endl;
  // // std::cin.get();
  // rReduceRank = rReduceNumRow - rReduceRankDeficiency;
  // rReduceColPerm = rReduceFactor.unPermute;
  // // print B^-1 * A_r for matlab analysis
  std::cout << "[" << std::endl;
  for (int i = 0; i < rReduceNumRow; ++i){
    for (int j = 0; j < rReduceNumCol; ++j){
      std::cout << printBinvAr[i][j] << " ";
    }
    std::cout << "; " << std::endl;
  }
  std::cout << "]" << std::endl;
}

void HQPrimal::buildGeneralQR(){
  HighsTimer& timer = workHMO.timer_;
  HighsSimplexInfo& simplex_info = workHMO.simplex_info_;
  spSubMatRowPerm.assign(workHMO.lp_.numRow_, -1);
  spSubMatRowUnperm.assign(workHMO.lp_.numRow_, -1);
  spSubMatColPerm.assign(workHMO.lp_.numCol_, -1);
  spSubMatColUnperm.assign(workHMO.lp_.numCol_, -1);
  for (int i = 0; i < workHMO.lp_.numResiduals_; ++i){
    computeReduceResiduals(workHMO.lp_.residuals_[i]);
  }
  std::cout << "Compute Res Col Clock: " << timer.clock_time[simplex_info.clock_[FtranPreClock]] << std::endl;
  std::cin.get();
  timer.start(timer.presolve_qr_clock);
  spSubMat.resize(rowPerms + 1, colPerms + 1);
  spSubMat.setFromTriplets(spTriplets.begin(), spTriplets.end());
  QR.analyzePattern(spSubMat);
  QR.factorize(spSubMat);
  timer.stop(timer.presolve_qr_clock);
  std::cout << "Compute QR Factor Clock: " << timer.clock_time[timer.presolve_qr_clock] << std::endl;
  std::cin.get();
  spSubMatQRPerm = QR.colsPermutation().indices().data();
  spSubMatRank = QR.rank(), spSubMatCols = QR.cols();
}

void HQPrimal::buildDependencyGraph(int slack){
  depInd.push_back(slack);
}

void HQPrimal::computeReduceResiduals(int iCol){
  HighsSimplexInfo& simplex_info = workHMO.simplex_info_;
  HighsTimer& timer = workHMO.timer_;
  HighsBasis& workBasis = workHMO.basis_;
  HighsLp& workLp = workHMO.lp_;
  bool degen;
  int i, iRow;
  double tol0, tol1, tolLB = 1e-6;
  std::fill(zeros.begin(), zeros.end(), 1);
  col_aq.clear();
  col_aq.packFlag = true;
  workHMO.matrix_.collect_aj(col_aq, iCol, 1);
  timer.start(simplex_info.clock_[FtranPreClock]);
  workHMO.factor_.ftran(col_aq, analysis->col_aq_density);
  timer.stop(simplex_info.clock_[FtranPreClock]);
  for (i = 0; i < col_aq.count; ++i){
    // if (rColPair[iCol]) break;
    iRow = col_aq.index[i];
    degen = (workLp.degenSlack_[iRow]);
    tol1 = std::fabs(col_aq.array[iRow] - 0);
    if (degen && tol1 > tolLB)
      zeros[iRow] = 0;
  }
  for (i = 0; i < col_aq.count; ++i){
    iRow = col_aq.index[i];
    degen = (workLp.degenSlack_[iRow]);
    tol0 = std::fabs(col_aq.array[iRow] - 0);
    if (degen && tol0 > tolLB){
      if (spSubMatRowPerm[iRow] == -1){
        spSubMatRowPerm[iRow] = ++rowPerms;
        spSubMatRowUnperm[rowPerms] = iRow;
      }
      if (spSubMatColPerm[iCol] == -1){
        spSubMatColPerm[iCol] = ++colPerms;
        spSubMatColUnperm[colPerms] = iCol;
      }
      rReduceAindex.push_back(spSubMatRowPerm[iRow]);
      rReduceAvalue.push_back(col_aq.array[iRow]);
      printBinvAr[spSubMatRowPerm[iRow]][spSubMatColPerm[iCol]] = col_aq.array[iRow];
      ++rReduceNnz;
      // spTriplets.push_back(Eigen::Triplet<double>(spSubMatRowPerm[iRow], spSubMatColPerm[iCol], col_aq.array[iRow]));
    }
  }
  rReduceAstart.push_back(rReduceNnz);
  rSwapBasis.push_back(spSubMatColPerm[iCol]);
  for (iRow = 0; iRow < workHMO.lp_.numS_; ++iRow)
    if (zeros[iRow] && workHMO.lp_.degenSlack_[iRow])
      rZeroIndex.push_back(iRow);
  rZeroStart.push_back(rZeroIndex.size());
}

void HQPrimal::buildZeroMatching(){
  int iCol, iRow, i, j;
  for (i = 0; i < workHMO.lp_.numResiduals_; ++i){
    iCol = i + workHMO.lp_.numX_;
    for (j = rZeroStart[i]; j < rZeroStart[i + 1]; ++j){
      iRow = rZeroIndex[j];
      if (degenRowPair[iRow] > -1) continue;
      rColPair[iCol] = iRow;
      degenRowPair[iRow] = iCol;
      numPairs++;
      break;
    }
  }
}

void HQPrimal::buildMaximumMatching(){
  hK.setUp(workHMO.lp_.numResiduals_, workHMO.lp_.numRow_ - workHMO.lp_.numResiduals_,
                  depInd, depStart);
  swap = hK.hopKarp();
}

void HQPrimal::buildNewBasis(){
  int residualStart = workHMO.lp_.numCol_ - workHMO.lp_.numResiduals_, res;
  for (int i = residualStart; i < rColPairing.size(); ++i){
    if (rColPairing[i] > -1){
      columnOut = rColPairing[i];
      columnIn = slackPairing[columnOut];
      workHMO.basis_.row_status[columnOut] = HighsBasisStatus::NONBASIC;
      workHMO.basis_.col_status[columnIn] = HighsBasisStatus::BASIC;
      workHMO.lp_.colLower_[columnIn] = -HIGHS_CONST_INF;
      workHMO.lp_.colUpper_[columnIn] = HIGHS_CONST_INF;
    }
  }
  workHMO.simplex_lp_status_.has_basis = false;
  // workHMO.basis_.col_status.assign(workHMO.basis_.col_status.size(), HighsBasisStatus::NONBASIC);
  transition(workHMO);
}

void HQPrimal::buildNewQRBasis(){
  HighsSimplexInfo& simplex_info = workHMO.simplex_info_;
  HighsTimer& timer = workHMO.timer_;
  int iColPerm, iRowPerm, iColUnperm, iRowUnperm;
  int iRank, numResiduals = workHMO.lp_.numResiduals_;
  int numCol = workHMO.lp_.numCol_ - workHMO.lp_.numResiduals_;
  HighsBasisStatus basic = HighsBasisStatus::BASIC;
  HighsBasisStatus lower = HighsBasisStatus::LOWER;
  rColSwapped.resize(numResiduals);
  for (iRank = 0; iRank < spSubMatRank; ++iRank){
    iColPerm = spSubMatQRPerm[iRank];
    iColUnperm = spSubMatColUnperm[iColPerm];
    iRowPerm = iRank;
    iRowUnperm = spSubMatRowUnperm[iRowPerm]; 
    workHMO.basis_.col_status[iColUnperm] = basic;
    workHMO.basis_.row_status[iRowUnperm] = lower;
    workHMO.lp_.colLower_[iColUnperm] = -HIGHS_CONST_INF;
    workHMO.lp_.colUpper_[iColUnperm] = HIGHS_CONST_INF;
    rColSwapped[iColUnperm - numCol] = 1;
  }
  workHMO.simplex_lp_status_.valid = false;
  workHMO.simplex_lp_status_.has_basis = false;
  transition(workHMO);
  std::cout << "Total LU Factor Clock Before Unfold: " << timer.clock_time[simplex_info.clock_[InvertClock]] << std::endl;
  std::cin.get();
}

void HQPrimal::buildNewLUBasis(){
  HighsSimplexInfo& simplex_info = workHMO.simplex_info_;
  HighsTimer& timer = workHMO.timer_;
  int iColPerm, iRowPerm, iRow, iColUnperm, iRowUnperm;
  int iSwap, numResiduals = workHMO.lp_.numResiduals_;
  int numCol = workHMO.lp_.numX_;
  HighsBasisStatus basic = HighsBasisStatus::BASIC;
  HighsBasisStatus lower = HighsBasisStatus::LOWER;
  rColSwapped.resize(numResiduals);
  // for (iRank = 0; iRank < rReduceRank; ++iRank){
  //   iColPerm = rReduceColPerm[iRank];
  //   iColUnperm = spSubMatColUnperm[iColPerm];
  //   iRowPerm = iRank;
  //   iRowUnperm = spSubMatRowUnperm[iRowPerm]; 
  //   workHMO.basis_.col_status[iColUnperm] = basic;
  //   workHMO.basis_.row_status[iRowUnperm] = lower;
  //   workHMO.lp_.colLower_[iColUnperm] = -HIGHS_CONST_INF;
  //   workHMO.lp_.colUpper_[iColUnperm] = HIGHS_CONST_INF;
  //   rColSwapped[iColUnperm - numCol] = 1;
  // }
  iSwap = 0;
  iRow = 0;
  while (iSwap < numPairs){
    iColUnperm = degenRowPair[iRow++];
    if (iColUnperm == -1) continue;
    iRowUnperm = rColPair[iColUnperm];
    workHMO.basis_.col_status[iColUnperm] = basic;
    workHMO.basis_.row_status[iRowUnperm] = lower;
    workHMO.lp_.colLower_[iColUnperm] = -HIGHS_CONST_INF;
    workHMO.lp_.colUpper_[iColUnperm] = HIGHS_CONST_INF;
    rColSwapped[iColUnperm - numCol] = 1;
    ++iSwap;
  }
  workHMO.simplex_lp_status_.valid = false;
  workHMO.simplex_lp_status_.has_basis = false;
  transition(workHMO);
}

void HQPrimal::appendLeftOverResiduals(){
  HighsLp& workLp = workHMO.lp_;
  int resStart = workLp.numCol_ - workLp.numResiduals_;
  int iCol;
  workLp.numResiduals_ = 0;
  workLp.residuals_.clear();
  for (iCol = resStart; iCol < workLp.numCol_; ++iCol){
    if (!rColPaired[iCol]){
      workLp.residuals_.push_back(iCol);
      workLp.numResiduals_++;
    }
  }
}

void HQPrimal::buildInitialFactor(){
  HighsSimplexInfo& simplex_info = workHMO.simplex_info_;
  HighsTimer& timer = workHMO.timer_;
  double init;
  timer.start(simplex_info.clock_[InvertClock]);
    init = timer.readRunHighsClock();
    int rankDeficiency = compute_factor(workHMO);
    invertTime += timer.readRunHighsClock() - init;
    timer.stop(simplex_info.clock_[InvertClock]);
    if (rankDeficiency) {
      throw runtime_error("Primal reInvert: singular-basis-matrix");
    }
    simplex_info.update_count = 0;
}

void HQPrimal::updateSimplexLp(){
  HighsLp &simplexLp = workHMO.simplex_lp_;
  HighsLp &modelLp = workHMO.lp_;
  SimplexBasis &simplexBasis = workHMO.simplex_basis_;
  // Shrink the simplex lp
  simplexLp.numCol_ = modelLp.numCol_ - numDrop;
  simplexLp.numRow_ = modelLp.numRow_ - numDrop;
  simplexLp.colLower_.resize(simplexLp.numCol_);
  simplexLp.colUpper_.resize(simplexLp.numCol_);
  simplexLp.rowLower_.resize(simplexLp.numRow_);
  simplexLp.rowUpper_.resize(simplexLp.numRow_);
  simplexLp.colCost_.resize(simplexLp.numCol_);
  std::vector<int> tempStart(1);
  std::vector<int> tempIndex;
  std::vector<int> tempValue;
  int iCol, iRow, rowVal, j, numEl = 0;
  for (iCol = 0; iCol < simplexLp.numCol_; ++i){
    for (j = simplexLp.Astart_[iCol]; j < simplexLp.Astart_[iCol + 1]; ++j){
      iRow = simplexLp.Aindex_[j];
      rowVal = simplexLp.Avalue_[j];
      if (iRow >= simplexLp.numRow_){
        tempIndex.push_back(iRow);
        tempValue.push_back(rowVal);
      }
    }
    tempStart.push_back(tempValue.size());
  }
  simplexLp.Astart_ = tempStart;
  simplexLp.Aindex_ = tempIndex;
  simplexLp.Avalue_ = tempValue;
  // Shrink simplex basis
  simplexBasis.basicIndex_.resize(simplexLp.numRow_);
  int sCol;
  int iPut = simplexLp.numCol_;
  for (sCol = solver_num_col; sCol < simplexLp.numRow_; ++sCol){
    simplexBasis.nonbasicFlag_[iPut] = simplexBasis.nonbasicFlag_[sCol];
    simplexBasis.nonbasicMove_[iPut++] = simplexBasis.nonbasicMove_[sCol];
  }
}

void HQPrimal::updateHighsBasis(){
  int iCol, iRow;
  std::fill(workHMO.basis_.col_status.begin(), workHMO.basis_.col_status.end(), HighsBasisStatus::NONBASIC);
  std::fill(workHMO.basis_.row_status.begin(), workHMO.basis_.row_status.end(), HighsBasisStatus::NONBASIC);
  for (int iRow = 0; iRow < solver_num_row; ++iRow){
    int basicCol = workHMO.simplex_basis_.basicIndex_[iRow];
    if (basicCol < workHMO.lp_.numCol_)
      workHMO.basis_.col_status[basicCol] = HighsBasisStatus::BASIC;
    else if (basicCol < solver_num_col)
      continue;
    else if (basicCol < solver_num_col + workHMO.lp_.numRow_)
      workHMO.basis_.row_status[basicCol - solver_num_col] = HighsBasisStatus::BASIC;
    else
      continue;
  }
  int numX = 0, numR = 0, numS = 0;
  for (iRow = 0; iRow < solver_num_row; ++iRow){
    int basicCol = workHMO.simplex_basis_.basicIndex_[iRow];
    if (basicCol < workHMO.lp_.numX_)
      numX++;
    else if (basicCol < solver_num_col)
      numR++;
    else numS++;
  }
  std::cout << "num x basic: " << numX << std::endl;
  std::cout << "num r basic: " << numR << std::endl;
  std::cout << "num slack basic: " << numS << std::endl;
  // std::cin.get();
}

void HQPrimal::updateSolver(){
  HighsSimplexInfo& simplex_info = workHMO.simplex_info_;
  HighsSimplexLpStatus& simplex_lp_status = workHMO.simplex_lp_status_;
  HighsTimer& timer = workHMO.timer_;
  printf("HQPrimal::solvePhase3 --- REMOVING BASIC R COLS/ROWS\n");
  printf("R COLS/ROWS REMOVED THUS FAR: %d\n", workHMO.lp_.unfoldIter);
  // When starting a new phase the (updated) primal objective function
  // value isn't known. Indicate this so that when the value
  // computed from scratch in build() isn't checked against the the
  // updated value
  simplex_lp_status.has_primal_objective_value = 0;
  // Set invertHint so that it's assigned when first tested
  invertHint = INVERT_HINT_NO;
  // Set solvePhase=2 so it's set if solvePhase2() is called directly
  solvePhase = 3;
  // Set up local copies of model dimensions
  solver_num_col = workHMO.simplex_lp_.numCol_;
  solver_num_row = workHMO.simplex_lp_.numRow_;
  solver_num_tot = solver_num_col + solver_num_row;

  analysis = &workHMO.simplex_analysis_;

  // Setup update limits
  simplex_info.update_limit =
      min(100 + solver_num_row / 100,
          1000);  // TODO: Consider allowing the dual limit to be used
  simplex_info.update_count = 0;

  // Setup local vectors
  col_aq.setup(solver_num_row);
  row_ep.setup(solver_num_row);
  row_ap.setup(solver_num_col);

  // Setup local dep graph vectors
  depStart.resize(workHMO.lp_.numResiduals_ + 1);

  ph1SorterR.reserve(solver_num_row);
  ph1SorterT.reserve(solver_num_row);

#ifdef HiGHSDEV
  printf("HQPrimal::solvePhase3 - WARNING: Setting analysis->col_aq_density = 0\n");
  printf("HQPrimal::solvePhase3 - WARNING: Setting analysis->row_ep_density = 0\n");
#endif
  analysis->col_aq_density = 0;
  analysis->row_ep_density = 0;

  devexReset();

  no_free_columns = true;
  for (int iCol = 0; iCol < solver_num_tot; iCol++) {
    if (highs_isInfinity(-workHMO.simplex_info_.workLower_[iCol])) {
      if (highs_isInfinity(workHMO.simplex_info_.workUpper_[iCol])) {
        // Free column
        no_free_columns = false;
        break;
      }
    }
  }
  #ifdef HiGHSDEV
    if (no_free_columns) {
      printf("Model has no free columns\n");
    } else {
      printf("Model has free columns\n");
    }
  #endif
}

void HQPrimal::updateAMatrix(){
  int iRow, iCol, idx, nnz = 0;
  std::vector<int> tempStart, tempIndex;
  std::vector<double> tempValue;
  tempStart.push_back(0);
  for (iCol = 0; iCol < workHMO.lp_.numCol_; ++iCol){
    for (idx = workHMO.lp_.Astart_[iCol]; idx < workHMO.lp_.Astart_[iCol + 1]; ++idx){
      iRow = workHMO.lp_.Aindex_[idx];
      if (iRow < workHMO.lp_.numRow_){
        ++nnz;
        tempIndex.push_back(workHMO.lp_.Aindex_[idx]);
        tempValue.push_back(workHMO.lp_.Avalue_[idx]);
      }
    }
    tempStart.push_back(nnz);
  }
  workHMO.lp_.Astart_ = tempStart;
  workHMO.lp_.Aindex_ = tempIndex;
  workHMO.lp_.Avalue_ = tempValue;
}

void HQPrimal::unfold() {
  ++workHMO.lp_.masterIter;
  double current_run_time = 0;
  HighsSimplexInfo& simplex_info = workHMO.simplex_info_;
  HighsTimer& timer = workHMO.timer_;
  bool run_highs_clock_already_running = timer.runningRunHighsClock();
  if (!run_highs_clock_already_running) timer.startRunHighsClock();
  // skipRow.assign(workHMO.lp_.numRow_,false);
  // rColPaired.assign(workHMO.lp_.numCol_, false);
  // rColPairing.assign(workHMO.lp_.numCol_, -1);
  // rInNeighborhood.assign(workHMO.lp_.numCol_, false);
  // slackPaired.assign(workHMO.lp_.numRow_, false);
  // slackPairing.assign(workHMO.lp_.numRow_, -1);
  // buildGeneralQR();
  // buildNewQRBasis();
  // buildReduceLUFactor();
  // buildZeroMatching();
  // buildNewLUBasis();
  // rColSwapped.assign(workHMO.lp_.numResiduals_, 0);
  // HighsSimplexInterface workInterface(workHMO);
  originalNumCol = workHMO.lp_.numCol_;
  originalNumRow = workHMO.lp_.numRow_;
  primalRebuild();
  int idx = workHMO.lp_.numCol_ - workHMO.lp_.numResiduals_;
  int offset = 0;
  int cnt = 0;
  int tPivCnt = 0;
  numDrop = 0;
  int update_limit = min(100 + workHMO.lp_.numRow_ / 100,
          1000);
  // update_limit = 1200;
  double chooseRowTime = 0;
  int fromCol, toCol, fromRow, toRow;
  bool fromColBool = false;
  ////////////////////// Orbital Crossover //////////////////////
  for (int i = workHMO.lp_.numResiduals_ - 1; i >= 0; --i){
  // for (int i = 0; i < workHMO.lp_.numResiduals_; ++i){
    // columnIn = workHMO.lp_.residuals_[i];
    // // Use this when we do linker swapping before main unfold
    // offset = workHMO.lp_.residuals_[i] - workHMO.lp_.numX_;
    // if (workHMO.lp_.basicResiduals_[offset]) continue;
    ++cnt;
    ++tPivCnt;
    ++numDrop;
    current_run_time = timer.readRunHighsClock();
    // degeneratePivot = false;
    if (current_run_time + workHMO.lp_.foldSolveTime > workHMO.options_.time_limit){
      workHMO.scaled_model_status_ = HighsModelStatus::REACHED_TIME_LIMIT;
      return;
    }
    ++workHMO.lp_.unfoldIter;
    timer.start(simplex_info.clock_[ChuzcPrimalClock]);
    columnIn = workHMO.lp_.residuals_[i];
    timer.stop(simplex_info.clock_[ChuzcPrimalClock]);
    if (!fromColBool){
      fromCol = columnIn;
      fromColBool = true;
      fromRow = workHMO.lp_.numS_ + i;
    }
    // Edit work and bounds for residual columns
    workHMO.simplex_info_.workCost_[columnIn] = 1;
    workHMO.simplex_info_.workUpper_[columnIn] = +HIGHS_CONST_INF;
    workHMO.simplex_info_.workLower_[columnIn] = -HIGHS_CONST_INF;
    workHMO.simplex_info_.workValue_[columnIn] = 0;
    workHMO.simplex_basis_.nonbasicMove_[columnIn] = 1;
    // Highs choose row function nothing changed from normal simplex
    double init = timer.readRunHighsClock();
    primalChooseRow();
    chooseRowTime += timer.readRunHighsClock() - init;
    // std::cout << "rowOut: " << rowOut << std::endl;
    // Highs primal update function nothing change from normal simplex
    primalUpdate();
    workHMO.simplex_info_.workCost_[columnIn] = 0;
    workHMO.lp_.colUpper_[columnIn] = +HIGHS_CONST_INF;
    workHMO.lp_.colLower_[columnIn] = -HIGHS_CONST_INF;
    // // Check for numerical instability in basis inverse
    // if (ivHint == INVERT_HINT_POSSIBLY_SINGULAR_BASIS){
    //   workHMO.lp_.numCol_ -= tPivCnt;
    //   workHMO.lp_.numRow_ -= tPivCnt;
    //   workHMO.basis_.numCol_ -= tPivCnt;
    //   workHMO.basis_.numRow_ -= tPivCnt;
    //   workHMO.simplex_lp_status_.valid = false;
    //   workHMO.simplex_lp_status_.has_basis = false;
    //   updateSolver();
    //   updateHighsBasis();
    //   updateAMatrix();
    //   transition(workHMO);
    //   primalRebuild(); 
    //   ivHint = 0;
    //   cnt = 0;
    //   tPivCnt = 0;
    //   chooseRowTime = 0;
    // }
  //   // Check if we need to refactor LU factorization based on pivots done
    if (cnt > update_limit){
      workHMO.lp_.numCol_ -= tPivCnt;
      workHMO.lp_.numRow_ -= tPivCnt;
      workHMO.basis_.numCol_ -= tPivCnt;
      workHMO.basis_.numRow_ -= tPivCnt;
      workHMO.simplex_lp_status_.valid = false;
      workHMO.simplex_lp_status_.has_basis = false;
      updateSolver();
      updateHighsBasis();
      updateAMatrix();
      transition(workHMO);
      primalRebuild();
      cnt = 0;
      tPivCnt = 0;
      // toCol = columnIn;
      // toRow = workHMO.lp_.numS_ + i;
      // workInterface.deleteCols(fromCol, toCol);
      // workInterface.deleteRows(fromRow, toRow);
      // fromColBool = false;
    }
  }
  ////////////////////// Orbital Crossover Done //////////////////////
  // Final Rebuild
  primalRebuild();
  // Clean up if we have dual infeas after all r pivots
  // std::sort(workHMO.simplex_basis_.basicIndex_.begin(), workHMO.simplex_basis_.basicIndex_.end());
  // for (int i = 0; i < workHMO.simplex_basis_.basicIndex_.size(); ++i){
  //   std::cout << workHMO.simplex_basis_.basicIndex_[i] << std::endl;
  // }
  if (workHMO.scaled_solution_params_.num_dual_infeasibilities > 0) solvePhase2();
}

void HQPrimal::swapDegenerate(){
  columnIn = -1;
  rowOut = -1;
  int numSwaps = 0;
  for (;;){
    primalChooseR();
    if (columnIn == workHMO.lp_.numCol_) break;
    // if (rSwapped[columnIn]) continue;
    primalChooseDegenerateSlack();
    if (rowOut == -1) continue;
    primalUpdate();
    if (++numSwaps == workHMO.lp_.numRow_) break;
  }
  primalRebuild();
}

void HQPrimal::primalChooseR(){
  ++columnIn;
}

void HQPrimal::primalChooseDegenerateSlack(){
  // std::vector<int>& As = workHMO.lp_.Astart_;
  // std::vector<int>& Ai = workHMO.lp_.Aindex_;
  // std::vector<double>& Av = workHMO.lp_.Avalue_;
  // int idx, iRow;
  double tol, tolUB = 1e-6;
  // for (idx = As[columnIn]; idx < As[columnIn + 1]; ++idx){
  //   tol = std::fabs(Av[idx] - 0);
  //   iRow = Ai[idx];
  //   if (tol > tolUB && !sSwapped[iRow]){
  //     rowOut = iRow;
  //     return;
  //   }
  // }
  rowOut = -1;
  // return;
  // Compute pivot column
  int iRow, iVal;
  HighsTimer& timer = workHMO.timer_;
  HighsSimplexInfo& simplex_info = workHMO.simplex_info_;
  timer.start(simplex_info.clock_[FtranClock]);
  double init = timer.readRunHighsClock();
  col_aq.clear();
  col_aq.packFlag = true;
  workHMO.matrix_.collect_aj(col_aq, columnIn, 1);
  workHMO.factor_.ftran(col_aq, analysis->col_aq_density);
  fTranTime += timer.readRunHighsClock() - init;
  timer.stop(simplex_info.clock_[FtranClock]);
  const double local_col_aq_density = (double)col_aq.count / solver_num_row;
  analysis->updateOperationResultDensity(local_col_aq_density, analysis->col_aq_density);
  // Pick Row
  for (int i = 0; i < col_aq.count; ++i){
    iRow = col_aq.index[i];
    iVal = col_aq.array[iRow];
    tol = std::fabs(iVal - 0);
    if (tol > tolUB && !sSwapped[iRow]){
      rowOut = iRow;
      sSwapped[iRow] = 1;
      return;
    }
  }
  return;
}

void HQPrimal::primalChooseColumn() {
  HighsRandom& random = workHMO.random_;
  HighsSimplexInfo& simplex_info = workHMO.simplex_info_;
  HighsTimer& timer = workHMO.timer_;
  const int* jFlag = &workHMO.simplex_basis_.nonbasicFlag_[0];
  const int* jMove = &workHMO.simplex_basis_.nonbasicMove_[0];
  double* workDual = &workHMO.simplex_info_.workDual_[0];
  const double* workLower = &workHMO.simplex_info_.workLower_[0];
  const double* workUpper = &workHMO.simplex_info_.workUpper_[0];
  const double dualTolerance = workHMO.scaled_solution_params_.dual_feasibility_tolerance;

  timer.start(simplex_info.clock_[ChuzcPrimalClock]);
  columnIn = -1;
  double bestInfeas = 0;
  if (no_free_columns) {
    const int numSection = 1;
    int startSection = random.integer() % numSection;
    int deltaCol = (solver_num_tot + numSection - 1) / numSection;
    int fromCol = startSection * deltaCol;
    int toCol = min(fromCol + deltaCol, solver_num_tot);
    int numPass = 1;
    //    printf("\nstartSection = %1d; deltaCol = %d\n", startSection,
    //    deltaCol);
    for (;;) {
      //      printf("CHUZC: %1d [%6d, %6d] %6d\n", numPass, fromCol, toCol,
      //      solver_num_tot);
      for (int iCol = fromCol; iCol < toCol; iCol++) {
        // Then look at dual infeasible
        if (jMove[iCol] * workDual[iCol] < -dualTolerance) {
          if (bestInfeas * devex_weight[iCol] < fabs(workDual[iCol])) {
            bestInfeas = fabs(workDual[iCol]) / devex_weight[iCol];
            columnIn = iCol;
          }
        }
      }
      if (columnIn >= 0 || numPass == numSection) {
        //	printf("Break from CHUZC after %d passes\n", numPass);
        break;
      }
      if (toCol == solver_num_tot) {
        fromCol = 0;
        toCol = deltaCol;
      } else {
        fromCol = toCol;
        toCol = min(fromCol + deltaCol, solver_num_tot);
      }
      numPass++;
    }
  } else {
    for (int iCol = 0; iCol < solver_num_tot; iCol++) {
      if (jFlag[iCol] && fabs(workDual[iCol]) > dualTolerance) {
        // Always take free
        // TODO: if we found free,
        // Then deal with it in dual phase 1
        if (workLower[iCol] == -HIGHS_CONST_INF &&
            workUpper[iCol] == HIGHS_CONST_INF) {
          columnIn = iCol;
          break;
        }
        // Then look at dual infeasible
        if (jMove[iCol] * workDual[iCol] < -dualTolerance) {
          if (bestInfeas * devex_weight[iCol]  < fabs(workDual[iCol])) {
            bestInfeas = fabs(workDual[iCol]) / devex_weight[iCol];
            columnIn = iCol;
          }
        }
      }
    }
  }
  timer.stop(simplex_info.clock_[ChuzcPrimalClock]);
}

void HQPrimal::primalChooseRow() {
  HighsSimplexInfo& simplex_info = workHMO.simplex_info_;
  HighsTimer& timer = workHMO.timer_;
  const double* baseLower = &workHMO.simplex_info_.baseLower_[0];
  const double* baseUpper = &workHMO.simplex_info_.baseUpper_[0];
  double* baseValue = &workHMO.simplex_info_.baseValue_[0];
  const double primalTolerance =
      workHMO.scaled_solution_params_.primal_feasibility_tolerance;

  // Compute pivot column
  timer.start(simplex_info.clock_[FtranClock]);
  double init = timer.readRunHighsClock();
  col_aq.clear();
  col_aq.packFlag = true;
  workHMO.matrix_.collect_aj(col_aq, columnIn, 1);
#ifdef HiGHSDEV
  if (simplex_info.analyse_iterations)
    analysis->operationRecordBefore(ANALYSIS_OPERATION_TYPE_FTRAN, col_aq, analysis->col_aq_density);
#endif
  workHMO.factor_.ftran(col_aq, analysis->col_aq_density);
  fTranTime += timer.readRunHighsClock() - init;
  timer.stop(simplex_info.clock_[FtranClock]);
#ifdef HiGHSDEV
  if (simplex_info.analyse_iterations)
    analysis->operationRecordAfter(ANALYSIS_OPERATION_TYPE_FTRAN, col_aq);
#endif

  const double local_col_aq_density = (double)col_aq.count / solver_num_row;
  analysis->updateOperationResultDensity(local_col_aq_density, analysis->col_aq_density);

  const bool check_dual = false;
  if (check_dual) {
    const double* workCost = &workHMO.simplex_info_.workCost_[0];
    const double* workDual = &workHMO.simplex_info_.workDual_[0];
    const int* basicIndex = &workHMO.simplex_basis_.basicIndex_[0];
    double check_dual_value = workCost[columnIn];
    for (int i = 0; i < col_aq.count; i++) {
      int row = col_aq.index[i];
      int col = basicIndex[row];
      double value = col_aq.array[row];
      double cost = workCost[col];
      check_dual_value -= value * cost;
      //    printf("Entry %2d: [%2d, %12g] Cost = %12g; check_dual_value =
      //    %12g\n", i, row, value, cost, check_dual_value);
    }
    thetaDual = workDual[columnIn];
    double dual_error =
        fabs(check_dual_value - thetaDual) / max(1.0, fabs(thetaDual));
    if (dual_error > 1e-8)
      printf("Checking dual: updated = %12g; direct = %12g; error = %12g\n",
             thetaDual, check_dual_value, dual_error);
  }

  timer.start(simplex_info.clock_[Chuzr1Clock]);
  init = timer.readRunHighsClock();
  // Initialize
  rowOut = -1;

  // Choose row pass 1
  double alphaTol = workHMO.simplex_info_.update_count < 10
                        ? 1e-9
                        : workHMO.simplex_info_.update_count < 20 ? 1e-8 : 1e-7;
  const int* jMove = &workHMO.simplex_basis_.nonbasicMove_[0];
  int moveIn = jMove[columnIn];
  if (moveIn == 0) {
    // If there's still free in the N
    // We would report not-solved
    // Need to handle free
  }
  double relaxTheta = 1e100;
  double relaxSpace;
  for (int i = 0; i < col_aq.count; i++) {
    int index = col_aq.index[i];
    // if (skipRow[index]) continue;
    alpha = col_aq.array[index] * moveIn;
    if (alpha > alphaTol) {
      relaxSpace = baseValue[index] - baseLower[index] + primalTolerance;
      if (relaxSpace < relaxTheta * alpha) relaxTheta = relaxSpace / alpha;
    } else if (alpha < -alphaTol) {
      relaxSpace = baseValue[index] - baseUpper[index] - primalTolerance;
      if (relaxSpace > relaxTheta * alpha) relaxTheta = relaxSpace / alpha;
    }
  }
  Chuzc1Time += timer.readRunHighsClock() - init;
  timer.stop(simplex_info.clock_[Chuzr1Clock]);

  timer.start(simplex_info.clock_[Chuzr2Clock]);
  init = timer.readRunHighsClock();
  double bestAlpha = 0;
  for (int i = 0; i < col_aq.count; i++) {
    int index = col_aq.index[i];
    // if (skipRow[index]) continue;
    alpha = col_aq.array[index] * moveIn;
    if (alpha > alphaTol) {
      // Positive pivotal column entry
      double tightSpace = baseValue[index] - baseLower[index];
      if (tightSpace < relaxTheta * alpha) {
        if (bestAlpha < alpha) {
          bestAlpha = alpha;
          rowOut = index;
        }
      }
    } else if (alpha < -alphaTol) {
      // Negative pivotal column entry
      double tightSpace = baseValue[index] - baseUpper[index];
      if (tightSpace > relaxTheta * alpha) {
        if (bestAlpha < -alpha) {
          bestAlpha = -alpha;
          rowOut = index;
        }
      }
    }
  }
  degeneratePivot = std::fabs(relaxTheta) < 1e-6;
  // if (degeneratePivot) skipRow[rowOut] = true;
  Chuzc2Time += timer.readRunHighsClock() - init;
  timer.stop(simplex_info.clock_[Chuzr2Clock]);
  // std::cout << "row out: " << rowOut << std::endl;
}

void HQPrimal::primalUpdate() {
  
  HighsTimer& timer = workHMO.timer_;
  int* jMove = &workHMO.simplex_basis_.nonbasicMove_[0];
  double* workDual = &workHMO.simplex_info_.workDual_[0];
  const double* workLower = &workHMO.simplex_info_.workLower_[0];
  const double* workUpper = &workHMO.simplex_info_.workUpper_[0];
  const double* baseLower = &workHMO.simplex_info_.baseLower_[0];
  const double* baseUpper = &workHMO.simplex_info_.baseUpper_[0];
  double* workValue = &workHMO.simplex_info_.workValue_[0];
  double* baseValue = &workHMO.simplex_info_.baseValue_[0];
  const double primalTolerance =
      workHMO.scaled_solution_params_.primal_feasibility_tolerance;
  HighsSimplexInfo& simplex_info = workHMO.simplex_info_;

  // Compute thetaPrimal
  int moveIn = jMove[columnIn];
  //  int
  columnOut = workHMO.simplex_basis_.basicIndex_[rowOut];
  //  double
  alpha = col_aq.array[rowOut];
  //  double
  thetaPrimal = 0;
  if (alpha * moveIn > 0) {
    // Lower bound
    thetaPrimal = (baseValue[rowOut] - baseLower[rowOut]) / alpha;
  } else {
    // Upper bound
    thetaPrimal = (baseValue[rowOut] - baseUpper[rowOut]) / alpha;
  }

  // 1. Make sure it is inside bounds or just flip bound
  double lowerIn = workLower[columnIn];
  double upperIn = workUpper[columnIn];
  double valueIn = workValue[columnIn] + thetaPrimal;
  bool flipped = false;
  if (jMove[columnIn] == 1) {
    if (valueIn > upperIn + primalTolerance) {
      // Flip to upper
      workValue[columnIn] = upperIn;
      thetaPrimal = upperIn - lowerIn;
      flipped = true;
      jMove[columnIn] = -1;
    }
  } else if (jMove[columnIn] == -1) {
    if (valueIn < lowerIn - primalTolerance) {
      // Flip to lower
      workValue[columnIn] = lowerIn;
      thetaPrimal = lowerIn - upperIn;
      flipped = true;
      jMove[columnIn] = 1;
    }
  }

  timer.start(simplex_info.clock_[UpdatePrimalClock]);
  double init = timer.readRunHighsClock();
  for (int i = 0; i < col_aq.count; i++) {
    int index = col_aq.index[i];
    baseValue[index] -= thetaPrimal * col_aq.array[index];
  }
  updatePrimalTime += timer.readRunHighsClock() - init;
  timer.stop(simplex_info.clock_[UpdatePrimalClock]);

  simplex_info.updated_primal_objective_value +=
      workDual[columnIn] * thetaPrimal;

  computePrimalInfeasible(workHMO);

  // If flipped, then no need touch the pivots
  if (flipped) {
    rowOut = -1;
    numericalTrouble = 0;
    thetaDual = workDual[columnIn];
    iterationAnalysis();
    num_flip_since_rebuild++;
    return;
  }

  // Pivot in
  int sourceOut = alpha * moveIn > 0 ? -1 : 1;
  update_pivots(workHMO, columnIn, rowOut, sourceOut);

  baseValue[rowOut] = valueIn;

  timer.start(simplex_info.clock_[CollectPrIfsClock]);
  init = timer.readRunHighsClock();
  // Check for any possible infeasible
  for (int iRow = 0; iRow < solver_num_row; iRow++) {
    if (baseValue[iRow] < baseLower[iRow] - primalTolerance) {
      invertHint = INVERT_HINT_PRIMAL_INFEASIBLE_IN_PRIMAL_SIMPLEX;
    } else if (baseValue[iRow] > baseUpper[iRow] + primalTolerance) {
      invertHint = INVERT_HINT_PRIMAL_INFEASIBLE_IN_PRIMAL_SIMPLEX;
    }
  }
  collectPrimalInfsTime += timer.readRunHighsClock() - init;
  timer.stop(simplex_info.clock_[CollectPrIfsClock]);

  // 2. Now we can update the dual

  timer.start(simplex_info.clock_[BtranClock]);
  init = timer.readRunHighsClock();
  row_ep.clear();
  row_ap.clear();
  row_ep.count = 1;
  row_ep.index[0] = rowOut;
  row_ep.array[rowOut] = 1;
  row_ep.packFlag = true;
#ifdef HiGHSDEV
  if (simplex_info.analyse_iterations)
    analysis->operationRecordBefore(ANALYSIS_OPERATION_TYPE_BTRAN_EP, row_ep, analysis->row_ep_density);
#endif
  workHMO.factor_.btran(row_ep, analysis->row_ep_density);
#ifdef HiGHSDEV
  if (simplex_info.analyse_iterations)
    analysis->operationRecordAfter(ANALYSIS_OPERATION_TYPE_BTRAN_EP, row_ep);
#endif
  bTranTime += timer.readRunHighsClock() - init;
  timer.stop(simplex_info.clock_[BtranClock]);
  //
  // PRICE
  //
  computeTableauRowFromPiP(workHMO, row_ep, row_ap);
  /*
#ifdef HiGHSDEV
  if (simplex_info.analyse_iterations) {
    analysis->operationRecordBefore(ANALYSIS_OPERATION_TYPE_PRICE_AP, row_ep, analysis->row_ap_density);
    analysis->num_row_price++;
  }
#endif
  timer.start(simplex_info.clock_[PriceClock]);
  workHMO.matrix_.priceByRowSparseResult(row_ap, row_ep);
  timer.stop(simplex_info.clock_[PriceClock]);
#ifdef HiGHSDEV
  if (simplex_info.analyse_iterations)
    analysis->operationRecordAfter(ANALYSIS_OPERATION_TYPE_PRICE_AP, row_ep);
#endif

  const double local_row_ep_density = (double)row_ep.count / solver_num_row;
  analysis->updateOperationResultDensity(local_row_ep_density, analysis->row_ep_density);
  */
  timer.start(simplex_info.clock_[UpdateDualClock]);
  init = timer.readRunHighsClock();
  //  double
  thetaDual = workDual[columnIn] / alpha;
  for (int i = 0; i < row_ap.count; i++) {
    int iCol = row_ap.index[i];
    workDual[iCol] -= thetaDual * row_ap.array[iCol];
  }
  for (int i = 0; i < row_ep.count; i++) {
    int iGet = row_ep.index[i];
    int iCol = iGet + solver_num_col;
    workDual[iCol] -= thetaDual * row_ep.array[iGet];
  }
  updateDualTime += timer.readRunHighsClock() - init;
  timer.stop(simplex_info.clock_[UpdateDualClock]);

  /* Update the devex weight */
  devexUpdate();

  // After dual update in primal simplex the dual objective value is not known
  workHMO.simplex_lp_status_.has_dual_objective_value = false;

  // updateVerify for primal
  numericalTrouble = 0;
  double aCol = fabs(alpha);
  double alphaRow;
  if (columnIn < workHMO.simplex_lp_.numCol_) {
    alphaRow = row_ap.array[columnIn];
  } else {
    alphaRow = row_ep.array[rowOut];
  }
  double aRow = fabs(alphaRow);
  double aDiff = fabs(aCol - aRow);
  numericalTrouble = aDiff / min(aCol, aRow);
  if (numericalTrouble > 1e-7 && solvePhase != 4)
    printf("NumericalTrouble - Reinverting\n");
  // Reinvert if the relative difference is large enough, and updates have been
  //performed
  if (numericalTrouble > 1e-7 && workHMO.simplex_info_.update_count > 0){
    invertHint = INVERT_HINT_POSSIBLY_SINGULAR_BASIS;
    ivHint = INVERT_HINT_POSSIBLY_SINGULAR_BASIS;
    // std::cout << "invertHint: " << (int)invertHint << std::endl;
    // // std::cout << "update_count: " << workHMO.simplex_info_.update_count << std::endl;
    // std::cin.get();
  }
  // Dual for the pivot
  workDual[columnIn] = 0;
  workDual[columnOut] = -thetaDual;

  // Update workHMO.factor_ basis
  init = timer.readRunHighsClock();
  update_factor(workHMO, &col_aq, &row_ep, &rowOut, &invertHint);
  // std::cout << "invertHint: " << (int)invertHint << std::endl;
  // std::cin.get();
  update_matrix(workHMO, columnIn, columnOut);
  if (simplex_info.update_count >= simplex_info.update_limit) {
    invertHint = INVERT_HINT_UPDATE_LIMIT_REACHED;
  }
  updateFactorTime += timer.readRunHighsClock() - init;
  // std::cout << "INVERT_HINT_UPDATE_LIMIT_REACHED: " << (int)INVERT_HINT_UPDATE_LIMIT_REACHED << std::endl;
  // std::cout << "UPDATE LIMIT: " << (int)simplex_info.update_limit << std::endl;
  // std::cout << "UPDATE COUNT: " << (int)simplex_info.update_count << std::endl;
  // std::cout << "INVERT HINT: " << (int)invertHint << std::endl;
  // std::cin.get();
  // Move this to Simplex class once it's created
  // simplex_method.record_pivots(columnIn, columnOut, alpha);
  workHMO.scaled_solution_params_.simplex_iteration_count++;

  /* Reset the devex when there are too many errors */
  if(num_bad_devex_weight > 3) {
    devexReset();
  }

  // Report on the iteration
  iterationAnalysis();
}

/* Compute the reduced cost for primal phase 1 with artificial cost. */
void HQPrimal::phase1ComputeDual() {
  /* Alias to problem size, tolerance and work arrays */
  const int nRow = workHMO.lp_.numRow_;
  const int nCol = workHMO.lp_.numCol_;
  const double dFeasTol = workHMO.scaled_solution_params_.primal_feasibility_tolerance;
  const double *baseLower = &workHMO.simplex_info_.baseLower_[0];
  const double *baseUpper = &workHMO.simplex_info_.baseUpper_[0];
  const double *baseValue = &workHMO.simplex_info_.baseValue_[0];

  HighsSimplexInfo& simplex_info = workHMO.simplex_info_;
  HighsTimer& timer = workHMO.timer_;
  timer.start(simplex_info.clock_[BtranClock]);
  /* Setup artificial cost and compute pi with BTran */
  HVector buffer;
  buffer.setup(nRow);
  buffer.clear();
  for (int iRow = 0; iRow < nRow; iRow++) {
    buffer.index[iRow] = iRow;
    if (baseValue[iRow] <  baseLower[iRow] - dFeasTol) {
      buffer.array[iRow] = -1.0;
    } else if (baseValue[iRow] > baseUpper[iRow] + dFeasTol) {
      buffer.array[iRow] = 1.0;
    } else {
      buffer.array[iRow] = 0.0;
    }
  }
#ifdef HiGHSDEV
  if (simplex_info.analyse_iterations)
    analysis->operationRecordBefore(ANALYSIS_OPERATION_TYPE_BTRAN_EP, buffer, analysis->row_ep_density);
#endif
  workHMO.factor_.btran(buffer, 1);
#ifdef HiGHSDEV
  if (simplex_info.analyse_iterations)
    analysis->operationRecordAfter(ANALYSIS_OPERATION_TYPE_BTRAN_EP, buffer);
#endif
  timer.stop(simplex_info.clock_[BtranClock]);

  timer.start(simplex_info.clock_[PriceClock]);
  /* The compute the phase 1 reduced cost for all variables by SpMV */
  HVector bufferLong;
  bufferLong.setup(nCol);
  bufferLong.clear();
#ifdef HiGHSDEV
    if (simplex_info.analyse_iterations) {
      analysis->operationRecordBefore(ANALYSIS_OPERATION_TYPE_PRICE_AP, buffer, 0.0);
      analysis->num_col_price++;
    }
#endif
  workHMO.matrix_.priceByColumn(bufferLong, buffer);
#ifdef HiGHSDEV
  if (simplex_info.analyse_iterations)
    analysis->operationRecordAfter(ANALYSIS_OPERATION_TYPE_PRICE_AP, row_ap);
#endif
  timer.stop(simplex_info.clock_[PriceClock]);

  const int* nbFlag = &workHMO.simplex_basis_.nonbasicFlag_[0];
  double *workDual = &workHMO.simplex_info_.workDual_[0];
  for (int iSeq = 0; iSeq < nCol + nRow; iSeq++) {
    workDual[iSeq] = 0.0;
  }
  for (int iSeq = 0; iSeq < nCol; iSeq++) {
    if (nbFlag[iSeq])
      workDual[iSeq] = -bufferLong.array[iSeq];
  }
  for (int iRow = 0, iSeq = nCol; iRow < nRow; iRow++, iSeq++) {
    if (nbFlag[iSeq])
      workDual[iSeq] = -buffer.array[iRow];
  }

  /* Recompute number of dual infeasible variables with the phase 1 cost */
  computeDualInfeasible(workHMO);
}

/* Choose a pivot column for the phase 1 primal simplex method */
void HQPrimal::phase1ChooseColumn() {
  const int nSeq = workHMO.lp_.numCol_ + workHMO.lp_.numRow_;
  const int* nbMove = &workHMO.simplex_basis_.nonbasicMove_[0];
  const double* workDual = &workHMO.simplex_info_.workDual_[0];
  const double dDualTol = workHMO.scaled_solution_params_.dual_feasibility_tolerance;
  HighsSimplexInfo& simplex_info = workHMO.simplex_info_;
  HighsTimer& timer = workHMO.timer_;
  timer.start(simplex_info.clock_[ChuzcPrimalClock]);
  double dBestScore = 0;
  columnIn = -1;
  for (int iSeq = 0; iSeq < nSeq; iSeq++) {
    double dMyDual = nbMove[iSeq] * workDual[iSeq];
    double dMyScore = dMyDual / devex_weight[iSeq];
    if (dMyDual < -dDualTol && dMyScore < dBestScore) {
      dBestScore = dMyScore;
      columnIn = iSeq;
    }
  }
  timer.stop(simplex_info.clock_[ChuzcPrimalClock]);
}

/* Choose a pivot row for the phase 1 primal simplex method */
void HQPrimal::phase1ChooseRow() {
  /* Alias to work arrays */
  const double dFeasTol = workHMO.scaled_solution_params_.primal_feasibility_tolerance;
  const double* baseLower = &workHMO.simplex_info_.baseLower_[0];
  const double* baseUpper = &workHMO.simplex_info_.baseUpper_[0];
  const double* baseValue = &workHMO.simplex_info_.baseValue_[0];

  /* Compute the transformed pivot column and update its density */
  HighsSimplexInfo& simplex_info = workHMO.simplex_info_;
  HighsTimer& timer = workHMO.timer_;
  timer.start(simplex_info.clock_[FtranClock]);
  col_aq.clear();
  col_aq.packFlag = true;
  workHMO.matrix_.collect_aj(col_aq, columnIn, 1);
#ifdef HiGHSDEV
  if (simplex_info.analyse_iterations)
    analysis->operationRecordBefore(ANALYSIS_OPERATION_TYPE_FTRAN, col_aq, analysis->col_aq_density);
#endif
  workHMO.factor_.ftran(col_aq, analysis->col_aq_density);
  timer.stop(simplex_info.clock_[FtranClock]);
#ifdef HiGHSDEV
  if (simplex_info.analyse_iterations)
    analysis->operationRecordAfter(ANALYSIS_OPERATION_TYPE_FTRAN, col_aq);
#endif

  const double local_col_aq_density = (double)col_aq.count / solver_num_row;
  analysis->updateOperationResultDensity(local_col_aq_density, analysis->col_aq_density);

  /* Compute the reduced cost for the pivot column and compare it with the kept value */
  double dCompDual = 0.0;
  for (int i = 0; i < col_aq.count; i++) {
    int iRow = col_aq.index[i];
    if (baseValue[iRow] < baseLower[iRow] - dFeasTol) {
      dCompDual -= col_aq.array[iRow] * -1.0;
    } else if (baseValue[iRow] > baseUpper[iRow] + dFeasTol) {
      dCompDual -= col_aq.array[iRow] * +1.0;
    }
  }
  if (fabs(workHMO.simplex_info_.workDual_[columnIn] - dCompDual) > (fabs(dCompDual) + 1.0) * 1e-9) {
    printf("==> Phase 1 reduced cost. Updated %g, Computed %g\n", workHMO.simplex_info_.workDual_[columnIn], dCompDual);
  }

  timer.start(simplex_info.clock_[Chuzr1Clock]);
  /* Collect phase 1 theta lists */
  int nRow = workHMO.lp_.numRow_;
  const int iMoveIn = workHMO.simplex_basis_.nonbasicMove_[columnIn];
  const double dPivotTol = workHMO.simplex_info_.update_count < 10 ? 1e-9 :
                           workHMO.simplex_info_.update_count < 20 ? 1e-8 : 1e-7;
  ph1SorterR.clear();
  ph1SorterT.clear();
  for (int i = 0; i < col_aq.count; i++) {
    int iRow = col_aq.index[i];
    double dAlpha = col_aq.array[iRow] * iMoveIn;

    /* When the basic variable x[i] decrease */
    if (dAlpha > +dPivotTol) {
      /* Whether it can become feasible by going below its upper bound */
      if (baseValue[iRow] > baseUpper[iRow] + dFeasTol) {
        double dFeasTheta = (baseValue[iRow] - baseUpper[iRow] - dFeasTol) / dAlpha;
        ph1SorterR.push_back(std::make_pair(dFeasTheta, iRow));
        ph1SorterT.push_back(std::make_pair(dFeasTheta, iRow));
      }
      /* Whether it can become infeasible (again) by going below its lower bound */
      if (baseValue[iRow] > baseLower[iRow] - dFeasTol && baseLower[iRow] > -HIGHS_CONST_INF) {
        double dRelaxTheta = (baseValue[iRow] - baseLower[iRow] + dFeasTol) / dAlpha;
        double dTightTheta = (baseValue[iRow] - baseLower[iRow]) / dAlpha;
        ph1SorterR.push_back(std::make_pair(dRelaxTheta, iRow - nRow));
        ph1SorterT.push_back(std::make_pair(dTightTheta, iRow - nRow));
      }
    }

    /* When the basic variable x[i] increase */
    if (dAlpha < -dPivotTol) {
      /* Whether it can become feasible by going above its lower bound */
      if (baseValue[iRow] < baseLower[iRow] - dFeasTol) {
        double dFeasTheta = (baseValue[iRow] - baseLower[iRow] + dFeasTol) / dAlpha;
        ph1SorterR.push_back(std::make_pair(dFeasTheta, iRow - nRow));
        ph1SorterT.push_back(std::make_pair(dFeasTheta, iRow - nRow));
      }

      /* Whether it can become infeasible (again) by going above its upper bound */
      if (baseValue[iRow] < baseUpper[iRow] + dFeasTol && baseUpper[iRow] < +HIGHS_CONST_INF) {
        double dRelaxTheta = (baseValue[iRow] - baseUpper[iRow] - dFeasTol) / dAlpha;
        double dTightTheta = (baseValue[iRow] - baseUpper[iRow]) / dAlpha;
        ph1SorterR.push_back(std::make_pair(dRelaxTheta, iRow));
        ph1SorterT.push_back(std::make_pair(dTightTheta, iRow));
      }
    }
  }

  timer.stop(simplex_info.clock_[Chuzr1Clock]);
  /* When there is no candidates at all, we can leave it here */
  if (ph1SorterR.empty()) {
    rowOut = -1;
    columnOut = -1;
    return;
  }

  /*
   * Now sort the relaxed theta to find the final break point.
   * TODO: Consider partial sort.
   *       Or heapify [O(n)] and then pop k points [kO(log(n))].
   */
  timer.start(simplex_info.clock_[Chuzr2Clock]);
  std::sort(ph1SorterR.begin(), ph1SorterR.end());
  double dMaxTheta = ph1SorterR.at(0).first;
  double dGradient = fabs(workHMO.simplex_info_.workDual_[columnIn]);
  for (unsigned int i = 0; i < ph1SorterR.size(); i++) {
    double dMyTheta = ph1SorterR.at(i).first;
    int index = ph1SorterR.at(i).second;
    int iRow = index >= 0 ? index : index + nRow;
    dGradient -= fabs(col_aq.array[iRow]);
    /* Stop when the gradient start to decrease */
    if (dGradient <= 0) {
      break;
    }
    dMaxTheta = dMyTheta;
  }

  /* Find out the biggest possible alpha for pivot */
  std::sort(ph1SorterT.begin(), ph1SorterT.end());
  double dMaxAlpha = 0.0;
  unsigned int iLast = ph1SorterT.size();
  for (unsigned int i = 0; i < ph1SorterT.size(); i++) {
    double dMyTheta = ph1SorterT.at(i).first;
    int index = ph1SorterT.at(i).second;
    int iRow = index >= 0 ? index : index + nRow;
    double dAbsAlpha = fabs(col_aq.array[iRow]);
    /* Stop when the theta is too large */
    if (dMyTheta > dMaxTheta) {
      iLast = i;
      break;
    }
    /* Update the maximal possible alpha */
    if (dMaxAlpha < dAbsAlpha) {
      dMaxAlpha = dAbsAlpha;
    }
  }

  /* Finally choose a pivot with good enough alpha, backwardly */
  rowOut = -1;
  columnOut = -1;
  phase1OutBnd = 0;
  for (int i = iLast - 1; i >= 0; i--) {
    int index = ph1SorterT.at(i).second;
    int iRow = index >= 0 ? index : index + nRow;
    double dAbsAlpha = fabs(col_aq.array[iRow]);
    if (dAbsAlpha > dMaxAlpha * 0.1) {
      rowOut = iRow;
      phase1OutBnd = index > 0 ? 1 : -1;
      break;
    }
  }
  if(rowOut != -1) {
    columnOut = workHMO.simplex_basis_.basicIndex_[rowOut];
  }
  timer.stop(simplex_info.clock_[Chuzr2Clock]);
}

/* Update the primal and dual solutions for the primal phase 1 */
void HQPrimal::phase1Update() {
  /* Alias to bounds arrays */
  const double* workLower = &workHMO.simplex_info_.workLower_[0];
  const double* workUpper = &workHMO.simplex_info_.workUpper_[0];
  const double* baseLower = &workHMO.simplex_info_.baseLower_[0];
  const double* baseUpper = &workHMO.simplex_info_.baseUpper_[0];
  double* workValue = &workHMO.simplex_info_.workValue_[0];
  double* baseValue = &workHMO.simplex_info_.baseValue_[0];
  const int iMoveIn = workHMO.simplex_basis_.nonbasicMove_[columnIn];
  const double dFeasTol = workHMO.scaled_solution_params_.primal_feasibility_tolerance;
  HighsSimplexInfo& simplex_info = workHMO.simplex_info_;
  HighsTimer& timer = workHMO.timer_;

  /* Compute the primal theta and see if we should have do bound flip instead */
  alpha = col_aq.array[rowOut];
  thetaPrimal = 0.0;
  if(phase1OutBnd == 1) {
    thetaPrimal = (baseValue[rowOut] - baseUpper[rowOut]) / alpha;
  } else {
    thetaPrimal = (baseValue[rowOut] - baseLower[rowOut]) / alpha;
  }

  double lowerIn = workLower[columnIn];
  double upperIn = workUpper[columnIn];
  double valueIn = workValue[columnIn] + thetaPrimal;
  int ifFlip = 0;
  if (iMoveIn == +1 && valueIn > upperIn + dFeasTol) {
    workValue[columnIn] = upperIn;
    thetaPrimal = upperIn - lowerIn;
    ifFlip = 1;
    workHMO.simplex_basis_.nonbasicMove_[columnIn] = -1;
  }
  if (iMoveIn == -1 && valueIn < lowerIn - dFeasTol) {
    workValue[columnIn] = lowerIn;
    thetaPrimal = lowerIn - upperIn;
    ifFlip = 1;
    workHMO.simplex_basis_.nonbasicMove_[columnIn] = +1;
  }

  /* Update for the flip case */
  if(ifFlip) {
    /* Recompute things on flip */
    if (invertHint == 0) {
      timer.start(simplex_info.clock_[ComputePrimalClock]);
      compute_primal(workHMO);
      timer.stop(simplex_info.clock_[ComputePrimalClock]);
      computePrimalInfeasible(workHMO);
      if (workHMO.scaled_solution_params_.num_primal_infeasibilities > 0) {
        isPrimalPhase1 = 1;
	timer.start(simplex_info.clock_[ComputeDualClock]);
        phase1ComputeDual();
	timer.stop(simplex_info.clock_[ComputeDualClock]);
      } else {
        invertHint = INVERT_HINT_UPDATE_LIMIT_REACHED;
      }
    }
    return;
  }

  /* Compute BTran for update LU */
  timer.start(simplex_info.clock_[BtranClock]);
  row_ep.clear();
  row_ep.count = 1;
  row_ep.index[0] = rowOut;
  row_ep.array[rowOut] = 1;
  row_ep.packFlag = true;
#ifdef HiGHSDEV
  if (simplex_info.analyse_iterations)
    analysis->operationRecordBefore(ANALYSIS_OPERATION_TYPE_BTRAN_EP, row_ep, analysis->row_ep_density);
#endif
  workHMO.factor_.btran(row_ep, analysis->row_ep_density);
#ifdef HiGHSDEV
  if (simplex_info.analyse_iterations)
    analysis->operationRecordAfter(ANALYSIS_OPERATION_TYPE_BTRAN_EP, row_ep);
#endif
  timer.stop(simplex_info.clock_[BtranClock]);

  const double local_row_ep_density = (double)row_ep.count / solver_num_row;
  analysis->updateOperationResultDensity(local_row_ep_density, analysis->row_ep_density);

  /* Compute the whole pivot row for updating the devex weight */
#ifdef HiGHSDEV
  if (simplex_info.analyse_iterations) {
    analysis->operationRecordBefore(ANALYSIS_OPERATION_TYPE_PRICE_AP, row_ep, analysis->row_ap_density);
    analysis->num_row_price++;
  }
#endif
  timer.start(simplex_info.clock_[PriceClock]);
  row_ap.clear();
  workHMO.matrix_.priceByRowSparseResult(row_ap, row_ep);
  timer.stop(simplex_info.clock_[PriceClock]);
#ifdef HiGHSDEV
  if (simplex_info.analyse_iterations)
    analysis->operationRecordAfter(ANALYSIS_OPERATION_TYPE_PRICE_AP, row_ep);
#endif

  /* Update the devex weight */
  devexUpdate();

   /* Update other things */
  update_pivots(workHMO, columnIn, rowOut, phase1OutBnd);
  update_factor(workHMO, &col_aq, &row_ep, &rowOut, &invertHint);
  update_matrix(workHMO, columnIn, columnOut);
  if (workHMO.simplex_info_.update_count >= workHMO.simplex_info_.update_limit) {
    invertHint = INVERT_HINT_UPDATE_LIMIT_REACHED;
  }

  /* Recompute dual and primal */
  if (invertHint == 0) {
    timer.start(simplex_info.clock_[ComputePrimalClock]);
    compute_primal(workHMO);
    timer.stop(simplex_info.clock_[ComputePrimalClock]);
    computePrimalInfeasible(workHMO);

    if (workHMO.scaled_solution_params_.num_primal_infeasibilities > 0) {
      isPrimalPhase1 = 1;
      timer.start(simplex_info.clock_[ComputeDualClock]);
      phase1ComputeDual();
      timer.stop(simplex_info.clock_[ComputeDualClock]);
    } else {
      invertHint = INVERT_HINT_UPDATE_LIMIT_REACHED;
    }
  }

  /* Reset the devex framework when necessary */
  if(num_bad_devex_weight > 3) {
    devexReset();
  }


  // Move this to Simplex class once it's created
  // simplex_method.record_pivots(columnIn, columnOut, alpha);
  workHMO.scaled_solution_params_.simplex_iteration_count++;
}

/* Reset the devex weight */
void HQPrimal::devexReset() {
  const int nSeq = workHMO.lp_.numCol_ + workHMO.lp_.numRow_;
  devex_weight.assign(nSeq, 1.0);
  devex_index.assign(nSeq, 0);
  for (int iSeq = 0; iSeq < nSeq; iSeq++) {
    const int nonbasicFlag = workHMO.simplex_basis_.nonbasicFlag_[iSeq];
    devex_index[iSeq] = nonbasicFlag*nonbasicFlag;
  }
  num_devex_iterations = 0;
  num_bad_devex_weight = 0;
}

void HQPrimal::devexUpdate() {
  /* Compute the pivot weight from the reference set */
  HighsSimplexInfo& simplex_info = workHMO.simplex_info_;
  HighsTimer& timer = workHMO.timer_;
  timer.start(simplex_info.clock_[DevexUpdateWeightClock]);
  double dPivotWeight = 0.0;
  for (int i = 0; i < col_aq.count; i++) {
    int iRow = col_aq.index[i];
    int iSeq = workHMO.simplex_basis_.basicIndex_[iRow];
    double dAlpha = devex_index[iSeq] * col_aq.array[iRow];
    dPivotWeight += dAlpha * dAlpha;
  }
  dPivotWeight += devex_index[columnIn] * 1.0;
  dPivotWeight = sqrt(dPivotWeight);

  /* Check if the saved weight is too large */
  if (devex_weight[columnIn] > 3.0 * dPivotWeight) {
    num_bad_devex_weight++;
  }

  /* Update the devex weight for all */
  double dPivot = col_aq.array[rowOut];
  dPivotWeight /= fabs(dPivot);

  for (int i = 0; i < row_ap.count; i++) {
    int iSeq = row_ap.index[i];
    double alpha = row_ap.array[iSeq];
    double devex = dPivotWeight * fabs(alpha);
    devex += devex_index[iSeq] * 1.0;
    if (devex_weight[iSeq] < devex) {
      devex_weight[iSeq] = devex;
    }
  }
  for (int i = 0; i < row_ep.count; i++) {
    int iPtr = row_ep.index[i];
    int iSeq = row_ep.index[i] + solver_num_col;
    double alpha = row_ep.array[iPtr];
    double devex = dPivotWeight * fabs(alpha);
    devex += devex_index[iSeq] * 1.0;
    if (devex_weight[iSeq] < devex) {
      devex_weight[iSeq] = devex;
    }
  }

  /* Update devex weight for the pivots */
  devex_weight[columnOut] = max(1.0, dPivotWeight);
  devex_weight[columnIn] = 1.0;
  num_devex_iterations++;
  timer.stop(simplex_info.clock_[DevexUpdateWeightClock]);
}

void HQPrimal::iterationAnalysisData() {
  HighsSolutionParams& scaled_solution_params = workHMO.scaled_solution_params_;
  HighsSimplexInfo& simplex_info = workHMO.simplex_info_;
  analysis->simplex_strategy = SIMPLEX_STRATEGY_PRIMAL;
  analysis->edge_weight_mode = DualEdgeWeightMode::DEVEX;
  analysis->solve_phase = solvePhase;
  analysis->simplex_iteration_count = scaled_solution_params.simplex_iteration_count;
  analysis->devex_iteration_count = num_devex_iterations;
  analysis->pivotal_row_index = rowOut;
  analysis->leaving_variable = columnOut;
  analysis->entering_variable = columnIn;
  analysis->invert_hint = invertHint;
  analysis->freelist_size = 0;
  analysis->reduced_rhs_value = 0;
  analysis->reduced_cost_value = 0;
  analysis->edge_weight = 0;
  analysis->primal_delta = 0;
  analysis->primal_step = thetaPrimal;
  analysis->dual_step = thetaDual;
  analysis->pivot_value_from_column = alpha;
  analysis->pivot_value_from_row = alpha;//Row;
  analysis->numerical_trouble = numericalTrouble;
  analysis->objective_value = simplex_info.updated_primal_objective_value;
  analysis->num_primal_infeasibilities = scaled_solution_params.num_primal_infeasibilities;
  analysis->num_dual_infeasibilities = scaled_solution_params.num_dual_infeasibilities;
  analysis->sum_primal_infeasibilities = scaled_solution_params.sum_primal_infeasibilities;
  analysis->sum_dual_infeasibilities = scaled_solution_params.sum_dual_infeasibilities;
#ifdef HiGHSDEV
  analysis->basis_condition = simplex_info.invert_condition;
#endif
  if ((analysis->edge_weight_mode == DualEdgeWeightMode::DEVEX) &&
      (num_devex_iterations == 0)) analysis->num_devex_framework++;
}

void HQPrimal::iterationAnalysis() {
  // Possibly report on the iteration
  iterationAnalysisData();
  analysis->iterationReport();

#ifdef HiGHSDEV
  analysis->iterationRecord();
#endif
}

void HQPrimal::reportRebuild(const int rebuild_invert_hint) {
  iterationAnalysisData();
  analysis->invert_hint = rebuild_invert_hint;
  analysis->invertReport();
}

