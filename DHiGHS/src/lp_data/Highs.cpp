/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License

              */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/Highs.cpp
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier

 */
#include "Highs.h"

#include <algorithm>
#include <iostream>
#include <memory>
#include <sstream>
#include <filesystem>
#include <iostream>

#include "HConfig.h"
#include "io/Filereader.h"
#include "io/HighsIO.h"
#include "io/LoadOptions.h"
#include "lp_data/HighsLpUtils.h"
#include "lp_data/HighsModelUtils.h"
#include "simplex/HApp.h"
#include "simplex/HighsSimplexInterface.h"
#include "equitable/HighsEquitable.h"

#ifdef OPENMP
#include "omp.h"
#endif

// until add_row_.. functions are moved to HighsLpUtils.h
#include "simplex/HSimplex.h"

#ifdef IPX_ON
#include "ipm/IpxWrapper.h"
#else
#include "ipm/IpxWrapperEmpty.h"
#endif

Highs::Highs() {
  hmos_.clear();
  HighsModelObject* hmo = new HighsModelObject(lp_, options_, timer_);
  hmos_.push_back(*hmo);
}

HighsStatus Highs::setHighsOptionValue(const std::string& option,
                                       const bool value) {
  if (setOptionValue(options_.logfile, option, options_.records, value) == OptionStatus::OK)
    return HighsStatus::OK;
  return HighsStatus::Error;
}

HighsStatus Highs::setHighsOptionValue(const std::string& option,
                                       const int value) {
  if (setOptionValue(options_.logfile, option, options_.records, value) == OptionStatus::OK)
    return HighsStatus::OK;
  return HighsStatus::Error;
}

HighsStatus Highs::setHighsOptionValue(const std::string& option,
                                       const double value) {
  if (setOptionValue(options_.logfile, option, options_.records, value) == OptionStatus::OK)
    return HighsStatus::OK;
  return HighsStatus::Error;
}

HighsStatus Highs::setHighsOptionValue(const std::string& option,
                                       const std::string value) {
  if (setOptionValue(options_.logfile, option, options_.records, value) == OptionStatus::OK)
    return HighsStatus::OK;
  return HighsStatus::Error;
}

HighsStatus Highs::setHighsOptionValue(const std::string& option,
                                       const char* value) {
  if (setOptionValue(options_.logfile, option, options_.records, value) == OptionStatus::OK)
    return HighsStatus::OK;
  return HighsStatus::Error;
}

HighsStatus Highs::setHighsLogfile(FILE* logfile) {
  options_.logfile = logfile;
  return HighsStatus::OK;
}

HighsStatus Highs::setHighsOutput(FILE* output) {
  options_.output = output;
  return HighsStatus::OK;
}

HighsStatus Highs::readHighsOptions(const std::string filename) {
  if (filename.size() <= 0) {
    HighsLogMessage(options_.logfile, HighsMessageType::WARNING,
                    "Empty file name so not reading options");
    return HighsStatus::Warning;
  }
  options_.options_file = filename;
  if (!loadOptionsFromFile(options_)) return HighsStatus::Error;
  return HighsStatus::OK;
}

HighsStatus Highs::passHighsOptions(const HighsOptions& options) {
  if (passOptions(options_.logfile, options, options_) == OptionStatus::OK){
    return HighsStatus::OK;
  }
  return HighsStatus::Error;
}

HighsStatus Highs::getHighsOptionValue(const std::string& option, bool& value) {
  if (getOptionValue(options_.logfile, option, options_.records, value) == OptionStatus::OK)
    return HighsStatus::OK;
  return HighsStatus::Error;
}

HighsStatus Highs::getHighsOptionValue(const std::string& option, int& value) {
  if (getOptionValue(options_.logfile, option, options_.records, value) == OptionStatus::OK)
    return HighsStatus::OK;
  return HighsStatus::Error;
}

HighsStatus Highs::getHighsOptionValue(const std::string& option,
                                       double& value) {
  if (getOptionValue(options_.logfile, option, options_.records, value) == OptionStatus::OK)
    return HighsStatus::OK;
  return HighsStatus::Error;
}

HighsStatus Highs::getHighsOptionValue(const std::string& option,
                                       std::string& value) {
  if (getOptionValue(options_.logfile, option, options_.records, value) == OptionStatus::OK)
    return HighsStatus::OK;
  return HighsStatus::Error;
}

HighsStatus Highs::writeHighsOptions(const std::string filename,
				     const bool report_only_non_default_values) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  HighsLp lp = this->lp_;
  FILE* file;
  bool html;
  call_status = openWriteFile(filename, "writeHighsOptions", file, html);
  return_status = interpretCallStatus(call_status, return_status, "openWriteFile");
  if (return_status == HighsStatus::Error) return return_status;

  call_status = writeOptionsToFile(file, options_.records,
				   report_only_non_default_values, html);
  return_status = interpretCallStatus(call_status, return_status, "writeOptionsToFile");
  return return_status;
}

const HighsOptions& Highs::getHighsOptions() const { return options_; }

const HighsInfo& Highs::getHighsInfo() const { return info_; }

HighsStatus Highs::getHighsInfoValue(const std::string& info, int& value) {
  if (getInfoValue(options_, info, info_.records, value) == InfoStatus::OK)
    return HighsStatus::OK;
  return HighsStatus::Error;
}

HighsStatus Highs::getHighsInfoValue(const std::string& info, double& value) const {
  if (getInfoValue(options_, info, info_.records, value) == InfoStatus::OK)
    return HighsStatus::OK;
  return HighsStatus::Error;
}

HighsStatus Highs::writeHighsInfo(const std::string filename) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  HighsLp lp = this->lp_;
  FILE* file;
  bool html;
  call_status = openWriteFile(filename, "writeHighsInfo", file, html);
  return_status = interpretCallStatus(call_status, return_status, "openWriteFile");
  if (return_status == HighsStatus::Error) return return_status;

  call_status = writeInfoToFile(file, info_.records, html);
  return_status = interpretCallStatus(call_status, return_status, "writeInfoToFile");
  return return_status;
}

HighsStatus Highs::passModel(const HighsLp& lp) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  // Copy the LP to the internal LP
  lp_ = lp;
  // Check validity of the LP, normalising its values (by default).
  call_status = assessLp(lp_, options_);
  return_status = interpretCallStatus(call_status, return_status, "assessLp");
  if (return_status == HighsStatus::Error) return return_status;
  hmos_.clear();
  hmos_.push_back(HighsModelObject(lp_, options_, timer_));
  return return_status;
}

HighsStatus Highs::readModel(const std::string filename) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  Filereader* reader = Filereader::getFilereader(filename.c_str());
  HighsLp model;
  this->options_.model_file = filename;

  FilereaderRetcode call_code = reader->readModelFromFile(this->options_, model);
  if (call_code != FilereaderRetcode::OK) {
    call_status = HighsStatus::Error;
    return_status = interpretCallStatus(call_status, return_status, "readModelFromFile");
    if (return_status == HighsStatus::Error) return return_status;
  }
  call_status = this->passModel(model);
  return_status = interpretCallStatus(call_status, return_status, "passModel");
  return return_status;
}

HighsStatus Highs::writeModel(const std::string filename) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  HighsLp model = this->lp_;

  if (filename == "") {
    // Empty file name: report model on stdout
    reportLp(options_, model, 2);
    return_status = HighsStatus::OK;
  } else {
    Filereader* writer = Filereader::getFilereader(filename.c_str());
    call_status = writer->writeModelToFile(options_, filename.c_str(), model);
    return_status = interpretCallStatus(call_status, return_status, "writeModelToFile");
  }
  return return_status;
}

// Checks the options calls presolve and postsolve if needed. Solvers are called
// with runLpSolver(..)
HighsStatus Highs::run() {
#ifdef OPENMP
  omp_max_threads = omp_get_max_threads();
  assert(omp_max_threads>0);
#ifdef HiGHSDEV
  if (omp_max_threads<=0)
    printf("WARNING: omp_get_max_threads() returns %d\n", omp_max_threads);
  printf("Running with %d OMP thread(s)\n", omp_max_threads);
#endif
  if (omp_max_threads < options_.highs_max_threads)
    HighsLogMessage(options_.logfile, HighsMessageType::WARNING,
		    "Number of OMP threads available = %d < %d = Max number of HiGHS threads requested: Parallel performance will be less than anticipated",
		    omp_max_threads, options_.highs_max_threads);
#endif
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
    /*
  if (options_.message_level >= 0) {
    printf("\n!! Actually solving an LP with %d cols, %d rows", lp_.numCol_, lp_.numRow_);
    if (lp_.numCol_) printf(" and %d nonzeros", lp_.Astart_[lp_.numCol_]);
    printf(":basis.valid_ = %d: basis_.valid_ = %d: simplex_lp_status_.has_basis = %d!!\n\n",
	   basis_.valid_,
	   hmos_[0].basis_.valid_,
	   hmos_[0].simplex_lp_status_.has_basis);
    if (basis_.valid_ != hmos_[0].basis_.valid_) {
      printf("NB %d = basis_.valid_ != hmos_[0].basis_.valid_ = %d\n", basis_.valid_, hmos_[0].basis_.valid_);
    }
  }
    */
  // If running as hsol, reset any changed options
  if (options_.run_as_hsol) setHsolOptions(options_);
  // Initialise the HiGHS model status values
  hmos_[0].scaled_model_status_ = HighsModelStatus::NOTSET;
  hmos_[0].unscaled_model_status_ = HighsModelStatus::NOTSET;
  model_status_ = hmos_[0].scaled_model_status_;
  scaled_model_status_ = hmos_[0].unscaled_model_status_;

#ifdef HIGHSDEV
  // Shouldn't have to check validity of the LP since this is done when it is
  // loaded or modified
  //  bool normalise = true;
  call_status = assessLp(lp_, options_);  //, normalise);
  assert(call_status == HighsStatus::OK);
  return_status = interpretCallStatus(call_status, return_status, "assessLp");
  if (return_status == HighsStatus::Error) return return_status;
#endif

  if (options_.icrash) {
    ICrashStrategy strategy = ICrashStrategy::kICA;
    bool strategy_ok = parseICrashStrategy(options_.icrash_strategy, strategy);
    if (!strategy_ok) {
      HighsPrintMessage(options_.output, options_.message_level, ML_ALWAYS,
			"ICrash error: unknown strategy.\n");
      return HighsStatus::Error;
    }
    ICrashOptions icrash_options{
        options_.icrash_dualize,
        strategy,
        options_.icrash_starting_weight,
        options_.icrash_iterations,
        options_.icrash_approximate_minimization_iterations,
        options_.icrash_exact,
	options_.icrash_breakpoints,
	options_.logfile,
	options_.output,
        options_.message_level};

    // todo: timing. some strange compile issue.
    HighsStatus icrash_status = callICrash(lp_, icrash_options, icrash_info_);
    return icrash_status;
  }

  // Return immediately if the LP has no columns
  if (!lp_.numCol_) {
    hmos_[0].unscaled_model_status_ = HighsModelStatus::MODEL_EMPTY;
    model_status_ = hmos_[0].unscaled_model_status_;
    return highsStatusFromHighsModelStatus(hmos_[0].unscaled_model_status_);
  }

  HighsSetIO(options_);
#ifdef HiGHSDEV
  if (checkOptions(options_.logfile, options_.records) != OptionStatus::OK) return HighsStatus::Error;
#endif
  HighsPrintMessage(options_.output, options_.message_level, ML_VERBOSE,
		    "Solving %s", lp_.model_name_.c_str());

  // Running as LP solver: start the HiGHS clock unless it's already running
  bool run_highs_clock_already_running = timer_.runningRunHighsClock();
  if (!run_highs_clock_already_running) timer_.startRunHighsClock();
  // Record the initial time and zero the overall iteration count
  double initial_time = timer_.readRunHighsClock();
  int postsolve_iteration_count = 0;
  // Define identifiers to refer to the HMO of the original LP
  // (0) and the HMO created when using presolve (1)
  const int original_hmo = 0;
  const int presolve_hmo = 1;
  // Keep track of the hmo that is the most recently solved. By default it's the
  // original LP
  int solved_hmo = original_hmo;
  // Initial solve. Presolve, choose solver (simplex, ipx), postsolve.
  //  printf("\nHighs::run() 1: basis_.valid_ = %d\n", basis_.valid_);
  //  fflush(stdout);
  // writeModel("../../DOKSmps/Original.lp");
  originalLp = lp_;
  // Count total number of r pivots that could be done
  int possibleRPivots = 0;
  // Count total number of actual r pivots done
  int actualRPivots = 0;
  // Count total pivots (ALP and all R pivots)
  int totalPivots = 0;
  // Count number of residuals swapped with pre-solve factorization
  int rSwapPivots = 0;
  if (options_.aggregate == on_string && options_.solver != ipm_sym_string){
    // Compute equitable partition
    timer_.start(timer_.equipart_clock);
    partitionLp();
    possibleRPivots = originalLp.numCol_ - part_->ncsplits;
    timer_.stop(timer_.equipart_clock);
    // Initial Aggregate
    timer_.start(timer_.fold_clock);
    foldLp();
    timer_.stop(timer_.fold_clock);
    // Run dual solver on ALP(A,b,c,P^0)
    passModel(*alp_);
    timer_.start(timer_.alp_solve_clock);
    call_status = runLpSolver(hmos_[solved_hmo], "Solving ALP");
    timer_.stop(timer_.alp_solve_clock);
    return_status = interpretCallStatus(call_status, return_status, "runLpSolver");
    totalPivots += hmos_[original_hmo].scaled_solution_params_.simplex_iteration_count;
    // Record basis and solution
    alpSolution_ = hmos_[original_hmo].solution_;
    alpBasis_ = hmos_[original_hmo].basis_;
    // // refine partition
    // timer_.start(timer_.equipart_clock);
    // refinePartitionFinal();
    // timer_.stop(timer_.equipart_clock);
    // // std::cout << "Refine done" << std::endl;
    // // std::cin.get();
    // // lift to next elp
    // timer_.start(timer_.lift_clock);
    // liftLpExtended();
    // // liftBasis();
    // timer_.stop(timer_.lift_clock);
    // int swaps = -1;
    // // Pass new elp
    // passModel(*alp_);
    // // writeModel("../../DOKSmps/Extended.lp");
    // setBasis(*startBasis_);
    // // std::cout << "set basis done" << std::endl;
    // // std::cin.get();
    // hmos_[solved_hmo].basis_ = basis_;
    // options_.simplex_strategy = SIMPLEX_STRATEGY_UNFOLD;
    // // std::string itercnt = std::string(cnt);
    // timer_.start(timer_.elp_solve_clock);
    // call_status = runLpSolver(hmos_[solved_hmo], "Solving ELP");
    // timer_.stop(timer_.elp_solve_clock);
    // return_status = interpretCallStatus(call_status, return_status, "runLpSolver"); 
    // actualRPivots += hmos_[original_hmo].scaled_solution_params_.simplex_iteration_count;
    // rSwapPivots = possibleRPivots - actualRPivots;
    // totalPivots += hmos_[original_hmo].scaled_solution_params_.simplex_iteration_count;
    // if (hmos_[original_hmo].readyForHighs){
    //   HighsModelObject* postOChmo = new HighsModelObject(originalLp, options_, timer_);
    //   int postOC_hmo = 1;
    //   hmos_.push_back(*postOChmo);
    //   trimOCSolution(hmos_[original_hmo], hmos_[postOC_hmo]);
    //   options_.solver = ipm_string;
    //   hmos_[postOC_hmo].readyForHighs = 1;
    //   runLpSolver(hmos_[postOC_hmo], "Cleaning OC with HighsCrossover");
    // }
    //////////////////////// Iterative Lifting /////////////////////////////////////////////
    // Start lift loop
    int cnt = 1;
    while (!discrete){
      // refine partition
      timer_.start(timer_.equipart_clock);
      refinePartition();
      timer_.stop(timer_.equipart_clock);
      // lift to next elp
      timer_.start(timer_.lift_clock);
      liftLpExtended();
      liftBasis();
      timer_.stop(timer_.lift_clock);
      // Pass new elp
      passModel(*elp_);
      setBasis(*startBasis_);
      hmos_[solved_hmo].basis_ = basis_;
      options_.simplex_strategy = SIMPLEX_STRATEGY_UNFOLD;
      std::string itercnt = std::to_string(cnt);
      std::string fName = "../../DOKSmps/ELP_" + itercnt + ".lp";
      writeModel(fName);
      timer_.start(timer_.elp_solve_clock);
      call_status = runLpSolver(hmos_[solved_hmo], "Solving ELP");
      timer_.stop(timer_.elp_solve_clock);
      return_status = interpretCallStatus(call_status, return_status, "runLpSolver"); 
      actualRPivots += hmos_[original_hmo].scaled_solution_params_.simplex_iteration_count;
      totalPivots += hmos_[original_hmo].scaled_solution_params_.simplex_iteration_count;
      alpSolution_ = hmos_[original_hmo].solution_;
      alpBasis_ = hmos_[original_hmo].basis_;
      ++cnt;
      if (hmos_[original_hmo].readyForHighs){
        HighsModelObject* postOChmo = new HighsModelObject(originalLp, options_, timer_);
        int postOC_hmo = 1;
        hmos_.push_back(*postOChmo);
        trimOCSolution(hmos_[original_hmo], hmos_[postOC_hmo]);
        options_.solver = ipm_string;
        hmos_[postOC_hmo].readyForHighs = 1;
        runLpSolver(hmos_[postOC_hmo], "Cleaning OC with HighsCrossover");
      }
    }
    // // Lift to elp
    // timer_.start(timer_.lift_clock);
    // liftLp(alpBasis_, alpSolution_);
    // timer_.stop(timer_.lift_clock);
    // // Do ELP UNFOLD
    // passModel(*elp_);
    // setBasis(*elpBasis_);
    // hmos_[solved_hmo].lp_.model_name_ = options_.model_file.c_str();
    // hmos_[solved_hmo].basis_ = basis_;
    // options_.simplex_strategy = SIMPLEX_STRATEGY_UNFOLD;
    // timer_.start(timer_.elp_solve_clock);
    // call_status = runLpSolver(hmos_[solved_hmo], "Solving ELP");
    // timer_.stop(timer_.elp_solve_clock);
    if (return_status == HighsStatus::Error) return return_status;
  }
  else if (options_.aggregate == on_string && options_.solver == ipm_sym_string){
    // Compute equitable parititon from original lp
    timer_.start(timer_.equipart_clock);
    partitionLp();
    timer_.stop(timer_.equipart_clock);
    // Fold original lp
    timer_.start(timer_.fold_clock);
    foldLp();
    timer_.stop(timer_.fold_clock);
    // Do alp solve 
    passModel(*alp_);
    timer_.start(timer_.alp_solve_clock);
    call_status = runLpSolver(hmos_[solved_hmo], "Solving ALP");
    timer_.stop(timer_.alp_solve_clock);
    return_status = interpretCallStatus(call_status, return_status, "runLpSolver");
  }
  else if (!basis_.valid_ && options_.presolve != off_string) {
    // cout << "here" << endl;
    // cin.get();
    // No HiGHS basis so consider presolve
    hmos_[original_hmo].scaled_model_status_ = HighsModelStatus::NOTSET;
    // Presolve. runPresolve handles the level of presolving (0 = don't
    // presolve)p
    timer_.start(timer_.presolve_clock);
    PresolveInfo presolve_info(options_.presolve, lp_);
    HighsPresolveStatus presolve_status = runPresolve(presolve_info);
    timer_.stop(timer_.presolve_clock);
    //    printf("\nHighs::run() 2: presolve status = %d\n",
    //    (int)presolve_status);fflush(stdout);
    // Run solver.
    switch (presolve_status) {
      case HighsPresolveStatus::NotPresolved: {
        hmos_[solved_hmo].lp_.lp_name_ = "Original LP";
        call_status = runLpSolver(hmos_[solved_hmo], "Not presolved: solving the LP");
	return_status = interpretCallStatus(call_status, return_status, "runLpSolver");
	if (return_status == HighsStatus::Error) return return_status;
        break;
      }
      case HighsPresolveStatus::NotReduced: {
        hmos_[solved_hmo].lp_.lp_name_ = "Unreduced LP";
        call_status = runLpSolver(hmos_[solved_hmo], "Problem not reduced by presolve: solving the LP");
	return_status = interpretCallStatus(call_status, return_status, "runLpSolver");
	if (return_status == HighsStatus::Error) return return_status;
        break;
      }
      case HighsPresolveStatus::Reduced: {
        HighsLp& reduced_lp = presolve_info.getReducedProblem();
        // Validate the reduced LP
        assert(assessLp(reduced_lp, options_) == HighsStatus::OK);
        // Add reduced lp object to vector of HighsModelObject,
        // so the last one in lp_ is the presolved one.

        hmos_.push_back(HighsModelObject(reduced_lp, options_, timer_));
        // Log the presolve reductions
        logPresolveReductions(hmos_[original_hmo].options_,
			      hmos_[original_hmo].lp_,
			      hmos_[presolve_hmo].lp_);
        // Record the HMO to be solved
        solved_hmo = presolve_hmo;
        hmos_[solved_hmo].lp_.lp_name_ = "Presolved LP";
	call_status = runLpSolver(hmos_[solved_hmo], "Solving the presolved LP");
	return_status = interpretCallStatus(call_status, return_status, "runLpSolver");
	if (return_status == HighsStatus::Error) return return_status;
        break;
      }
      case HighsPresolveStatus::ReducedToEmpty: {
        hmos_[0].unscaled_model_status_ = HighsModelStatus::OPTIMAL;
        // Proceed to postsolve.
        break;
      }
        //	printf("\nHighs::run() 3: presolve status = %d\n",
        //(int)presolve_status);fflush(stdout);
      case HighsPresolveStatus::Infeasible:
      case HighsPresolveStatus::Unbounded: {
        if (presolve_status == HighsPresolveStatus::Infeasible) {
          hmos_[original_hmo].unscaled_model_status_ =
              HighsModelStatus::PRIMAL_INFEASIBLE;
        } else {
          hmos_[original_hmo].unscaled_model_status_ =
              HighsModelStatus::PRIMAL_UNBOUNDED;
        }
        HighsLogMessage(options_.logfile, HighsMessageType::INFO,
			"Problem status detected on presolve: %s",
			highsModelStatusToString(hmos_[original_hmo].unscaled_model_status_).c_str());

        // Report this way for the moment. May modify after merge with
        // OSIinterface branch which has new way of setting up a
        // HighsModelObject and can support multiple calls to run(). Stop and
        // read the HiGHS clock, then work out time for this call
        if (!run_highs_clock_already_running) timer_.stopRunHighsClock();

	model_status_ = hmos_[original_hmo].unscaled_model_status_;
        return HighsStatus::OK;
      }
      default: {
        // case HighsPresolveStatus::Error
        HighsPrintMessage(options_.output, options_.message_level, ML_ALWAYS,
			  "Presolve failed.");
        if (!run_highs_clock_already_running) timer_.stopRunHighsClock();
	hmos_[original_hmo].unscaled_model_status_ = HighsModelStatus::PRESOLVE_ERROR;
	model_status_ = hmos_[original_hmo].unscaled_model_status_;
        return HighsStatus::Error;
      }
    }
    // Postsolve. Does nothing if there were no reductions during presolve.
    if (hmos_[solved_hmo].scaled_model_status_ == HighsModelStatus::OPTIMAL) {
      if (presolve_status == HighsPresolveStatus::Reduced ||
          presolve_status == HighsPresolveStatus::ReducedToEmpty) {
        // If presolve is nontrivial, extract the optimal solution
        // and basis for the presolved problem in order to generate
        // the solution and basis for postsolve to use to generate a
        // solution(?) and basis that is, hopefully, optimal. This is
        // confirmed or corrected by hot-starting the simplex solver
        presolve_info.reduced_solution_ = hmos_[solved_hmo].solution_;
        presolve_info.presolve_[0].setBasisInfo(
            hmos_[solved_hmo].basis_.col_status,
            hmos_[solved_hmo].basis_.row_status);
        // Run postsolve
        timer_.start(timer_.postsolve_clock);
        HighsPostsolveStatus postsolve_status = runPostsolve(presolve_info);
        timer_.stop(timer_.postsolve_clock);
        if (postsolve_status == HighsPostsolveStatus::SolutionRecovered) {
          HighsPrintMessage(options_.output, options_.message_level, ML_VERBOSE,
			    "Postsolve finished.");
	  //
          // Now hot-start the simplex solver for the original_hmo:
	  //
	  // The original model hasn't been solved, so set up its solution parameters
	  resetModelStatusAndSolutionParams(hmos_[original_hmo]);
	  // Set solution and its status
          hmos_[original_hmo].solution_ = presolve_info.recovered_solution_;
	  //
	  // Set basis and its status
          hmos_[original_hmo].basis_.col_status =
              presolve_info.presolve_[0].getColStatus();
          hmos_[original_hmo].basis_.row_status =
              presolve_info.presolve_[0].getRowStatus();
          hmos_[original_hmo].basis_.valid_ = true;
	  analyseHighsBasicSolution(options_.logfile,
				    hmos_[original_hmo],
				    "after returning from postsolve");
          // Now hot-start the simplex solver for the original_hmo
          solved_hmo = original_hmo;
          // Save the options to allow the best simplex strategy to
          // be used
          HighsOptions& options = hmos_[solved_hmo].options_;
          HighsOptions save_options = options;
	  const bool full_logging = false;
	  if (full_logging) options.message_level = ML_ALWAYS;
	  // Force the use of simplex to clean up if IPM has been used
	  // to solve the presolved problem
	  if (options.solver == ipm_string) options.solver = simplex_string;
          options.simplex_strategy = SIMPLEX_STRATEGY_CHOOSE;
	  // Ensure that the parallel solver isn't used
	  options.highs_min_threads = 1;
	  options.highs_max_threads = 1;
          hmos_[solved_hmo].lp_.lp_name_ = "Postsolve LP";
	  int iteration_count0 = hmos_[solved_hmo].unscaled_solution_params_.simplex_iteration_count;
	  call_status = runLpSolver(hmos_[solved_hmo], "Solving the original LP from the solution after postsolve");
	  return_status = interpretCallStatus(call_status, return_status, "runLpSolver");
          // Recover the options
          options = save_options;
	  if (return_status == HighsStatus::Error) return return_status;
	  int iteration_count1 = hmos_[solved_hmo].unscaled_solution_params_.simplex_iteration_count;
	  postsolve_iteration_count = iteration_count1 - iteration_count0;
        }
      }
    } else {
      // Optimal solution of presolved problem has not been found
      // The original model inherits the solved model's status
      hmos_[original_hmo].unscaled_model_status_ = hmos_[solved_hmo].unscaled_model_status_;
    }
  }
  else {
    // The problem has been solved before so we ignore presolve/postsolve/ipx.
    solved_hmo = original_hmo;
    // hmos_[solved_hmo].lp_.lp_name_ = "Aggregate LP";
    // hmos_[solved_hmo].basis_ = basis_;
    call_status = runLpSolver(hmos_[solved_hmo], "Unfolding the current aggregate");
    return_status = interpretCallStatus(call_status, return_status, "runLpSolver");
    if (return_status == HighsStatus::Error) return return_status;
  }
  // else if (reduced problem failed to solve) {
  //   todo: handle case when presolved problem failed to solve. Try to solve
  //   again with no presolve.
  // }

  //   assert(solved_hmo == original_hmo);
  // solved_hmo will be original_hmo unless the presolved LP is found to be infeasible or unbounded

  if (!getHighsModelStatusAndInfo(solved_hmo)) return HighsStatus::Error;

  // Copy HMO solution/basis to HiGHS solution/basis: this resizes solution_ and basis_
  // The HiGHS solution and basis have to come from the original_hmo
  // for them to have the right dimension.
  solution_ = hmos_[original_hmo].solution_;
  basis_ = hmos_[original_hmo].basis_;
  // tableau_ = hmos_[original_hmo].tableau_;
  // Report times
  if (hmos_[original_hmo].report_model_operations_clock) {
    std::vector<int> clockList{timer_.presolve_clock, timer_.solve_clock,
                               timer_.postsolve_clock};
    timer_.report("ModelOperations", clockList);
  }
  // Stop and read the HiGHS clock, then work out time for this call
  if (!run_highs_clock_already_running) timer_.stopRunHighsClock();
  // Record all revelavant times for OC and symmetry reducs
  sTimes = new struct solveTimeInfo;
  reducs = new struct symmetryReductionInfo;
  if (options_.aggregate == on_string){
    // Times
    sTimes->saucyTime = timer_.clock_time[timer_.saucy_clock];
    sTimes->foldTime = timer_.clock_time[timer_.fold_clock];
    sTimes->alpSolveTime = timer_.clock_time[timer_.elp_solve_clock];
    sTimes->liftTime = timer_.clock_time[timer_.lift_clock];
    sTimes->elpSolveTime = hmos_[solved_hmo].timer_.clock_time[hmos_[solved_hmo].timer_.orbital_crossover_clock];
    sTimes->solveTime = timer_.clock_time[timer_.solve_clock];
    sTimes->runTime = timer_.clock_time[timer_.run_highs_clock];
    // Reductions
    // reducs = (struct symmetryReductionInfo*)calloc((1), sizeof(struct symmetryReductionInfo));
    // Original numCol, numRow, numNnz
    int oNCol = originalLp.numCol_;
    int oNRow = originalLp.numRow_;
    int oNnz = originalLp.nnz_;
    // ALP numCol, numRwo, numNnz
    int alpNCol = elp_->numCol_;
    int alpNRow = elp_->numRow_;
    int alpNnz = elp_->nnz_;
    reducs->colReductions = (double)(oNCol - alpNCol)/oNCol * 100;
    reducs->rowReductions = (double)(oNRow - alpNRow)/oNRow * 100;
    reducs->nnzReductions = (double)(oNnz - alpNnz)/oNnz * 100;
  }
  else if (options_.solver == ipm_string){
    sTimes->solveTime = getCrossoverTime();
  }
  else{
    sTimes->solveTime = timer_.clock_time[timer_.solve_clock];
    sTimes->runTime = timer_.clock_time[timer_.run_highs_clock];
  }
  runTime = timer_.clock_time[timer_.run_highs_clock];

  // double lp_solve_final_time = timer_.readRunHighsClock();
  HighsPrintMessage(options_.output, options_.message_level, ML_MINIMAL,
		    "Postsolve            : %d\n", postsolve_iteration_count);
  HighsPrintMessage(options_.output, options_.message_level, ML_MINIMAL,
		    "Time                 : %0.3g\n", runTime);
  // HighsPrintMessage(options_.output, options_.message_level, ML_MINIMAL,
  //       "Possible R Piots : %0.3g\n", possibleRPivots);
  // HighsPrintMessage(options_.output, options_.message_level, ML_MINIMAL,
  //       "Actual R Piots   : %0.3g\n", actualRPivots);     
  std::cout << "Possible R Pivots    : " << possibleRPivots << std::endl;
  std::cout << "Actual R Pivots      : " << actualRPivots << std::endl;
  std::cout << "Presolve R Swaps     : " << rSwapPivots << std::endl;
  std::cout << "Total Simplex Pivots : " << totalPivots << std::endl;
  // totUnfoldTime_ += (lp_solve_final_time - initial_time);
  // totIter_ += hmos_[solved_hmo].unscaled_solution_params_.simplex_iteration_count;
  /* Determin whether to write time info to output file and reduction
  info to output file */
  if (options_.time_file != "") writeTimes(options_.time_file, solved_hmo);
  delete sTimes;
  if (options_.reduction_file != "") writeReductions(options_.reduction_file);
  delete reducs;
  // Assess success according to the scaled model status, unless
  // something worse has happened earlier
  call_status = highsStatusFromHighsModelStatus(scaled_model_status_);
  return_status = interpretCallStatus(call_status, return_status);
  return return_status;
}

// Time Orbital Crossover solution for highs crossover for clean up if needed.
void Highs::trimOCSolution(HighsModelObject& ocModel, HighsModelObject& model){
  // Original Model
  int numCol = model.lp_.numCol_;
  int numRow = model.lp_.numRow_;
  std::vector<double>& col_value = model.solution_.col_value;
  std::vector<double>& col_dual = model.solution_.col_dual;
  std::vector<double>& row_value = model.solution_.row_value;
  std::vector<double>& row_dual = model.solution_.row_dual;
  std::vector<HighsBasisStatus>& col_status = model.basis_.col_status;
  std::vector<HighsBasisStatus>& row_status = model.basis_.row_status;
  std::vector<double>& x = model.interior_.x;
  std::vector<double>& y = model.interior_.y;
  std::vector<double>& z = model.interior_.z;
  std::vector<double>& slack = model.interior_.slack;
  x.resize(numCol);
  y.resize(numRow);
  z.resize(numCol + numRow);
  slack.resize(numRow);
  col_value.resize(numCol);
  col_dual.resize(numCol);
  row_value.resize(numRow);
  row_dual.resize(numRow);
  col_status.resize(numCol);
  row_status.resize(numRow);
  int& basisNumCol = model.basis_.numCol_;
  int& basisNumRow = model.basis_.numRow_;
  basisNumCol = numCol;
  basisNumRow = numRow;
  // OC model after OC
  int numX = ocModel.lp_.numX_;
  int numS = ocModel.lp_.numS_;
  std::vector<int>& colrep = agg_.colrep;
  std::vector<int>& rowrep = agg_.rowrep;
  std::vector<double>& rowUpper = ocModel.lp_.rowUpper_;
  std::vector<double>& rowLower = ocModel.lp_.rowLower_;
  std::vector<double>& colCost = ocModel.lp_.colCost_;
  std::vector<double>& occol_value = ocModel.solution_.col_value;
  std::vector<double>& occol_dual = ocModel.solution_.col_dual;
  std::vector<double>& ocrow_value = ocModel.solution_.row_value;
  std::vector<double>& ocrow_dual = ocModel.solution_.row_dual;
  std::vector<HighsBasisStatus>& occol_status = ocModel.basis_.col_status;
  std::vector<HighsBasisStatus>& ocrow_status = ocModel.basis_.row_status;
  // Fill in original model basis and solution by trimming OC solution/basis
  for (int iCol = 0; iCol < numX; ++iCol){
    int rep = colrep[iCol];
    col_value[rep] = occol_value[iCol];
    col_dual[rep] = occol_dual[iCol];
    col_status[rep] = occol_status[iCol];
    x[rep] = occol_value[iCol];
    z[rep] = occol_dual[iCol]; 
  }
  for (int iRow; iRow < numS; ++iRow){
    int rep = rowrep[iRow] - numX;
    row_value[rep] = ocrow_value[iRow];
    row_dual[rep] = ocrow_dual[iRow];
    row_status[rep] = ocrow_status[iRow];
    y[rep] = - ocrow_dual[iRow];
    double low = std::fabs(rowLower[iRow]);
    double high = std::fabs(rowUpper[iRow]);
    double rhs = std::min(low, high);
    slack[rep] = rhs - ocrow_value[iRow];
    // z[rep + numX] = - ocrow_dual[iRow];
    z[rep + numX] = 0;
  }
}

// Compute equitable partition
void Highs::partitionLp(){
  // ep_.setUp(lp, options_.model_file.c_str());
  discrete = nep_.allocatePartition(&originalLp);
  part_ = nep_.getPartition();
}

void Highs::refinePartition(){
  discrete = nep_.isolate();
  part_ = nep_.getPartition();
}

void Highs::refinePartitionFinal(){
  nep_.runToDiscrete();
  part_ = nep_.getPartition();
}

// Fold Lp
void Highs::foldLp(){
  agg_.allocate(&originalLp, part_);
  alp_ = agg_.getLp();
}

// Lift ALP to LP
void Highs::liftLpFinal(){
  agg_.buildLp(part_, &alpBasis_, &alpSolution_, true, false);
  elp_ = agg_.getLp();
}

// Lift ALP to ELP
void Highs::liftLpExtended(){
  agg_.buildLp(part_, &alpBasis_, &alpSolution_, false, false);
  elp_ = agg_.getLp();
  slp_ = agg_.getSubLp();
  startBasis_ = agg_.getBasis();
  slpSymBasis_ = agg_.getSubBasis();
}

void Highs::liftLpExtendedFinal(){
  agg_.buildLp(part_, &alpBasis_, &alpSolution_, true, true);
  elp_ = agg_.getLp();
  slp_ = agg_.getSubLp();
  startBasis_ = agg_.getBasis();
  slpSymBasis_ = agg_.getSubBasis();
}

void Highs::liftBasis(){
  // agg_.buildBasis(false, false);
  startBasis_ = agg_.getBasis();
}

void Highs::liftBasisFinal(){
  agg_.buildBasis(true, false);
  startBasis_ = agg_.getBasis();
}

void Highs::liftBasisExtendedFinal(){
  agg_.buildBasis(true, true);
  startBasis_ = agg_.getBasis();
}

void Highs::liftSolutionFinal(){
  agg_.buildSolution(false, false);
  lpSymSolution_ = agg_.getSolution();
}

void Highs::liftSolutionExtended(){
  agg_.buildSolution(false, true);
  lpSymSolution_ = agg_.getSolution();
}

void Highs::liftSolutionExtendedFinal(){
  agg_.buildSolution(true, true);
  lpSymSolution_ = agg_.getSolution();
}

void Highs::ipxBasisToHighsBasis(){
  vbasis = getVbasis();
  cbasis = getCbasis();
  int numX_ = elp_->numX_;
}

int Highs::swapInRColumns(){
  std::vector<int>& unPerm = agg_.getUnPerm();
  int i, iRow, iCol, swaps = 0;
  for (i = 0; i < slp_->numRow_; ++i){
    if (alpBasis_.row_status[i] != HighsBasisStatus::BASIC){
      iRow = unPerm[i];
      startBasis_->row_status[iRow] = alpBasis_.row_status[i];
    }
  }
  for (i = 0; i < slp_->numCol_; ++i){
    if (alpBasis_.col_status[i] == HighsBasisStatus::BASIC){
      iCol = i + elp_->numX_;
      startBasis_->col_status[iCol] = HighsBasisStatus::BASIC;
      elp_->basicResiduals_[i] = 1;
      swaps++;
    }
  }
  return swaps;
}

const HighsLp& Highs::getLp() const { return lp_; }

const HighsSolution& Highs::getSolution() const { return solution_; }

const ICrashInfo& Highs::getICrashInfo() const { return icrash_info_; }

const HighsBasis& Highs::getBasis() const { return basis_; }

const HighsTableau& Highs::getTableau() const { return tableau_; }

const HighsModelStatus& Highs::getModelStatus(const bool scaled_model) const {
  if (scaled_model) {
    return scaled_model_status_;
  } else {
    return model_status_;
  }
}

HighsStatus Highs::getBasicVariables(int* basic_variables) {
  if (hmos_.size() == 0) return HighsStatus::Error;
  if (!hmos_[0].simplex_lp_status_.has_basis) {
    HighsLogMessage(options_.logfile, HighsMessageType::ERROR,
                    "No basis available in getBasicVariables");
    return HighsStatus::Error;
  }
  int numRow = hmos_[0].lp_.numRow_;
  int numCol = hmos_[0].lp_.numCol_;
  if (numRow != hmos_[0].simplex_lp_.numRow_) {
    HighsLogMessage(options_.logfile, HighsMessageType::ERROR,
        "Model LP and simplex LP row dimension difference (%d-%d=%d", numRow,
        hmos_[0].simplex_lp_.numRow_, numRow - hmos_[0].simplex_lp_.numRow_);
    return HighsStatus::Error;
  }
  for (int row = 0; row < numRow; row++) {
    int var = hmos_[0].simplex_basis_.basicIndex_[row];
    if (var < numCol) {
      basic_variables[row] = var;
    } else {
      basic_variables[row] = -(1 + var - numCol);
    }
  }
  return HighsStatus::OK;
}

HighsStatus Highs::getBasisInverseRow(const int row, double* row_vector,
                                      int* row_num_nz, int* row_indices) {
  if (hmos_.size() == 0) return HighsStatus::Error;
  int numRow = hmos_[0].lp_.numRow_;
  if (row < 0 || row >= numRow) {
    HighsLogMessage(options_.logfile, HighsMessageType::ERROR,
                    "Row index %d out of range [0, %d] in getBasisInverseRow",
                    row, numRow - 1);
    return HighsStatus::Error;
  }
  if (!hmos_[0].simplex_lp_status_.has_invert) {
    HighsLogMessage(options_.logfile, HighsMessageType::ERROR,
                    "No invertible representation for getBasisInverseRow");
    return HighsStatus::Error;
  }
  // Compute a row i of the inverse of the basis matrix by solving B^Tx=e_i
  vector<double> rhs;
  rhs.assign(numRow, 0);
  rhs[row] = 1;
  HighsSimplexInterface simplex_interface(hmos_[0]);
  simplex_interface.basisSolve(rhs, row_vector, row_num_nz, row_indices, true);
  return HighsStatus::OK;
}

HighsStatus Highs::getBasisInverseCol(const int col, double* col_vector,
                                      int* col_num_nz, int* col_indices) {
  if (hmos_.size() == 0) return HighsStatus::Error;
  int numRow = hmos_[0].lp_.numRow_;
  if (col < 0 || col >= numRow) {
    HighsLogMessage(options_.logfile, HighsMessageType::ERROR,
        "Column index %d out of range [0, %d] in getBasisInverseCol", col,
        numRow - 1);
    return HighsStatus::Error;
  }
  if (!hmos_[0].simplex_lp_status_.has_invert) {
    HighsLogMessage(options_.logfile, HighsMessageType::ERROR,
                    "No invertible representation for getBasisInverseCol");
    return HighsStatus::Error;
  }
  // Compute a col i of the inverse of the basis matrix by solving Bx=e_i
  vector<double> rhs;
  rhs.assign(numRow, 0);
  rhs[col] = 1;
  HighsSimplexInterface simplex_interface(hmos_[0]);
  simplex_interface.basisSolve(rhs, col_vector, col_num_nz, col_indices, false);
  return HighsStatus::OK;
}

HighsStatus Highs::getBasisSolve(const double* Xrhs, double* solution_vector,
                                 int* solution_num_nz, int* solution_indices) {
  if (hmos_.size() == 0) return HighsStatus::Error;
  if (!hmos_[0].simplex_lp_status_.has_invert) {
    HighsLogMessage(options_.logfile, HighsMessageType::ERROR,
                    "No invertible representation for getBasisSolve");
    return HighsStatus::Error;
  }
  int numRow = hmos_[0].lp_.numRow_;
  vector<double> rhs;
  rhs.assign(numRow, 0);
  for (int row = 0; row < numRow; row++) rhs[row] = Xrhs[row];
  HighsSimplexInterface simplex_interface(hmos_[0]);
  simplex_interface.basisSolve(rhs, solution_vector, solution_num_nz,
                               solution_indices, false);
  return HighsStatus::OK;
}

HighsStatus Highs::getBasisTransposeSolve(const double* Xrhs,
                                          double* solution_vector,
                                          int* solution_num_nz,
                                          int* solution_indices) {
  if (hmos_.size() == 0) return HighsStatus::Error;
  if (!hmos_[0].simplex_lp_status_.has_invert) {
    HighsLogMessage(options_.logfile, HighsMessageType::ERROR,
                    "No invertible representation for getBasisTransposeSolve");
    return HighsStatus::Error;
  }
  int numRow = hmos_[0].lp_.numRow_;
  vector<double> rhs;
  rhs.assign(numRow, 0);
  for (int row = 0; row < numRow; row++) rhs[row] = Xrhs[row];
  HighsSimplexInterface simplex_interface(hmos_[0]);
  simplex_interface.basisSolve(rhs, solution_vector, solution_num_nz,
                               solution_indices, true);
  return HighsStatus::OK;
}

HighsStatus Highs::getReducedRow(const int row, double* row_vector,
                                 int* row_num_nz, int* row_indices) {
  if (hmos_.size() == 0) return HighsStatus::Error;
  if (row < 0 || row >= hmos_[0].lp_.numRow_) {
    HighsLogMessage(options_.logfile, HighsMessageType::ERROR,
                    "Row index %d out of range [0, %d] in getReducedRow", row,
                    hmos_[0].lp_.numRow_ - 1);
    return HighsStatus::Error;
  }
  if (!hmos_[0].simplex_lp_status_.has_invert) {
    HighsLogMessage(options_.logfile, HighsMessageType::ERROR,
                    "No invertible representation for getReducedRow");
    return HighsStatus::Error;
  }
  HighsLp& lp = hmos_[0].lp_;
  int numRow = lp.numRow_;
  vector<double> rhs;
  vector<double> col_vector;
  vector<int> col_indices;
  int col_num_nz;
  rhs.assign(numRow, 0);
  rhs[row] = 1;
  col_vector.resize(numRow, 0);
  col_indices.resize(numRow, 0);
  HighsSimplexInterface simplex_interface(hmos_[0]);
  // Form B^{-T}e_{row}
  simplex_interface.basisSolve(rhs, &col_vector[0], &col_num_nz,
                               &col_indices[0], true);
  bool return_indices = row_num_nz != NULL;
  if (return_indices) *row_num_nz = 0;
  for (int col = 0; col < lp.numCol_; col++) {
    double value = 0;
    for (int el = lp.Astart_[col]; el < lp.Astart_[col + 1]; el++) {
      int row = lp.Aindex_[el];
      value += lp.Avalue_[el] * col_vector[row];
    }
    row_vector[col] = 0;
    if (fabs(value) > HIGHS_CONST_TINY) {
      if (return_indices) row_indices[(*row_num_nz)++] = col;
      row_vector[col] = value;
    }
  }
  return HighsStatus::OK;
}

HighsStatus Highs::getReducedColumn(const int col, double* col_vector,
                                    int* col_num_nz, int* col_indices) {
  if (hmos_.size() == 0) return HighsStatus::Error;
  if (col < 0 || col >= hmos_[0].lp_.numCol_) {
    HighsLogMessage(options_.logfile, HighsMessageType::ERROR,
                    "Column index %d out of range [0, %d] in getReducedColumn",
                    col, hmos_[0].lp_.numCol_ - 1);
    return HighsStatus::Error;
  }
  if (!hmos_[0].simplex_lp_status_.has_invert) {
    HighsLogMessage(options_.logfile, HighsMessageType::ERROR,
                    "No invertible representation for getReducedColumn");
    return HighsStatus::Error;
  }
  HighsLp& lp = hmos_[0].lp_;
  int numRow = lp.numRow_;
  vector<double> rhs;
  rhs.assign(numRow, 0);
  for (int el = lp.Astart_[col]; el < lp.Astart_[col + 1]; el++)
    rhs[lp.Aindex_[el]] = lp.Avalue_[el];
  HighsSimplexInterface simplex_interface(hmos_[0]);
  simplex_interface.basisSolve(rhs, col_vector, col_num_nz, col_indices, false);
  return HighsStatus::OK;
}

HighsStatus Highs::setSolution(const HighsSolution& solution) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  underDevelopmentLogMessage("setSolution");
  // Check if solution is valid.
  assert((int)solution_.col_value.size() != 0 ||
         (int)solution_.col_value.size() != lp_.numCol_);
  assert((int)solution.col_dual.size() == 0 ||
         (int)solution.col_dual.size() == lp_.numCol_);
  assert((int)solution.row_dual.size() == 0 ||
         (int)solution.row_dual.size() == lp_.numRow_);

  if (solution.col_value.size()) solution_.col_value = solution.col_value;
  if (solution.col_dual.size()) solution_.col_dual = solution.col_dual;
  if (solution.row_dual.size()) solution_.row_dual = solution.row_dual;

  if (solution.col_value.size() > 0) {
    call_status = calculateRowValues(lp_, solution_);
    return_status = interpretCallStatus(call_status, return_status, "calculateRowValues");
    if (return_status == HighsStatus::Error) return return_status;
  }
  if (solution.row_dual.size() > 0) {
    call_status = calculateColDuals(lp_, solution_);
    return_status = interpretCallStatus(call_status, return_status, "calculateColDuals");
    if (return_status == HighsStatus::Error) return return_status;
  }
  return return_status;
}

HighsStatus Highs::setBasis(const HighsBasis& basis) {
  underDevelopmentLogMessage("setBasis");
  if (!basisOk(options_.logfile, lp_, basis)) {
    HighsLogMessage(options_.logfile, HighsMessageType::ERROR,
		    "setBasis: invalid basis");
    return HighsStatus::Error;
  }
  basis_ = basis;
  basis_.valid_ = true;
  return HighsStatus::OK;
}

bool Highs::addRow(const double lower_bound, const double upper_bound,
                   const int num_new_nz, const int* indices,
                   const double* values) {
  int starts = 0;
  return addRows(1, &lower_bound, &upper_bound, num_new_nz, &starts, indices,
                 values);
}

bool Highs::addRows(const int num_new_row, const double* lower_bounds,
                    const double* upper_bounds, const int num_new_nz,
                    const int* starts, const int* indices,
                    const double* values) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  underDevelopmentLogMessage("addRows");
  // Check that there is a HighsModelObject
  if (!haveHmo("addRows")) return false;
  HighsSimplexInterface interface(hmos_[0]);
  call_status = interface.addRows(num_new_row, lower_bounds, upper_bounds,
				  num_new_nz, starts, indices, values);
  return_status = interpretCallStatus(call_status, return_status, "addRows");
  if (return_status == HighsStatus::Error) return false;
  if (!updateHighsSolutionBasis()) return false;
  return return_status != HighsStatus::Error;
}

bool Highs::addCol(const double cost, const double lower_bound,
                   const double upper_bound, const int num_new_nz,
                   const int* indices, const double* values) {
  int starts = 0;
  return addCols(1, &cost, &lower_bound, &upper_bound, num_new_nz, &starts,
                 indices, values);
}

bool Highs::addCols(const int num_new_col, const double* costs,
                    const double* lower_bounds, const double* upper_bounds,
                    const int num_new_nz, const int* starts, const int* indices,
                    const double* values) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  underDevelopmentLogMessage("addCols");
  if (!haveHmo("addCols")) return false;
  HighsSimplexInterface interface(hmos_[0]);
  call_status = interface.addCols(num_new_col, costs, lower_bounds, upper_bounds,
				  num_new_nz, starts, indices, values);
  return_status = interpretCallStatus(call_status, return_status, "addCols");
  if (return_status == HighsStatus::Error) return false;
  if (!updateHighsSolutionBasis()) return false;
  return return_status != HighsStatus::Error;
}

bool Highs::changeObjectiveSense(const int sense) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  underDevelopmentLogMessage("changeObjectiveSense");
  if (!haveHmo("changeObjectiveSense")) return false;
  HighsSimplexInterface interface(hmos_[0]);
  call_status = interface.changeObjectiveSense(sense);
  return_status = interpretCallStatus(call_status, return_status, "changeObjectiveSense");
  if (return_status == HighsStatus::Error) return false;
  return return_status != HighsStatus::Error;
}

bool Highs::changeColCost(const int col, const double cost) {
  return changeColsCost(1, &col, &cost);
}

bool Highs::changeColsCost(const int num_set_entries, const int* set,
                           const double* cost) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  underDevelopmentLogMessage("changeColsCost");
  if (!haveHmo("changeColsCost")) return false;
  HighsSimplexInterface interface(hmos_[0]);
  call_status = interface.changeCosts(num_set_entries, set, cost);
  return_status = interpretCallStatus(call_status, return_status, "changeCosts");
  if (return_status == HighsStatus::Error) return false;
  return return_status != HighsStatus::Error;
}

bool Highs::changeColsCost(const int* mask, const double* cost) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  underDevelopmentLogMessage("changeColsCost");
  if (!haveHmo("changeColsCost")) return false;
  HighsSimplexInterface interface(hmos_[0]);
  call_status = interface.changeCosts(mask, cost);
  return_status = interpretCallStatus(call_status, return_status, "changeCosts");
  if (return_status == HighsStatus::Error) return false;
  return return_status != HighsStatus::Error;
}

bool Highs::changeColBounds(const int col, const double lower,
                            const double upper) {
  return changeColsBounds(1, &col, &lower, &upper);
}

bool Highs::changeColsBounds(const int num_set_entries, const int* set,
                             const double* lower, const double* upper) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  underDevelopmentLogMessage("changeColsBounds");
  if (!haveHmo("changeColsBounds")) return false;
  HighsSimplexInterface interface(hmos_[0]);
  call_status = interface.changeColBounds(num_set_entries, set, lower, upper);
  return_status = interpretCallStatus(call_status, return_status, "changeColBounds");
  if (return_status == HighsStatus::Error) return false;
  return return_status != HighsStatus::Error;
}

bool Highs::changeColsBounds(const int from_col, const int to_col,
                             const double* lower, const double* upper) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  underDevelopmentLogMessage("changeColsBounds");
  if (!haveHmo("changeColsBounds")) return false;
  HighsSimplexInterface interface(hmos_[0]);
  call_status = interface.changeColBounds(from_col, to_col, lower, upper);
  return_status = interpretCallStatus(call_status, return_status, "changeColBounds");
  if (return_status == HighsStatus::Error) return false;
  return return_status != HighsStatus::Error;
}

bool Highs::changeColsBounds(const int* mask, const double* lower,
                             const double* upper) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  underDevelopmentLogMessage("changeColsBounds");
  if (!haveHmo("changeColsBounds")) return false;
  HighsSimplexInterface interface(hmos_[0]);
  call_status = interface.changeColBounds(mask, lower, upper);
  return_status = interpretCallStatus(call_status, return_status, "changeColBounds");
  if (return_status == HighsStatus::Error) return false;
  return return_status != HighsStatus::Error;
}

bool Highs::changeRowBounds(const int row, const double lower,
                            const double upper) {
  return changeRowsBounds(1, &row, &lower, &upper);
}

bool Highs::changeRowsBounds(const int num_set_entries, const int* set,
                             const double* lower, const double* upper) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  underDevelopmentLogMessage("changeRowsBounds");
  if (!haveHmo("changeRowsBounds")) return false;
  HighsSimplexInterface interface(hmos_[0]);
  call_status = interface.changeRowBounds(num_set_entries, set, lower, upper);
  return_status = interpretCallStatus(call_status, return_status, "changeRowBounds");
  if (return_status == HighsStatus::Error) return false;
  return return_status != HighsStatus::Error;
}

bool Highs::changeRowsBounds(const int* mask, const double* lower,
                             const double* upper) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  underDevelopmentLogMessage("changeRowsBounds");
  if (!haveHmo("changeRowsBounds")) return false;
  HighsSimplexInterface interface(hmos_[0]);
  call_status = interface.changeRowBounds(mask, lower, upper);
  return_status = interpretCallStatus(call_status, return_status, "changeRowBounds");
  if (return_status == HighsStatus::Error) return false;
  return return_status != HighsStatus::Error;
}

bool Highs::changeCoeff(const int row, const int col, const double value) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  underDevelopmentLogMessage("changeCoeff");
  if (!haveHmo("changeCoeff")) return false;
  HighsSimplexInterface interface(hmos_[0]);
  call_status = interface.changeCoefficient(row, col, value);
  return_status = interpretCallStatus(call_status, return_status, "changeCoefficient");
  if (return_status == HighsStatus::Error) return false;
  return return_status != HighsStatus::Error;
}

bool Highs::getCols(const int from_col, const int to_col, int& num_col,
                    double* costs, double* lower, double* upper, int& num_nz,
                    int* start, int* index, double* value) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  underDevelopmentLogMessage("getCols");
  if (!haveHmo("getCols")) return false;
  HighsSimplexInterface interface(hmos_[0]);
  call_status = interface.getCols(from_col, to_col, num_col, costs, lower,
				  upper, num_nz, start, index, value);
  return_status = interpretCallStatus(call_status, return_status, "getCols");
  if (return_status == HighsStatus::Error) return false;
  return return_status != HighsStatus::Error;
}

bool Highs::getCols(const int n, const int* set, int& num_col, double* costs,
                    double* lower, double* upper, int& num_nz, int* start,
                    int* index, double* value) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  underDevelopmentLogMessage("getCols");
  if (!haveHmo("getCols")) return false;
  HighsSimplexInterface interface(hmos_[0]);
  call_status = interface.getCols(n, set, num_col, costs, lower, upper,
				  num_nz, start, index, value);
  return_status = interpretCallStatus(call_status, return_status, "getCols");
  if (return_status == HighsStatus::Error) return false;
  return return_status != HighsStatus::Error;
}

bool Highs::getCols(const int* col_mask, int& num_col, double* costs,
                    double* lower, double* upper, int& num_nz, int* start,
                    int* index, double* value) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  underDevelopmentLogMessage("getCols");
  if (!haveHmo("getCols")) return false;
  HighsSimplexInterface interface(hmos_[0]);
  call_status = interface.getCols(col_mask, num_col, costs, lower, upper,
				  num_nz, start, index, value);
  return_status = interpretCallStatus(call_status, return_status, "getCols");
  if (return_status == HighsStatus::Error) return false;
  return return_status != HighsStatus::Error;
}

bool Highs::getRows(const int from_row, const int to_row, int& num_row,
                    double* lower, double* upper, int& num_nz, int* start,
                    int* index, double* value) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  underDevelopmentLogMessage("getRows");
  if (!haveHmo("getRows")) return false;
  HighsSimplexInterface interface(hmos_[0]);
  call_status = interface.getRows(from_row, to_row, num_row, lower, upper,
				  num_nz, start, index, value);
  return_status = interpretCallStatus(call_status, return_status, "getRows");
  if (return_status == HighsStatus::Error) return false;
  return return_status != HighsStatus::Error;
}

bool Highs::getRows(const int num_set_entries, const int* set, int& num_row,
                    double* lower, double* upper, int& num_nz, int* start,
                    int* index, double* value) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  underDevelopmentLogMessage("getRows");
  if (!haveHmo("getRows")) return false;
  HighsSimplexInterface interface(hmos_[0]);
  call_status = interface.getRows(num_set_entries, set, num_row, lower, upper,
				  num_nz, start, index, value);
  return_status = interpretCallStatus(call_status, return_status, "getRows");
  if (return_status == HighsStatus::Error) return false;
  return return_status != HighsStatus::Error;
}

bool Highs::getRows(const int* mask, int& num_row, double* lower, double* upper,
                    int& num_nz, int* start, int* index, double* value) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  underDevelopmentLogMessage("getRows");
  if (!haveHmo("getRows")) return false;
  HighsSimplexInterface interface(hmos_[0]);
  call_status = interface.getRows(mask, num_row, lower, upper, num_nz, start,
				  index, value);
  return_status = interpretCallStatus(call_status, return_status, "getRows");
  if (return_status == HighsStatus::Error) return false;
  return return_status != HighsStatus::Error;
}

bool Highs::getCoeff(const int row, const int col, double& value) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  underDevelopmentLogMessage("getCoeff");
  if (!haveHmo("getCoeff")) return false;
  HighsSimplexInterface interface(hmos_[0]);

  call_status = interface.getCoefficient(row, col, value);
  return_status = interpretCallStatus(call_status, return_status, "getCoefficient");
  if (return_status == HighsStatus::Error) return false;
  return return_status != HighsStatus::Error;
}

bool Highs::deleteCols(const int from_col, const int to_col) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  underDevelopmentLogMessage("deleteCols");
  if (!haveHmo("deleteCols")) return false;
  HighsSimplexInterface interface(hmos_[0]);
  call_status = interface.deleteCols(from_col, to_col);
  return_status = interpretCallStatus(call_status, return_status, "deleteCols");
  if (return_status == HighsStatus::Error) return false;
  if (!updateHighsSolutionBasis()) return false;
  return return_status != HighsStatus::Error;
}

bool Highs::deleteCols(const int num_set_entries, const int* set) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  underDevelopmentLogMessage("deleteCols");
  if (!haveHmo("deleteCols")) return false;
  HighsSimplexInterface interface(hmos_[0]);
  call_status = interface.deleteCols(num_set_entries, set);
  return_status = interpretCallStatus(call_status, return_status, "deleteCols");
  if (return_status == HighsStatus::Error) return false;
  if (!updateHighsSolutionBasis()) return false;
  return return_status != HighsStatus::Error;
}

bool Highs::deleteCols(int* mask) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  underDevelopmentLogMessage("deleteCols");
  if (!haveHmo("deleteCols")) return false;
  HighsSimplexInterface interface(hmos_[0]);
  call_status = interface.deleteCols(mask);
  return_status = interpretCallStatus(call_status, return_status, "deleteCols");
  if (return_status == HighsStatus::Error) return false;
  if (!updateHighsSolutionBasis()) return false;
  return return_status != HighsStatus::Error;
}

bool Highs::deleteRows(const int from_row, const int to_row) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  underDevelopmentLogMessage("deleteRows");
  if (!haveHmo("deleteRows")) return false;
  HighsSimplexInterface interface(hmos_[0]);
  call_status = interface.deleteRows(from_row, to_row);
  return_status = interpretCallStatus(call_status, return_status, "deleteRows");
  if (return_status == HighsStatus::Error) return false;
  if (!updateHighsSolutionBasis()) return false;
  return return_status != HighsStatus::Error;
}

bool Highs::deleteRows(const int num_set_entries, const int* set) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  underDevelopmentLogMessage("deleteRows");
  if (!haveHmo("deleteRows")) return false;
  HighsSimplexInterface interface(hmos_[0]);
  call_status = interface.deleteRows(num_set_entries, set);
  return_status = interpretCallStatus(call_status, return_status, "deleteRows");
  if (return_status == HighsStatus::Error) return false;
  if (!updateHighsSolutionBasis()) return false;
  return return_status != HighsStatus::Error;
}

bool Highs::deleteRows(int* mask) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  underDevelopmentLogMessage("deleteRows");
  if (!haveHmo("deleteRows")) return false;
  HighsSimplexInterface interface(hmos_[0]);
  call_status = interface.deleteRows(mask);
  return_status = interpretCallStatus(call_status, return_status, "deleteRows");
  if (return_status == HighsStatus::Error) return false;
  if (!updateHighsSolutionBasis()) return false;
  return return_status != HighsStatus::Error;
}

HighsStatus Highs::clearSolver() {
  underDevelopmentLogMessage("clearSolver");
  basis_.valid_ = false;
  return HighsStatus::OK;
}

#ifdef HiGHSDEV
void Highs::reportModelStatusSolutionBasis(const std::string message,
					   const HighsModelStatus model_status,
					   const HighsLp &lp,
					   const HighsSolution &solution,
					   const HighsBasis &basis) {
  printf("\n%s\nModelStatus = %s; LP(%d, %d); solution (%d, %d; %d, %d); basis %d (%d, %d)\n\n",
	 message.c_str(), utilHighsModelStatusToString(model_status).c_str(), lp.numCol_, lp.numRow_,
	 (int)solution.col_value.size(), (int)solution.row_value.size(), (int)solution.col_dual.size(), (int)solution.row_dual.size(),
	 basis.valid_, (int)basis.col_status.size(), (int)basis.row_status.size());
}
#endif

std::string Highs::highsModelStatusToString(const HighsModelStatus model_status) const {
  return utilHighsModelStatusToString(model_status);
}

std::string Highs::highsPrimalDualStatusToString(const int primal_dual_status) {
  return utilPrimalDualStatusToString(primal_dual_status);
}

// Private methods
HighsPresolveStatus Highs::runPresolve(PresolveInfo& info) {
  if (options_.presolve == off_string) return HighsPresolveStatus::NotPresolved;

  if (info.lp_ == nullptr) return HighsPresolveStatus::NullError;

  if (info.presolve_.size() == 0) return HighsPresolveStatus::NotReduced;

  info.presolve_[0].load(*(info.lp_));

  // Initialize a new presolve class instance for the LP given in presolve info
  return info.presolve_[0].presolve();
}

HighsPostsolveStatus Highs::runPostsolve(PresolveInfo& info) {
  if (info.presolve_.size() != 0) {
    bool solution_ok =
        isSolutionConsistent(info.getReducedProblem(), info.reduced_solution_);
    if (!solution_ok)
      return HighsPostsolveStatus::ReducedSolutionDimenionsError;

    // todo: error handling + see todo in run()
    info.presolve_[0].postsolve(info.reduced_solution_,
                                info.recovered_solution_);

    return HighsPostsolveStatus::SolutionRecovered;
  } else {
    return HighsPostsolveStatus::NoPostsolve;
  }
}

// The method below runs simplex or ipx solver on the lp.
HighsStatus Highs::runLpSolver(HighsModelObject& model, const string message) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  // Reset unscaled and scaled model status and solution params - except for iteration counts
  resetModelStatusAndSolutionParams(model);
  HighsLogMessage(options_.logfile, HighsMessageType::INFO,
		  message.c_str());
#ifdef HIGHSDEV
  // Shouldn't have to check validity of the LP since this is done when it is
  // loaded or modified
  //  bool normalise = true;
  call_status = assessLp(lp_, options_);
  assert(call_status == HighsStatus::OK);
  return_status = interpretCallStatus(call_status, return_status, "assessLp");
  if (return_status == HighsStatus::Error) return return_status;
#endif
  if (!model.lp_.numRow_) {
    // Unconstrained LP so solve directly
    call_status = solveUnconstrainedLp(model);
    return_status = interpretCallStatus(call_status, return_status, "solveUnconstrainedLp");
    if (return_status == HighsStatus::Error) return return_status;
  } 
  // else if (options_.aggregate == on_string && 
  //            options_.solver == ipm_sym_string){
  //   // Use Simplex
  //   call_status = solveLpSimplex(model);
  //   return_status = interpretCallStatus(call_status, return_status, "solveLpSimplex");
  //   if (return_status == HighsStatus::Error) return return_status;

  //   if (!isSolutionConsistent(model.lp_, model.solution_)) {
  //     HighsLogMessage(options_.logfile, HighsMessageType::ERROR,
	// 	      "Inconsistent solution returned from solver");
  //     return HighsStatus::Error;
  //   }
  //   // Record basis and solution
  //   alpSolution_ = hmos_[0].solution_;
  //   alpBasis_ = hmos_[0].basis_;
  //   liftLpExtendedFinal();
  //   liftBasisExtendedFinal();
  //   liftSolutionExtendedFinal();
  //   // liftLpFinal();
  //   // liftBasisFinal();
  //   // liftSolutionFinal();
  //   HighsPrintMessage(options_.output, options_.message_level, ML_ALWAYS,
	// 	      "Starting Crossover From Symmetric Solution...\n");
  //   passModel(*elp_);
  //   call_status = CrossoverFromSymmetricSolution(model.lp_, options_,
	// 		     *startBasis_, *lpSymSolution_,
	// 		     model.unscaled_model_status_,
	// 		     model.unscaled_solution_params_);
  // } 
  // else if (options_.aggregate == on_string &&
  //            options_.solver == ipm_string){
  //   // Use Simplex
  //   call_status = solveLpSimplex(model);
  //   return_status = interpretCallStatus(call_status, return_status, "solveLpSimplex");
  //   if (return_status == HighsStatus::Error) return return_status;

  //   if (!isSolutionConsistent(model.lp_, model.solution_)) {
  //     HighsLogMessage(options_.logfile, HighsMessageType::ERROR,
	// 	      "Inconsistent solution returned from solver");
  //     return HighsStatus::Error;
  //   }
  //   // Record basis and solution
  //   alpSolution_ = hmos_[0].solution_;
  //   alpBasis_ = hmos_[0].basis_;
  //   liftLpFinal();
  //   // liftBasis();
  //   liftSolutionFinal();
  //   HighsPrintMessage(options_.output, options_.message_level, ML_ALWAYS,
	// 	      "Starting Crossover From Symmetric Solution...\n");
  //   passModel(*elp_);
  //   call_status = CrossoverFromSymmetricSolution(model.lp_, options_,
	// 		     *startBasis_, *lpSymSolution_,
	// 		     model.unscaled_model_status_,
	// 		     model.unscaled_solution_params_);
  // } 
  else if (options_.solver == ipm_string) {
    // Use IPM
    if (!model.readyForHighs){
#ifdef IPX_ON
    HighsPrintMessage(options_.output, options_.message_level, ML_ALWAYS,
		      "Starting IPX...\n");
    call_status = solveLpIpx(model.lp_, options_,
			     model.basis_, model.solution_,
			     model.unscaled_model_status_,
			     model.unscaled_solution_params_);
    return_status = interpretCallStatus(call_status, return_status, "solveLpIpx");
    if (return_status == HighsStatus::Error) return return_status;
    // Set the scaled model status and solution params for completeness
    model.scaled_model_status_ = model.unscaled_model_status_;
    model.scaled_solution_params_ = model.unscaled_solution_params_;
#else
    HighsLogMessage(options_.logfile, HighsMessageType::ERROR,
		    "Model cannot be solved with IPM");
    return HighsStatus::Error;
#endif
    }
    else{
      #ifdef IPX_ON
    HighsPrintMessage(options_.output, options_.message_level, ML_ALWAYS,
		      "Starting IPX...\n");
    call_status = BuildOCFinalBasis(model.lp_, options_,
			     model.basis_, model.solution_,
			     model.unscaled_model_status_,
			     model.unscaled_solution_params_,
           model.interior_);
    return_status = interpretCallStatus(call_status, return_status, "Clean up OC Interior Basis");
    if (return_status == HighsStatus::Error) return return_status;
    // Set the scaled model status and solution params for completeness
    model.scaled_model_status_ = model.unscaled_model_status_;
    model.scaled_solution_params_ = model.unscaled_solution_params_;
#else
    HighsLogMessage(options_.logfile, HighsMessageType::ERROR,
		    "Model cannot be solved with IPM");
    return HighsStatus::Error;
#endif
    }
  } else {
    // Use Simplex
    call_status = solveLpSimplex(model);
    return_status = interpretCallStatus(call_status, return_status, "solveLpSimplex");
    if (return_status == HighsStatus::Error) return return_status;

    if (!isSolutionConsistent(model.lp_, model.solution_)) {
      HighsLogMessage(options_.logfile, HighsMessageType::ERROR,
		      "Inconsistent solution returned from solver");
      return HighsStatus::Error;
    }
  }
  call_status = analyseHighsBasicSolution(options_.logfile,
					  model.lp_, model.basis_, model.solution_,
					  model.unscaled_model_status_,
					  model.unscaled_solution_params_,
					  message);
  return_status = interpretCallStatus(call_status, return_status, "analyseHighsBasicSolution");
  return return_status;
}

HighsStatus Highs::writeSolution(const std::string filename, const bool pretty) const {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  HighsLp lp = this->lp_;
  HighsBasis basis = this->basis_;
  HighsSolution solution = this->solution_;
  FILE* file;
  bool html;
  call_status = openWriteFile(filename, "writeSolution", file, html);
  return_status = interpretCallStatus(call_status, return_status, "openWriteFile");
  if (return_status == HighsStatus::Error) return return_status;

  writeSolutionToFile(file, lp, basis, solution, pretty);
  return HighsStatus::OK;
}

HighsStatus Highs::writeTimes(const std::string filename, int solvedHmo){
  if (!options_.aggregate.compare(on_string)) writeTimesToFile(filename, sTimes, on_string,
                                                        options_.model_file, hmos_[solvedHmo].simplex_info_.primal_objective_value,
                                                        hmos_[solvedHmo].simplex_info_.dual_objective_value);
  else if (!options_.solver.compare(ipm_string)) writeTimesToFile(filename, sTimes, ipm_string,
                                                        options_.model_file, hmos_[solvedHmo].simplex_info_.primal_objective_value,
                                                        hmos_[solvedHmo].simplex_info_.dual_objective_value);
  else writeTimesToFile(filename, sTimes, off_string,
                        options_.model_file, hmos_[solvedHmo].simplex_info_.primal_objective_value,
                        hmos_[solvedHmo].simplex_info_.dual_objective_value);
  return HighsStatus::OK;
}

HighsStatus Highs::writeReductions(const std::string filename){
  if (!options_.aggregate.compare(on_string)) writeReductionsToFile(filename, reducs, options_.model_file);
  else HighsLogMessage(options_.logfile, HighsMessageType::INFO,
                        "Cannot write reductions stats without using Orbital Crossover");
  return HighsStatus::OK;
}

bool Highs::updateHighsSolutionBasis() {
  if (!haveHmo("updateHighsSolutionBasis")) return false;
  solution_.col_value.resize(lp_.numCol_);
  solution_.row_value.resize(lp_.numRow_);
  solution_.col_dual.resize(lp_.numCol_);
  solution_.row_dual.resize(lp_.numRow_);
  hmos_[0].solution_.col_value.resize(lp_.numCol_);
  hmos_[0].solution_.row_value.resize(lp_.numRow_);
  hmos_[0].solution_.col_dual.resize(lp_.numCol_);
  hmos_[0].solution_.row_dual.resize(lp_.numRow_);

  if (hmos_[0].basis_.valid_) {
    basis_ = hmos_[0].basis_;
  } else {
    basis_.valid_ = false;
    basis_.col_status.resize(lp_.numCol_);
    basis_.row_status.resize(lp_.numRow_);
  }
  return true;
}

bool Highs::getHighsModelStatusAndInfo(const int solved_hmo) {
  if (!haveHmo("getHighsModelStatusAndInfo")) return false;

  model_status_ = hmos_[solved_hmo].unscaled_model_status_;
  scaled_model_status_ = hmos_[solved_hmo].scaled_model_status_;

  HighsSolutionParams& solution_params = hmos_[solved_hmo].unscaled_solution_params_;

  // Get the total simplex IPM and crossover iteration counts over all HMO
  info_.simplex_iteration_count = 0;
  info_.ipm_iteration_count = 0;
  info_.crossover_iteration_count = 0;
  int hmos_size = hmos_.size();
  for (int k = 0; k < hmos_size; k++) {
    info_.simplex_iteration_count += hmos_[k].unscaled_solution_params_.simplex_iteration_count;
    info_.ipm_iteration_count += hmos_[k].unscaled_solution_params_.ipm_iteration_count;
    info_.crossover_iteration_count += hmos_[k].unscaled_solution_params_.crossover_iteration_count;
  }
  info_.primal_status = solution_params.primal_status;
  info_.dual_status = solution_params.dual_status;
  info_.objective_function_value = solution_params.objective_function_value;
  info_.num_primal_infeasibilities = solution_params.num_primal_infeasibilities;
  info_.max_primal_infeasibility = solution_params.max_primal_infeasibility;
  info_.sum_primal_infeasibilities = solution_params.sum_primal_infeasibilities;
  info_.num_dual_infeasibilities = solution_params.num_dual_infeasibilities;
  info_.max_dual_infeasibility = solution_params.max_dual_infeasibility;
  info_.sum_dual_infeasibilities = solution_params.sum_dual_infeasibilities;
  return true;
}

HighsStatus Highs::
openWriteFile(const string filename, const string method_name, FILE*& file, bool& html) const {
  html = false;
  if (filename == "") {
    // Empty file name: use stdout
    file = stdout;
  } else {
    file = fopen(filename.c_str(), "w");
    if (file == 0) {
      HighsLogMessage(options_.logfile, HighsMessageType::ERROR,
		      "Cannot open writeable file \"%s\" in %s",
		      filename.c_str(), method_name.c_str());
      return HighsStatus::Error;
    }
    const char* dot = strrchr(filename.c_str(), '.');
    if (dot && dot != filename) html = strcmp(dot + 1, "html") == 0;
  }
  return HighsStatus::OK;
}

bool Highs::haveHmo(const string method_name) {
  bool have_hmo = hmos_.size() > 0;
  assert(have_hmo);
#ifdef HiGHSDEV
  if (!have_hmo)
    HighsLogMessage(options_.logfile, HighsMessageType::ERROR,
		    "Method %s called without any HighsModelObject",
		    method_name.c_str());
#endif
  return have_hmo;
}

void Highs::underDevelopmentLogMessage(const string method_name) {
  HighsLogMessage(options_.logfile, HighsMessageType::WARNING,
      "Method %s is still under development and behaviour may be unpredictable",
      method_name.c_str());
}
