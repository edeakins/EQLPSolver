
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2019 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file ../app/RunHighs.cpp
 * @brief HiGHS main
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#include "Highs.h"
#include "HighsIO.h"
#include "HighsMipSolver.h"
#include "HighsOptions.h"
#include "HighsRuntimeOptions.h"
#include "HighsTimer.h"
#include "LoadProblem.h"
#include "Aggregate.h"

void printHighsVersionCopyright(FILE* output, const int message_level,
                                const char* message = nullptr);
void reportLpStatsOrError(FILE* output, int message_level,
                          const HighsStatus read_status, const HighsLp& lp);
void reportSolvedLpStats(FILE* output, int message_level,
                         const HighsStatus run_status, const Highs& highs);
HighsStatus callLpSolver(const HighsOptions& options, HighsLp& lp,
                         FILE* output, int message_level, bool run_quiet);
HighsStatus callMipSolver(const HighsOptions& options, const HighsLp& lp,
                          FILE* output, int message_level, bool run_quiet);

int main(int argc, char** argv) {
  // std::ofstream resultsFile("results.txt", std::ios_base::app);
  // std::clock_t start;
  // start = std::clock();
  // double duration;
  printHighsVersionCopyright(stdout, ML_ALWAYS);

  // Load user options.
  HighsOptions options;
  bool options_ok = loadOptions(argc, argv, options);
  if (!options_ok) return 0;

  // Set message level.
  FILE* output;
  int message_level;
  output = options.output;
  message_level = options.message_level;

  bool run_quiet = false;  // true;//
  if (run_quiet) {
    HighsPrintMessage(output, message_level, ML_ALWAYS,
                      "In main: running highs.run() quietly\n");
  }

  output = options.output;
  message_level = options.message_level;

  // Load problem.
  HighsLp lp;
  HighsStatus read_status = loadLpFromFile(options, lp);
  reportLpStatsOrError(output, message_level, read_status, lp);
  if (read_status == HighsStatus::Error) return (int)HighsStatus::Error;

  // Run LP or MIP solver.
  HighsStatus run_status = HighsStatus::Error;
  // If no integrality constraints shrink member.
  bool mip = false;
  for  (unsigned int i=0; i < lp.integrality_.size(); i++) {
    if (lp.integrality_[i]) {
      mip = true;
      break;
    }
  }
  if (!mip) {
    run_status = callLpSolver(options, lp, output, message_level, run_quiet);
  } else {
    run_status = callMipSolver(options, lp, output, message_level, run_quiet);
  }
  // duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
  // // std::cout << "time: " << duration << std::endl;
  // resultsFile << options.model_file + " | " + std::to_string(duration) << + "\n";
  // resultsFile.close();
  return (int)run_status;
}

void printHighsVersionCopyright(FILE* output, const int message_level,
                                const char* message) {
  HighsPrintMessage(output, message_level, ML_ALWAYS,
                    "Running HiGHS %d.%d.%d [date: %s, git hash: %s]\n",
                    HIGHS_VERSION_MAJOR, HIGHS_VERSION_MINOR,
                    HIGHS_VERSION_PATCH, HIGHS_COMPILATION_DATE, HIGHS_GITHASH);
  HighsPrintMessage(output, message_level, ML_ALWAYS,
                    "Copyright (c) 2019 ERGO-Code under MIT licence terms\n\n");
#ifdef HiGHSDEV
  // Report on preprocessing macros
  if (message != nullptr) {
    HighsPrintMessage(output, message_level, ML_ALWAYS, "In %s\n", message);
  }
#ifdef OPENMP
  HighsPrintMessage(output, message_level, ML_ALWAYS,
                    "OPENMP           is     defined\n");
#else
  HighsPrintMessage(output, message_level, ML_ALWAYS,
                    "OPENMP           is not defined\n");
#endif

#ifdef SCIP_DEV
  HighsPrintMessage(output, message_level, ML_ALWAYS,
                    "SCIP_DEV         is     defined\n");
#else
  HighsPrintMessage(output, message_level, ML_ALWAYS,
                    "SCIP_DEV         is not defined\n");
#endif

#ifdef HiGHSDEV
  HighsPrintMessage(output, message_level, ML_ALWAYS,
                    "HiGHSDEV         is     defined\n");
#else
  HighsPrintMessage(output, message_level, ML_ALWAYS,
                    "HiGHSDEV         is not defined\n");
#endif
  HighsPrintMessage(output, message_level, ML_ALWAYS,
                    "Built with CMAKE_BUILD_TYPE=%s\n", CMAKE_BUILD_TYPE);
#endif
}

void reportLpStatsOrError(FILE* output, int message_level,
                          const HighsStatus read_status, const HighsLp& lp) {
  if (read_status == HighsStatus::Error) {
    HighsPrintMessage(output, message_level, ML_ALWAYS, "Error loading file\n");
  } else {
    HighsPrintMessage(output, message_level, ML_ALWAYS, "LP       : %s\n",
                      lp.model_name_.c_str());
    HighsPrintMessage(output, message_level, ML_ALWAYS, "Rows     : %d\n",
                      lp.numRow_);
    HighsPrintMessage(output, message_level, ML_ALWAYS, "Cols     : %d\n",
                      lp.numCol_);
    HighsPrintMessage(output, message_level, ML_ALWAYS, "Nonzeros : %d\n",
                      lp.Avalue_.size());
    if (lp.numInt_)
      HighsPrintMessage(output, message_level, ML_ALWAYS, "Integer  : %d\n",
                        lp.numInt_);
    HighsPrintMessage(output, message_level, ML_ALWAYS, "\n");
  }
}

void reportSolvedLpStats(FILE* output, int message_level,
                         const HighsStatus run_status, const Highs& highs) {
  if (run_status == HighsStatus::Error) {
    std::string statusname = HighsStatusToString(run_status);
    HighsPrintMessage(output, message_level, ML_ALWAYS, "HiGHS status: %s\n",
                      statusname.c_str());
  } else {
    HighsPrintMessage(output, message_level, ML_ALWAYS, "\n");
    HighsModelStatus model_status = highs.getModelStatus();
    HighsModelStatus scaled_model_status = highs.getModelStatus(true);
    HighsInfo highs_info = highs.getHighsInfo();
    if (model_status != scaled_model_status) {
      if (scaled_model_status == HighsModelStatus::OPTIMAL) {
        // The scaled model has been solved to optimality, but not the
        // unscaled model, flag this up, but report the scaled model
        // status
        HighsPrintMessage(output, message_level, ML_ALWAYS,
                          "Primal infeasibility: %10.3e (%d)\n",
                          highs_info.max_primal_infeasibility,
                          highs_info.num_primal_infeasibilities);
        HighsPrintMessage(output, message_level, ML_ALWAYS,
                          "Dual   infeasibility: %10.3e (%d)\n",
                          highs_info.max_dual_infeasibility,
                          highs_info.num_dual_infeasibilities);
        model_status = scaled_model_status;
      }
    }
    HighsPrintMessage(output, message_level, ML_ALWAYS,
                      "Model   status      : %s\n",
                      highs.highsModelStatusToString(model_status).c_str());
    /*
    HighsPrintMessage(output, message_level, ML_ALWAYS,
                      "Primal  status      : %s\n",
                      highs.highsPrimalDualStatusToString(highs_info.primal_status).c_str());
    HighsPrintMessage(output, message_level, ML_ALWAYS,
                      "Dual    status      : %s\n",
                      highs.highsPrimalDualStatusToString(highs_info.dual_status).c_str());
    */
    HighsPrintMessage(output, message_level, ML_ALWAYS,
                      "Simplex   iterations: %d\n",
                      highs_info.simplex_iteration_count);
    if (highs_info.ipm_iteration_count)
      HighsPrintMessage(output, message_level, ML_ALWAYS,
                        "IPM       iterations: %d\n",
                        highs_info.ipm_iteration_count);
    if (highs_info.crossover_iteration_count)
      HighsPrintMessage(output, message_level, ML_ALWAYS,
                        "Crossover iterations: %d\n",
                        highs_info.crossover_iteration_count);
    if (model_status == HighsModelStatus::OPTIMAL) {
      double objective_function_value;
      highs.getHighsInfoValue("objective_function_value",
                              objective_function_value);
      HighsPrintMessage(output, message_level, ML_ALWAYS,
                        "Objective value     : %13.6e\n",
                        objective_function_value);
    }

    // Possibly write the solution to a file
    const HighsOptions& options = highs.getHighsOptions();
    if (options.write_solution_to_file)
      highs.writeSolution(options.solution_file, options.write_solution_pretty);
  }

  /*
  highs.writeSolution("", true);
  highs.writeSolution("", false);
  highs.writeHighsInfo("");
  highs.writeHighsInfo("HighsInfo.html");
  */
}

// void foldAndUnfoldIter(const HighsLp& lp){
//   HighsLp alp;
//   HighsAggregate aggregate;
//   HighsSolution solution;
//   HighsBasis basis;
// }

HighsStatus callLpSolver(const HighsOptions& options, HighsLp& lp,
  		         FILE* output, int message_level, bool run_quiet) {
  // New options for aggregate models (work around for const input)
  double time;
  HighsTimer timer;
  HighsOptions alpOpt;
  Highs highs;
  alpOpt.presolve = string("off");
  alpOpt.simplex_scale_strategy = 0;
  alpOpt.model_file = string("No File, Reduced Model");
  alpOpt.solver = string("simplex");
  alpOpt.parallel = string("off");
  // initialize partitioning tool/class
  bool run_highs_clock_already_running = timer.runningRunHighsClock();
  if (!run_highs_clock_already_running) timer.startRunHighsClock();
  double initial_time = timer.readRunHighsClock();
  //std::cout << "refine" << std::endl;
  HighsEquitable ep;
  ep.setup(lp);
  //std::cout << "refine done" << std::endl;
  highs.totPartTime_ += timer.readRunHighsClock() - initial_time;
  // Basis and solution to store from unfold interations
  HighsBasis basis;
  HighsSolution solution;
  HighsTableau tableau;
  // Use aggregator to obtain an aggreated lp and return it
  initial_time = timer.readRunHighsClock();
  //std::cout << "fold start" << std::endl;
  HighsAggregate lpFolder(lp, ep, solution, basis, tableau, false);
  //std::cout << "fold done" << std::endl;
  highs.totFoldTime_ += timer.readRunHighsClock() - initial_time;
  // initial_time = timer.readRunHighsClock();
  HighsLp& alp = lpFolder.getAlp();
  // cout << "fold time: " << lp_folder_time - initial_time << endl;
  // cin.get();
  // Solve LP case.
  HighsStatus return_status = highs.passHighsOptions(alpOpt);
  if (return_status != HighsStatus::OK) {
    if (return_status == HighsStatus::Warning) {
#ifdef HiGHSDEV
      HighsPrintMessage(output, message_level, ML_ALWAYS,
                        "HighsStatus::Warning return from passHighsOptions\n");
#endif
    } else {
      HighsPrintMessage(output, message_level, ML_ALWAYS,
                        "In main: fail return from passHighsOptions\n");
      return return_status;
    }
  }

  if (run_quiet) {
    highs.setHighsLogfile(NULL);
    highs.setHighsOutput(NULL);
  }

  HighsStatus init_status = highs.passModel(alp);
  if (init_status != HighsStatus::OK) {
    if (init_status == HighsStatus::Warning) {
#ifdef HiGHSDEV
      HighsPrintMessage(output, message_level, ML_ALWAYS,
                        "HighsStatus::Warning return setting HighsLp\n");
#endif
    } else {
      HighsPrintMessage(output, message_level, ML_ALWAYS,
                        "Error setting HighsLp\n");
      return HighsStatus::Error;
    }
  }

  HighsStatus write_status;
  //write_status = highs.writeModel("initial.mps");
  /*
  HighsStatus write_status;
  write_status = highs.writeModel("write.mps");
  if (write_status != HighsStatus::OK) {
    if (write_status == HighsStatus::Warning) {
#ifdef HiGHSDEV
      HighsPrintMessage(output, message_level, ML_ALWAYS,
                        "HighsStatus::Warning return from highs.writeModel\n");
#endif
    } else {
      HighsPrintMessage(output, message_level, ML_ALWAYS,
                        "Error return from highs.writeModel\n");
    }
  }
  */

  // Write all the options to an options file
  // highs.writeHighsOptions("Highs.set", false);
  // Write all the options as HTML
  // highs.writeHighsOptions("Highs.html", false);
  // Possibly report options settings
  highs.writeHighsOptions("");  //, false);

  if (run_quiet)
    HighsPrintMessage(output, message_level, ML_ALWAYS,
                      "Before calling highs.run()\n");

  // Run HiGHS.

  HighsStatus run_status = highs.run();
  basis = highs.getBasis();
  solution = highs.getSolution();
  // tableau = highs.getTableau();

  if (run_quiet)
    HighsPrintMessage(output, message_level, ML_ALWAYS,
                      "After calling highs.run()\n");

  reportSolvedLpStats(output, message_level, run_status, highs);
  //highs.totUnfoldTime_ += timer.readRunHighsClock() - initial_time;
  while(!ep.isPartitionDiscrete()){
    //highs = Highs();
    // std::cout << "refine" << std::endl;
    run_highs_clock_already_running = timer.runningRunHighsClock();
    if (!run_highs_clock_already_running) timer.startRunHighsClock();
    initial_time = timer.readRunHighsClock();
    ep.refine();
    highs.totPartTime_ += timer.readRunHighsClock() - initial_time;
    // std::cout << "refine done" << std::endl;
    // std::cin.get();
    // try{
    //   // run_highs_clock_already_running = timer.runningRunHighsClock();
    //   // if (!run_highs_clock_already_running) timer.startRunHighsClock();
    //   // double initial_time = timer.readRunHighsClock();
    // std::cout << "fold" << std::endl;
    initial_time = timer.readRunHighsClock();
    lpFolder = HighsAggregate(lp, ep, solution, basis, tableau, true);
    highs.totFoldTime_ += timer.readRunHighsClock() - initial_time;
    // std::cout << "fold done"  << std::endl;
    // std::cin.get();
    //   // double lp_folder_time = timer.readRunHighsClock();
    //   // time += lp_folder_time - initial_time;
    //   // cout << "fold time: " << lp_folder_time - initial_time << endl;
    //   // cin.get();
      
    // }
    // catch(exception& e){
    //   cout << e.what() << endl;
    // }
    initial_time = timer.readRunHighsClock();
    basis = HighsBasis();
    solution = HighsSolution();
    alpOpt = HighsOptions();
    alpOpt.presolve = string("off");
    alpOpt.simplex_scale_strategy = 0;
    alpOpt.model_file = string("No File, Reduced Model");
    alpOpt.solver = string("simplex");
    alpOpt.parallel = string("off");
    alpOpt.simplex_strategy = SIMPLEX_STRATEGY_UNFOLD;
    return_status = highs.passHighsOptions(alpOpt);
    if (return_status != HighsStatus::OK) {
      if (return_status == HighsStatus::Warning) {
        HighsPrintMessage(output, message_level, ML_ALWAYS,
                        "HighsStatus::Warning return from passHighsOptions\n");
      }
      else {
        HighsPrintMessage(output, message_level, ML_ALWAYS,
                        "In main: fail return from passHighsOptions\n");
        return return_status;
      }
    }
    HighsLp& alp = lpFolder.getAlp();
    HighsBasis& alpBasis = lpFolder.getAlpBasis();
  //   // cout << "Start Basis" << endl;
  //   // for (int j = 0; j < alpBasis.col_status.size(); ++j){
  //   //   cout << "col: " << j << " basis is " << (int)alpBasis.col_status[j] << endl;
  //   // }
  //   // for (int j = 0; j < alpBasis.row_status.size(); ++j){
  //   //   cout << "row: " << j << " basis is " << (int)alpBasis.row_status[j] << endl;
  //   // }
  //   // cout << "\n" << endl;
    HighsStatus init_status = highs.passModel(alp);
    HighsStatus write_status;
    //write_status = highs.writeModel("write.mps");
    HighsStatus basisStatus = highs.setBasis(alpBasis);
    if (init_status != HighsStatus::OK) {
      if (init_status == HighsStatus::Warning) {
        HighsPrintMessage(output, message_level, ML_ALWAYS,
                        "HighsStatus::Warning return setting HighsLp\n");
      }
      else {
        HighsPrintMessage(output, message_level, ML_ALWAYS,
                        "Error setting HighsLp\n");
        return HighsStatus::Error;
      }
    }
    run_status = highs.run();
    basis = highs.getBasis();
  //     // cout << "Finish Basis \n" << endl;
  //     // for (int j = 0; j < basis.col_status.size(); ++j){
  //     //   cout << "col: " << j << " basis is " << (int)basis.col_status[j] << endl;
  //     // }
  //     // for (int j = 0; j < basis.row_status.size(); ++j){
  //     //   cout << "row: " << j << " basis is " << (int)basis.row_status[j] << endl;
  //     // }
  //     // cin.get();
    solution = highs.getSolution();
    highs.totUnfoldTime_ += timer.readRunHighsClock() - initial_time;
    // tableau = highs.getTableau();
  //   // cout << "\n DOKS UNFOLD SOLUTION:\n " << endl;
  //   // for (int i = 0; i < solution.col_value.size(); ++i){
  //   //   cout << "var_" << i << " = " << solution.col_value[i] << endl;
  //   // }
  }
  std::ofstream resultsFile("results.txt", std::ios_base::app);
  std::cout << "Partition time: " << highs.totPartTime_ << std::endl;
  std::cout << "Folding time: " << highs.totFoldTime_ << std::endl;
  std::cout << "Unfolding time: " << highs.totUnfoldTime_ << std::endl;
  double totTime = highs.totUnfoldTime_ + highs.totFoldTime_ + highs.totPartTime_;
  std::cout << "Total OC time: " << totTime << std::endl;
  resultsFile << lp.model_name_ + ".mps " << std::to_string(highs.totPartTime_) + " " << std::to_string(highs.totFoldTime_) + " " << std::to_string(highs.totUnfoldTime_) + " " << std::to_string(totTime) << + "\n";
  // // cout << "\n DOKS UNFOLD SOLUTION:\n " << endl;
  // // for (int i = 0; i < solution.col_value.size(); ++i){
  // //   cout << "var_" << i << " = " << solution.col_value[i] << endl;
  // // }
  // // cout << "\nUnfold iterations: " << highs.totIter_ << endl;
  // cout << "\nFold time total: " << time << endl;
  return run_status;
}

HighsStatus callMipSolver(const HighsOptions& options, const HighsLp& lp,
                          FILE* output, int message_level, bool run_quiet) {
  HighsMipSolver solver(options, lp);
  HighsMipStatus status = solver.runMipSolver();
  switch (status) {
    case HighsMipStatus::kOptimal:
      return HighsStatus::OK;
    default:
      break;
  }
  return HighsStatus::Error;
}
