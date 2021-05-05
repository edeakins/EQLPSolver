
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
  if (lp.numCol_ > 80000 || lp.numRow_ > 80000){
    HighsStatus run_status;
    const char *fileName = "DHiGHS_mittleman_timings.csv";
    std::ofstream resultsFile(fileName, std::ios_base::app);
    std::ifstream in(fileName);
    std::string name = options.model_file.c_str();
    std::string pTime = std::to_string(0);
    std::string fTime = std::to_string(0);
    std::string uTime = std::to_string(0);
    std::string tTime = std::to_string(0);
    std::string objval = std::to_string(0);
    std::string big = "to big";
    name.erase(0,33);
    std::string outP = name + "," + pTime + "," + fTime + "," + uTime + "," + tTime + "," + objval + "," + big + "\n";
    if (in.peek() == std::ifstream::traits_type::eof()){
      std::string column0 = "Instance";
      std::string column1 = "Partition Time";
      std::string column2 = "Fold Time";
      std::string column3 = "Unfold Time";
      std::string column4 = "Total Time";
      std::string column5 = "Objective";
      std::string outCols = column0 + "," + column1 + "," + column2 + "," + column3 + "," + column4 + "," + column5 + "\n";
      resultsFile << outCols;
    }
    resultsFile << outP;
    resultsFile.close();
    return run_status;
  }
  double time;
  // Timer, options, highs solver, write_status, run_status
  HighsTimer timer;
  HighsOptions alpOpt;
  Highs highs;
  HighsStatus run_status;
  HighsStatus write_status;
  HighsStatus return_status;
  HighsStatus init_status;
  HighsStatus basis_status;
  // HighsLp class to contain aggregates and HighsBasis class to contain bases
  HighsLp* alp;
  HighsBasis* alpBasis;
  // Basis and solution to store from unfold interations
  HighsBasis basis;
  HighsSolution solution;
  init_status = highs.passModel(lp);
  write_status = highs.writeModel("Original.lp");
  // Set options
  alpOpt.presolve = string("off");
  alpOpt.simplex_scale_strategy = 0;
  alpOpt.model_file = string("No File, Reduced Model");
  // alpOpt.solver = string("simplex");
  alpOpt.parallel = string("off");
  alpOpt.time_limit = (double)3600;
  // initialize partitioning tool/class
  bool run_highs_clock_already_running = timer.runningRunHighsClock();
  if (!run_highs_clock_already_running) timer.startRunHighsClock();
  double initial_time = timer.readRunHighsClock();
  HighsEquitable ep(lp);
  struct eq_part *partitions;
  partitions = ep.refine();
  int numRefinements = ep.getNumRefinements();
  highs.totPartTime_ += timer.readRunHighsClock() - initial_time;
  // Use aggregator to get first aggregate lp (level 0)
  initial_time = timer.readRunHighsClock();
  HighsAggregate lpFolder(lp, partitions, solution, basis, numRefinements);
  highs.totFoldTime_ += timer.readRunHighsClock() - initial_time;
  alp = lpFolder.getAlp();
  alpBasis = lpFolder.getBasis();
  // Intitial solve of level 0 aggregate
  return_status = highs.passHighsOptions(alpOpt);
  init_status = highs.passModel(*alp);
  write_status = highs.writeModel("level_0.lp");
  initial_time = timer.readRunHighsClock();
  run_status = highs.run(); 
  highs.totUnfoldTime_ += timer.readRunHighsClock() - initial_time; // Add this timer to highs
  basis = highs.getBasis();
  solution = highs.getSolution();
  // cout << "\n DOKS UNFOLD SOLUTION:\n " << endl;
  // for (int i = 0; i < solution.col_value.size(); ++i){
  //   cout << "var_" << i << " = " << solution.col_value[i] << endl;
  // }
  // Start loop for level 1+ aggregates
  alpOpt.simplex_strategy = SIMPLEX_STRATEGY_UNFOLD;
  return_status = highs.passHighsOptions(alpOpt);
  // std::string mitName = options.model_file.c_str();
  // mitName.erase(0,24);
  // mitName.erase(mitName.length()-4, mitName.length());
  // std::string totName = "../../../mittleman_runs_for_jim/DHiGHS_" + mitName + ".csv";
  // char fileName[totName.length() + 1];
  // std::strcpy(fileName, totName.c_str()); 
  // std::ofstream resultsFile(fileName, std::ios_base::app);
  // std::string column0 = "Instance";
  // std::string column1 = "Master Iteration";
  // std::string column2 = "r Pivots";
  // std::string column3 = "Fold Time";
  // std::string column4 = "Initial Invert Time";
  // std::string column5 = "Pivoting Time";
  // std::string outCols = column0 + "," + column1 + "," + column2 + "," + column3 + "," + column4 + "," + column5 + "\n";
  // resultsFile << outCols;
  for (int i = 1; i < numRefinements; ++i){
    // double foldTime = timer.readRunHighsClock();
    initial_time = timer.readRunHighsClock();
    int links = lpFolder.update(solution, basis);
    if (!links) break;
    if (links == 1) continue;
    // foldTime = timer.readRunHighsClock() - foldTime;
    highs.totFoldTime_ += timer.readRunHighsClock() - initial_time;
    alp = lpFolder.getAlp();
    alpBasis = lpFolder.getBasis();
    init_status = highs.passModel(*alp);
    basis_status = highs.setBasis(*alpBasis);
    std::string name = "level_" + std::to_string(i) + ".lp";
    // if (i == 1 or i == 2)
    write_status = highs.writeModel(name);
    initial_time = timer.readRunHighsClock();
    run_status = highs.run(); 
    highs.totUnfoldTime_ += timer.readRunHighsClock() - initial_time; // Add this timer to highs
    basis = highs.getBasis();
    solution = highs.getSolution();
    std::cout << "Basis for Linkers" << std::endl;
    for (int i = alp->numCol_ - alp->numLinkers_; i < alp->numCol_; ++i)
      if ((int)basis.col_status[i] != 1)
        std::cout << i << std::endl;
  }
    
      // for (int i = 0; i < solution.col_value.size(); ++i)
      //   std::cout << "x_" << i << " = " << solution.col_value[i] << std::endl;
      
    
    // std::string iTime = std::to_string(highs.getLp().invertTime);
    // std::string pTime = std::to_string(highs.getLp().pivotTime);
    // std::string uTime = std::to_string(foldTime);
    // std::string mIter = std::to_string(i);
    // std::string rIter = std::to_string(highs.getLp().unfoldIter);
    // std::string mName = options.model_file.c_str(); 
    // mName.erase(0,24);
    // mName.erase(mName.length() - 4, mName.length());
    // std::string outN = mName + "," + mIter + "," + rIter + "," + uTime + "," + iTime + "," + pTime + "\n";
    // resultsFile << outN;
  // resultsFile.close();
  double obj = 0;
  for (int i = 0; i < solution.col_value.size(); ++i)
    obj += solution.col_value[i] * alp->colCost_[i];
  // cout << "Partition Time: " << highs.totPartTime_ << endl;
  // cout << "Fold Time: " << highs.totFoldTime_ << endl;
  // cout << "Unfold Time: " << highs.totUnfoldTime_ << endl;
  // cout << "Total Time: " << highs.totPartTime_ + highs.totFoldTime_ + highs.totUnfoldTime_ << endl;
  // Report timings
  const char *fileName = "DHiGHS_mittleman_timings.csv";
  std::ofstream resultsFile(fileName, std::ios_base::app);
  std::ifstream in(fileName);
  std::string name = options.model_file.c_str();
  std::string pTime = std::to_string(highs.totPartTime_);
  std::string fTime = std::to_string(highs.totFoldTime_);
  std::string uTime = std::to_string(highs.totUnfoldTime_);
  std::string tTime = std::to_string(highs.totPartTime_ + highs.totFoldTime_ + highs.totUnfoldTime_);
  std::string objval = std::to_string(obj);
  // name.erase(0,33);
  name.erase(0,42);
  std::string outP = name + "," + pTime + "," + fTime + "," + uTime + "," + tTime + "," + objval + "\n";
  if (in.peek() == std::ifstream::traits_type::eof()){
    std::string column0 = "Instance";
    std::string column1 = "Partition Time";
    std::string column2 = "Fold Time";
    std::string column3 = "Unfold Time";
    std::string column4 = "Total Time";
    std::string column5 = "Objective";
    std::string outCols = column0 + "," + column1 + "," + column2 + "," + column3 + "," + column4 + "," + column5 + "\n";
    resultsFile << outCols;
  }
  resultsFile << outP;
  resultsFile.close(); 
  
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
