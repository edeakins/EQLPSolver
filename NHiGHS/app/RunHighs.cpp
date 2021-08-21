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
/**@file ../app/RunHighs.cpp
 * @brief HiGHS main
 */
#include "Highs.h"
#include "lp_data/HighsRuntimeOptions.h"

void printHighsVersionCopyright(const HighsLogOptions& log_options);
void reportModelStatsOrError(const HighsLogOptions& log_options,
                             const HighsStatus read_status,
                             const HighsModel& model);
void reportSolvedLpStats(const HighsLogOptions& log_options,
                         const HighsStatus run_status, Highs& highs);

int main(int argc, char** argv) {
  // Load user options.
  HighsOptions options;
  std::string model_file;
  printHighsVersionCopyright(options.log_options);

  bool options_ok = loadOptions(argc, argv, options, model_file);
  if (!options_ok) return 0;
  Highs highs;
  HighsLp* presolvedLp;
  // Load the model from model_file
  HighsStatus read_status = highs.readModel(model_file);
  reportModelStatsOrError(options.log_options, read_status, highs.getModel());
  if (read_status == HighsStatus::kError)
    return 1;  // todo: change to read error
  //
  // Pass the option settings to HiGHS
  highs.passOptions(options);
  highs.presolve();
  // Start Aggregation
  presolvedLp = &highs.getPresolvedLp();
  HighsOCEquitablePartition OCEquitable(presolvedLp);
  OCEquitable.runToDiscrete();
  // HighsOCAggregate OCAggregator;
  // Solve the model
  HighsStatus run_status = highs.run();
  //
  // Report solution stats
  reportSolvedLpStats(options.log_options, run_status, highs);

  // Possibly write the solution to a file
  if (options.write_solution_to_file)
    highs.writeSolution(options.solution_file, options.write_solution_pretty);

  return (int)run_status;
}

void printHighsVersionCopyright(const HighsLogOptions& log_options) {
  highsLogUser(log_options, HighsLogType::kInfo,
               "Running HiGHS %" HIGHSINT_FORMAT ".%" HIGHSINT_FORMAT
               ".%" HIGHSINT_FORMAT " [date: %s, git hash: %s]\n",
               HIGHS_VERSION_MAJOR, HIGHS_VERSION_MINOR, HIGHS_VERSION_PATCH,
               HIGHS_COMPILATION_DATE, HIGHS_GITHASH);
  highsLogUser(log_options, HighsLogType::kInfo,
               "Copyright (c) 2021 ERGO-Code under MIT licence terms\n\n");
}

void reportModelStatsOrError(const HighsLogOptions& log_options,
                             const HighsStatus read_status,
                             const HighsModel& model) {
  const HighsLp& lp = model.lp_;
  const HighsHessian& hessian = model.hessian_;
  if (read_status == HighsStatus::kError) {
    highsLogUser(log_options, HighsLogType::kInfo, "Error loading file\n");
  } else {
    HighsInt num_int = 0;
    for (HighsUInt i = 0; i < lp.integrality_.size(); i++)
      if (lp.integrality_[i] != HighsVarType::kContinuous) num_int++;
    std::string problem_type;
    if (hessian.dim_) {
      if (num_int) {
        problem_type = "MIQP";
      } else {
        problem_type = "QP  ";
      }
    } else {
      if (num_int) {
        problem_type = "MIP ";
      } else {
        problem_type = "LP  ";
      }
    }
    const HighsInt a_num_nz = lp.Astart_[lp.numCol_];
    HighsInt q_num_nz = 0;
    if (hessian.dim_) q_num_nz = hessian.q_start_[lp.numCol_];
    highsLogUser(log_options, HighsLogType::kInfo, "%4s      : %s\n",
                 problem_type.c_str(), lp.model_name_.c_str());
    highsLogUser(log_options, HighsLogType::kInfo,
                 "Rows      : %" HIGHSINT_FORMAT "\n", lp.numRow_);
    highsLogUser(log_options, HighsLogType::kInfo,
                 "Cols      : %" HIGHSINT_FORMAT "\n", lp.numCol_);
    if (q_num_nz) {
      highsLogUser(log_options, HighsLogType::kInfo,
                   "Matrix Nz : %" HIGHSINT_FORMAT "\n", a_num_nz);
      highsLogUser(log_options, HighsLogType::kInfo,
                   "Hessian Nz: %" HIGHSINT_FORMAT "\n", q_num_nz);
    } else {
      highsLogUser(log_options, HighsLogType::kInfo,
                   "Nonzeros  : %" HIGHSINT_FORMAT "\n", a_num_nz);
    }
    if (num_int)
      highsLogUser(log_options, HighsLogType::kInfo,
                   "Integer  : %" HIGHSINT_FORMAT "\n", num_int);
  }
}

void reportSolvedLpStats(const HighsLogOptions& log_options,
                         const HighsStatus run_status, Highs& highs) {
  if (run_status == HighsStatus::kError) {
    std::string statusname = HighsStatusToString(run_status);
    highsLogUser(log_options, HighsLogType::kInfo, "HiGHS status: %s\n",
                 statusname.c_str());
  } else {
    highsLogUser(log_options, HighsLogType::kInfo, "\n");
    HighsModelStatus model_status = highs.getModelStatus();
    HighsModelStatus scaled_model_status = highs.getModelStatus(true);
    HighsInfo highs_info = highs.getInfo();
    if (model_status != scaled_model_status) {
      if (scaled_model_status == HighsModelStatus::kOptimal) {
        // The scaled model has been solved to optimality, but not the
        // unscaled model, flag this up, but report the scaled model
        // status
        highsLogUser(log_options, HighsLogType::kInfo,
                     "Primal infeasibility: %10.3e (%" HIGHSINT_FORMAT ")\n",
                     highs_info.max_primal_infeasibility,
                     highs_info.num_primal_infeasibilities);
        highsLogUser(log_options, HighsLogType::kInfo,
                     "Dual   infeasibility: %10.3e (%" HIGHSINT_FORMAT ")\n",
                     highs_info.max_dual_infeasibility,
                     highs_info.num_dual_infeasibilities);
        model_status = scaled_model_status;
      }
    }
    highsLogUser(log_options, HighsLogType::kInfo, "Model   status      : %s\n",
                 highs.modelStatusToString(model_status).c_str());
    highsLogUser(log_options, HighsLogType::kInfo, "Primal  status      : %s\n",
                 highs.solutionStatusToString(highs_info.primal_solution_status)
                     .c_str());
    highsLogUser(
        log_options, HighsLogType::kInfo, "Dual    status      : %s\n",
        highs.solutionStatusToString(highs_info.dual_solution_status).c_str());
    highsLogUser(log_options, HighsLogType::kInfo,
                 "Simplex   iterations: %" HIGHSINT_FORMAT "\n",
                 highs_info.simplex_iteration_count);
    if (highs_info.qp_iteration_count)
      highsLogUser(log_options, HighsLogType::kInfo,
                   "QP ASM    iterations: %" HIGHSINT_FORMAT "\n",
                   highs_info.qp_iteration_count);
    if (highs_info.ipm_iteration_count)
      highsLogUser(log_options, HighsLogType::kInfo,
                   "IPM       iterations: %" HIGHSINT_FORMAT "\n",
                   highs_info.ipm_iteration_count);
    if (highs_info.crossover_iteration_count)
      highsLogUser(log_options, HighsLogType::kInfo,
                   "Crossover iterations: %" HIGHSINT_FORMAT "\n",
                   highs_info.crossover_iteration_count);
    if (model_status == HighsModelStatus::kOptimal) {
      double objective_function_value;
      highs.getInfoValue("objective_function_value", objective_function_value);
      highsLogUser(log_options, HighsLogType::kInfo,
                   "Objective value     : %17.10e\n", objective_function_value);
    }
    double run_time = highs.getRunTime();
    highsLogUser(log_options, HighsLogType::kInfo,
                 "HiGHS run time      : %13.2f\n", run_time);
  }
}
