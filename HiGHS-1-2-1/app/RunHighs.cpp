/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2022 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/*    Authors: Julian Hall, Ivet Galabova, Leona Gottwald and Michael    */
/*    Feldmeier                                                          */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file ../app/RunHighs.cpp
 * @brief HiGHS main
 */
#include <bits/stdc++.h>
#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <sys/types.h>
#include "Highs.h"
#include "lp_data/HighsRuntimeOptions.h"

void printHighsVersionCopyright(const HighsLogOptions& log_options);
void reportModelStatsOrError(const HighsLogOptions& log_options,
                             const HighsStatus read_status,
                             const HighsModel& model);

int main(int argc, char** argv) {
  // Create the Highs instance
  Highs highs;
  const HighsOptions& options = highs.getOptions();
  const HighsLogOptions& log_options = options.log_options;

  // Load user options
  std::string model_file;
  std::string orbit_file;
  std::string basis_file;
  std::string orbit_link_file;
  HighsOptions loaded_options;
  // Set "HiGHS.log" as the default log_file for the app so that
  // log_file has this value if it isn't set in the file
  loaded_options.log_file = "HiGHS.log";
  // When loading the options file, any messages are reported using
  // the default HighsLogOptions
  if (!loadOptions(log_options, argc, argv, loaded_options, model_file,
                   orbit_file, basis_file, orbit_link_file))
    return (int)HighsStatus::kError;
  // Open the app log file - unless output_flag is false, to avoid
  // creating an empty file. It does nothing if its name is "".
  if (loaded_options.output_flag) highs.openLogFile(loaded_options.log_file);

  // Pass the option settings to HiGHS. Only error-checking produces
  // output, but values are checked in loadOptions, so it's safe to
  // call this first so that printHighsVersionCopyright uses reporting
  // settings defined in any options file.
  highs.passOptions(loaded_options);

  printHighsVersionCopyright(log_options);

  // Load the model from model_file
  HighsStatus read_status = highs.readModel(model_file);
  HighsInt original_num_col;
  if (loaded_options.solver == kOrbitalCutGenerationString)
    original_num_col = highs.getLp().num_col_;

  // Load the orbit from the orbit_file
  read_status = highs.readOrbits(orbit_file);
  // // Load the basis from basis_file
  // read_status = highs.readBasis(basis_file);
  // Load in the orbital partition linkers or orbital crossover
  read_status = highs.readOrbitLinkers(orbit_link_file);
  reportModelStatsOrError(log_options, read_status, highs.getModel());
  if (read_status == HighsStatus::kError) return (int)read_status;

  // Solve the model
  HighsStatus run_status = highs.run();
  if (loaded_options.time_file != ""){
    highs.populateTimesInInfo();
    mkdir("./Timings", 0777);
    std::stringstream report_time_file;
    report_time_file << "./Timings/" << loaded_options.time_file << ".csv";
    highs.writeTimes(report_time_file.str().c_str());
  }
  if (run_status == HighsStatus::kError) return (int)run_status;
  if (loaded_options.solver == kOrbitalCutGenerationString){
    // Write the basis to be used in python for setting next lp basis
    highs.writeBasis("../../debugBuild/python_input_basis.txt");
    // Grab the nonbasic variables for gomory cut generation
    HighsLp lp = highs.getLp();
    HighsInt lp_num_col = lp.num_col_;
    HighsInt lp_num_row = lp.num_row_;
    HighsInt lp_nz = lp.a_matrix_.start_.at(lp_num_col);
    std::vector<HighsInt> nonbasic_cols; nonbasic_cols.resize(lp_num_col + lp_num_row);
    std::vector<HighsInt> basic_cols; basic_cols.resize(lp_num_row);
    std::vector<double> col_solution = highs.getSolution().col_value;
    std::vector<double> rhs;
    highs.getNonbasicVariables(nonbasic_cols.data());
    highs.getBasicVariables(basic_cols.data());
    // Grab the submatrix Ar_sub from Ar matrix pertaining to fractional
    // rhs's after solve
    HighsInt num_row_get = 0; 
    std::vector<HighsInt> row_set;
    HighsInt num_row_got; 
    std::vector<double> lower(lp_num_row);      
    std::vector<double> upper(lp_num_row);      
    HighsInt num_row_got_nz;
    std::vector<HighsInt> frac_row_start(lp_num_row + 1);
    std::vector<HighsInt> frac_row_index(lp_nz);
    std::vector<double> frac_row_value(lp_nz);
    for (HighsInt i_row = 0; i_row < lp_num_row; ++i_row){
      HighsInt frac_col = basic_cols.at(i_row);
      double frac_val = col_solution.at(frac_col);
      if (std::abs(std::round(frac_val) - frac_val) > kHighsTiny){
        num_row_get++;
        row_set.push_back(i_row);
        double lb = lp.row_lower_.at(i_row);
        double ub = lp.row_upper_.at(i_row);
        if (lb == -kHighsInf) rhs.push_back(ub);
        else rhs.push_back(lb);
      }
    }
    highs.getRows(num_row_get, row_set.data(), num_row_got, lower.data(),
                  upper.data(), num_row_got_nz, frac_row_start.data(), frac_row_index.data(),
                  frac_row_value.data());
    // Grab reduced rows, generate cuts, lift cuts, write cuts
    std::ofstream cut_file("~/LP/MIPSymmetryCuts/cuts.txt");
    HighsInt frac_matrix_map = 0;
    for (auto i_row : row_set){
      std::string gomory_sense = ">=";
      HighsInt basic_col = basic_cols.at(i_row);
      double b_value = col_solution.at(basic_col);
      std::vector<double> cut(lp_num_col + lp_num_row, 0);
      HighsInt num_row_nz;
      std::vector<double> row_vals; row_vals.resize(lp.num_col_);
      std::vector<HighsInt> row_inds; row_inds.resize(lp.num_row_);
      highs.getReducedRow(0, row_vals.data(), &num_row_nz, row_inds.data());
      for (HighsInt j = 0; j < num_row_nz; ++j){
        HighsInt i_col = row_inds.at(j);
        if (nonbasic_cols.at(i_col)) cut.at(i_col) = (row_vals.at(i_col) - 
                                                      std::floor(row_vals.at(i_col)));
      }
      HighsInt i_col = i_row + lp_num_col;
      HighsInt nonbasic_slack = nonbasic_cols.at(i_col);
      double lb = lp.row_lower_.at(i_row);
      double ub = lp.row_upper_.at(i_row);
      double slack_coeff;
      if (lb == -kHighsInf) slack_coeff = 1.0;
      else slack_coeff = -1.0;
      if (nonbasic_slack){ 
        cut.at(i_col) = slack_coeff;
        std::vector<double> inv_row_vals; inv_row_vals.resize(lp.num_col_);
        std::vector<HighsInt> inv_row_inds; inv_row_inds.resize(lp.num_row_);
        HighsInt inv_num_row_nz;
        highs.getBasisInverseRow(i_row, inv_row_vals.data(), &inv_num_row_nz, inv_row_inds.data());
        highs.getSlackCoeff(inv_row_vals, inv_row_inds, i_row, cut.at(i_col));
        cut.at(i_col) -= std::floor(cut.at(i_col));
      }
      b_value -= std::floor(b_value);
      // Now we need to substitute in variables from non standard form lp 
      // to get rid of nonbasic slack vars if they exist
      if (nonbasic_slack){
        HighsInt row_start = frac_row_start.at(frac_matrix_map);
        HighsInt row_end = frac_row_start.at(frac_matrix_map + 1);
        std::vector<HighsInt> sub_ind(frac_row_index.begin() + row_start, 
                                      frac_row_index.begin() + row_end);
        std::vector<double> sub_val(frac_row_value.begin() + row_start, 
                                    frac_row_value.begin() + row_end);
        for (auto& val : sub_val){
          double move_to_rhs = slack_coeff * (0.0 - val);
          val = move_to_rhs;
        }     
        rhs.at(frac_matrix_map) *= slack_coeff;
        HighsInt sub_map = 0;
        for (auto& j_col : sub_ind){
          cut.at(j_col) += cut.at(i_col) * sub_val.at(sub_map);
          sub_map++;
        }
        b_value += -1.0 * cut.at(i_col) * rhs.at(frac_matrix_map);
        cut.at(i_col) = 0.0;
      }
      // Lift the cut and write it
      std::vector<HighsInt> lift_cut_ind;
      std::vector<double> lift_cut_val;
      std::vector<double> lift_cut(original_num_col);
      highs.getLiftedCut(cut, b_value, lp_num_col, lift_cut_ind,  lift_cut_val, lift_cut);
      for (HighsInt i_col = 0; i_col < original_num_col; ++i_col)
        cut_file << lift_cut.at(i_col) << ",";
      cut_file << ">=" << "," << b_value << "\n";
      frac_matrix_map++;
    }
    // // Grab the reduced inverse row and multiply times slack column to obtain reduced
    // // slack column
    // std::vector<double> inv_row_vals; inv_row_vals.resize(lp.num_col_);
    // std::vector<HighsInt> inv_row_inds; inv_row_inds.resize(lp.num_row_);
    // HighsInt inv_num_row_nz;
    // HighsInt slack_idx = 0;
    // double slack_coeff;
    // highs.getBasisInverseRow(0, inv_row_vals.data(), &inv_num_row_nz, inv_row_inds.data());
    // highs.getSlackCoeff(inv_row_vals, inv_row_inds, slack_idx, slack_coeff);
    // std::cout << "we shall see" << std::endl;
  }

  // Possibly compute the ranging information
  if (options.ranging == kHighsOnString) highs.getRanging();

  // Possibly write the solution to a file
  if (options.write_solution_to_file)
    highs.writeSolution(options.solution_file, options.write_solution_style);

  // Possibly write the model to a file
  if (options.write_model_to_file) {
    HighsStatus write_model_status = highs.writeModel(options.write_model_file);
    if (write_model_status == HighsStatus::kError)
      return (int)write_model_status;  // todo: change to write model error
  }
  return (int)run_status;
}

void printHighsVersionCopyright(const HighsLogOptions& log_options) {
  highsLogUser(log_options, HighsLogType::kInfo,
               "Running HiGHS %d.%d.%d [date: %s, git hash: %s]\n",
               (int)HIGHS_VERSION_MAJOR, (int)HIGHS_VERSION_MINOR,
               (int)HIGHS_VERSION_PATCH, HIGHS_COMPILATION_DATE, HIGHS_GITHASH);
  highsLogUser(log_options, HighsLogType::kInfo,
               "Copyright (c) 2022 ERGO-Code under MIT licence terms\n");
}

void reportModelStatsOrError(const HighsLogOptions& log_options,
                             const HighsStatus read_status,
                             const HighsModel& model) {
  const HighsLp& lp = model.lp_;
  const HighsHessian& hessian = model.hessian_;
  if (read_status == HighsStatus::kError) {
    highsLogUser(log_options, HighsLogType::kInfo, "Error loading file\n");
  } else {
    HighsInt num_integer = 0;
    HighsInt num_semi_continuous = 0;
    HighsInt num_semi_integer = 0;
    for (HighsUInt i = 0; i < lp.integrality_.size(); i++) {
      switch (lp.integrality_[i]) {
        case HighsVarType::kInteger:
          num_integer++;
          break;
        case HighsVarType::kSemiContinuous:
          num_semi_continuous++;
          break;
        case HighsVarType::kSemiInteger:
          num_semi_integer++;
          break;
        default:
          break;
      }
    }
    std::string problem_type;
    const bool non_continuous =
        num_integer + num_semi_continuous + num_semi_integer;
    if (hessian.dim_) {
      if (non_continuous) {
        problem_type = "MIQP";
      } else {
        problem_type = "QP  ";
      }
    } else {
      if (non_continuous) {
        problem_type = "MIP ";
      } else {
        problem_type = "LP  ";
      }
    }
    const HighsInt a_num_nz = lp.a_matrix_.numNz();
    HighsInt q_num_nz = hessian.numNz();
    if (*log_options.log_dev_level) {
      highsLogDev(log_options, HighsLogType::kInfo, "%4s      : %s\n",
                  problem_type.c_str(), lp.model_name_.c_str());
      highsLogDev(log_options, HighsLogType::kInfo,
                  "Rows      : %" HIGHSINT_FORMAT "\n", lp.num_row_);
      highsLogDev(log_options, HighsLogType::kInfo,
                  "Cols      : %" HIGHSINT_FORMAT "\n", lp.num_col_);
      if (q_num_nz) {
        highsLogDev(log_options, HighsLogType::kInfo,
                    "Matrix Nz : %" HIGHSINT_FORMAT "\n", a_num_nz);
        highsLogDev(log_options, HighsLogType::kInfo,
                    "Hessian Nz: %" HIGHSINT_FORMAT "\n", q_num_nz);
      } else {
        highsLogDev(log_options, HighsLogType::kInfo,
                    "Nonzeros  : %" HIGHSINT_FORMAT "\n", a_num_nz);
      }
      if (num_integer)
        highsLogDev(log_options, HighsLogType::kInfo,
                    "Integer  : %" HIGHSINT_FORMAT "\n", num_integer);
      if (num_semi_continuous)
        highsLogDev(log_options, HighsLogType::kInfo,
                    "SemiConts: %" HIGHSINT_FORMAT "\n", num_semi_continuous);
      if (num_semi_integer)
        highsLogDev(log_options, HighsLogType::kInfo,
                    "SemiInt  : %" HIGHSINT_FORMAT "\n", num_semi_integer);
    } else {
      highsLogUser(log_options, HighsLogType::kInfo, "%s",
                   problem_type.c_str());
      if (lp.model_name_.length())
        highsLogUser(log_options, HighsLogType::kInfo, " %s",
                     lp.model_name_.c_str());
      highsLogUser(log_options, HighsLogType::kInfo,
                   " has %" HIGHSINT_FORMAT " rows; %" HIGHSINT_FORMAT " cols",
                   lp.num_row_, lp.num_col_);
      if (q_num_nz) {
        highsLogUser(log_options, HighsLogType::kInfo,
                     "; %" HIGHSINT_FORMAT " matrix nonzeros", a_num_nz);
        highsLogUser(log_options, HighsLogType::kInfo,
                     "; %" HIGHSINT_FORMAT " Hessian nonzeros", q_num_nz);
      } else {
        highsLogUser(log_options, HighsLogType::kInfo,
                     "; %" HIGHSINT_FORMAT " nonzeros", a_num_nz);
      }
      if (num_integer)
        highsLogUser(log_options, HighsLogType::kInfo,
                     "; %" HIGHSINT_FORMAT " integer variables", num_integer);
      if (num_semi_continuous)
        highsLogUser(log_options, HighsLogType::kInfo,
                     "; %" HIGHSINT_FORMAT " semi-continuous variables",
                     num_semi_continuous);
      if (num_semi_integer)
        highsLogUser(log_options, HighsLogType::kInfo,
                     "; %" HIGHSINT_FORMAT " semi-integer variables",
                     num_semi_integer);
      highsLogUser(log_options, HighsLogType::kInfo, "\n");
    }
  }
}
