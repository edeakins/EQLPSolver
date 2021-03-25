/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HighsLp.h
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef LP_DATA_HIGHS_LP_H_
#define LP_DATA_HIGHS_LP_H_

#include <cassert>
#include <iostream>
#include <string>
#include <vector>

#include "HConfig.h"
#include "lp_data/HConst.h"        // For HiGHS strategy options
#include "simplex/SimplexConst.h"  // For simplex strategy options

enum class LpAction {
  DUALISE = 0,
  PERMUTE,
  SCALE,
  NEW_COSTS,
  NEW_BOUNDS,
  NEW_BASIS,
  NEW_COLS,
  NEW_ROWS,
  DEL_COLS,
  DEL_ROWS,
  DEL_ROWS_BASIS_OK
};

enum class HighsModelStatus {
  NOTSET = 0,
    LOAD_ERROR,
    MODEL_ERROR,
    MODEL_EMPTY,
    PRESOLVE_ERROR,
    SOLVE_ERROR,
    POSTSOLVE_ERROR,
    PRIMAL_INFEASIBLE,
    PRIMAL_UNBOUNDED,
    OPTIMAL,
    REACHED_DUAL_OBJECTIVE_VALUE_UPPER_BOUND,
    REACHED_TIME_LIMIT,
    REACHED_ITERATION_LIMIT
    };

class HighsLp {
 public:
  // ~HighsLp(){
  //   printf("\ndtor called\n");
  // }
  // Specifically for Deakins (nice to have to clear things out)
  void clear(){
    numCol_ = 0;
    numRow_ = 0;
    numInt_ = 0;
    nnz_ = 0;
    Astart_.resize(0);
    Aindex_.resize(0);
    Avalue_.resize(0);
    colCost_.resize(0);
    colLower_.resize(0);
    colUpper_.resize(0);
    rowLower_.resize(0);
    rowUpper_.resize(0);
    sense_ = 1;
    offset_ = 0;
    model_name_ = "";
    lp_name_ = "";
    row_names_.resize(0);
    col_names_.resize(0);
    integrality_.resize(0);
  }
  // Iteration Data
  int masterIter = 0;
  // Model data
  int numCol_ = 0;
  // int addNumCol_ = 0;
  int numRow_ = 0;
  // int addNumRow_ = 0;
  int numInt_ = 0;
  int nnz_ = 0;
  // int realNumCol_ = 0;
  // int realNumRow_ = 0;
  int numXCol_;
  int numSCol_;
  int numRCol_;

  std::vector<int> Astart_;
  std::vector<int> addARstart_;
  std::vector<int> Aindex_;
  std::vector<int> addARindex_;
  std::vector<double> Avalue_;
  std::vector<double> addARvalue_;
  std::vector<double> colCost_;
  std::vector<double> colLower_;
  std::vector<double> colUpper_;
  std::vector<double> rowLower_;
  std::vector<double> rowUpper_;
  std::vector<double> linkLower_;
  std::vector<double> linkUpper_;
  std::vector<int> linkers;
  std::vector<int> artificialVariables;
  std::vector<bool> activeColorHistory;
  std::vector<int> cell;
  std::vector<int> labels;
  std::vector<int> fronts;


  // sense 1 = minimize, -1 = maximize
  int sense_ = 1;
  double offset_ = 0;
  int numLinkers;

  std::string model_name_ = "";
  std::string lp_name_ = "";

  std::vector<std::string> row_names_;
  std::vector<std::string> col_names_;

  std::vector<int> integrality_;
  std::vector<int> rowColor;

  bool operator==(const HighsLp& lp) {
    if (numCol_ != lp.numCol_ || numRow_ != lp.numRow_ || nnz_ != lp.nnz_ ||
        sense_ != lp.sense_ || offset_ != lp.offset_ ||
        model_name_ != lp.model_name_)
      return false;

    if (row_names_ != lp.row_names_ || col_names_ != lp.col_names_)
      return false;

    if (colCost_ != lp.colCost_) return false;

    if (colUpper_ != lp.colUpper_ || colLower_ != lp.colLower_ ||
        rowUpper_ != lp.rowUpper_ || rowLower_ != lp.rowLower_)
      return false;

    if (Astart_ != lp.Astart_ || Aindex_ != lp.Aindex_ || Avalue_ != lp.Avalue_)
      return false;

    return true;
  }

  int numRealRows;
  int numRealCols;
};

// Cost, column and row scaling factors
struct HighsScale {
  bool is_scaled_;
  double cost_;
  std::vector<double> col_;
  std::vector<double> row_;
};

struct SimplexBasis {
  // The basis for the simplex method consists of basicIndex,
  // nonbasicFlag and nonbasicMove. If HighsSimplexLpStatus has_basis
  // is true then it is assumed that basicIndex_ and nonbasicFlag_ are
  // self-consistent and correpond to the dimensions of an associated
  // HighsLp, but the basis matrix B is not necessarily nonsingular.
  std::vector<int> basicIndex_;
  std::vector<int> nonbasicFlag_;
  std::vector<int> nonbasicMove_;
};

struct HighsSimplexLpStatus {
  // Status of LP solved by the simplex method and its data
  bool valid = false;
  bool is_dualised = false;
  bool is_permuted = false;
  bool scaling_tried = false;
  bool has_basis = false;  // The simplex LP has a valid simplex basis
  bool has_matrix_col_wise = false;  // The HMatrix column-wise matrix is valid
  bool has_matrix_row_wise = false;  // The HMatrix row-wise matrix is valid
  bool has_factor_arrays =
      false;  // Has the arrays for the representation of B^{-1}
  bool has_dual_steepest_edge_weights = false;  // The DSE weights are known
  bool has_nonbasic_dual_values = false;  // The nonbasic dual values are known
  bool has_basic_primal_values = false;   // The basic primal values are known
  bool has_invert =
      false;  // The representation of B^{-1} corresponds to the current basis
  bool has_fresh_invert = false;  // The representation of B^{-1} corresponds to
                                  // the current basis and is fresh
  bool has_fresh_rebuild = false;  // The data are fresh from rebuild
  bool has_dual_objective_value =
      false;  // The dual objective function value is known
  bool has_primal_objective_value =
      false;  // The dual objective function value is known
  SimplexSolutionStatus solution_status =
      SimplexSolutionStatus::UNSET;  // The solution status is UNSET
};

struct HighsSimplexInfo {
  bool initialised = false;
  // Simplex information regarding primal and dual solution, objective
  // and iteration counts for this Highs Model Object. This is
  // information which should be retained from one run to the next in
  // order to provide hot starts.
  //
  // Part of working model which are assigned and populated as much as
  // possible when a model is being defined

  // workCost: Originally just costs from the model but, in solve(), may
  // be perturbed or set to alternative values in Phase I??
  //
  // workDual: Values of the dual variables corresponding to
  // workCost. Latter not known until solve() is called since B^{-1}
  // is required to compute them. Knowledge of them is indicated by
  // has_nonbasic_dual_values
  //
  // workShift: WTF
  //
  std::vector<double> workCost_;
  std::vector<double> workDual_;
  std::vector<double> workShift_;

  // workLower/workUpper: Originally just lower (upper) bounds from
  // the model but, in solve(), may be perturbed or set to
  // alternative values in Phase I??
  //
  // workRange: Distance between lower and upper bounds
  //
  // workValue: Values of the nonbasic variables corresponding to
  // workLower/workUpper and the basis. Always known.
  //
  std::vector<double> workLower_;
  std::vector<double> workUpper_;
  std::vector<double> workRange_;
  std::vector<double> workValue_;

  // baseLower/baseUpper/baseValue: Lower and upper bounds on the
  // basic variables and their values. Latter not known until solve()
  // is called since B^{-1} is required to compute them. Knowledge of
  // them is indicated by has_basic_primal_values
  //
  std::vector<double> baseLower_;
  std::vector<double> baseUpper_;
  std::vector<double> baseValue_;
  //
  // Vectors of random reals for column cost perturbation, a random
  // permutation of all indices for CHUZR and a random permutation of
  // column indices for permuting the columns
  std::vector<double> numTotRandomValue_;
  std::vector<int> numTotPermutation_;
  std::vector<int> numColPermutation_;

  std::vector<int> devex_index_;

  // Values of iClock for simplex timing clocks
  std::vector<int> clock_;
  //
  // Options from HighsOptions for the simplex solver
  int simplex_strategy;
  int dual_edge_weight_strategy;
  int primal_edge_weight_strategy;
  int price_strategy;

  //  double primal_feasibility_tolerance;
  //  double dual_feasibility_tolerance;
  bool perturb_costs;
  int update_limit;
  //  int iteration_limit;
  //  double dual_objective_value_upper_bound;

  // Internal options - can't be changed externally
  bool store_squared_primal_infeasibility = false;
  bool allow_primal_flips_for_dual_feasibility = true;
#ifndef HiGHSDEV
  bool analyse_lp_solution = false; //true;// 
#else  
  bool analyse_lp_solution = true;
  // Options for reporting timing
  bool report_simplex_inner_clock = false;
  bool report_simplex_outer_clock = false;
  bool report_simplex_phases_clock = false;
  bool report_HFactor_clock = false;
  // Option for analysing the LP simplex iterations, INVERT time and rebuild
  // time
  bool analyse_lp = false;
  bool analyse_iterations = false;
  bool analyse_invert_form = false;
  bool analyse_invert_condition = false;
  bool analyse_invert_time = false;
  bool analyse_rebuild_time = false;
#endif
  // Simplex runtime information
  int costs_perturbed = 0;
  // Records of cumulative iteration counts - updated at the end of a phase
  int dual_phase1_iteration_count = 0;
  int dual_phase2_iteration_count = 0;
  int primal_phase1_iteration_count = 0;
  int primal_phase2_iteration_count = 0;

  int min_threads = 1;
  int num_threads = 1;
  int max_threads = HIGHS_THREAD_LIMIT;
  
  // Cutoff for PAMI
  double pami_cutoff = 0.95;

  // Info on PAMI iterations
  int multi_iteration = 0;

  // Number of UPDATE operations performed - should be zeroed when INVERT is
  // performed
  int update_count;
  // Value of dual objective - only set when computed from scratch in dual
  // rebuild()
  double dual_objective_value;
  // Value of primal objective - only set when computed from scratch in primal
  // rebuild()
  double primal_objective_value;

  // Value of dual objective that is updated in dual simplex solver
  double updated_dual_objective_value;
  // Value of primal objective that is updated in primal simplex solver
  double updated_primal_objective_value;
  // Number of logical variables in the basis
  int num_basic_logicals;

#ifdef HiGHSDEV
  // Analysis of INVERT
  int num_invert = 0;
  // Analysis of INVERT form
  int num_kernel = 0;
  int num_major_kernel = 0;
  const double major_kernel_relative_dim_threshhold = 0.1;
  double max_kernel_dim = 0;
  double sum_kernel_dim = 0;
  double running_average_kernel_dim = 0;
  double sum_invert_fill_factor = 0;
  double sum_kernel_fill_factor = 0;
  double sum_major_kernel_fill_factor = 0;
  double running_average_invert_fill_factor = 1;
  double running_average_kernel_fill_factor = 1;
  double running_average_major_kernel_fill_factor = 1;

  int total_inverts;
  double total_invert_time;
  double invert_condition = 1;
#endif

  /*
#ifdef HiGHSDEV
  // Move this to Simplex class once it's created
  vector<int> historyColumnIn;
  vector<int> historyColumnOut;
  vector<double> historyAlpha;
#endif
  */
};

struct HighsSolutionParams {
  // Input to solution analysis method
  double primal_feasibility_tolerance;
  double dual_feasibility_tolerance;
  int simplex_iteration_count = 0;
  int ipm_iteration_count = 0;
  int crossover_iteration_count = 0;
  int primal_status = PrimalDualStatus::STATUS_NOTSET;
  int dual_status = PrimalDualStatus::STATUS_NOTSET;
  // Output from solution analysis method
  double objective_function_value;
  int num_primal_infeasibilities;
  double sum_primal_infeasibilities;
  double max_primal_infeasibility;
  int num_dual_infeasibilities;
  double sum_dual_infeasibilities;
  double max_dual_infeasibility;
};

struct HighsSolution {
  void clear(){
    col_value.resize(0);
    col_dual.resize(0);
    row_value.resize(0);
    row_dual.resize(0);
  }
  std::vector<double> col_value;
  std::vector<double> col_dual;
  std::vector<double> row_value;
  std::vector<double> row_dual;
};

// This is a csc collection of the previous iterations reduced A matrix
// that we use with the EP to construct the current LP
// Struct to contain tableau info
struct HighsTableau{
  int nnz;
  int numCol;
  int numRow;
  int numXCol;
  int numSCol;
  int numRCol;
  std::vector<int> Aindex;
  std::vector<int> Astart;
  std::vector<double> Avalue;
  std::vector<double> basicValue;
  std::vector<int> basicIndex;
  std::vector<int> reps;
};

// To be the basis representation given back to the user. Values of
// HighsBasisStatus are defined in HConst.h
struct HighsBasis {
  int numCol_;
  int numRow_;
  bool valid_ = false;
  std::vector<HighsBasisStatus> col_status;
  std::vector<HighsBasisStatus> row_status;
};

struct HighsPreviousSolutionInfo{
  HighsBasis basis;
  HighsSolution solution;
};

struct HighsRanging {
  std::vector<double> colCostRangeUpValue_;
  std::vector<double> colCostRangeUpObjective_;
  std::vector<int> colCostRangeUpInCol_;
  std::vector<int> colCostRangeUpOutCol_;
  std::vector<double> colCostRangeDnValue_;
  std::vector<double> colCostRangeDnObjective_;
  std::vector<int> colCostRangeDnInCol_;
  std::vector<int> colCostRangeDnOutCol_;
  std::vector<double> rowBoundRangeUpValue_;
  std::vector<double> rowBoundRangeUpObjective_;
  std::vector<int> rowBoundRangeUpInCol_;
  std::vector<int> rowBoundRangeUpOutCol_;
  std::vector<double> rowBoundRangeDnValue_;
  std::vector<double> rowBoundRangeDnObjective_;
  std::vector<int> rowBoundRangeDnInCol_;
  std::vector<int> rowBoundRangeDnOutCol_;
};

// Make sure the dimensions of solution are the same as numRow_ and numCol_.
bool isSolutionConsistent(const HighsLp& lp, const HighsSolution& solution);

// If debug this method terminates the program when the status is not OK. If
// standard build it only prints a message.
// void checkStatus(HighsStatus status);

#endif
