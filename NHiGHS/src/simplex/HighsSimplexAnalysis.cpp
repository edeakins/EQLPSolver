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
/**@file simplex/HighsSimplexAnalysis.cpp
 * @brief
 */
#include <cmath>
//#include <cstdio>
#include <iomanip>

#include "HConfig.h"
#include "simplex/FactorTimer.h"
#include "simplex/HEkkDebug.h"
#include "simplex/HFactor.h"
#include "simplex/HighsSimplexAnalysis.h"
#include "simplex/SimplexTimer.h"

#ifdef OPENMP
#include "omp.h"
#endif

void HighsSimplexAnalysis::setup(const std::string lp_name, const HighsLp& lp,
                                 const HighsOptions& options,
                                 const HighsInt simplex_iteration_count_) {
  // Copy Problem size
  numRow = lp.numRow_;
  numCol = lp.numCol_;
  numTot = numRow + numCol;
  model_name_ = lp.model_name_;
  lp_name_ = lp_name;
  // Set up analysis logic short-cuts
  analyse_lp_data = kHighsAnalysisLevelModelData & options.highs_analysis_level;
  analyse_simplex_data =
      kHighsAnalysisLevelSolverData & options.highs_analysis_level;
  analyse_simplex_time =
      kHighsAnalysisLevelSolverTime & options.highs_analysis_level;
  analyse_factor_data =
      kHighsAnalysisLevelNlaData & options.highs_analysis_level;
  analyse_factor_time =
      kHighsAnalysisLevelNlaTime & options.highs_analysis_level;

  // Set up the thread clocks
  HighsInt omp_max_threads = 1;
#ifdef OPENMP
  omp_max_threads = omp_get_max_threads();
#endif
  if (analyse_simplex_time) {
    for (HighsInt i = 0; i < omp_max_threads; i++) {
      HighsTimerClock clock(timer_reference);
      thread_simplex_clocks.push_back(clock);
    }
  }
  if (analyse_factor_time) {
    for (HighsInt i = 0; i < omp_max_threads; i++) {
      HighsTimerClock clock(timer_reference);
      thread_factor_clocks.push_back(clock);
    }
    pointer_serial_factor_clocks = &thread_factor_clocks[0];
  } else {
    pointer_serial_factor_clocks = NULL;
  }

  // Copy tolerances from options
  allow_dual_steepest_edge_to_devex_switch =
      options.simplex_dual_edge_weight_strategy ==
      kSimplexDualEdgeWeightStrategyChoose;
  dual_steepest_edge_weight_log_error_threshold =
      options.dual_steepest_edge_weight_log_error_threshold;
  //
  AnIterIt0 = simplex_iteration_count_;
  AnIterCostlyDseFq = 0;
  AnIterNumCostlyDseIt = 0;
  // Copy messaging parameter from options
  messaging(options.log_options);
  // Initialise the densities
  col_aq_density = 0;
  row_ep_density = 0;
  row_ap_density = 0;
  row_DSE_density = 0;
  col_basic_feasibility_change_density = 0;
  row_basic_feasibility_change_density = 0;
  col_BFRT_density = 0;
  primal_col_density = 0;
  // Set the row_dual_density to 1 since it's assumed all costs are at
  // least perturbed from zero, if not initially nonzero
  dual_col_density = 1;
  // Set up the data structures for scatter data
  tran_stage.resize(NUM_TRAN_STAGE_TYPE);
  tran_stage[TRAN_STAGE_FTRAN_LOWER].name_ = "FTRAN lower";
  tran_stage[TRAN_STAGE_FTRAN_UPPER_FT].name_ = "FTRAN upper FT";
  tran_stage[TRAN_STAGE_FTRAN_UPPER].name_ = "FTRAN upper";
  tran_stage[TRAN_STAGE_BTRAN_UPPER].name_ = "BTRAN upper";
  tran_stage[TRAN_STAGE_BTRAN_UPPER_FT].name_ = "BTRAN upper FT";
  tran_stage[TRAN_STAGE_BTRAN_LOWER].name_ = "BTRAN lower";
  for (HighsInt tran_stage_type = 0; tran_stage_type < NUM_TRAN_STAGE_TYPE;
       tran_stage_type++) {
    TranStageAnalysis& stage = tran_stage[tran_stage_type];
    initialiseScatterData(20, stage.rhs_density_);
    stage.num_decision_ = 0;
    stage.num_wrong_original_sparse_decision_ = 0;
    stage.num_wrong_original_hyper_decision_ = 0;
    stage.num_wrong_new_sparse_decision_ = 0;
    stage.num_wrong_new_hyper_decision_ = 0;
  }
  original_start_density_tolerance.resize(NUM_TRAN_STAGE_TYPE);
  new_start_density_tolerance.resize(NUM_TRAN_STAGE_TYPE);
  historical_density_tolerance.resize(NUM_TRAN_STAGE_TYPE);
  predicted_density_tolerance.resize(NUM_TRAN_STAGE_TYPE);

  for (HighsInt tran_stage_type = 0; tran_stage_type < NUM_TRAN_STAGE_TYPE;
       tran_stage_type++) {
    original_start_density_tolerance[tran_stage_type] = 0.05;
    new_start_density_tolerance[tran_stage_type] = 0.05;
  }
  historical_density_tolerance[TRAN_STAGE_FTRAN_LOWER] = 0.15;
  historical_density_tolerance[TRAN_STAGE_FTRAN_UPPER] = 0.10;
  historical_density_tolerance[TRAN_STAGE_BTRAN_UPPER] = 0.10;
  historical_density_tolerance[TRAN_STAGE_BTRAN_LOWER] = 0.15;
  predicted_density_tolerance[TRAN_STAGE_FTRAN_LOWER] = 0.10;
  predicted_density_tolerance[TRAN_STAGE_FTRAN_UPPER] = 0.10;
  predicted_density_tolerance[TRAN_STAGE_BTRAN_UPPER] = 0.10;
  predicted_density_tolerance[TRAN_STAGE_BTRAN_LOWER] = 0.10;

  // Initialise the measures used to analyse accuracy of steepest edge weights
  //
  const HighsInt dual_edge_weight_strategy =
      options.simplex_dual_edge_weight_strategy;
  if (dual_edge_weight_strategy == kSimplexDualEdgeWeightStrategyChoose ||
      dual_edge_weight_strategy == kSimplexDualEdgeWeightStrategySteepestEdge ||
      dual_edge_weight_strategy ==
          kSimplexDualEdgeWeightStrategySteepestEdgeUnitInitial) {
    // Initialise the measures used to analyse accuracy of steepest edge weights
    num_dual_steepest_edge_weight_check = 0;
    num_dual_steepest_edge_weight_reject = 0;
    num_wrong_low_dual_steepest_edge_weight = 0;
    num_wrong_high_dual_steepest_edge_weight = 0;
    average_frequency_low_dual_steepest_edge_weight = 0;
    average_frequency_high_dual_steepest_edge_weight = 0;
    average_log_low_dual_steepest_edge_weight_error = 0;
    average_log_high_dual_steepest_edge_weight_error = 0;
    max_average_frequency_low_dual_steepest_edge_weight = 0;
    max_average_frequency_high_dual_steepest_edge_weight = 0;
    max_sum_average_frequency_extreme_dual_steepest_edge_weight = 0;
    max_average_log_low_dual_steepest_edge_weight_error = 0;
    max_average_log_high_dual_steepest_edge_weight_error = 0;
    max_sum_average_log_extreme_dual_steepest_edge_weight_error = 0;
  }
  num_devex_framework = 0;

  num_iteration_report_since_last_header = -1;
  num_invert_report_since_last_header = -1;
  // Set values used to skip an iteration report when column (row)
  // choice has found nothing in primal (dual) simplex
  entering_variable = -1;
  pivotal_row_index = -1;

  // Set following averages to illegal values so that first average is
  // set equal to first value
  average_num_threads = -1;
  average_fraction_of_possible_minor_iterations_performed = -1;
  sum_multi_chosen = 0;
  sum_multi_finished = 0;

  if (analyse_simplex_data) {
    AnIterPrevIt = simplex_iteration_count_;

    AnIterOpRec* AnIter;
    AnIter = &AnIterOp[ANALYSIS_OPERATION_TYPE_BTRAN_FULL];
    AnIter->AnIterOpName = "BTRAN Full";
    AnIter = &AnIterOp[ANALYSIS_OPERATION_TYPE_PRICE_FULL];
    AnIter->AnIterOpName = "PRICE Full";
    AnIter = &AnIterOp[ANALYSIS_OPERATION_TYPE_BTRAN_BASIC_FEASIBILITY_CHANGE];
    AnIter->AnIterOpName = "BTRAN BcFsCg";
    AnIter = &AnIterOp[ANALYSIS_OPERATION_TYPE_PRICE_BASIC_FEASIBILITY_CHANGE];
    AnIter->AnIterOpName = "PRICE BcFsCg";
    AnIter = &AnIterOp[ANALYSIS_OPERATION_TYPE_BTRAN_EP];
    AnIter->AnIterOpName = "BTRAN e_p";
    AnIter = &AnIterOp[ANALYSIS_OPERATION_TYPE_PRICE_AP];
    AnIter->AnIterOpName = "PRICE a_p";
    AnIter = &AnIterOp[ANALYSIS_OPERATION_TYPE_FTRAN];
    AnIter->AnIterOpName = "FTRAN";
    AnIter = &AnIterOp[ANALYSIS_OPERATION_TYPE_FTRAN_BFRT];
    AnIter->AnIterOpName = "FTRAN BFRT";
    AnIter = &AnIterOp[ANALYSIS_OPERATION_TYPE_FTRAN_DSE];
    AnIter->AnIterOpName = "FTRAN DSE";
    for (HighsInt k = 0; k < NUM_ANALYSIS_OPERATION_TYPE; k++) {
      AnIter = &AnIterOp[k];
      if ((k == ANALYSIS_OPERATION_TYPE_PRICE_AP) ||
          (k == ANALYSIS_OPERATION_TYPE_PRICE_BASIC_FEASIBILITY_CHANGE) ||
          (k == ANALYSIS_OPERATION_TYPE_PRICE_FULL)) {
        AnIter->AnIterOpHyperCANCEL = 1.0;
        AnIter->AnIterOpHyperTRAN = 1.0;
        AnIter->AnIterOpRsDim = numCol;
      } else {
        if ((k == ANALYSIS_OPERATION_TYPE_BTRAN_EP) ||
            (k == ANALYSIS_OPERATION_TYPE_BTRAN_BASIC_FEASIBILITY_CHANGE) ||
            (k == ANALYSIS_OPERATION_TYPE_BTRAN_FULL)) {
          AnIter->AnIterOpHyperCANCEL = kHyperCancel;
          AnIter->AnIterOpHyperTRAN = kHyperBtranU;
        } else {
          AnIter->AnIterOpHyperCANCEL = kHyperCancel;
          AnIter->AnIterOpHyperTRAN = kHyperFtranL;
        }
        AnIter->AnIterOpRsDim = numRow;
      }
      AnIter->AnIterOpNumCa = 0;
      AnIter->AnIterOpNumHyperOp = 0;
      AnIter->AnIterOpNumHyperRs = 0;
      AnIter->AnIterOpSumLog10RsDensity = 0;
      initialiseValueDistribution("", "density ", 1e-8, 1.0, 10.0,
                                  AnIter->AnIterOp_density);
    }
    HighsInt last_rebuild_reason = kRebuildReasonCount - 1;
    for (HighsInt k = 1; k <= last_rebuild_reason; k++) AnIterNumInvert[k] = 0;
    num_col_price = 0;
    num_row_price = 0;
    num_row_price_with_switch = 0;
    HighsInt last_dual_edge_weight_mode =
        (HighsInt)DualEdgeWeightMode::kSteepestEdge;
    for (HighsInt k = 0; k <= last_dual_edge_weight_mode; k++)
      AnIterNumEdWtIt[k] = 0;
    AnIterTraceNumRec = 0;
    AnIterTraceIterDl = 1;
    AnIterTraceRec* lcAnIter = &AnIterTrace[0];
    lcAnIter->AnIterTraceIter = AnIterIt0;
    lcAnIter->AnIterTraceTime = timer_->getWallTime();
    initialiseValueDistribution("Primal step summary", "", 1e-16, 1e16, 10.0,
                                primal_step_distribution);
    initialiseValueDistribution("Dual step summary", "", 1e-16, 1e16, 10.0,
                                dual_step_distribution);
    initialiseValueDistribution("Simplex pivot summary", "", 1e-8, 1e16, 10.0,
                                simplex_pivot_distribution);
    initialiseValueDistribution("Factor pivot threshold summary", "",
                                kMinPivotThreshold, kMaxPivotThreshold,
                                kPivotThresholdChangeFactor,
                                factor_pivot_threshold_distribution);
    initialiseValueDistribution("Numerical trouble summary", "", 1e-16, 1.0,
                                10.0, numerical_trouble_distribution);
    initialiseValueDistribution("", "1 ", 1e-16, 1e16, 10.0,
                                cost_perturbation1_distribution);
    initialiseValueDistribution("", "2 ", 1e-16, 1e16, 10.0,
                                cost_perturbation2_distribution);
    initialiseValueDistribution("FTRAN upper sparse summary - before", "", 1e-8,
                                1.0, 10.0, before_ftran_upper_sparse_density);
    initialiseValueDistribution("FTRAN upper sparse summary - after", "", 1e-8,
                                1.0, 10.0, before_ftran_upper_hyper_density);
    initialiseValueDistribution("FTRAN upper hyper-sparse summary - before", "",
                                1e-8, 1.0, 10.0, ftran_upper_sparse_density);
    initialiseValueDistribution("FTRAN upper hyper-sparse summary - after", "",
                                1e-8, 1.0, 10.0, ftran_upper_hyper_density);
    initialiseValueDistribution("Cleanup dual change summary", "", 1e-16, 1e16,
                                10.0, cleanup_dual_change_distribution);
    initialiseValueDistribution("Cleanup primal change summary", "", 1e-16,
                                1e16, 10.0, cleanup_primal_step_distribution);
    initialiseValueDistribution("Cleanup primal step summary", "", 1e-16, 1e16,
                                10.0, cleanup_primal_change_distribution);
    initialiseValueDistribution("Cleanup dual step summary", "", 1e-16, 1e16,
                                10.0, cleanup_dual_step_distribution);
  }

  if (analyse_simplex_time) {
    SimplexTimer simplex_timer;
    for (HighsTimerClock& clock : thread_simplex_clocks)
      simplex_timer.initialiseSimplexClocks(clock);
  }
  if (analyse_factor_time) {
    FactorTimer factor_timer;
    for (HighsTimerClock& clock : thread_factor_clocks)
      factor_timer.initialiseFactorClocks(clock);
  }
}

void HighsSimplexAnalysis::messaging(const HighsLogOptions& log_options_) {
  log_options = log_options_;
}

void HighsSimplexAnalysis::updateOperationResultDensity(
    const double local_density, double& density) {
  density = (1 - kRunningAverageMultiplier) * density +
            kRunningAverageMultiplier * local_density;
}

void HighsSimplexAnalysis::iterationReport() {
  if ((HighsInt)kIterationReportLogType > *log_options.log_dev_level) return;
  const bool header = (num_iteration_report_since_last_header < 0) ||
                      (num_iteration_report_since_last_header > 49);
  if (header) {
    iterationReport(header);
    num_iteration_report_since_last_header = 0;
  }
  iterationReport(false);
}

void HighsSimplexAnalysis::invertReport() {
  const bool header = (num_invert_report_since_last_header < 0) ||
                      (num_invert_report_since_last_header > 49) ||
                      (num_iteration_report_since_last_header >= 0);
  if (header) {
    invertReport(header);
    num_invert_report_since_last_header = 0;
  }
  invertReport(false);
  // Force an iteration report header if this is an INVERT report without an
  // rebuild_reason
  if (!rebuild_reason) num_iteration_report_since_last_header = -1;
}

void HighsSimplexAnalysis::invertReport(const bool header) {
  analysis_log = std::unique_ptr<std::stringstream>(new std::stringstream());
  reportAlgorithmPhaseIterationObjective(header);
  if (analyse_simplex_data) {
    if (simplex_strategy == kSimplexStrategyDualMulti) {
      // Report on threads and PAMI
      reportThreads(header);
      reportMulti(header);
    }
    reportDensity(header);
    reportInvert(header);
    //  reportCondition(header);
  }
  reportInfeasibility(header);
  highsLogUser(log_options, HighsLogType::kInfo, "%s\n",
               analysis_log->str().c_str());
  if (!header) num_invert_report_since_last_header++;
}

void HighsSimplexAnalysis::dualSteepestEdgeWeightError(
    const double computed_edge_weight, const double updated_edge_weight) {
  const double kWeightErrorThreshold = 4.0;
  const bool accept_weight =
      updated_edge_weight >= kAcceptDseWeightThreshold * computed_edge_weight;
  HighsInt low_weight_error = 0;
  HighsInt high_weight_error = 0;
  double weight_error;
  string error_type = "  OK";
  num_dual_steepest_edge_weight_check++;
  if (!accept_weight) num_dual_steepest_edge_weight_reject++;
  if (updated_edge_weight < computed_edge_weight) {
    // Updated weight is low
    weight_error = computed_edge_weight / updated_edge_weight;
    if (weight_error > kWeightErrorThreshold) {
      low_weight_error = 1;
      error_type = " Low";
    }
    average_log_low_dual_steepest_edge_weight_error =
        0.99 * average_log_low_dual_steepest_edge_weight_error +
        0.01 * log(weight_error);
  } else {
    // Updated weight is correct or high
    weight_error = updated_edge_weight / computed_edge_weight;
    if (weight_error > kWeightErrorThreshold) {
      high_weight_error = 1;
      error_type = "High";
    }
    average_log_high_dual_steepest_edge_weight_error =
        0.99 * average_log_high_dual_steepest_edge_weight_error +
        0.01 * log(weight_error);
  }
  average_frequency_low_dual_steepest_edge_weight =
      0.99 * average_frequency_low_dual_steepest_edge_weight +
      0.01 * low_weight_error;
  average_frequency_high_dual_steepest_edge_weight =
      0.99 * average_frequency_high_dual_steepest_edge_weight +
      0.01 * high_weight_error;
  max_average_frequency_low_dual_steepest_edge_weight =
      max(max_average_frequency_low_dual_steepest_edge_weight,
          average_frequency_low_dual_steepest_edge_weight);
  max_average_frequency_high_dual_steepest_edge_weight =
      max(max_average_frequency_high_dual_steepest_edge_weight,
          average_frequency_high_dual_steepest_edge_weight);
  max_sum_average_frequency_extreme_dual_steepest_edge_weight =
      max(max_sum_average_frequency_extreme_dual_steepest_edge_weight,
          average_frequency_low_dual_steepest_edge_weight +
              average_frequency_high_dual_steepest_edge_weight);
  max_average_log_low_dual_steepest_edge_weight_error =
      max(max_average_log_low_dual_steepest_edge_weight_error,
          average_log_low_dual_steepest_edge_weight_error);
  max_average_log_high_dual_steepest_edge_weight_error =
      max(max_average_log_high_dual_steepest_edge_weight_error,
          average_log_high_dual_steepest_edge_weight_error);
  max_sum_average_log_extreme_dual_steepest_edge_weight_error =
      max(max_sum_average_log_extreme_dual_steepest_edge_weight_error,
          average_log_low_dual_steepest_edge_weight_error +
              average_log_high_dual_steepest_edge_weight_error);
  if (analyse_simplex_data) {
    const bool report_weight_error = false;
    if (report_weight_error && weight_error > 0.5 * kWeightErrorThreshold) {
      printf(
          "DSE Wt Ck |%8" HIGHSINT_FORMAT "| OK = %1d (%4" HIGHSINT_FORMAT
          " / %6" HIGHSINT_FORMAT
          ") (c %10.4g, u %10.4g, er %10.4g "
          "- "
          "%s): Low (Fq %10.4g, Er %10.4g); High (Fq%10.4g, Er%10.4g) | %10.4g "
          "%10.4g %10.4g %10.4g %10.4g %10.4g\n",
          simplex_iteration_count, accept_weight,
          num_dual_steepest_edge_weight_check,
          num_dual_steepest_edge_weight_reject, computed_edge_weight,
          updated_edge_weight, weight_error, error_type.c_str(),
          average_frequency_low_dual_steepest_edge_weight,
          average_log_low_dual_steepest_edge_weight_error,
          average_frequency_high_dual_steepest_edge_weight,
          average_log_high_dual_steepest_edge_weight_error,
          max_average_frequency_low_dual_steepest_edge_weight,
          max_average_frequency_high_dual_steepest_edge_weight,
          max_sum_average_frequency_extreme_dual_steepest_edge_weight,
          max_average_log_low_dual_steepest_edge_weight_error,
          max_average_log_high_dual_steepest_edge_weight_error,
          max_sum_average_log_extreme_dual_steepest_edge_weight_error);
    }
  }
}

bool HighsSimplexAnalysis::switchToDevex() {
  bool switch_to_devex = false;
  // Firstly consider switching on the basis of NLA cost
  const double kAnIterCostlyDseMeasureLimit = 1000.0;
  const double kAnIterCostlyDseMnDensity = 0.01;
  const double kAnIterFracNumTotItBfSw = 0.1;
  const double kAnIterFracNumCostlyDseItbfSw = 0.05;
  double AnIterCostlyDseMeasureDen;
  AnIterCostlyDseMeasureDen =
      max(max(row_ep_density, col_aq_density), row_ap_density);
  if (AnIterCostlyDseMeasureDen > 0) {
    AnIterCostlyDseMeasure = row_DSE_density / AnIterCostlyDseMeasureDen;
    AnIterCostlyDseMeasure = AnIterCostlyDseMeasure * AnIterCostlyDseMeasure;
  } else {
    AnIterCostlyDseMeasure = 0;
  }
  bool CostlyDseIt = AnIterCostlyDseMeasure > kAnIterCostlyDseMeasureLimit &&
                     row_DSE_density > kAnIterCostlyDseMnDensity;
  AnIterCostlyDseFq = (1 - kRunningAverageMultiplier) * AnIterCostlyDseFq;
  if (CostlyDseIt) {
    AnIterNumCostlyDseIt++;
    AnIterCostlyDseFq += kRunningAverageMultiplier * 1.0;
    HighsInt lcNumIter = simplex_iteration_count - AnIterIt0;
    // Switch to Devex if at least 5% of the (at least) 0.1NumTot iterations
    // have been costly
    switch_to_devex =
        allow_dual_steepest_edge_to_devex_switch &&
        (AnIterNumCostlyDseIt > lcNumIter * kAnIterFracNumCostlyDseItbfSw) &&
        (lcNumIter > kAnIterFracNumTotItBfSw * numTot);
    if (switch_to_devex) {
      highsLogUser(log_options, HighsLogType::kInfo,
                   "Switch from DSE to Devex after %" HIGHSINT_FORMAT
                   " costly DSE iterations of %" HIGHSINT_FORMAT
                   " with "
                   "densities C_Aq = %11.4g; R_Ep = %11.4g; R_Ap = "
                   "%11.4g; DSE = %11.4g\n",
                   AnIterNumCostlyDseIt, lcNumIter, col_aq_density,
                   row_ep_density, row_ap_density, row_DSE_density);
    }
  }
  if (!switch_to_devex) {
    // Secondly consider switching on the basis of weight accuracy
    double dse_weight_error_measure =
        average_log_low_dual_steepest_edge_weight_error +
        average_log_high_dual_steepest_edge_weight_error;
    double dse_weight_error_threshold =
        dual_steepest_edge_weight_log_error_threshold;
    switch_to_devex = allow_dual_steepest_edge_to_devex_switch &&
                      dse_weight_error_measure > dse_weight_error_threshold;
    if (switch_to_devex) {
      highsLogUser(log_options, HighsLogType::kInfo,
                   "Switch from DSE to Devex with log error measure of %g > "
                   "%g = threshold",
                   dse_weight_error_measure, dse_weight_error_threshold);
    }
  }
  return switch_to_devex;
}

bool HighsSimplexAnalysis::predictEndDensity(const HighsInt tran_stage_type,
                                             const double start_density,
                                             double& end_density) {
  return predictFromScatterData(tran_stage[tran_stage_type].rhs_density_,
                                start_density, end_density);
}

void HighsSimplexAnalysis::afterTranStage(
    const HighsInt tran_stage_type, const double start_density,
    const double end_density, const double historical_density,
    const double predicted_end_density,
    const bool use_solve_sparse_original_HFactor_logic,
    const bool use_solve_sparse_new_HFactor_logic) {
  TranStageAnalysis& stage = tran_stage[tran_stage_type];
  const double rp = false;
  const double kMaxHyperDensity = 0.1;

  if (predicted_end_density > 0) {
    stage.num_decision_++;
    if (end_density <= kMaxHyperDensity) {
      // Should have done hyper-sparse TRAN
      if (use_solve_sparse_original_HFactor_logic) {
        // Original logic makes wrong decision to use sparse TRAN
        if (rp) {
          printf("Original: Wrong sparse: ");
          const double start_density_tolerance =
              original_start_density_tolerance[tran_stage_type];
          const double this_historical_density_tolerance =
              historical_density_tolerance[tran_stage_type];
          if (start_density > start_density_tolerance) {
            printf("(start = %10.4g >  %4.2f)  or ", start_density,
                   start_density_tolerance);
          } else {
            printf(" start = %10.4g              ", start_density);
          }
          if (historical_density > this_historical_density_tolerance) {
            printf("(historical = %10.4g  > %4.2f); ", historical_density,
                   this_historical_density_tolerance);
          } else {
            printf(" historical = %10.4g           ", historical_density);
          }
          printf("end = %10.4g", end_density);
          if (end_density < 0.1 * historical_density) printf(" !! OG");
          printf("\n");
        }
        stage.num_wrong_original_sparse_decision_++;
      }
      if (use_solve_sparse_new_HFactor_logic) {
        // New logic makes wrong decision to use sparse TRAN
        if (rp) {
          printf("New     : Wrong sparse: ");
          const double start_density_tolerance =
              original_start_density_tolerance[tran_stage_type];
          const double end_density_tolerance =
              predicted_density_tolerance[tran_stage_type];
          if (start_density > start_density_tolerance) {
            printf("(start = %10.4g >  %4.2f)  or ", start_density,
                   start_density_tolerance);
          } else {
            printf(" start = %10.4g                       ", start_density);
          }
          if (predicted_end_density > end_density_tolerance) {
            printf("( predicted = %10.4g  > %4.2f); ", predicted_end_density,
                   end_density_tolerance);
          } else {
            printf("  predicted = %10.4g           ", predicted_end_density);
          }
          printf("end = %10.4g", end_density);
          if (end_density < 0.1 * predicted_end_density) printf(" !! NW");
          printf("\n");
        }
        stage.num_wrong_new_sparse_decision_++;
      }
    } else {
      // Should have done sparse TRAN
      if (!use_solve_sparse_original_HFactor_logic) {
        // Original logic makes wrong decision to use hyper TRAN
        if (rp) {
          printf(
              "Original: Wrong  hyper: (start = %10.4g <= %4.2f) and "
              "(historical = %10.4g <= %4.2f); end = %10.4g",
              start_density, original_start_density_tolerance[tran_stage_type],
              historical_density, historical_density_tolerance[tran_stage_type],
              end_density);
          if (end_density > 10.0 * historical_density) printf(" !! OG");
          printf("\n");
        }
        stage.num_wrong_original_hyper_decision_++;
      }
      if (!use_solve_sparse_new_HFactor_logic) {
        // New logic makes wrong decision to use hyper TRAN
        if (rp) {
          printf(
              "New     : Wrong  hyper: (start = %10.4g <= %4.2f) and ( "
              "predicted = %10.4g <= %4.2f); end = %10.4g",
              start_density, new_start_density_tolerance[tran_stage_type],
              predicted_end_density,
              predicted_density_tolerance[tran_stage_type], end_density);
          if (end_density > 10.0 * predicted_end_density) printf(" !! NW");
          printf("\n");
        }
        stage.num_wrong_new_hyper_decision_++;
      }
    }
  }
  updateScatterData(start_density, end_density, stage.rhs_density_);
  regressScatterData(stage.rhs_density_);
}

void HighsSimplexAnalysis::simplexTimerStart(const HighsInt simplex_clock,
                                             const HighsInt thread_id) {
  if (!analyse_simplex_time) return;
  // assert(analyse_simplex_time);
  thread_simplex_clocks[thread_id].timer_.start(
      thread_simplex_clocks[thread_id].clock_[simplex_clock]);
}

void HighsSimplexAnalysis::simplexTimerStop(const HighsInt simplex_clock,
                                            const HighsInt thread_id) {
  if (!analyse_simplex_time) return;
  // assert(analyse_simplex_time);
  thread_simplex_clocks[thread_id].timer_.stop(
      thread_simplex_clocks[thread_id].clock_[simplex_clock]);
}

bool HighsSimplexAnalysis::simplexTimerRunning(const HighsInt simplex_clock,
                                               const HighsInt thread_id) {
  if (!analyse_simplex_time) return false;
  // assert(analyse_simplex_time);
  return thread_simplex_clocks[thread_id].timer_.clock_start
             [thread_simplex_clocks[thread_id].clock_[simplex_clock]] < 0;
}

HighsInt HighsSimplexAnalysis::simplexTimerNumCall(const HighsInt simplex_clock,
                                                   const HighsInt thread_id) {
  if (!analyse_simplex_time) return -1;
  // assert(analyse_simplex_time);
  return thread_simplex_clocks[thread_id]
      .timer_
      .clock_num_call[thread_simplex_clocks[thread_id].clock_[simplex_clock]];
}

double HighsSimplexAnalysis::simplexTimerRead(const HighsInt simplex_clock,
                                              const HighsInt thread_id) {
  if (!analyse_simplex_time) return -1.0;
  // assert(analyse_simplex_time);
  return thread_simplex_clocks[thread_id].timer_.read(
      thread_simplex_clocks[thread_id].clock_[simplex_clock]);
}

HighsTimerClock* HighsSimplexAnalysis::getThreadFactorTimerClockPointer() {
  HighsTimerClock* factor_timer_clock_pointer = NULL;
  if (analyse_factor_time) {
    HighsInt thread_id = 0;
#ifdef OPENMP
    thread_id = omp_get_thread_num();
#endif
    factor_timer_clock_pointer = &thread_factor_clocks[thread_id];
  }
  return factor_timer_clock_pointer;
}

void HighsSimplexAnalysis::iterationRecord() {
  assert(analyse_simplex_data);
  HighsInt AnIterCuIt = simplex_iteration_count;
  if (rebuild_reason > 0) AnIterNumInvert[rebuild_reason]++;
  if (AnIterCuIt > AnIterPrevIt)
    AnIterNumEdWtIt[(HighsInt)edge_weight_mode] += (AnIterCuIt - AnIterPrevIt);

  AnIterTraceRec& lcAnIter = AnIterTrace[AnIterTraceNumRec];
  //  if (simplex_iteration_count ==
  //  AnIterTraceIterRec[AnIterTraceNumRec]+AnIterTraceIterDl) {
  if (simplex_iteration_count == lcAnIter.AnIterTraceIter + AnIterTraceIterDl) {
    if (AnIterTraceNumRec == kAnIterTraceMaxNumRec) {
      for (HighsInt rec = 1; rec <= kAnIterTraceMaxNumRec / 2; rec++)
        AnIterTrace[rec] = AnIterTrace[2 * rec];
      AnIterTraceNumRec = AnIterTraceNumRec / 2;
      AnIterTraceIterDl = AnIterTraceIterDl * 2;
    } else {
      AnIterTraceNumRec++;
      AnIterTraceRec& lcAnIter = AnIterTrace[AnIterTraceNumRec];
      lcAnIter.AnIterTraceIter = simplex_iteration_count;
      lcAnIter.AnIterTraceTime = timer_->getWallTime();
      if (average_fraction_of_possible_minor_iterations_performed > 0) {
        lcAnIter.AnIterTraceMulti =
            average_fraction_of_possible_minor_iterations_performed;
      } else {
        lcAnIter.AnIterTraceMulti = 0;
      }
      lcAnIter.AnIterTraceDensity[ANALYSIS_OPERATION_TYPE_FTRAN] =
          col_aq_density;
      lcAnIter.AnIterTraceDensity[ANALYSIS_OPERATION_TYPE_BTRAN_EP] =
          row_ep_density;
      lcAnIter.AnIterTraceDensity[ANALYSIS_OPERATION_TYPE_PRICE_AP] =
          row_ap_density;
      lcAnIter.AnIterTraceDensity[ANALYSIS_OPERATION_TYPE_FTRAN_BFRT] =
          col_aq_density;
      if (edge_weight_mode == DualEdgeWeightMode::kSteepestEdge) {
        lcAnIter.AnIterTraceDensity[ANALYSIS_OPERATION_TYPE_FTRAN_DSE] =
            row_DSE_density;
        lcAnIter.AnIterTraceCostlyDse = AnIterCostlyDseMeasure;
      } else {
        lcAnIter.AnIterTraceDensity[ANALYSIS_OPERATION_TYPE_FTRAN_DSE] = 0;
        lcAnIter.AnIterTraceCostlyDse = 0;
      }
      lcAnIter.AnIterTrace_dual_edge_weight_mode = (HighsInt)edge_weight_mode;
    }
  }
  AnIterPrevIt = AnIterCuIt;
  updateValueDistribution(primal_step, cleanup_primal_step_distribution);
  updateValueDistribution(dual_step, cleanup_dual_step_distribution);
  updateValueDistribution(primal_step, primal_step_distribution);
  updateValueDistribution(dual_step, dual_step_distribution);
  updateValueDistribution(pivot_value_from_column, simplex_pivot_distribution);
  updateValueDistribution(factor_pivot_threshold,
                          factor_pivot_threshold_distribution);
  // Only update the distribution of legal values for
  // numerical_trouble. Illegal values are set in PAMI since it's not
  // known in minor iterations
  if (numerical_trouble >= 0)
    updateValueDistribution(numerical_trouble, numerical_trouble_distribution);
}

void HighsSimplexAnalysis::iterationRecordMajor() {
  assert(analyse_simplex_data);
  sum_multi_chosen += multi_chosen;
  sum_multi_finished += multi_finished;
  assert(multi_chosen > 0);
  const double fraction_of_possible_minor_iterations_performed =
      1.0 * multi_finished / multi_chosen;
  if (average_fraction_of_possible_minor_iterations_performed < 0) {
    average_fraction_of_possible_minor_iterations_performed =
        fraction_of_possible_minor_iterations_performed;
  } else {
    average_fraction_of_possible_minor_iterations_performed =
        kRunningAverageMultiplier *
            fraction_of_possible_minor_iterations_performed +
        (1 - kRunningAverageMultiplier) *
            average_fraction_of_possible_minor_iterations_performed;
  }
  if (average_num_threads < 0) {
    average_num_threads = num_threads;
  } else {
    average_num_threads = kRunningAverageMultiplier * num_threads +
                          (1 - kRunningAverageMultiplier) * average_num_threads;
  }
}

void HighsSimplexAnalysis::operationRecordBefore(
    const HighsInt operation_type, const HVector& vector,
    const double historical_density) {
  assert(analyse_simplex_data);
  operationRecordBefore(operation_type, vector.count, historical_density);
}

void HighsSimplexAnalysis::operationRecordBefore(
    const HighsInt operation_type, const HighsInt current_count,
    const double historical_density) {
  double current_density = 1.0 * current_count / numRow;
  AnIterOpRec& AnIter = AnIterOp[operation_type];
  AnIter.AnIterOpNumCa++;
  if (current_density <= AnIter.AnIterOpHyperCANCEL &&
      historical_density <= AnIter.AnIterOpHyperTRAN)
    AnIter.AnIterOpNumHyperOp++;
}

void HighsSimplexAnalysis::operationRecordAfter(const HighsInt operation_type,
                                                const HVector& vector) {
  assert(analyse_simplex_data);
  operationRecordAfter(operation_type, vector.count);
}

void HighsSimplexAnalysis::operationRecordAfter(const HighsInt operation_type,
                                                const HighsInt result_count) {
  AnIterOpRec& AnIter = AnIterOp[operation_type];
  const double result_density = 1.0 * result_count / AnIter.AnIterOpRsDim;
  if (result_density <= kHyperResult) AnIter.AnIterOpNumHyperRs++;
  if (result_density > 0) {
    AnIter.AnIterOpSumLog10RsDensity += log(result_density) / log(10.0);
  } else {
    /*
    // TODO Investigate these zero norms
    double vectorNorm = 0;

    for (HighsInt index = 0; index < AnIter.AnIterOpRsDim; index++) {
      double vectorValue = vector.array[index];
      vectorNorm += vectorValue * vectorValue;
    }
    vectorNorm = sqrt(vectorNorm);
    printf("Strange: operation %s has result density = %g: ||vector|| = %g\n",
    AnIter.AnIterOpName.c_str(), result_density, vectorNorm);
    */
  }
  updateValueDistribution(result_density, AnIter.AnIterOp_density);
}

void HighsSimplexAnalysis::summaryReport() {
  assert(analyse_simplex_data);
  HighsInt AnIterNumIter = simplex_iteration_count - AnIterIt0;
  if (AnIterNumIter <= 0) return;
  printf("\nAnalysis of %" HIGHSINT_FORMAT " iterations (%" HIGHSINT_FORMAT
         " to %" HIGHSINT_FORMAT ")\n",
         AnIterNumIter, AnIterIt0 + 1, simplex_iteration_count);
  if (AnIterNumIter <= 0) return;
  HighsInt lc_EdWtNumIter;
  lc_EdWtNumIter = AnIterNumEdWtIt[(HighsInt)DualEdgeWeightMode::kSteepestEdge];
  if (lc_EdWtNumIter > 0)
    printf("DSE for %12" HIGHSINT_FORMAT " (%3" HIGHSINT_FORMAT
           "%%) iterations\n",
           lc_EdWtNumIter, (100 * lc_EdWtNumIter) / AnIterNumIter);
  lc_EdWtNumIter = AnIterNumEdWtIt[(HighsInt)DualEdgeWeightMode::kDevex];
  if (lc_EdWtNumIter > 0)
    printf("Dvx for %12" HIGHSINT_FORMAT " (%3" HIGHSINT_FORMAT
           "%%) iterations\n",
           lc_EdWtNumIter, (100 * lc_EdWtNumIter) / AnIterNumIter);
  lc_EdWtNumIter = AnIterNumEdWtIt[(HighsInt)DualEdgeWeightMode::kDantzig];
  if (lc_EdWtNumIter > 0)
    printf("Dan for %12" HIGHSINT_FORMAT " (%3" HIGHSINT_FORMAT
           "%%) iterations\n",
           lc_EdWtNumIter, (100 * lc_EdWtNumIter) / AnIterNumIter);
  for (HighsInt k = 0; k < NUM_ANALYSIS_OPERATION_TYPE; k++) {
    AnIterOpRec& AnIter = AnIterOp[k];
    HighsInt lcNumCa = AnIter.AnIterOpNumCa;
    printf("\n%-10s performed %" HIGHSINT_FORMAT " times\n",
           AnIter.AnIterOpName.c_str(), AnIter.AnIterOpNumCa);
    if (lcNumCa > 0) {
      HighsInt lcHyperOp = AnIter.AnIterOpNumHyperOp;
      HighsInt lcHyperRs = AnIter.AnIterOpNumHyperRs;
      HighsInt pctHyperOp = (100 * lcHyperOp) / lcNumCa;
      HighsInt pctHyperRs = (100 * lcHyperRs) / lcNumCa;
      double lcRsDensity =
          pow(10.0, AnIter.AnIterOpSumLog10RsDensity / lcNumCa);
      HighsInt lcAnIterOpRsDim = AnIter.AnIterOpRsDim;
      HighsInt lcNumNNz = lcRsDensity * lcAnIterOpRsDim;
      printf("%12" HIGHSINT_FORMAT
             " hyper-sparse operations (%3" HIGHSINT_FORMAT "%%)\n",
             lcHyperOp, pctHyperOp);
      printf("%12" HIGHSINT_FORMAT
             " hyper-sparse results    (%3" HIGHSINT_FORMAT "%%)\n",
             lcHyperRs, pctHyperRs);
      printf("%12g density of result (%" HIGHSINT_FORMAT " / %" HIGHSINT_FORMAT
             " nonzeros)\n",
             lcRsDensity, lcNumNNz, lcAnIterOpRsDim);
      logValueDistribution(log_options, AnIter.AnIterOp_density,
                           AnIter.AnIterOpRsDim);
    }
  }
  HighsInt NumInvert = 0;

  HighsInt last_rebuild_reason = kRebuildReasonCount - 1;
  for (HighsInt k = 1; k <= last_rebuild_reason; k++)
    NumInvert += AnIterNumInvert[k];
  if (NumInvert > 0) {
    HighsInt lcNumInvert = 0;
    printf("\nInvert    performed %" HIGHSINT_FORMAT
           " times: average frequency = %" HIGHSINT_FORMAT "\n",
           NumInvert, AnIterNumIter / NumInvert);
    lcNumInvert = AnIterNumInvert[kRebuildReasonUpdateLimitReached];
    if (lcNumInvert > 0)
      printf("%12" HIGHSINT_FORMAT " (%3" HIGHSINT_FORMAT
             "%%) Invert operations due to update limit reached\n",
             lcNumInvert, (100 * lcNumInvert) / NumInvert);
    lcNumInvert = AnIterNumInvert[kRebuildReasonSyntheticClockSaysInvert];
    if (lcNumInvert > 0)
      printf("%12" HIGHSINT_FORMAT " (%3" HIGHSINT_FORMAT
             "%%) Invert operations due to pseudo-clock\n",
             lcNumInvert, (100 * lcNumInvert) / NumInvert);
    lcNumInvert = AnIterNumInvert[kRebuildReasonPossiblyOptimal];
    if (lcNumInvert > 0)
      printf("%12" HIGHSINT_FORMAT " (%3" HIGHSINT_FORMAT
             "%%) Invert operations due to possibly optimal\n",
             lcNumInvert, (100 * lcNumInvert) / NumInvert);
    lcNumInvert = AnIterNumInvert[kRebuildReasonPossiblyPrimalUnbounded];
    if (lcNumInvert > 0)
      printf("%12" HIGHSINT_FORMAT " (%3" HIGHSINT_FORMAT
             "%%) Invert operations due to possibly primal unbounded\n",
             lcNumInvert, (100 * lcNumInvert) / NumInvert);
    lcNumInvert = AnIterNumInvert[kRebuildReasonPossiblyDualUnbounded];
    if (lcNumInvert > 0)
      printf("%12" HIGHSINT_FORMAT " (%3" HIGHSINT_FORMAT
             "%%) Invert operations due to possibly dual unbounded\n",
             lcNumInvert, (100 * lcNumInvert) / NumInvert);
    lcNumInvert = AnIterNumInvert[kRebuildReasonPossiblySingularBasis];
    if (lcNumInvert > 0)
      printf("%12" HIGHSINT_FORMAT " (%3" HIGHSINT_FORMAT
             "%%) Invert operations due to possibly singular basis\n",
             lcNumInvert, (100 * lcNumInvert) / NumInvert);
    lcNumInvert =
        AnIterNumInvert[kRebuildReasonPrimalInfeasibleInPrimalSimplex];
    if (lcNumInvert > 0)
      printf("%12" HIGHSINT_FORMAT " (%3" HIGHSINT_FORMAT
             "%%) Invert operations due to primal infeasible in primal "
             "simplex\n",
             lcNumInvert, (100 * lcNumInvert) / NumInvert);
  }
  HighsInt suPrice = num_col_price + num_row_price + num_row_price_with_switch;
  if (suPrice > 0) {
    printf("\n%12" HIGHSINT_FORMAT " Price operations:\n", suPrice);
    printf("%12" HIGHSINT_FORMAT " Col Price      (%3" HIGHSINT_FORMAT "%%)\n",
           num_col_price, (100 * num_col_price) / suPrice);
    printf("%12" HIGHSINT_FORMAT " Row Price      (%3" HIGHSINT_FORMAT "%%)\n",
           num_row_price, (100 * num_row_price) / suPrice);
    printf("%12" HIGHSINT_FORMAT " Row PriceWSw   (%3" HIGHSINT_FORMAT "%%)\n",
           num_row_price_with_switch,
           (100 * num_row_price_with_switch / suPrice));
  }
  printf("\n%12" HIGHSINT_FORMAT " (%3" HIGHSINT_FORMAT
         "%%) costly DSE        iterations\n",
         AnIterNumCostlyDseIt, (100 * AnIterNumCostlyDseIt) / AnIterNumIter);

  // Look for any Devex data to summarise
  if (num_devex_framework) {
    printf("\nDevex summary\n");
    printf("%12" HIGHSINT_FORMAT " Devex frameworks\n", num_devex_framework);
    printf("%12" HIGHSINT_FORMAT " average number of iterations\n",
           AnIterNumEdWtIt[(HighsInt)DualEdgeWeightMode::kDevex] /
               num_devex_framework);
  }
  // Look for any PAMI data to summarise
  if (sum_multi_chosen > 0) {
    const HighsInt pct_minor_iterations_performed =
        (100 * sum_multi_finished) / sum_multi_chosen;
    printf("\nPAMI summary: for average of %0.1g threads \n",
           average_num_threads);
    printf("%12" HIGHSINT_FORMAT " Major iterations\n", multi_iteration_count);
    printf("%12" HIGHSINT_FORMAT " Minor iterations\n", sum_multi_finished);
    printf("%12" HIGHSINT_FORMAT
           " Total rows chosen: performed %3" HIGHSINT_FORMAT
           "%% of possible minor "
           "iterations\n\n",
           sum_multi_chosen, pct_minor_iterations_performed);
  }

  highsLogDev(log_options, HighsLogType::kInfo,
              "\nCost perturbation summary\n");
  logValueDistribution(log_options, cost_perturbation1_distribution);
  logValueDistribution(log_options, cost_perturbation2_distribution);

  logValueDistribution(log_options, before_ftran_upper_sparse_density, numRow);
  logValueDistribution(log_options, ftran_upper_sparse_density, numRow);
  logValueDistribution(log_options, before_ftran_upper_hyper_density, numRow);
  logValueDistribution(log_options, ftran_upper_hyper_density, numRow);
  logValueDistribution(log_options, primal_step_distribution);
  logValueDistribution(log_options, dual_step_distribution);
  logValueDistribution(log_options, simplex_pivot_distribution);
  logValueDistribution(log_options, factor_pivot_threshold_distribution);
  logValueDistribution(log_options, numerical_trouble_distribution);
  logValueDistribution(log_options, cleanup_dual_change_distribution);
  logValueDistribution(log_options, cleanup_primal_step_distribution);
  logValueDistribution(log_options, cleanup_dual_step_distribution);
  logValueDistribution(log_options, cleanup_primal_change_distribution);

  if (AnIterTraceIterDl >= 100) {
    // Possibly (usually) add a temporary record for the final
    // iterations: may end up with one more than
    // kAnIterTraceMaxNumRec records, so ensure that there is
    // enough space in the arrays
    //
    const bool add_extra_record =
        simplex_iteration_count >
        AnIterTrace[AnIterTraceNumRec].AnIterTraceIter;
    if (add_extra_record) {
      AnIterTraceNumRec++;
      AnIterTraceRec& lcAnIter = AnIterTrace[AnIterTraceNumRec];
      lcAnIter.AnIterTraceIter = simplex_iteration_count;
      lcAnIter.AnIterTraceTime = timer_->getWallTime();
      if (average_fraction_of_possible_minor_iterations_performed > 0) {
        lcAnIter.AnIterTraceMulti =
            average_fraction_of_possible_minor_iterations_performed;
      } else {
        lcAnIter.AnIterTraceMulti = 0;
      }
      lcAnIter.AnIterTraceDensity[ANALYSIS_OPERATION_TYPE_FTRAN] =
          col_aq_density;
      lcAnIter.AnIterTraceDensity[ANALYSIS_OPERATION_TYPE_BTRAN_EP] =
          row_ep_density;
      lcAnIter.AnIterTraceDensity[ANALYSIS_OPERATION_TYPE_PRICE_AP] =
          row_ap_density;
      lcAnIter.AnIterTraceDensity[ANALYSIS_OPERATION_TYPE_FTRAN_BFRT] =
          col_aq_density;
      if (edge_weight_mode == DualEdgeWeightMode::kSteepestEdge) {
        lcAnIter.AnIterTraceDensity[ANALYSIS_OPERATION_TYPE_FTRAN_DSE] =
            row_DSE_density;
        lcAnIter.AnIterTraceCostlyDse = AnIterCostlyDseMeasure;
      } else {
        lcAnIter.AnIterTraceDensity[ANALYSIS_OPERATION_TYPE_FTRAN_DSE] = 0;
        lcAnIter.AnIterTraceCostlyDse = 0;
      }
      lcAnIter.AnIterTrace_dual_edge_weight_mode = (HighsInt)edge_weight_mode;
    }
    // Determine whether the Multi and steepest edge columns should be reported
    double su_multi_values = 0;
    double su_dse_values = 0;
    for (HighsInt rec = 1; rec <= AnIterTraceNumRec; rec++) {
      AnIterTraceRec& lcAnIter = AnIterTrace[rec];
      su_multi_values += fabs(lcAnIter.AnIterTraceMulti);
      su_dse_values +=
          fabs(lcAnIter.AnIterTraceDensity[ANALYSIS_OPERATION_TYPE_FTRAN_DSE]);
    }
    const bool report_multi = su_multi_values > 0;
    const bool rp_dual_steepest_edge = su_dse_values > 0;
    printf("\n Iteration speed analysis\n");
    AnIterTraceRec& lcAnIter = AnIterTrace[0];
    HighsInt fmIter = lcAnIter.AnIterTraceIter;
    double fmTime = lcAnIter.AnIterTraceTime;
    printf("        Iter (      FmIter:      ToIter)      Time      Iter/sec ");
    if (report_multi) printf("| PAMI ");
    printf("| C_Aq R_Ep R_Ap ");
    if (rp_dual_steepest_edge) printf(" DSE ");
    printf("| EdWt ");
    if (rp_dual_steepest_edge) {
      printf("| CostlyDse\n");
    } else {
      printf("\n");
    }

    for (HighsInt rec = 1; rec <= AnIterTraceNumRec; rec++) {
      AnIterTraceRec& lcAnIter = AnIterTrace[rec];
      HighsInt toIter = lcAnIter.AnIterTraceIter;
      double toTime = lcAnIter.AnIterTraceTime;
      HighsInt dlIter = toIter - fmIter;
      if (rec < AnIterTraceNumRec && dlIter != AnIterTraceIterDl)
        printf("STRANGE: %" HIGHSINT_FORMAT
               " = dlIter != AnIterTraceIterDl = %" HIGHSINT_FORMAT "\n",
               dlIter, AnIterTraceIterDl);
      double dlTime = toTime - fmTime;
      HighsInt iterSpeed = 0;
      if (dlTime > 0) iterSpeed = dlIter / dlTime;
      HighsInt lc_dual_edge_weight_mode =
          lcAnIter.AnIterTrace_dual_edge_weight_mode;
      std::string str_dual_edge_weight_mode;
      if (lc_dual_edge_weight_mode ==
          (HighsInt)DualEdgeWeightMode::kSteepestEdge)
        str_dual_edge_weight_mode = "DSE";
      else if (lc_dual_edge_weight_mode == (HighsInt)DualEdgeWeightMode::kDevex)
        str_dual_edge_weight_mode = "Dvx";
      else if (lc_dual_edge_weight_mode ==
               (HighsInt)DualEdgeWeightMode::kDantzig)
        str_dual_edge_weight_mode = "Dan";
      else
        str_dual_edge_weight_mode = "XXX";
      printf("%12" HIGHSINT_FORMAT " (%12" HIGHSINT_FORMAT
             ":%12" HIGHSINT_FORMAT ") %9.4f  %12" HIGHSINT_FORMAT " ",
             dlIter, fmIter, toIter, dlTime, iterSpeed);
      if (report_multi) {
        const HighsInt pct = (100 * lcAnIter.AnIterTraceMulti);
        printf("|  %3" HIGHSINT_FORMAT " ", pct);
      }
      printf("|");
      printOneDensity(
          lcAnIter.AnIterTraceDensity[ANALYSIS_OPERATION_TYPE_FTRAN]);
      printOneDensity(
          lcAnIter.AnIterTraceDensity[ANALYSIS_OPERATION_TYPE_BTRAN_EP]);
      printOneDensity(
          lcAnIter.AnIterTraceDensity[ANALYSIS_OPERATION_TYPE_PRICE_AP]);
      double use_row_DSE_density;
      if (rp_dual_steepest_edge) {
        if (lc_dual_edge_weight_mode ==
            (HighsInt)DualEdgeWeightMode::kSteepestEdge) {
          use_row_DSE_density =
              lcAnIter.AnIterTraceDensity[ANALYSIS_OPERATION_TYPE_FTRAN_DSE];
        } else {
          use_row_DSE_density = 0;
        }
        printOneDensity(use_row_DSE_density);
      }
      printf(" |  %3s ", str_dual_edge_weight_mode.c_str());
      if (rp_dual_steepest_edge) {
        double use_costly_dse;
        printf("|     ");
        if (lc_dual_edge_weight_mode ==
            (HighsInt)DualEdgeWeightMode::kSteepestEdge) {
          use_costly_dse = lcAnIter.AnIterTraceCostlyDse;
        } else {
          use_costly_dse = 0;
        }
        printOneDensity(use_costly_dse);
      }
      printf("\n");
      fmIter = toIter;
      fmTime = toTime;
    }
    printf("\n");
    // Remove any temporary record added for the final iterations
    if (add_extra_record) AnIterTraceNumRec--;
  }
}

void HighsSimplexAnalysis::summaryReportFactor() {
  for (HighsInt tran_stage_type = 0; tran_stage_type < NUM_TRAN_STAGE_TYPE;
       tran_stage_type++) {
    TranStageAnalysis& stage = tran_stage[tran_stage_type];
    //    printScatterData(stage.name_, stage.rhs_density_);
    printScatterDataRegressionComparison(stage.name_, stage.rhs_density_);
    if (!stage.num_decision_) return;
    printf("Of %10" HIGHSINT_FORMAT
           " Sps/Hyper decisions made using regression:\n",
           stage.num_decision_);
    printf("   %10" HIGHSINT_FORMAT " wrong sparseTRAN; %10" HIGHSINT_FORMAT
           " wrong hyperTRAN: using original "
           "logic\n",
           stage.num_wrong_original_sparse_decision_,
           stage.num_wrong_original_hyper_decision_);
    printf("   %10" HIGHSINT_FORMAT " wrong sparseTRAN; %10" HIGHSINT_FORMAT
           " wrong hyperTRAN: using new      "
           "logic\n",
           stage.num_wrong_new_sparse_decision_,
           stage.num_wrong_new_hyper_decision_);
  }
}

void HighsSimplexAnalysis::reportSimplexTimer() {
  assert(analyse_simplex_time);
  SimplexTimer simplex_timer;
  simplex_timer.reportSimplexInnerClock(thread_simplex_clocks[0]);
}

void HighsSimplexAnalysis::reportFactorTimer() {
  assert(analyse_factor_time);
  FactorTimer factor_timer;
  HighsInt omp_max_threads = 1;
#ifdef OPENMP
  omp_max_threads = omp_get_max_threads();
#endif
  for (HighsInt i = 0; i < omp_max_threads; i++) {
    //  for (HighsTimerClock clock : thread_factor_clocks) {
    printf("reportFactorTimer: HFactor clocks for OMP thread %" HIGHSINT_FORMAT
           " / %" HIGHSINT_FORMAT "\n",
           i, omp_max_threads - 1);
    factor_timer.reportFactorClock(thread_factor_clocks[i]);
  }
  if (omp_max_threads > 1) {
    HighsTimer& timer = thread_factor_clocks[0].timer_;
    HighsTimerClock all_factor_clocks(timer);
    vector<HighsInt>& clock = all_factor_clocks.clock_;
    factor_timer.initialiseFactorClocks(all_factor_clocks);
    for (HighsInt i = 0; i < omp_max_threads; i++) {
      vector<HighsInt>& thread_clock = thread_factor_clocks[i].clock_;
      for (HighsInt clock_id = 0; clock_id < FactorNumClock; clock_id++) {
        HighsInt all_factor_iClock = clock[clock_id];
        HighsInt thread_factor_iClock = thread_clock[clock_id];
        timer.clock_num_call[all_factor_iClock] +=
            timer.clock_num_call[thread_factor_iClock];
        timer.clock_time[all_factor_iClock] +=
            timer.clock_time[thread_factor_iClock];
      }
    }
    printf("reportFactorTimer: HFactor clocks for all %" HIGHSINT_FORMAT
           " threads\n",
           omp_max_threads);
    factor_timer.reportFactorClock(all_factor_clocks);
  }
}

void HighsSimplexAnalysis::updateInvertFormData(const HFactor& factor) {
  assert(analyse_factor_data);
  const bool report_kernel = false;
  num_invert++;
  assert(factor.basis_matrix_num_el);
  double invert_fill_factor =
      ((1.0 * factor.invert_num_el) / factor.basis_matrix_num_el);
  if (report_kernel) printf("INVERT fill = %6.2f", invert_fill_factor);
  sum_invert_fill_factor += invert_fill_factor;
  running_average_invert_fill_factor =
      0.95 * running_average_invert_fill_factor + 0.05 * invert_fill_factor;

  double kernel_relative_dim = (1.0 * factor.kernel_dim) / numRow;
  if (report_kernel) printf("; kernel dim = %11.4g", kernel_relative_dim);
  if (factor.kernel_dim) {
    num_kernel++;
    max_kernel_dim = max(kernel_relative_dim, max_kernel_dim);
    sum_kernel_dim += kernel_relative_dim;
    running_average_kernel_dim =
        0.95 * running_average_kernel_dim + 0.05 * kernel_relative_dim;

    HighsInt kernel_invert_num_el =
        factor.invert_num_el -
        (factor.basis_matrix_num_el - factor.kernel_num_el);
    assert(factor.kernel_num_el);
    double kernel_fill_factor =
        (1.0 * kernel_invert_num_el) / factor.kernel_num_el;
    sum_kernel_fill_factor += kernel_fill_factor;
    running_average_kernel_fill_factor =
        0.95 * running_average_kernel_fill_factor + 0.05 * kernel_fill_factor;
    if (report_kernel) printf("; fill = %6.2f", kernel_fill_factor);
    const double kMajorKernelRelativeDimThreshold = 0.1;
    if (kernel_relative_dim > kMajorKernelRelativeDimThreshold) {
      num_major_kernel++;
      sum_major_kernel_fill_factor += kernel_fill_factor;
      running_average_major_kernel_fill_factor =
          0.95 * running_average_major_kernel_fill_factor +
          0.05 * kernel_fill_factor;
    }
  }
  if (report_kernel) printf("\n");
}

void HighsSimplexAnalysis::reportInvertFormData() {
  assert(analyse_factor_data);
  printf("grep_kernel,%s,%s,%" HIGHSINT_FORMAT ",%" HIGHSINT_FORMAT
         ",%" HIGHSINT_FORMAT ",",
         model_name_.c_str(), lp_name_.c_str(), num_invert, num_kernel,
         num_major_kernel);
  if (num_kernel) printf("%g", sum_kernel_dim / num_kernel);
  printf(",%g,%g,", running_average_kernel_dim, max_kernel_dim);
  if (num_invert) printf("Fill-in,%g", sum_invert_fill_factor / num_invert);
  printf(",");
  if (num_kernel) printf("%g", sum_kernel_fill_factor / num_kernel);
  printf(",");
  if (num_major_kernel)
    printf("%g", sum_major_kernel_fill_factor / num_major_kernel);
  printf(",%g,%g,%g\n", running_average_invert_fill_factor,
         running_average_kernel_fill_factor,
         running_average_major_kernel_fill_factor);
}

void HighsSimplexAnalysis::iterationReport(const bool header) {
  analysis_log = std::unique_ptr<std::stringstream>(new std::stringstream());
  if (!header) {
    if (dualAlgorithm()) {
      if (pivotal_row_index < 0) return;
    } else {
      if (entering_variable < 0) return;
    }
  }
  reportAlgorithmPhaseIterationObjective(header);
  if (analyse_simplex_data) {
    reportDensity(header);
    reportIterationData(header);
  }
  highsLogDev(log_options, kIterationReportLogType, "%s\n",
              analysis_log->str().c_str());
  if (!header) num_iteration_report_since_last_header++;
}

void HighsSimplexAnalysis::reportAlgorithmPhaseIterationObjective(
    const bool header) {
  if (header) {
    *analysis_log << "       Iteration        Objective    ";
  } else {
    std::string algorithm_name;
    if (dualAlgorithm()) {
      algorithm_name = "Du";
    } else {
      algorithm_name = "Pr";
    }
    *analysis_log << highsFormatToString(
        "%2sPh%1" HIGHSINT_FORMAT " %10" HIGHSINT_FORMAT " %20.10e",
        algorithm_name.c_str(), solve_phase, simplex_iteration_count,
        objective_value);
  }
}

void HighsSimplexAnalysis::reportInfeasibility(const bool header) {
  if (header) {
    *analysis_log << " Infeasibilities num(sum)";
  } else {
    if (solve_phase == 1) {
      *analysis_log << highsFormatToString(" Ph1: %" HIGHSINT_FORMAT "(%g)",
                                           num_primal_infeasibility,
                                           sum_primal_infeasibility);
    } else {
      *analysis_log << highsFormatToString(" Pr: %" HIGHSINT_FORMAT "(%g)",
                                           num_primal_infeasibility,
                                           sum_primal_infeasibility);
    }
    if (sum_dual_infeasibility > 0) {
      *analysis_log << highsFormatToString("; Du: %" HIGHSINT_FORMAT "(%g)",
                                           num_dual_infeasibility,
                                           sum_dual_infeasibility);
    }
  }
}

void HighsSimplexAnalysis::reportThreads(const bool header) {
  assert(analyse_simplex_data);
  if (header) {
    *analysis_log << highsFormatToString("  Threads");
  } else if (num_threads > 0) {
    *analysis_log << highsFormatToString(
        " %2" HIGHSINT_FORMAT "|%2" HIGHSINT_FORMAT "|%2" HIGHSINT_FORMAT "",
        min_threads, num_threads, max_threads);
  } else {
    *analysis_log << highsFormatToString("   |  |  ");
  }
}

void HighsSimplexAnalysis::reportMulti(const bool header) {
  assert(analyse_simplex_data);
  if (header) {
    *analysis_log << highsFormatToString("  Multi");
  } else if (average_fraction_of_possible_minor_iterations_performed >= 0) {
    *analysis_log << highsFormatToString(
        "   %3" HIGHSINT_FORMAT "%%",
        (HighsInt)(100 *
                   average_fraction_of_possible_minor_iterations_performed));
  } else {
    *analysis_log << highsFormatToString("       ");
  }
}

void HighsSimplexAnalysis::reportOneDensity(const double density) {
  assert(analyse_simplex_data);
  const HighsInt log_10_density = intLog10(density);
  if (log_10_density > -99) {
    *analysis_log << highsFormatToString(" %4" HIGHSINT_FORMAT "",
                                         log_10_density);
  } else {
    *analysis_log << highsFormatToString("     ");
  }
}

void HighsSimplexAnalysis::printOneDensity(const double density) {
  assert(analyse_simplex_data);
  const HighsInt log_10_density = intLog10(density);
  if (log_10_density > -99) {
    printf(" %4" HIGHSINT_FORMAT "", log_10_density);
  } else {
    printf("     ");
  }
}

void HighsSimplexAnalysis::reportDensity(const bool header) {
  assert(analyse_simplex_data);
  const bool rp_dual_steepest_edge =
      edge_weight_mode == DualEdgeWeightMode::kSteepestEdge;
  if (header) {
    *analysis_log << highsFormatToString(" C_Aq R_Ep R_Ap");
    if (rp_dual_steepest_edge) {
      *analysis_log << highsFormatToString("  DSE");
    } else {
      *analysis_log << highsFormatToString("     ");
    }
  } else {
    reportOneDensity(col_aq_density);
    reportOneDensity(row_ep_density);
    reportOneDensity(row_ap_density);
    double use_row_DSE_density;
    if (rp_dual_steepest_edge) {
      use_row_DSE_density = row_DSE_density;
    } else {
      use_row_DSE_density = 0;
    }
    reportOneDensity(use_row_DSE_density);
  }
}

void HighsSimplexAnalysis::reportInvert(const bool header) {
  if (header) {
    *analysis_log << highsFormatToString(" Inv");
  } else {
    *analysis_log << highsFormatToString("  %2" HIGHSINT_FORMAT "",
                                         rebuild_reason);
  }
}
/*
void HighsSimplexAnalysis::reportCondition(const bool header) {
  assert(analyse_simplex_data);
  if (header) {
    *analysis_log << highsFormatToString("       k(B)");
  } else {
    *analysis_log << highsFormatToString(" %10.4g",
                      basis_condition);
  }
}
*/

// Primal:
// * primal_delta - 0
// * dual_step    - theta_dual
// * primal_step  - theta_primal
//
// Dual:
// * primal_delta - delta_primal
// * dual_step    - theta_dual
// * primal_step  - theta_primal
void HighsSimplexAnalysis::reportIterationData(const bool header) {
  if (header) {
    *analysis_log << highsFormatToString(
        "     EnC     LvC     LvR        ThDu        ThPr        "
        "DlPr       NumCk          Aa");
  } else if (pivotal_row_index >= 0) {
    *analysis_log << highsFormatToString(
        " %7" HIGHSINT_FORMAT " %7" HIGHSINT_FORMAT " %7" HIGHSINT_FORMAT
        " %11.4g %11.4g %11.4g %11.4g %11.4g",
        entering_variable, leaving_variable, pivotal_row_index, dual_step,
        primal_step, primal_delta, numerical_trouble, pivot_value_from_column);
  } else {
    *analysis_log << highsFormatToString(
        " %7" HIGHSINT_FORMAT " %7" HIGHSINT_FORMAT " %7" HIGHSINT_FORMAT
        " %11.4g %11.4g                                    ",
        entering_variable, leaving_variable, pivotal_row_index, dual_step,
        primal_step);
  }
}

HighsInt HighsSimplexAnalysis::intLog10(const double v) {
  HighsInt intLog10V = -99;
  if (v > 0) intLog10V = log(v) / log(10.0);
  return intLog10V;
}

bool HighsSimplexAnalysis::dualAlgorithm() {
  return (simplex_strategy == kSimplexStrategyDual ||
          simplex_strategy == kSimplexStrategyDualTasks ||
          simplex_strategy == kSimplexStrategyDualMulti);
}
