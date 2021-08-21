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
/**@file lp_data/HEkkDebug.cpp
 * @brief
 */

#include "simplex/HEkkDebug.h"

#include <cassert>
#include <cmath>
#include <string>

#include "lp_data/HighsDebug.h"
#include "lp_data/HighsModelUtils.h"

using std::abs;
using std::max;

const double ok_feasibility_difference = 1e-3;

const double large_basic_dual = 1e-12;
const double excessive_basic_dual = sqrt(large_basic_dual);

const double large_residual_error = 1e-12;
const double excessive_residual_error = sqrt(large_residual_error);

const double updated_dual_small_relative_error = 1e-12;
const double updated_dual_large_relative_error =
    sqrt(updated_dual_small_relative_error);
const double updated_dual_small_absolute_error = 1e-6;
const double updated_dual_large_absolute_error =
    sqrt(updated_dual_small_absolute_error);

HighsDebugStatus ekkDebugSimplex(const std::string message,
                                 const HEkk& ekk_instance,
                                 const SimplexAlgorithm algorithm,
                                 const HighsInt phase, const bool initialise) {
  if (ekk_instance.options_.highs_debug_level < kHighsDebugLevelCheap)
    return HighsDebugStatus::kNotChecked;
  static double max_max_basic_dual;
  static double max_max_primal_residual;
  static double max_max_dual_residual;
  if (initialise) {
    max_max_basic_dual = 0;
    max_max_primal_residual = 0;
    max_max_dual_residual = 0;
    return ::HighsDebugStatus::kOk;
  }
  HighsDebugStatus return_status = HighsDebugStatus::kOk;
  const HighsLp& lp = ekk_instance.lp_;
  const HighsSimplexInfo& info = ekk_instance.info_;
  const SimplexBasis& basis = ekk_instance.basis_;
  const HighsOptions& options = ekk_instance.options_;

  const HighsInt num_col = lp.numCol_;
  const HighsInt num_row = lp.numRow_;
  const HighsInt num_tot = num_col + num_row;
  const HighsInt iteration_count = ekk_instance.iteration_count_;
  std::string value_adjective;
  HighsLogType report_level;

  // Check the nonbasic flags are all kNonbasicFlagTrue or kNonbasicFlagFalse
  for (HighsInt iVar = 0; iVar < num_tot; iVar++) {
    HighsInt flag = basis.nonbasicFlag_[iVar];
    bool flag_error = flag != kNonbasicFlagTrue && flag != kNonbasicFlagFalse;
    if (flag_error) {
      highsLogUser(options.log_options, HighsLogType::kError,
                   "ekkDebugSimplex - %s: Iteration %" HIGHSINT_FORMAT
                   " Variable %" HIGHSINT_FORMAT
                   " has "
                   "nonbasic flag = %" HIGHSINT_FORMAT "\n",
                   message.c_str(), iteration_count, iVar, flag);
      assert(!flag_error);
      return HighsDebugStatus::kLogicalError;
    }
  }
  const double primal_feasibility_tolerance =
      options.primal_feasibility_tolerance;
  const double dual_feasibility_tolerance = options.dual_feasibility_tolerance;
  HighsInt num_dual_infeasibility = 0;
  double max_dual_infeasibility = 0;
  double sum_dual_infeasibility = 0;
  HighsInt num_primal_infeasibility = 0;
  double max_primal_infeasibility = 0;
  double sum_primal_infeasibility = 0;
  // Check the nonbasic variables
  for (HighsInt iVar = 0; iVar < num_tot; iVar++) {
    if (basis.nonbasicFlag_[iVar] == kNonbasicFlagFalse) continue;
    // For nonbasic variables, check that they are on a bound (or free
    // at 0 with correct nonbasic move. Determine dual infeasibilities
    double dual = info.workDual_[iVar];
    double lower = info.workLower_[iVar];
    double upper = info.workUpper_[iVar];
    double value = info.workValue_[iVar];
    double primal_error = 0;
    double dual_infeasibility = 0;
    HighsInt move;
    if (lower == upper) {
      primal_error = abs(lower - value);
      move = kNonbasicMoveZe;
    } else if (value == lower) {
      move = kNonbasicMoveUp;
      dual_infeasibility = max(-dual, 0.);
    } else if (value == upper) {
      move = kNonbasicMoveDn;
      dual_infeasibility = max(dual, 0.);
    } else {
      // If not fixed or at a bound, only valid status can be zero and free
      primal_error = abs(value);
      move = kNonbasicMoveZe;
      dual_infeasibility = abs(dual);
    }
    if (dual_infeasibility > 0) {
      if (dual_infeasibility > dual_feasibility_tolerance) {
        num_dual_infeasibility++;
      }
      max_dual_infeasibility = max(dual_infeasibility, max_dual_infeasibility);
      sum_dual_infeasibility += dual_infeasibility;
    }
    if (primal_error) {
      highsLogUser(options.log_options, HighsLogType::kError,
                   "ekkDebugSimplex - %s: Iteration %" HIGHSINT_FORMAT
                   " Nonbasic variable %" HIGHSINT_FORMAT
                   " "
                   "has primal error "
                   "= %g for [%g, %g, %g]\n",
                   message.c_str(), iteration_count, iVar, primal_error, lower,
                   value, upper);
      assert(!primal_error);
      return HighsDebugStatus::kLogicalError;
    }
    bool move_error = move != basis.nonbasicMove_[iVar];
    if (move_error) {
      highsLogUser(options.log_options, HighsLogType::kError,
                   "ekkDebugSimplex - %s: Iteration %" HIGHSINT_FORMAT
                   " Nonbasic variable %" HIGHSINT_FORMAT
                   " "
                   "has move error "
                   "[%" HIGHSINT_FORMAT " <> %" HIGHSINT_FORMAT
                   "] for [%g, %g, %g]\n",
                   message.c_str(), iteration_count, iVar, move,
                   basis.nonbasicMove_[iVar], lower, value, upper);
      assert(!move_error);
      return HighsDebugStatus::kLogicalError;
    }
  }
  // Check the basic variables
  double max_basic_dual = 0;
  const double base =
      info.primal_simplex_phase1_cost_perturbation_multiplier * 5e-7;
  for (HighsInt iRow = 0; iRow < num_row; iRow++) {
    HighsInt iVar = basis.basicIndex_[iRow];
    // For basic variables, check that the nonbasic flag isn't set,
    // that nonbasicMove is zero, that baseLower/Upper are correct,
    // that the dual is zero and, in primal phase 1, that the cost is
    // correct. Determine primal infeasibilities
    bool nonbasicFlag_error = basis.nonbasicFlag_[iVar] == kNonbasicFlagTrue;
    if (nonbasicFlag_error) {
      highsLogUser(options.log_options, HighsLogType::kError,
                   "ekkDebugSimplex - %s: Iteration %" HIGHSINT_FORMAT
                   " Basic variable %" HIGHSINT_FORMAT
                   " "
                   "has nonbasicFlag = %" HIGHSINT_FORMAT "\n",
                   message.c_str(), iteration_count, iVar,
                   basis.nonbasicFlag_[iVar]);
      assert(!nonbasicFlag_error);
      return HighsDebugStatus::kLogicalError;
    }
    bool nonbasicMove_error = basis.nonbasicMove_[iVar];
    if (nonbasicMove_error) {
      highsLogUser(options.log_options, HighsLogType::kError,
                   "ekkDebugSimplex - %s: Iteration %" HIGHSINT_FORMAT
                   " Basic variable %" HIGHSINT_FORMAT
                   " "
                   "has nonbasicMove = %" HIGHSINT_FORMAT "\n",
                   message.c_str(), iteration_count, iVar,
                   basis.nonbasicMove_[iVar]);
      assert(!nonbasicMove_error);
      return HighsDebugStatus::kLogicalError;
    }
    double workLower = info.workLower_[iVar];
    double workUpper = info.workUpper_[iVar];
    double cost = info.workCost_[iVar];
    double dual = info.workDual_[iVar];
    double lower = info.baseLower_[iRow];
    double upper = info.baseUpper_[iRow];
    double value = info.baseValue_[iRow];
    bool baseBound_error = workLower != lower || workUpper != upper;
    if (baseBound_error) {
      highsLogUser(options.log_options, HighsLogType::kError,
                   "ekkDebugSimplex - %s: Iteration %" HIGHSINT_FORMAT
                   " Basic variable %" HIGHSINT_FORMAT
                   " "
                   "(in row %" HIGHSINT_FORMAT
                   ") has "
                   "baseBound [%g, %g] and workBound [%g, %g]\n",
                   message.c_str(), iteration_count, iVar, iRow, lower, upper,
                   workLower, workUpper);
      assert(!baseBound_error);
      return HighsDebugStatus::kLogicalError;
    }
    max_basic_dual = max(abs(dual), max_basic_dual);
    HighsInt bound_violated = 0;
    if (value < lower - primal_feasibility_tolerance) {
      bound_violated = -1;
    } else if (value > upper + primal_feasibility_tolerance) {
      bound_violated = 1;
    }
    if (algorithm == SimplexAlgorithm::kPrimal && phase == 1) {
      double primal_phase1_cost = bound_violated;
      if (base) primal_phase1_cost *= 1 + base * info.numTotRandomValue_[iRow];
      bool primal_phase1_cost_error = abs(cost - primal_phase1_cost);
      if (primal_phase1_cost_error) {
        highsLogUser(options.log_options, HighsLogType::kError,
                     "ekkDebugSimplex - %s: Iteration %" HIGHSINT_FORMAT
                     " Basic variable %" HIGHSINT_FORMAT
                     " "
                     "(in row %" HIGHSINT_FORMAT
                     ") has "
                     "primal phase 1 cost %g for [%g, %g, %g]\n",
                     message.c_str(), iteration_count, iVar, iRow, cost, lower,
                     value, upper);
        assert(!primal_phase1_cost_error);
        return HighsDebugStatus::kLogicalError;
      }
    }
    if (!bound_violated) continue;
    // @primal_infeasibility calculation
    double primal_infeasibility = 0;
    if (bound_violated < 0) {
      primal_infeasibility = lower - value;
    } else {
      primal_infeasibility = value - upper;
    }
    if (primal_infeasibility > primal_feasibility_tolerance)
      num_primal_infeasibility++;
    max_primal_infeasibility =
        max(primal_infeasibility, max_primal_infeasibility);
    sum_primal_infeasibility += primal_infeasibility;
  }
  // Report on basic dual values
  if (max_basic_dual > excessive_basic_dual) {
    value_adjective = "Excessive";
    report_level = HighsLogType::kInfo;
    return_status = debugWorseStatus(HighsDebugStatus::kError, return_status);
  } else if (max_basic_dual > large_basic_dual) {
    value_adjective = "Large";
    report_level = HighsLogType::kDetailed;
    return_status = debugWorseStatus(HighsDebugStatus::kWarning, return_status);
  } else {
    value_adjective = "OK";
    report_level = HighsLogType::kVerbose;
    return_status = debugWorseStatus(HighsDebugStatus::kOk, return_status);
  }
  if (max_basic_dual > 2 * max_max_basic_dual) {
    highsLogDev(options.log_options, report_level,
                "ekkDebugSimplex - %s: Iteration %" HIGHSINT_FORMAT
                " %-9s max   basic dual = %9.4g\n",
                message.c_str(), iteration_count, value_adjective.c_str(),
                max_basic_dual);
    max_max_basic_dual = max_basic_dual;
  }
  // Check that the number, max and sums of primal and dual infeasibilities (if
  // known) are correct
  const HighsInt info_num_primal_infeasibility =
      ekk_instance.info_.num_primal_infeasibility;
  if (info_num_primal_infeasibility >= 0) {
    const bool illegal_num_primal_infeasibility =
        num_primal_infeasibility != info_num_primal_infeasibility;
    if (illegal_num_primal_infeasibility) {
      highsLogUser(options.log_options, HighsLogType::kError,
                   "ekkDebugSimplex - %s: Iteration %" HIGHSINT_FORMAT
                   " Should have %" HIGHSINT_FORMAT
                   " not "
                   "%" HIGHSINT_FORMAT " primal infeasibilities\n",
                   message.c_str(), iteration_count, num_primal_infeasibility,
                   info_num_primal_infeasibility);
      assert(!illegal_num_primal_infeasibility);
      return HighsDebugStatus::kLogicalError;
    }
  }
  const double info_max_primal_infeasibility =
      ekk_instance.info_.max_primal_infeasibility;
  if (info_max_primal_infeasibility >= 0) {
    const bool illegal_max_primal_infeasibility =
        abs(max_primal_infeasibility - info_max_primal_infeasibility) >
        ok_feasibility_difference;
    if (illegal_max_primal_infeasibility) {
      highsLogUser(options.log_options, HighsLogType::kError,
                   "ekkDebugSimplex - %s: Iteration %" HIGHSINT_FORMAT
                   " Should have %g not "
                   "%g max primal infeasibility\n",
                   message.c_str(), iteration_count, max_primal_infeasibility,
                   info_max_primal_infeasibility);
      assert(!illegal_max_primal_infeasibility);
      return HighsDebugStatus::kLogicalError;
    }
  }
  const double info_sum_primal_infeasibility =
      ekk_instance.info_.sum_primal_infeasibility;
  if (info_sum_primal_infeasibility >= 0) {
    const bool illegal_sum_primal_infeasibility =
        abs(sum_primal_infeasibility - info_sum_primal_infeasibility) >
        ok_feasibility_difference;
    if (illegal_sum_primal_infeasibility) {
      highsLogUser(options.log_options, HighsLogType::kError,
                   "ekkDebugSimplex - %s: Iteration %" HIGHSINT_FORMAT
                   " Should have %g not "
                   "%g sum primal infeasibilities\n",
                   message.c_str(), iteration_count, sum_primal_infeasibility,
                   info_sum_primal_infeasibility);
      assert(!illegal_sum_primal_infeasibility);
      return HighsDebugStatus::kLogicalError;
    }
  }
  const HighsInt info_num_dual_infeasibility =
      ekk_instance.info_.num_dual_infeasibility;
  if (info_num_dual_infeasibility >= 0) {
    const bool illegal_num_dual_infeasibility =
        num_dual_infeasibility != info_num_dual_infeasibility;
    if (illegal_num_dual_infeasibility) {
      highsLogUser(options.log_options, HighsLogType::kError,
                   "ekkDebugSimplex - %s: Iteration %" HIGHSINT_FORMAT
                   " Should have %" HIGHSINT_FORMAT
                   " not "
                   "%" HIGHSINT_FORMAT " dual infeasibilities\n",
                   message.c_str(), iteration_count, num_dual_infeasibility,
                   info_num_dual_infeasibility);
      assert(!illegal_num_dual_infeasibility);
      return HighsDebugStatus::kLogicalError;
    }
  }
  const double info_max_dual_infeasibility =
      ekk_instance.info_.max_dual_infeasibility;
  if (info_max_dual_infeasibility >= 0) {
    const bool illegal_max_dual_infeasibility =
        abs(max_dual_infeasibility - info_max_dual_infeasibility) >
        ok_feasibility_difference;
    if (illegal_max_dual_infeasibility) {
      highsLogUser(options.log_options, HighsLogType::kError,
                   "ekkDebugSimplex - %s: Iteration %" HIGHSINT_FORMAT
                   " Should have %g not "
                   "%g max dual infeasibility\n",
                   message.c_str(), iteration_count, max_dual_infeasibility,
                   info_max_dual_infeasibility);
      assert(!illegal_max_dual_infeasibility);
      return HighsDebugStatus::kLogicalError;
    }
  }
  const double info_sum_dual_infeasibility =
      ekk_instance.info_.sum_dual_infeasibility;
  if (info_sum_dual_infeasibility >= 0) {
    const bool illegal_sum_dual_infeasibility =
        abs(sum_dual_infeasibility - info_sum_dual_infeasibility) >
        ok_feasibility_difference;
    if (illegal_sum_dual_infeasibility) {
      highsLogUser(options.log_options, HighsLogType::kError,
                   "ekkDebugSimplex - %s: Iteration %" HIGHSINT_FORMAT
                   " Should have %g not "
                   "%g sum dual infeasibilities\n",
                   message.c_str(), iteration_count, sum_dual_infeasibility,
                   info_sum_dual_infeasibility);
      assert(!illegal_sum_dual_infeasibility);
      return HighsDebugStatus::kLogicalError;
    }
  }
  // Check any assumed feasibility
  bool require_primal_feasible_in_primal_simplex =
      algorithm == SimplexAlgorithm::kPrimal && (phase == 0 || phase == 2);
  bool require_primal_feasible_in_dual_simplex =
      algorithm == SimplexAlgorithm::kDual && phase == 0;
  bool require_primal_feasible = require_primal_feasible_in_primal_simplex ||
                                 require_primal_feasible_in_dual_simplex;
  bool illegal_primal_infeasibility =
      require_primal_feasible && num_primal_infeasibility > 0;
  if (illegal_primal_infeasibility) {
    // Should be primal feasible, but isn't
    highsLogUser(options.log_options, HighsLogType::kError,
                 "ekkDebugSimplex - %s: Iteration %" HIGHSINT_FORMAT
                 " Should be primal "
                 "feasible, but num / max "
                 "/ sum primal infeasibility is %" HIGHSINT_FORMAT
                 " / %g / %g\n",
                 message.c_str(), iteration_count, num_primal_infeasibility,
                 max_primal_infeasibility, sum_primal_infeasibility);
    assert(!illegal_primal_infeasibility);
    return HighsDebugStatus::kLogicalError;
  }
  bool require_dual_feasible_in_dual_simplex =
      algorithm == SimplexAlgorithm::kDual &&
      ekk_instance.status_.has_fresh_rebuild &&
      ekk_instance.info_.allow_cost_perturbation;

  bool illegal_dual_infeasibility =
      (require_dual_feasible_in_dual_simplex || phase == 0) &&
      num_dual_infeasibility > 0;
  if (illegal_dual_infeasibility) {
    // Dual simplex or optimal but has dual infeasibilities
    highsLogUser(options.log_options, HighsLogType::kError,
                 "ekkDebugSimplex - %s: Iteration %" HIGHSINT_FORMAT
                 " Should be dual "
                 "feasible, but num / max / "
                 "sum dual infeasibility is %" HIGHSINT_FORMAT
                 " / %g / %g; Phase = %" HIGHSINT_FORMAT "; status = %s\n",
                 message.c_str(), iteration_count, num_dual_infeasibility,
                 max_dual_infeasibility, sum_dual_infeasibility, phase,
                 utilModelStatusToString(ekk_instance.model_status_).c_str());
    assert(!illegal_dual_infeasibility);
    return HighsDebugStatus::kLogicalError;
  }

  if (ekk_instance.options_.highs_debug_level < kHighsDebugLevelCostly)
    return return_status;
  // Now determine the primal and dual residuals.
  //
  // This uses the primal values for the columns to determine row
  // activities that are checked against the primal values for the
  // rows. It uses the pi vector to determine column duals. The
  // entries of the pi vector are the negated duals for nonbasic rows,
  // and costs for basic rows. The latter are normally zero, but will
  // be nonzero if the constraint is violated in primal phase 1, or if
  // the row cost is a perturbed zero in dual simplex.
  vector<double> primal_value(num_tot);
  vector<double> dual_value(num_tot);
  for (HighsInt iVar = 0; iVar < num_tot; iVar++) {
    primal_value[iVar] = info.workValue_[iVar];
    dual_value[iVar] = info.workDual_[iVar];
  }
  for (HighsInt iRow = 0; iRow < num_row; iRow++) {
    HighsInt iVar = basis.basicIndex_[iRow];
    primal_value[iVar] = info.baseValue_[iRow];
    dual_value[iVar] = -info.workCost_[iVar];
  }
  // Accumulate primal_activities
  double max_dual_residual = 0;
  vector<double> primal_activity(num_row, 0);
  for (HighsInt iCol = 0; iCol < num_col; iCol++) {
    double dual = info.workCost_[iCol];
    double value = primal_value[iCol];
    for (HighsInt iEl = lp.Astart_[iCol]; iEl < lp.Astart_[iCol + 1]; iEl++) {
      HighsInt iRow = lp.Aindex_[iEl];
      HighsInt iVar = num_col + iRow;
      double Avalue = lp.Avalue_[iEl];
      primal_activity[iRow] += value * Avalue;
      dual += dual_value[iVar] * Avalue;
    }
    double dual_residual = abs(dual - info.workDual_[iCol]);
    max_dual_residual = max(dual_residual, max_dual_residual);
  }
  // Remember that simplex row values are the negated row activities
  double max_primal_residual = 0;
  for (HighsInt iRow = 0; iRow < num_row; iRow++) {
    HighsInt iVar = num_col + iRow;
    double primal_residual = abs(primal_activity[iRow] + primal_value[iVar]);
    max_primal_residual = max(primal_residual, max_primal_residual);
  }
  if (max_primal_residual > excessive_residual_error) {
    value_adjective = "Excessive";
    report_level = HighsLogType::kInfo;
    return_status = debugWorseStatus(HighsDebugStatus::kError, return_status);
  } else if (max_primal_residual > large_residual_error) {
    value_adjective = "Large";
    report_level = HighsLogType::kDetailed;
    return_status = debugWorseStatus(HighsDebugStatus::kWarning, return_status);
  } else {
    value_adjective = "OK";
    report_level = HighsLogType::kVerbose;
    return_status = debugWorseStatus(HighsDebugStatus::kOk, return_status);
  }
  if (max_primal_residual > 2 * max_max_primal_residual) {
    highsLogDev(options.log_options, report_level,
                "ekkDebugSimplex - %s: Iteration %" HIGHSINT_FORMAT
                " %-9s max primal residual = %9.4g\n",
                message.c_str(), iteration_count, value_adjective.c_str(),
                max_primal_residual);
    max_max_primal_residual = max_primal_residual;
  }
  if (max_dual_residual > excessive_residual_error) {
    value_adjective = "Excessive";
    report_level = HighsLogType::kInfo;
    return_status = debugWorseStatus(HighsDebugStatus::kError, return_status);
  } else if (max_dual_residual > large_residual_error) {
    value_adjective = "Large";
    report_level = HighsLogType::kDetailed;
    return_status = debugWorseStatus(HighsDebugStatus::kWarning, return_status);
  } else {
    value_adjective = "OK";
    report_level = HighsLogType::kVerbose;
    return_status = debugWorseStatus(HighsDebugStatus::kOk, return_status);
  }
  if (max_dual_residual > 2 * max_max_dual_residual) {
    highsLogDev(options.log_options, report_level,
                "ekkDebugSimplex - %s: Iteration %" HIGHSINT_FORMAT
                " %-9s max   dual residual = %9.4g\n",
                message.c_str(), iteration_count, value_adjective.c_str(),
                max_dual_residual);
    max_max_dual_residual = max_dual_residual;
  }
  return return_status;
}

HighsDebugStatus ekkDebugBasisCorrect(const HEkk& ekk_instance) {
  // Nontrivially expensive analysis of a Simplex basis, checking
  // consistency, and then correctness of nonbasicMove
  const HighsOptions& options = ekk_instance.options_;
  if (options.highs_debug_level < kHighsDebugLevelCheap)
    return HighsDebugStatus::kNotChecked;
  HighsDebugStatus return_status = HighsDebugStatus::kOk;
  const bool consistent =
      ekkDebugBasisConsistent(ekk_instance) != HighsDebugStatus::kLogicalError;
  if (!consistent) {
    highsLogUser(options.log_options, HighsLogType::kError,
                 "Supposed to be a Simplex basis, but not consistent\n");
    assert(consistent);
    return_status = HighsDebugStatus::kLogicalError;
  }
  if (options.highs_debug_level < kHighsDebugLevelCostly) return return_status;
  const bool correct_nonbasicMove =
      ekkDebugNonbasicMove(ekk_instance) != HighsDebugStatus::kLogicalError;
  if (!correct_nonbasicMove) {
    highsLogUser(
        options.log_options, HighsLogType::kError,
        "Supposed to be a Simplex basis, but nonbasicMove is incorrect\n");
    assert(correct_nonbasicMove);
    return_status = HighsDebugStatus::kLogicalError;
  }
  return return_status;
}

HighsDebugStatus ekkDebugNonbasicMove(const HEkk& ekk_instance) {
  // Non-trivially expensive check of NonbasicMove
  if (ekk_instance.options_.highs_debug_level < kHighsDebugLevelCostly)
    return HighsDebugStatus::kNotChecked;
  HighsDebugStatus return_status = HighsDebugStatus::kOk;
  const HighsOptions& options = ekk_instance.options_;
  const HighsLp& lp = ekk_instance.lp_;
  const SimplexBasis& basis = ekk_instance.basis_;
  HighsInt num_free_variable_move_errors = 0;
  HighsInt num_lower_bounded_variable_move_errors = 0;
  HighsInt num_upper_bounded_variable_move_errors = 0;
  HighsInt num_boxed_variable_move_errors = 0;
  HighsInt num_fixed_variable_move_errors = 0;
  const HighsInt numTot = lp.numCol_ + lp.numRow_;
  bool right_size = (HighsInt)basis.nonbasicMove_.size() == numTot;
  // Check consistency of nonbasicMove
  if (!right_size) {
    highsLogUser(options.log_options, HighsLogType::kError,
                 "nonbasicMove size error\n");
    assert(right_size);
    return_status = HighsDebugStatus::kLogicalError;
  }
  double lower;
  double upper;

  for (HighsInt iVar = 0; iVar < numTot; iVar++) {
    if (!basis.nonbasicFlag_[iVar]) continue;
    // Nonbasic variable
    if (iVar < lp.numCol_) {
      lower = lp.colLower_[iVar];
      upper = lp.colUpper_[iVar];
    } else {
      HighsInt iRow = iVar - lp.numCol_;
      lower = -lp.rowUpper_[iRow];
      upper = -lp.rowLower_[iRow];
    }

    if (highs_isInfinity(upper)) {
      if (highs_isInfinity(-lower)) {
        // Free
        if (basis.nonbasicMove_[iVar]) {
          num_free_variable_move_errors++;
        }
      } else {
        // Only lower bounded
        if (basis.nonbasicMove_[iVar] != kNonbasicMoveUp) {
          num_lower_bounded_variable_move_errors++;
        }
      }
    } else {
      if (highs_isInfinity(-lower)) {
        // Only upper bounded
        if (basis.nonbasicMove_[iVar] != kNonbasicMoveDn) {
          num_upper_bounded_variable_move_errors++;
        }
      } else {
        // Boxed or fixed
        if (lower != upper) {
          // Boxed
          if (!basis.nonbasicMove_[iVar]) {
            num_boxed_variable_move_errors++;
          }
        } else {
          // Fixed
          if (basis.nonbasicMove_[iVar]) {
            num_fixed_variable_move_errors++;
          }
        }
      }
    }
  }
  HighsInt num_errors =
      num_free_variable_move_errors + num_lower_bounded_variable_move_errors +
      num_upper_bounded_variable_move_errors + num_boxed_variable_move_errors +
      num_fixed_variable_move_errors;

  if (num_errors) {
    highsLogUser(
        options.log_options, HighsLogType::kError,
        "There are %" HIGHSINT_FORMAT " nonbasicMove errors: %" HIGHSINT_FORMAT
        " free; %" HIGHSINT_FORMAT " lower; %" HIGHSINT_FORMAT
        " upper; %" HIGHSINT_FORMAT
        " "
        "boxed; %" HIGHSINT_FORMAT " fixed\n",
        num_errors, num_free_variable_move_errors,
        num_lower_bounded_variable_move_errors,
        num_upper_bounded_variable_move_errors, num_boxed_variable_move_errors,
        num_fixed_variable_move_errors);
    assert(num_errors == 0);
    return_status = HighsDebugStatus::kLogicalError;
  }
  return return_status;
}

HighsDebugStatus ekkDebugBasisConsistent(const HEkk& ekk_instance) {
  // Cheap analysis of a Simplex basis, checking vector sizes, numbers
  // of basic/nonbasic variables and non-repetition of basic variables
  if (ekk_instance.options_.highs_debug_level < kHighsDebugLevelCheap)
    return HighsDebugStatus::kNotChecked;
  HighsDebugStatus return_status = HighsDebugStatus::kOk;
  const HighsOptions& options = ekk_instance.options_;
  const HighsLp& lp = ekk_instance.lp_;
  const SimplexBasis& basis = ekk_instance.basis_;
  // Check consistency of nonbasicFlag
  if (ekkDebugNonbasicFlagConsistent(ekk_instance) ==
      HighsDebugStatus::kLogicalError) {
    highsLogUser(options.log_options, HighsLogType::kError,
                 "nonbasicFlag inconsistent\n");
    return_status = HighsDebugStatus::kLogicalError;
  }
  const bool right_size = (HighsInt)basis.basicIndex_.size() == lp.numRow_;
  // Check consistency of basicIndex
  if (!right_size) {
    highsLogUser(options.log_options, HighsLogType::kError,
                 "basicIndex size error\n");
    assert(right_size);
    return_status = HighsDebugStatus::kLogicalError;
  }
  // Use localNonbasicFlag so that duplicate entries in basicIndex can
  // be spotted
  vector<int8_t> localNonbasicFlag = basis.nonbasicFlag_;
  for (HighsInt iRow = 0; iRow < lp.numRow_; iRow++) {
    HighsInt iCol = basis.basicIndex_[iRow];
    HighsInt flag = localNonbasicFlag[iCol];
    // Indicate that this column has been found in basicIndex
    localNonbasicFlag[iCol] = -1;
    if (flag) {
      // Nonzero value for localNonbasicFlag entry means that column is either
      if (flag == kNonbasicFlagTrue) {
        // Nonbasic...
        highsLogUser(options.log_options, HighsLogType::kError,
                     "Entry basicIndex_[%" HIGHSINT_FORMAT
                     "] = %" HIGHSINT_FORMAT " is not basic\n",
                     iRow, iCol);
      } else {
        // .. or is -1 since it has already been found in basicIndex
        highsLogUser(options.log_options, HighsLogType::kError,
                     "Entry basicIndex_[%" HIGHSINT_FORMAT
                     "] = %" HIGHSINT_FORMAT " is already basic\n",
                     iRow, iCol);
        assert(flag == -1);
      }
      assert(!flag);
      return_status = HighsDebugStatus::kLogicalError;
    }
  }
  return return_status;
}

HighsDebugStatus ekkDebugNonbasicFlagConsistent(const HEkk& ekk_instance) {
  if (ekk_instance.options_.highs_debug_level < kHighsDebugLevelCheap)
    return HighsDebugStatus::kNotChecked;
  HighsDebugStatus return_status = HighsDebugStatus::kOk;
  const HighsOptions& options = ekk_instance.options_;
  const HighsLp& lp = ekk_instance.lp_;
  const SimplexBasis& basis = ekk_instance.basis_;
  HighsInt numTot = lp.numCol_ + lp.numRow_;
  const bool right_size = (HighsInt)basis.nonbasicFlag_.size() == numTot;
  if (!right_size) {
    highsLogUser(options.log_options, HighsLogType::kError,
                 "nonbasicFlag size error\n");
    assert(right_size);
    return_status = HighsDebugStatus::kLogicalError;
  }
  HighsInt num_basic_variables = 0;
  for (HighsInt var = 0; var < numTot; var++) {
    if (basis.nonbasicFlag_[var] == kNonbasicFlagFalse) {
      num_basic_variables++;
    } else {
      assert(basis.nonbasicFlag_[var] == kNonbasicFlagTrue);
    }
  }
  bool right_num_basic_variables = num_basic_variables == lp.numRow_;
  if (!right_num_basic_variables) {
    highsLogUser(options.log_options, HighsLogType::kError,
                 "nonbasicFlag has %" HIGHSINT_FORMAT ", not %" HIGHSINT_FORMAT
                 " basic variables\n",
                 num_basic_variables, lp.numRow_);
    assert(right_num_basic_variables);
    return_status = HighsDebugStatus::kLogicalError;
  }
  return return_status;
}

HighsDebugStatus ekkDebugOkForSolve(const HEkk& ekk_instance,
                                    const SimplexAlgorithm algorithm,
                                    const HighsInt phase,
                                    const HighsModelStatus model_status) {
  if (ekk_instance.options_.highs_debug_level < kHighsDebugLevelCheap)
    return HighsDebugStatus::kNotChecked;
  const HighsDebugStatus return_status = HighsDebugStatus::kOk;
  const HighsLp& lp = ekk_instance.lp_;
  const HighsSimplexStatus& status = ekk_instance.status_;
  const SimplexBasis& basis = ekk_instance.basis_;
  const HighsOptions& options = ekk_instance.options_;
  bool ok;
  // Minimal check - just look at flags. This means we trust them!
  ok = status.has_basis && status.has_matrix && status.has_factor_arrays &&
       //       status.has_dual_steepest_edge_weights &&
       status.has_invert;
  if (!ok) {
    if (!status.has_basis)
      highsLogUser(options.log_options, HighsLogType::kError,
                   "Not OK to solve since status.has_basis = "
                   "%" HIGHSINT_FORMAT "\n",
                   status.has_basis);
    if (!status.has_matrix)
      highsLogUser(options.log_options, HighsLogType::kError,
                   "Not OK to solve since status.has_matrix = "
                   "%" HIGHSINT_FORMAT "\n",
                   status.has_matrix);
    if (!status.has_factor_arrays)
      highsLogUser(options.log_options, HighsLogType::kError,
                   "Not OK to solve since status.has_factor_arrays "
                   "= %" HIGHSINT_FORMAT "\n",
                   status.has_factor_arrays);
    if (!status.has_dual_steepest_edge_weights)
      highsLogUser(options.log_options, HighsLogType::kError,
                   "Not OK to solve since "
                   "status.has_dual_steepest_edge_weights = %" HIGHSINT_FORMAT
                   "\n",
                   status.has_dual_steepest_edge_weights);
    if (!status.has_invert)
      highsLogUser(options.log_options, HighsLogType::kError,
                   "Not OK to solve since status.has_invert = "
                   "%" HIGHSINT_FORMAT "\n",
                   status.has_invert);
  }
  if (ekk_instance.options_.highs_debug_level < kHighsDebugLevelCostly)
    return return_status;
  // Basis and data check
  if (ekkDebugBasisConsistent(ekk_instance) == HighsDebugStatus::kLogicalError)
    return HighsDebugStatus::kLogicalError;
  // Check work cost, lower, upper and range
  if (!ekkDebugWorkArraysOk(ekk_instance, algorithm, phase, model_status))
    return HighsDebugStatus::kLogicalError;
  const HighsInt numTot = lp.numCol_ + lp.numRow_;
  // Check nonbasic move against work cost, lower, upper and range
  for (HighsInt var = 0; var < numTot; ++var) {
    if (basis.nonbasicFlag_[var]) {
      // Nonbasic variable
      if (!ekkDebugOneNonbasicMoveVsWorkArraysOk(ekk_instance, var))
        return HighsDebugStatus::kLogicalError;
    }
  }
  return return_status;
}

// Methods below are not called externally

bool ekkDebugWorkArraysOk(const HEkk& ekk_instance,
                          const SimplexAlgorithm algorithm,
                          const HighsInt phase,
                          const HighsModelStatus model_status) {
  const HighsLp& lp = ekk_instance.lp_;
  const HighsSimplexInfo& info = ekk_instance.info_;
  const HighsOptions& options = ekk_instance.options_;
  bool ok = true;
  // Don't check dual simplex phase 1 bounds or perturbed bounds
  const bool dual_phase1 = algorithm == SimplexAlgorithm::kDual && phase == 1;
  const bool primal_phase1 =
      algorithm == SimplexAlgorithm::kPrimal && phase == 1;
  if (!(dual_phase1 || info.bounds_perturbed)) {
    for (HighsInt col = 0; col < lp.numCol_; ++col) {
      HighsInt var = col;
      if (!highs_isInfinity(-info.workLower_[var])) {
        double lp_lower = info.workLower_[var];
        ok = lp_lower == lp.colLower_[col];
        if (!ok) {
          highsLogUser(options.log_options, HighsLogType::kError,
                       "For col %" HIGHSINT_FORMAT
                       ", info.workLower_ should be %g but is %g\n",
                       col, lp.colLower_[col], lp_lower);
          return ok;
        }
      }
      if (!highs_isInfinity(info.workUpper_[var])) {
        double lp_upper = info.workUpper_[var];
        ok = lp_upper == lp.colUpper_[col];
        if (!ok) {
          highsLogUser(options.log_options, HighsLogType::kError,
                       "For col %" HIGHSINT_FORMAT
                       ", info.workUpper_ should be %g but is %g\n",
                       col, lp.colUpper_[col], lp_upper);
          return ok;
        }
      }
    }
    for (HighsInt row = 0; row < lp.numRow_; ++row) {
      HighsInt var = lp.numCol_ + row;
      if (!highs_isInfinity(-info.workLower_[var])) {
        double lp_lower = info.workLower_[var];
        ok = lp_lower == -lp.rowUpper_[row];
        if (!ok) {
          highsLogUser(options.log_options, HighsLogType::kError,
                       "For row %" HIGHSINT_FORMAT
                       ", info.workLower_ should be %g but is %g\n",
                       row, -lp.rowUpper_[row], lp_lower);
          return ok;
        }
      }
      if (!highs_isInfinity(info.workUpper_[var])) {
        double lp_upper = info.workUpper_[var];
        ok = lp_upper == -lp.rowLower_[row];
        if (!ok) {
          highsLogUser(options.log_options, HighsLogType::kError,
                       "For row %" HIGHSINT_FORMAT
                       ", info.workUpper_ should be %g but is %g\n",
                       row, -lp.rowLower_[row], lp_upper);
          return ok;
        }
      }
    }
    const HighsInt numTot = lp.numCol_ + lp.numRow_;
    for (HighsInt var = 0; var < numTot; ++var) {
      ok =
          info.workRange_[var] == (info.workUpper_[var] - info.workLower_[var]);
      if (!ok) {
        highsLogUser(options.log_options, HighsLogType::kError,
                     "For variable %" HIGHSINT_FORMAT
                     ", info.workRange_ should be %g = %g - %g "
                     "but is %g\n",
                     var, info.workUpper_[var] - info.workLower_[var],
                     info.workUpper_[var], info.workLower_[var],
                     info.workRange_[var]);
        return ok;
      }
    }
  }
  // Don't check costs against the LP, when using primal simplex in
  // primal phase 1, if the LP is primal infeasible, or if the costs
  // have been perturbed
  if (!(primal_phase1 || model_status == HighsModelStatus::kInfeasible ||
        info.costs_perturbed)) {
    for (HighsInt col = 0; col < lp.numCol_; ++col) {
      HighsInt var = col;
      double work_cost = info.workCost_[var];
      double ok_cost = (HighsInt)lp.sense_ * lp.colCost_[col];
      ok = work_cost == ok_cost;
      if (!ok) {
        highsLogUser(options.log_options, HighsLogType::kError,
                     "For col %" HIGHSINT_FORMAT
                     ", info.workCost_ should be %g but is %g\n",
                     col, ok_cost, info.workCost_[var]);
        return ok;
      }
    }
    for (HighsInt row = 0; row < lp.numRow_; ++row) {
      HighsInt var = lp.numCol_ + row;
      ok = info.workCost_[var] == 0.;
      if (!ok) {
        highsLogUser(options.log_options, HighsLogType::kError,
                     "For row %" HIGHSINT_FORMAT
                     ", info.workCost_ should be zero but is %g\n",
                     row, info.workCost_[var]);
        return ok;
      }
    }
  }
  // ok must be true if we reach here
  assert(ok);
  return ok;
}

bool ekkDebugOneNonbasicMoveVsWorkArraysOk(const HEkk& ekk_instance,
                                           const HighsInt var) {
  const HighsLp& lp = ekk_instance.lp_;
  const HighsSimplexInfo& info = ekk_instance.info_;
  const SimplexBasis& basis = ekk_instance.basis_;
  const HighsOptions& options = ekk_instance.options_;
  assert(var >= 0);
  assert(var < lp.numCol_ + lp.numRow_);
  // Make sure we're not checking a basic variable
  if (!basis.nonbasicFlag_[var]) return true;
  bool ok;
  if (!highs_isInfinity(-info.workLower_[var])) {
    if (!highs_isInfinity(info.workUpper_[var])) {
      // Finite lower and upper bounds so nonbasic move depends on whether they
      // are equal
      if (info.workLower_[var] == info.workUpper_[var]) {
        // Fixed variable
        ok = basis.nonbasicMove_[var] == kNonbasicMoveZe;
        if (!ok) {
          highsLogUser(options.log_options, HighsLogType::kError,
                       "Fixed variable %" HIGHSINT_FORMAT
                       " (lp.numCol_ = %" HIGHSINT_FORMAT
                       ") [%11g, %11g, "
                       "%11g] so nonbasic "
                       "move should be zero but is %" HIGHSINT_FORMAT "\n",
                       var, lp.numCol_, info.workLower_[var],
                       info.workValue_[var], info.workUpper_[var],
                       basis.nonbasicMove_[var]);
          return ok;
        }
        ok = info.workValue_[var] == info.workLower_[var];
        if (!ok) {
          highsLogUser(options.log_options, HighsLogType::kError,
                       "Fixed variable %" HIGHSINT_FORMAT
                       " (lp.numCol_ = %" HIGHSINT_FORMAT
                       ") so "
                       "info.work value should be %g but "
                       "is %g\n",
                       var, lp.numCol_, info.workLower_[var],
                       info.workValue_[var]);
          return ok;
        }
      } else {
        // Boxed variable
        ok = (basis.nonbasicMove_[var] == kNonbasicMoveUp) ||
             (basis.nonbasicMove_[var] == kNonbasicMoveDn);
        if (!ok) {
          highsLogUser(
              options.log_options, HighsLogType::kError,
              "Boxed variable %" HIGHSINT_FORMAT
              " (lp.numCol_ = %" HIGHSINT_FORMAT
              ") [%11g, %11g, "
              "%11g] range %g so "
              "nonbasic move should be up/down but is  %" HIGHSINT_FORMAT "\n",
              var, lp.numCol_, info.workLower_[var], info.workValue_[var],
              info.workUpper_[var], info.workUpper_[var] - info.workLower_[var],
              basis.nonbasicMove_[var]);
          return ok;
        }
        if (basis.nonbasicMove_[var] == kNonbasicMoveUp) {
          ok = info.workValue_[var] == info.workLower_[var];
          if (!ok) {
            highsLogUser(options.log_options, HighsLogType::kError,
                         "Boxed variable %" HIGHSINT_FORMAT
                         " (lp.numCol_ = %" HIGHSINT_FORMAT
                         ") with "
                         "kNonbasicMoveUp so work "
                         "value should be %g but is %g\n",
                         var, lp.numCol_, info.workLower_[var],
                         info.workValue_[var]);
            return ok;
          }
        } else {
          ok = info.workValue_[var] == info.workUpper_[var];
          if (!ok) {
            highsLogUser(options.log_options, HighsLogType::kError,
                         "Boxed variable %" HIGHSINT_FORMAT
                         " (lp.numCol_ = %" HIGHSINT_FORMAT
                         ") with "
                         "kNonbasicMoveDn so work "
                         "value should be %g but is %g\n",
                         var, lp.numCol_, info.workUpper_[var],
                         info.workValue_[var]);
            return ok;
          }
        }
      }
    } else {
      // Infinite upper bound
      ok = basis.nonbasicMove_[var] == kNonbasicMoveUp;
      if (!ok) {
        highsLogUser(options.log_options, HighsLogType::kError,
                     "Finite lower bound and infinite upper bound variable "
                     "%" HIGHSINT_FORMAT
                     " "
                     "(lp.numCol_ = "
                     "%" HIGHSINT_FORMAT
                     ") [%11g, %11g, %11g] so nonbasic move should be "
                     "up=%2" HIGHSINT_FORMAT
                     " but is  "
                     "%" HIGHSINT_FORMAT "\n",
                     var, lp.numCol_, info.workLower_[var],
                     info.workValue_[var], info.workUpper_[var],
                     kNonbasicMoveUp, basis.nonbasicMove_[var]);
        return ok;
      }
      ok = info.workValue_[var] == info.workLower_[var];
      if (!ok) {
        highsLogUser(
            options.log_options, HighsLogType::kError,
            "Finite lower bound and infinite upper bound variable "
            "%" HIGHSINT_FORMAT
            " "
            "(lp.numCol_ = "
            "%" HIGHSINT_FORMAT ") so work value should be %g but is %g\n",
            var, lp.numCol_, info.workLower_[var], info.workValue_[var]);
        return ok;
      }
    }
  } else {
    // Infinite lower bound
    if (!highs_isInfinity(info.workUpper_[var])) {
      ok = basis.nonbasicMove_[var] == kNonbasicMoveDn;
      if (!ok) {
        highsLogUser(
            options.log_options, HighsLogType::kError,
            "Finite upper bound and infinite lower bound variable "
            "%" HIGHSINT_FORMAT
            " "
            "(lp.numCol_ = "
            "%" HIGHSINT_FORMAT
            ") [%11g, %11g, %11g] so nonbasic move should be down but is  "
            "%" HIGHSINT_FORMAT "\n",
            var, lp.numCol_, info.workLower_[var], info.workValue_[var],
            info.workUpper_[var], basis.nonbasicMove_[var]);
        return ok;
      }
      ok = info.workValue_[var] == info.workUpper_[var];
      if (!ok) {
        highsLogUser(
            options.log_options, HighsLogType::kError,
            "Finite upper bound and infinite lower bound variable "
            "%" HIGHSINT_FORMAT
            " "
            "(lp.numCol_ = "
            "%" HIGHSINT_FORMAT ") so work value should be %g but is %g\n",
            var, lp.numCol_, info.workUpper_[var], info.workValue_[var]);
        return ok;
      }
    } else {
      // Infinite upper bound
      ok = basis.nonbasicMove_[var] == kNonbasicMoveZe;
      if (!ok) {
        highsLogUser(
            options.log_options, HighsLogType::kError,
            "Free variable %" HIGHSINT_FORMAT " (lp.numCol_ = %" HIGHSINT_FORMAT
            ") [%11g, %11g, %11g] "
            "so nonbasic "
            "move should be zero but is  %" HIGHSINT_FORMAT "\n",
            var, lp.numCol_, info.workLower_[var], info.workValue_[var],
            info.workUpper_[var], basis.nonbasicMove_[var]);
        return ok;
      }
      ok = info.workValue_[var] == 0.0;
      if (!ok) {
        highsLogUser(options.log_options, HighsLogType::kError,
                     "Free variable %" HIGHSINT_FORMAT
                     " (lp.numCol_ = %" HIGHSINT_FORMAT
                     ") so work value should "
                     "be zero but "
                     "is %g\n",
                     var, lp.numCol_, info.workValue_[var]);
        return ok;
      }
    }
  }
  // ok must be true if we reach here
  assert(ok);
  return ok;
}

void ekkDebugReportReinvertOnNumericalTrouble(
    const std::string method_name, const HEkk& ekk_instance,
    const double numerical_trouble_measure, const double alpha_from_col,
    const double alpha_from_row, const double numerical_trouble_tolerance,
    const bool reinvert) {
  if (ekk_instance.options_.highs_debug_level < kHighsDebugLevelCheap) return;
  const double abs_alpha_from_col = abs(alpha_from_col);
  const double abs_alpha_from_row = abs(alpha_from_row);
  const double abs_alpha_diff = abs(abs_alpha_from_col - abs_alpha_from_row);
  const HighsInt iteration_count = ekk_instance.iteration_count_;
  const HighsInt update_count = ekk_instance.info_.update_count;
  const std::string model_name = ekk_instance.lp_.model_name_;

  const bool numerical_trouble =
      numerical_trouble_measure > numerical_trouble_tolerance;
  const bool near_numerical_trouble =
      10 * numerical_trouble_measure > numerical_trouble_tolerance;

  const bool wrong_sign = alpha_from_col * alpha_from_row <= 0;
  if (!near_numerical_trouble && !wrong_sign) return;
  std::string adjective;
  if (numerical_trouble) {
    adjective = "       exceeds";
  } else if (near_numerical_trouble) {
    adjective = "almost exceeds";
  } else {
    adjective = "clearly satisfies";
  }
  highsLogUser(ekk_instance.options_.log_options, HighsLogType::kWarning,
               "%s (%s) [Iter %" HIGHSINT_FORMAT "; Update %" HIGHSINT_FORMAT
               "] Col: %11.4g; Row: %11.4g; Diff "
               "= %11.4g: Measure %11.4g %s %11.4g\n",
               method_name.c_str(), model_name.c_str(), iteration_count,
               update_count, abs_alpha_from_col, abs_alpha_from_row,
               abs_alpha_diff, numerical_trouble_measure, adjective.c_str(),
               numerical_trouble_tolerance);
  if (wrong_sign) {
    highsLogUser(ekk_instance.options_.log_options, HighsLogType::kWarning,
                 "   Incompatible signs for Col: %11.4g and Row: %11.4g\n",
                 alpha_from_col, alpha_from_row);
  }
  if ((numerical_trouble || wrong_sign) && !reinvert) {
    highsLogUser(ekk_instance.options_.log_options, HighsLogType::kWarning,
                 "   Numerical trouble or wrong sign and not reinverting\n");
  }
}

HighsDebugStatus ekkDebugUpdatedDual(const HighsOptions& options,
                                     const double updated_dual,
                                     const double computed_dual) {
  if (options.highs_debug_level < kHighsDebugLevelCheap)
    return HighsDebugStatus::kNotChecked;
  HighsDebugStatus return_status = HighsDebugStatus::kOk;
  std::string error_adjective;
  HighsLogType report_level;
  double updated_dual_absolute_error = abs(updated_dual - computed_dual);
  double updated_dual_relative_error =
      updated_dual_absolute_error / max(abs(computed_dual), 1.0);
  bool sign_error = updated_dual * computed_dual <= 0;
  bool at_least_small_error =
      sign_error ||
      updated_dual_absolute_error > updated_dual_small_absolute_error ||
      updated_dual_relative_error > updated_dual_small_relative_error;
  if (!at_least_small_error) return return_status;

  if (updated_dual_relative_error > updated_dual_large_relative_error ||
      updated_dual_absolute_error > updated_dual_large_absolute_error) {
    error_adjective = "Large";
    report_level = HighsLogType::kInfo;
    return_status = HighsDebugStatus::kLargeError;
  } else if (updated_dual_relative_error > updated_dual_small_relative_error ||
             updated_dual_absolute_error > updated_dual_small_absolute_error) {
    error_adjective = "Small";
    report_level = HighsLogType::kDetailed;
    return_status = HighsDebugStatus::kSmallError;
  } else {
    error_adjective = "OK";
    report_level = HighsLogType::kVerbose;
    return_status = HighsDebugStatus::kOk;
  }
  if (sign_error) {
    report_level = HighsLogType::kInfo;
    return_status = HighsDebugStatus::kLargeError;
  }
  highsLogDev(
      options.log_options, report_level,
      "UpdatedDual:  %-9s absolute (%9.4g) or relative (%9.4g) error in "
      "updated dual value",
      error_adjective.c_str(), updated_dual_absolute_error,
      updated_dual_relative_error);
  if (sign_error) {
    highsLogDev(options.log_options, report_level,
                ": Also sign error with (%9.4g, %9.4g)\n", updated_dual,
                computed_dual);
  } else {
    highsLogDev(options.log_options, report_level, "\n");
  }
  return return_status;
}

HighsDebugStatus ekkDebugNonbasicFreeColumnSet(
    const HEkk& ekk_instance, const HighsInt num_free_col,
    const HSet nonbasic_free_col_set) {
  const HighsOptions& options = ekk_instance.options_;
  if (options.highs_debug_level < kHighsDebugLevelCheap)
    return HighsDebugStatus::kNotChecked;
  const HighsLp& lp = ekk_instance.lp_;
  const HighsSimplexInfo& info = ekk_instance.info_;
  const SimplexBasis& basis = ekk_instance.basis_;
  HighsInt num_tot = lp.numCol_ + lp.numRow_;

  // Check the number of free columns
  HighsInt check_num_free_col = 0;
  for (HighsInt iVar = 0; iVar < num_tot; iVar++) {
    if (info.workLower_[iVar] <= -kHighsInf &&
        info.workUpper_[iVar] >= kHighsInf)
      check_num_free_col++;
  }
  if (check_num_free_col != num_free_col) {
    highsLogDev(options.log_options, HighsLogType::kInfo,
                "NonbasicFreeColumnData: Number of free columns should be "
                "%" HIGHSINT_FORMAT ", not %" HIGHSINT_FORMAT "\n",
                check_num_free_col, num_free_col);
    return HighsDebugStatus::kLogicalError;
  }
  if (!num_free_col) return HighsDebugStatus::kOk;
  // Debug HSet nonbasic_free_col
  bool nonbasic_free_col_ok = nonbasic_free_col_set.debug();
  if (!nonbasic_free_col_ok) {
    highsLogDev(options.log_options, HighsLogType::kInfo,
                "NonbasicFreeColumnData: HSet error\n");
    return HighsDebugStatus::kLogicalError;
  }

  // Check that we have the right number of nonbasic free columns
  const HighsInt& num_nonbasic_free_col = nonbasic_free_col_set.count();
  HighsInt check_num_nonbasic_free_col = 0;
  for (HighsInt iVar = 0; iVar < num_tot; iVar++) {
    bool nonbasic_free = basis.nonbasicFlag_[iVar] == kNonbasicFlagTrue &&
                         info.workLower_[iVar] <= -kHighsInf &&
                         info.workUpper_[iVar] >= kHighsInf;
    if (nonbasic_free) check_num_nonbasic_free_col++;
  }
  if (check_num_nonbasic_free_col != num_nonbasic_free_col) {
    highsLogDev(options.log_options, HighsLogType::kInfo,
                "NonbasicFreeColumnData: Set should have %" HIGHSINT_FORMAT
                " entries, not %" HIGHSINT_FORMAT "\n",
                check_num_nonbasic_free_col, num_nonbasic_free_col);
    return HighsDebugStatus::kLogicalError;
  }
  // Check that all in the set are nonbasic free columns
  const vector<HighsInt>& nonbasic_free_col_set_entry =
      nonbasic_free_col_set.entry();
  for (HighsInt ix = 0; ix < num_nonbasic_free_col; ix++) {
    HighsInt iVar = nonbasic_free_col_set_entry[ix];
    bool nonbasic_free = basis.nonbasicFlag_[iVar] == kNonbasicFlagTrue &&
                         info.workLower_[iVar] <= -kHighsInf &&
                         info.workUpper_[iVar] >= kHighsInf;
    if (!nonbasic_free) {
      highsLogDev(options.log_options, HighsLogType::kInfo,
                  "NonbasicFreeColumnData: Variable %" HIGHSINT_FORMAT
                  " in nonbasic free "
                  "set has nonbasicFlag = %" HIGHSINT_FORMAT
                  " and bounds [%g, %g]\n",
                  iVar, basis.nonbasicFlag_[iVar], info.workLower_[iVar],
                  info.workUpper_[iVar]);
      return HighsDebugStatus::kLogicalError;
    }
  }
  return HighsDebugStatus::kOk;
}

HighsDebugStatus ekkDebugRowMatrix(const HEkk& ekk_instance) {
  /*
  printf("Checking row-wise matrix\n");
  for (HighsInt row = 0; row < numRow; row++) {
    for (HighsInt el = ARstart[row]; el < AR_Nend[row]; el++) {
      HighsInt col = ARindex[el];
      if (!nonbasicFlag_[col]) {
        printf("Row-wise matrix error: col %" HIGHSINT_FORMAT ", (el = %"
  HIGHSINT_FORMAT " for row %" HIGHSINT_FORMAT ") is basic\n", col, el, row);
        return false;
      }
    }
    for (HighsInt el = AR_Nend[row]; el < ARstart[row + 1]; el++) {
      HighsInt col = ARindex[el];
      if (nonbasicFlag_[col]) {
        printf(
            "Row-wise matrix error: col %" HIGHSINT_FORMAT ", (el = %"
  HIGHSINT_FORMAT " for row %" HIGHSINT_FORMAT ") is nonbasic\n", col, el, row);
        return false;
      }
    }
  }
  */
  return HighsDebugStatus::kOk;
}
