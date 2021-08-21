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
/**@file lp_data/HighsUtils.cpp
 * @brief Class-independent utilities for HiGHS
 */
#include "lp_data/HighsLpUtils.h"

#include <algorithm>
#include <cassert>

#include "HConfig.h"
#include "io/Filereader.h"
#include "io/HMPSIO.h"
#include "io/HighsIO.h"
#include "lp_data/HighsModelUtils.h"
#include "lp_data/HighsStatus.h"
#include "util/HighsSort.h"
#include "util/HighsTimer.h"

HighsStatus assessLp(HighsLp& lp, const HighsOptions& options) {
  HighsStatus return_status = HighsStatus::kOk;
  HighsStatus call_status;
  // Assess the LP dimensions and vector sizes, returning on error
  call_status = assessLpDimensions(options, lp);
  return_status =
      interpretCallStatus(call_status, return_status, "assessLpDimensions");
  if (return_status == HighsStatus::kError) return return_status;

  // If the LP has no columns there is nothing left to test
  if (lp.numCol_ == 0) return HighsStatus::kOk;
  assert(lp.orientation_ == MatrixOrientation::kColwise);

  // From here, any LP has lp.numCol_ > 0 and lp.Astart_[lp.numCol_] exists (as
  // the number of nonzeros)
  assert(lp.numCol_ > 0);

  // Assess the LP column costs
  HighsIndexCollection index_collection;
  index_collection.dimension_ = lp.numCol_;
  index_collection.is_interval_ = true;
  index_collection.from_ = 0;
  index_collection.to_ = lp.numCol_ - 1;
  call_status = assessCosts(options, 0, index_collection, lp.colCost_,
                            options.infinite_cost);
  return_status =
      interpretCallStatus(call_status, return_status, "assessCosts");
  if (return_status == HighsStatus::kError) return return_status;
  // Assess the LP column bounds
  call_status = assessBounds(options, "Col", 0, index_collection, lp.colLower_,
                             lp.colUpper_, options.infinite_bound);
  return_status =
      interpretCallStatus(call_status, return_status, "assessBounds");
  if (return_status == HighsStatus::kError) return return_status;
  if (lp.numRow_) {
    // Assess the LP row bounds
    index_collection.dimension_ = lp.numRow_;
    index_collection.is_interval_ = true;
    index_collection.from_ = 0;
    index_collection.to_ = lp.numRow_ - 1;
    call_status =
        assessBounds(options, "Row", 0, index_collection, lp.rowLower_,
                     lp.rowUpper_, options.infinite_bound);
    return_status =
        interpretCallStatus(call_status, return_status, "assessBounds");
    if (return_status == HighsStatus::kError) return return_status;
    // Assess the LP matrix
    call_status =
        assessMatrix(options.log_options, "LP", lp.numRow_, lp.numCol_,
                     lp.Astart_, lp.Aindex_, lp.Avalue_,
                     options.small_matrix_value, options.large_matrix_value);
    return_status =
        interpretCallStatus(call_status, return_status, "assessMatrix");
    if (return_status == HighsStatus::kError) return return_status;
    HighsInt lp_num_nz = lp.Astart_[lp.numCol_];
    // If entries have been removed from the matrix, resize the index
    // and value vectors to prevent bug in presolve
    if ((HighsInt)lp.Aindex_.size() > lp_num_nz) lp.Aindex_.resize(lp_num_nz);
    if ((HighsInt)lp.Avalue_.size() > lp_num_nz) lp.Avalue_.resize(lp_num_nz);
  }
  if (return_status == HighsStatus::kError)
    return_status = HighsStatus::kError;
  else
    return_status = HighsStatus::kOk;
  if (return_status != HighsStatus::kOk)
    highsLogDev(options.log_options, HighsLogType::kInfo,
                "assessLp returns HighsStatus = %s\n",
                HighsStatusToString(return_status).c_str());
  return return_status;
}

HighsStatus assessLpDimensions(const HighsOptions& options, const HighsLp& lp) {
  HighsStatus return_status = HighsStatus::kOk;

  // Use error_found to track whether an error has been found in multiple tests
  bool error_found = false;

  // Don't expect the matrix_start_size to be legal if there are no columns
  bool check_matrix_start_size = lp.numCol_ > 0;

  // Assess column-related dimensions
  bool legal_num_col = lp.numCol_ >= 0;
  if (!legal_num_col) {
    highsLogUser(options.log_options, HighsLogType::kError,
                 "LP has illegal number of cols = %" HIGHSINT_FORMAT "\n",
                 lp.numCol_);
    error_found = true;
  } else {
    // Check the size of the column vectors
    HighsInt col_cost_size = lp.colCost_.size();
    HighsInt col_lower_size = lp.colLower_.size();
    HighsInt col_upper_size = lp.colUpper_.size();
    HighsInt matrix_start_size = lp.Astart_.size();
    bool legal_col_cost_size = col_cost_size >= lp.numCol_;
    bool legal_col_lower_size = col_lower_size >= lp.numCol_;
    bool legal_col_upper_size = col_lower_size >= lp.numCol_;

    if (!legal_col_cost_size) {
      highsLogUser(options.log_options, HighsLogType::kError,
                   "LP has illegal colCost size = %" HIGHSINT_FORMAT
                   " < %" HIGHSINT_FORMAT "\n",
                   col_cost_size, lp.numCol_);
      error_found = true;
    }
    if (!legal_col_lower_size) {
      highsLogUser(options.log_options, HighsLogType::kError,
                   "LP has illegal colLower size = %" HIGHSINT_FORMAT
                   " < %" HIGHSINT_FORMAT "\n",
                   col_lower_size, lp.numCol_);
      error_found = true;
    }
    if (!legal_col_upper_size) {
      highsLogUser(options.log_options, HighsLogType::kError,
                   "LP has illegal colUpper size = %" HIGHSINT_FORMAT
                   " < %" HIGHSINT_FORMAT "\n",
                   col_upper_size, lp.numCol_);
      error_found = true;
    }
  }

  // Assess row-related dimensions
  bool legal_num_row = lp.numRow_ >= 0;
  if (!legal_num_row) {
    highsLogUser(options.log_options, HighsLogType::kError,
                 "LP has illegal number of rows = %" HIGHSINT_FORMAT "\n",
                 lp.numRow_);
    error_found = true;
  } else {
    HighsInt row_lower_size = lp.rowLower_.size();
    HighsInt row_upper_size = lp.rowUpper_.size();
    bool legal_row_lower_size = row_lower_size >= lp.numRow_;
    bool legal_row_upper_size = row_lower_size >= lp.numRow_;
    if (!legal_row_lower_size) {
      highsLogUser(options.log_options, HighsLogType::kError,
                   "LP has illegal rowLower size = %" HIGHSINT_FORMAT
                   " < %" HIGHSINT_FORMAT "\n",
                   row_lower_size, lp.numRow_);
      error_found = true;
    }
    if (!legal_row_upper_size) {
      highsLogUser(options.log_options, HighsLogType::kError,
                   "LP has illegal rowUpper size = %" HIGHSINT_FORMAT
                   " < %" HIGHSINT_FORMAT "\n",
                   row_upper_size, lp.numRow_);
      error_found = true;
    }
  }

  // Assess matrix-related dimensions
  if (assessMatrixDimensions(options.log_options, "LP", lp.numCol_, lp.Astart_,
                             lp.Aindex_, lp.Avalue_) == HighsStatus::kError) {
    error_found = true;
  }
  assert(!error_found);
  if (error_found)
    return_status = HighsStatus::kError;
  else
    return_status = HighsStatus::kOk;

  return return_status;
}

HighsStatus assessCosts(const HighsOptions& options, const HighsInt ml_col_os,
                        const HighsIndexCollection& index_collection,
                        vector<double>& cost, const double infinite_cost) {
  HighsStatus return_status = HighsStatus::kOk;
  // Check parameters for technique and, if OK set the loop limits
  if (!assessIndexCollection(options.log_options, index_collection))
    return interpretCallStatus(HighsStatus::kError, return_status,
                               "assessIndexCollection");
  HighsInt from_k;
  HighsInt to_k;
  if (!limitsForIndexCollection(options.log_options, index_collection, from_k,
                                to_k))
    return interpretCallStatus(HighsStatus::kError, return_status,
                               "limitsForIndexCollection");
  if (from_k > to_k) return return_status;

  return_status = HighsStatus::kOk;
  bool error_found = false;
  // Work through the data to be assessed.
  //
  // Loop is k \in [from_k...to_k) covering the entries in the
  // interval, set or mask to be considered.
  //
  // For an interval or mask, these values of k are the columns to be
  // considered in a local sense, as well as the entries in the
  // cost data to be assessed
  //
  // For a set, these values of k are the indices in the set, from
  // which the columns to be considered in a local sense are
  // drawn. The entries in the cost data to be assessed correspond
  // to the values of k
  //
  // Adding the value of ml_col_os to local_col yields the value of
  // ml_col, being the column in a global (whole-model) sense. This is
  // necessary when assessing the costs of columns being added to a
  // model, since they are specified using an interval
  // [0...num_new_col) which must be offset by the current number of
  // columns in the model.
  //
  HighsInt local_col;
  HighsInt ml_col;
  HighsInt usr_col = -1;
  for (HighsInt k = from_k; k < to_k + 1; k++) {
    if (index_collection.is_interval_ || index_collection.is_mask_) {
      local_col = k;
    } else {
      local_col = index_collection.set_[k];
    }
    if (index_collection.is_interval_) {
      usr_col++;
    } else {
      usr_col = k;
    }
    ml_col = ml_col_os + local_col;
    if (index_collection.is_mask_ && !index_collection.mask_[local_col])
      continue;
    double abs_cost = fabs(cost[usr_col]);
    bool legal_cost = abs_cost < infinite_cost;
    if (!legal_cost) {
      error_found = !kHighsAllowInfiniteCosts;
      HighsLogType log_type = HighsLogType::kWarning;
      if (error_found) log_type = HighsLogType::kError;
      highsLogUser(options.log_options, log_type,
                   "Col  %12" HIGHSINT_FORMAT " has |cost| of %12g >= %12g\n",
                   ml_col, abs_cost, infinite_cost);
    }
  }
  if (error_found)
    return_status = HighsStatus::kError;
  else
    return_status = HighsStatus::kOk;

  return return_status;
}

HighsStatus assessBounds(const HighsOptions& options, const char* type,
                         const HighsInt ml_ix_os,
                         const HighsIndexCollection& index_collection,
                         vector<double>& lower, vector<double>& upper,
                         const double infinite_bound) {
  HighsStatus return_status = HighsStatus::kOk;
  // Check parameters for technique and, if OK set the loop limits
  if (!assessIndexCollection(options.log_options, index_collection))
    return interpretCallStatus(HighsStatus::kError, return_status,
                               "assessIndexCollection");
  HighsInt from_k;
  HighsInt to_k;
  if (!limitsForIndexCollection(options.log_options, index_collection, from_k,
                                to_k))
    return interpretCallStatus(HighsStatus::kError, return_status,
                               "limitsForIndexCollection");
  if (from_k > to_k) return HighsStatus::kOk;

  return_status = HighsStatus::kOk;
  bool error_found = false;
  bool warning_found = false;
  // Work through the data to be assessed.
  //
  // Loop is k \in [from_k...to_k) covering the entries in the
  // interval, set or mask to be considered.
  //
  // For an interval or mask, these values of k are the row/column
  // indices to be considered in a local sense, as well as the entries
  // in the lower and upper bound data to be assessed
  //
  // For a set, these values of k are the indices in the set, from
  // which the indices to be considered in a local sense are
  // drawn. The entries in the lower and
  // upper bound data to be assessed correspond to the values of
  // k.
  //
  // Adding the value of ml_ix_os to local_ix yields the value of
  // ml_ix, being the index in a global (whole-model) sense. This is
  // necessary when assessing the bounds of rows/columns being added
  // to a model, since they are specified using an interval
  // [0...num_new_row/col) which must be offset by the current number
  // of rows/columns (generically indices) in the model.
  //
  HighsInt num_infinite_lower_bound = 0;
  HighsInt num_infinite_upper_bound = 0;
  HighsInt local_ix;
  HighsInt ml_ix;
  HighsInt usr_ix = -1;
  for (HighsInt k = from_k; k < to_k + 1; k++) {
    if (index_collection.is_interval_ || index_collection.is_mask_) {
      local_ix = k;
    } else {
      local_ix = index_collection.set_[k];
    }
    if (index_collection.is_interval_) {
      usr_ix++;
    } else {
      usr_ix = k;
    }
    ml_ix = ml_ix_os + local_ix;
    if (index_collection.is_mask_ && !index_collection.mask_[local_ix])
      continue;

    if (!highs_isInfinity(-lower[usr_ix])) {
      // Check whether a finite lower bound will be treated as -Infinity
      bool infinite_lower_bound = lower[usr_ix] <= -infinite_bound;
      if (infinite_lower_bound) {
        lower[usr_ix] = -kHighsInf;
        num_infinite_lower_bound++;
      }
    }
    if (!highs_isInfinity(upper[usr_ix])) {
      // Check whether a finite upper bound will be treated as Infinity
      bool infinite_upper_bound = upper[usr_ix] >= infinite_bound;
      if (infinite_upper_bound) {
        upper[usr_ix] = kHighsInf;
        num_infinite_upper_bound++;
      }
    }
    // Check that the lower bound does not exceed the upper bound
    bool legalLowerUpperBound = lower[usr_ix] <= upper[usr_ix];
    if (!legalLowerUpperBound) {
      // Leave inconsistent bounds to be used to deduce infeasibility
      highsLogUser(options.log_options, HighsLogType::kWarning,
                   "%3s  %12" HIGHSINT_FORMAT
                   " has inconsistent bounds [%12g, %12g]\n",
                   type, ml_ix, lower[usr_ix], upper[usr_ix]);
      warning_found = true;
    }
    // Check that the lower bound is not as much as +Infinity
    bool legalLowerBound = lower[usr_ix] < infinite_bound;
    if (!legalLowerBound) {
      highsLogUser(options.log_options, HighsLogType::kError,
                   "%3s  %12" HIGHSINT_FORMAT
                   " has lower bound of %12g >= %12g\n",
                   type, ml_ix, lower[usr_ix], infinite_bound);
      error_found = true;
    }
    // Check that the upper bound is not as little as -Infinity
    bool legalUpperBound = upper[usr_ix] > -infinite_bound;
    if (!legalUpperBound) {
      highsLogUser(options.log_options, HighsLogType::kError,
                   "%3s  %12" HIGHSINT_FORMAT
                   " has upper bound of %12g <= %12g\n",
                   type, ml_ix, upper[usr_ix], -infinite_bound);
      error_found = true;
    }
  }
  if (num_infinite_lower_bound) {
    highsLogUser(options.log_options, HighsLogType::kInfo,
                 "%3ss:%12" HIGHSINT_FORMAT
                 " lower bounds exceeding %12g are treated as -Infinity\n",
                 type, num_infinite_lower_bound, -infinite_bound);
  }
  if (num_infinite_upper_bound) {
    highsLogUser(options.log_options, HighsLogType::kInfo,
                 "%3ss:%12" HIGHSINT_FORMAT
                 " upper bounds exceeding %12g are treated as +Infinity\n",
                 type, num_infinite_upper_bound, infinite_bound);
  }

  if (error_found)
    return_status = HighsStatus::kError;
  else if (warning_found)
    return_status = HighsStatus::kWarning;
  else
    return_status = HighsStatus::kOk;

  return return_status;
}

HighsStatus cleanBounds(const HighsOptions& options, HighsLp& lp) {
  double max_residual = 0;
  HighsInt num_change = 0;
  for (HighsInt iCol = 0; iCol < lp.numCol_; iCol++) {
    double residual = lp.colLower_[iCol] - lp.colUpper_[iCol];
    if (residual > options.primal_feasibility_tolerance) {
      highsLogUser(options.log_options, HighsLogType::kError,
                   "Column %" HIGHSINT_FORMAT
                   " has inconsistent bounds [%g, %g] (residual = "
                   "%g) after presolve\n",
                   iCol, lp.colLower_[iCol], lp.colUpper_[iCol], residual);
      return HighsStatus::kError;
    } else if (residual > 0) {
      num_change++;
      max_residual = std::max(residual, max_residual);
      double mid = 0.5 * (lp.colLower_[iCol] + lp.colUpper_[iCol]);
      lp.colLower_[iCol] = mid;
      lp.colUpper_[iCol] = mid;
    }
  }
  for (HighsInt iRow = 0; iRow < lp.numRow_; iRow++) {
    double residual = lp.rowLower_[iRow] - lp.rowUpper_[iRow];
    if (residual > options.primal_feasibility_tolerance) {
      highsLogUser(options.log_options, HighsLogType::kError,
                   "Row %" HIGHSINT_FORMAT
                   " has inconsistent bounds [%g, %g] (residual = %g) "
                   "after presolve\n",
                   iRow, lp.rowLower_[iRow], lp.rowUpper_[iRow], residual);
      return HighsStatus::kError;
    } else if (residual > 0) {
      num_change++;
      max_residual = std::max(residual, max_residual);
      double mid = 0.5 * (lp.rowLower_[iRow] + lp.rowUpper_[iRow]);
      lp.rowLower_[iRow] = mid;
      lp.rowUpper_[iRow] = mid;
    }
  }
  if (num_change) {
    highsLogUser(options.log_options, HighsLogType::kWarning,
                 "Resolved %" HIGHSINT_FORMAT
                 " inconsistent bounds (maximum residual = "
                 "%9.4g) after presolve\n",
                 num_change, max_residual);
    return HighsStatus::kWarning;
  }
  return HighsStatus::kOk;
}

HighsStatus applyScalingToLp(const HighsLogOptions& log_options, HighsLp& lp,
                             const HighsScale& scale) {
  if (!scale.is_scaled) return HighsStatus::kOk;
  if ((HighsInt)scale.col.size() < lp.numCol_) return HighsStatus::kError;
  if ((HighsInt)scale.row.size() < lp.numRow_) return HighsStatus::kError;
  bool scale_error = false;
  // Set up column and row index collections for scaling
  HighsIndexCollection all_cols;
  all_cols.is_interval_ = true;
  all_cols.dimension_ = lp.numCol_;
  all_cols.from_ = 0;
  all_cols.to_ = lp.numCol_ - 1;
  HighsIndexCollection all_rows;
  all_rows.is_interval_ = true;
  all_rows.dimension_ = lp.numRow_;
  all_rows.from_ = 0;
  all_rows.to_ = lp.numRow_ - 1;

  scale_error = applyScalingToLpColCost(log_options, lp, scale.col, all_cols) !=
                    HighsStatus::kOk ||
                scale_error;
  scale_error = applyScalingToLpColBounds(log_options, lp, scale.col,
                                          all_cols) != HighsStatus::kOk ||
                scale_error;
  scale_error = applyScalingToLpRowBounds(log_options, lp, scale.row,
                                          all_rows) != HighsStatus::kOk ||
                scale_error;
  scale_error = applyScalingToLpMatrix(log_options, lp, &scale.col[0],
                                       &scale.row[0], 0, lp.numCol_ - 1, 0,
                                       lp.numRow_ - 1) != HighsStatus::kOk ||
                scale_error;
  if (scale_error) return HighsStatus::kError;
  return HighsStatus::kOk;
}

HighsStatus applyScalingToLpColCost(
    const HighsLogOptions& log_options, HighsLp& lp,
    const vector<double>& colScale,
    const HighsIndexCollection& index_collection) {
  HighsStatus return_status = HighsStatus::kOk;
  // Check parameters for technique and, if OK set the loop limits
  if (!assessIndexCollection(log_options, index_collection))
    return interpretCallStatus(HighsStatus::kError, return_status,
                               "assessIndexCollection");

  HighsInt from_k;
  HighsInt to_k;
  if (!limitsForIndexCollection(log_options, index_collection, from_k, to_k))
    return interpretCallStatus(HighsStatus::kError, return_status,
                               "limitsForIndexCollection");
  if (from_k > to_k) return HighsStatus::kOk;

  const bool& interval = index_collection.is_interval_;
  const bool& mask = index_collection.is_mask_;
  const HighsInt* col_set = index_collection.set_;
  const HighsInt* col_mask = index_collection.mask_;

  HighsInt local_col;
  HighsInt ml_col;
  const HighsInt ml_col_os = 0;
  for (HighsInt k = from_k; k < to_k + 1; k++) {
    if (interval || mask) {
      local_col = k;
    } else {
      local_col = col_set[k];
    }
    ml_col = ml_col_os + local_col;
    if (mask && !col_mask[local_col]) continue;
    lp.colCost_[ml_col] *= colScale[ml_col];
  }

  return HighsStatus::kOk;
}

HighsStatus applyScalingToLpColBounds(
    const HighsLogOptions& log_options, HighsLp& lp,
    const vector<double>& colScale,
    const HighsIndexCollection& index_collection) {
  HighsStatus return_status = HighsStatus::kOk;
  // Check parameters for technique and, if OK set the loop limits
  if (!assessIndexCollection(log_options, index_collection))
    return interpretCallStatus(HighsStatus::kError, return_status,
                               "assessIndexCollection");

  HighsInt from_k;
  HighsInt to_k;
  if (!limitsForIndexCollection(log_options, index_collection, from_k, to_k))
    return interpretCallStatus(HighsStatus::kError, return_status,
                               "limitsForIndexCollection");
  if (from_k > to_k) return HighsStatus::kOk;

  const bool& interval = index_collection.is_interval_;
  const bool& mask = index_collection.is_mask_;
  const HighsInt* col_set = index_collection.set_;
  const HighsInt* col_mask = index_collection.mask_;

  HighsInt local_col;
  HighsInt ml_col;
  const HighsInt ml_col_os = 0;
  for (HighsInt k = from_k; k < to_k + 1; k++) {
    if (interval || mask) {
      local_col = k;
    } else {
      local_col = col_set[k];
    }
    ml_col = ml_col_os + local_col;
    if (mask && !col_mask[local_col]) continue;
    if (!highs_isInfinity(-lp.colLower_[ml_col]))
      lp.colLower_[ml_col] /= colScale[ml_col];
    if (!highs_isInfinity(lp.colUpper_[ml_col]))
      lp.colUpper_[ml_col] /= colScale[ml_col];
  }

  return HighsStatus::kOk;
}

HighsStatus applyScalingToLpRowBounds(
    const HighsLogOptions& log_options, HighsLp& lp,
    const vector<double>& rowScale,
    const HighsIndexCollection& index_collection) {
  HighsStatus return_status = HighsStatus::kOk;
  // Check parameters for technique and, if OK set the loop limits
  if (!assessIndexCollection(log_options, index_collection))
    return interpretCallStatus(HighsStatus::kError, return_status,
                               "assessIndexCollection");

  HighsInt from_k;
  HighsInt to_k;
  if (!limitsForIndexCollection(log_options, index_collection, from_k, to_k))
    return interpretCallStatus(HighsStatus::kError, return_status,
                               "limitsForIndexCollection");
  if (from_k > to_k) return HighsStatus::kOk;

  const bool& interval = index_collection.is_interval_;
  const bool& mask = index_collection.is_mask_;
  const HighsInt* row_set = index_collection.set_;
  const HighsInt* row_mask = index_collection.mask_;

  HighsInt local_row;
  HighsInt ml_row;
  const HighsInt ml_row_os = 0;
  for (HighsInt k = from_k; k < to_k + 1; k++) {
    if (interval || mask) {
      local_row = k;
    } else {
      local_row = row_set[k];
    }
    ml_row = ml_row_os + local_row;
    if (mask && !row_mask[local_row]) continue;
    if (!highs_isInfinity(-lp.rowLower_[ml_row]))
      lp.rowLower_[ml_row] *= rowScale[ml_row];
    if (!highs_isInfinity(lp.rowUpper_[ml_row]))
      lp.rowUpper_[ml_row] *= rowScale[ml_row];
  }

  return HighsStatus::kOk;
}

HighsStatus applyScalingToLpMatrix(
    const HighsLogOptions& log_options, HighsLp& lp, const double* colScale,
    const double* rowScale, const HighsInt from_col, const HighsInt to_col,
    const HighsInt from_row, const HighsInt to_row) {
  if (from_col < 0) return HighsStatus::kError;
  if (to_col >= lp.numCol_) return HighsStatus::kError;
  if (from_row < 0) return HighsStatus::kError;
  if (to_row >= lp.numRow_) return HighsStatus::kError;
  if (colScale != NULL) {
    if (rowScale != NULL) {
      for (HighsInt iCol = from_col; iCol <= to_col; iCol++) {
        for (HighsInt iEl = lp.Astart_[iCol]; iEl < lp.Astart_[iCol + 1];
             iEl++) {
          HighsInt iRow = lp.Aindex_[iEl];
          if (iRow < from_row || iRow > to_row) continue;
          lp.Avalue_[iEl] *= (colScale[iCol] * rowScale[iRow]);
        }
      }
    } else {
      // No row scaling
      for (HighsInt iCol = from_col; iCol <= to_col; iCol++) {
        for (HighsInt iEl = lp.Astart_[iCol]; iEl < lp.Astart_[iCol + 1];
             iEl++) {
          HighsInt iRow = lp.Aindex_[iEl];
          if (iRow < from_row || iRow > to_row) continue;
          lp.Avalue_[iEl] *= colScale[iCol];
        }
      }
    }
  } else {
    // No column scaling
    if (rowScale != NULL) {
      for (HighsInt iCol = from_col; iCol <= to_col; iCol++) {
        for (HighsInt iEl = lp.Astart_[iCol]; iEl < lp.Astart_[iCol + 1];
             iEl++) {
          HighsInt iRow = lp.Aindex_[iEl];
          if (iRow < from_row || iRow > to_row) continue;
          lp.Avalue_[iEl] *= rowScale[iRow];
        }
      }
    }
  }
  return HighsStatus::kOk;
}

void applyRowScalingToMatrix(const vector<double>& rowScale,
                             const HighsInt numCol,
                             const vector<HighsInt>& Astart,
                             const vector<HighsInt>& Aindex,
                             vector<double>& Avalue) {
  for (HighsInt iCol = 0; iCol < numCol; iCol++) {
    for (HighsInt el = Astart[iCol]; el < Astart[iCol + 1]; el++) {
      Avalue[el] *= rowScale[Aindex[el]];
    }
  }
}

void colScaleMatrix(const HighsInt max_scale_factor_exponent, double* colScale,
                    const HighsInt numCol, const vector<HighsInt>& Astart,
                    const vector<HighsInt>& Aindex, vector<double>& Avalue) {
  const double log2 = log(2.0);
  const double max_allow_scale = pow(2.0, max_scale_factor_exponent);
  const double min_allow_scale = 1 / max_allow_scale;

  const double min_allow_col_scale = min_allow_scale;
  const double max_allow_col_scale = max_allow_scale;

  for (HighsInt iCol = 0; iCol < numCol; iCol++) {
    double col_max_value = 0;
    for (HighsInt k = Astart[iCol]; k < Astart[iCol + 1]; k++)
      col_max_value = max(fabs(Avalue[k]), col_max_value);
    if (col_max_value) {
      double col_scale_value = 1 / col_max_value;
      // Convert the col scale factor to the nearest power of two, and
      // ensure that it is not excessively large or small
      col_scale_value = pow(2.0, floor(log(col_scale_value) / log2 + 0.5));
      col_scale_value =
          min(max(min_allow_col_scale, col_scale_value), max_allow_col_scale);
      colScale[iCol] = col_scale_value;
      // Scale the column
      for (HighsInt k = Astart[iCol]; k < Astart[iCol + 1]; k++)
        Avalue[k] *= colScale[iCol];
    } else {
      // Empty column
      colScale[iCol] = 1;
    }
  }
}

HighsStatus applyScalingToLpCol(const HighsLogOptions& log_options, HighsLp& lp,
                                const HighsInt col, const double colScale) {
  if (col < 0) return HighsStatus::kError;
  if (col >= lp.numCol_) return HighsStatus::kError;
  if (!colScale) return HighsStatus::kError;

  for (HighsInt el = lp.Astart_[col]; el < lp.Astart_[col + 1]; el++)
    lp.Avalue_[el] *= colScale;
  lp.colCost_[col] *= colScale;
  if (colScale > 0) {
    lp.colLower_[col] /= colScale;
    lp.colUpper_[col] /= colScale;
  } else {
    const double new_upper = lp.colLower_[col] / colScale;
    lp.colLower_[col] = lp.colUpper_[col] / colScale;
    lp.colUpper_[col] = new_upper;
  }
  return HighsStatus::kOk;
}

HighsStatus applyScalingToLpRow(const HighsLogOptions& log_options, HighsLp& lp,
                                const HighsInt row, const double rowScale) {
  if (row < 0) return HighsStatus::kError;
  if (row >= lp.numRow_) return HighsStatus::kError;
  if (!rowScale) return HighsStatus::kError;

  for (HighsInt col = 0; col < lp.numCol_; col++) {
    for (HighsInt el = lp.Astart_[col]; el < lp.Astart_[col + 1]; el++) {
      if (lp.Aindex_[el] == row) lp.Avalue_[el] *= rowScale;
    }
  }
  if (rowScale > 0) {
    lp.rowLower_[row] /= rowScale;
    lp.rowUpper_[row] /= rowScale;
  } else {
    const double new_upper = lp.rowLower_[row] / rowScale;
    lp.rowLower_[row] = lp.rowUpper_[row] / rowScale;
    lp.rowUpper_[row] = new_upper;
  }
  return HighsStatus::kOk;
}

HighsStatus appendColsToLpVectors(HighsLp& lp, const HighsInt num_new_col,
                                  const vector<double>& colCost,
                                  const vector<double>& colLower,
                                  const vector<double>& colUpper) {
  if (num_new_col < 0) return HighsStatus::kError;
  if (num_new_col == 0) return HighsStatus::kOk;
  HighsInt new_num_col = lp.numCol_ + num_new_col;
  lp.colCost_.resize(new_num_col);
  lp.colLower_.resize(new_num_col);
  lp.colUpper_.resize(new_num_col);
  bool have_names = lp.col_names_.size();
  if (have_names) lp.col_names_.resize(new_num_col);
  for (HighsInt new_col = 0; new_col < num_new_col; new_col++) {
    HighsInt iCol = lp.numCol_ + new_col;
    lp.colCost_[iCol] = colCost[new_col];
    lp.colLower_[iCol] = colLower[new_col];
    lp.colUpper_[iCol] = colUpper[new_col];
    // Cannot guarantee to create unique names, so name is blank
    if (have_names) lp.col_names_[iCol] = "";
  }
  return HighsStatus::kOk;
}

HighsStatus appendRowsToLpVectors(HighsLp& lp, const HighsInt num_new_row,
                                  const vector<double>& rowLower,
                                  const vector<double>& rowUpper) {
  if (num_new_row < 0) return HighsStatus::kError;
  if (num_new_row == 0) return HighsStatus::kOk;
  HighsInt new_num_row = lp.numRow_ + num_new_row;
  lp.rowLower_.resize(new_num_row);
  lp.rowUpper_.resize(new_num_row);
  bool have_names = lp.row_names_.size();
  if (have_names) lp.row_names_.resize(new_num_row);

  for (HighsInt new_row = 0; new_row < num_new_row; new_row++) {
    HighsInt iRow = lp.numRow_ + new_row;
    lp.rowLower_[iRow] = rowLower[new_row];
    lp.rowUpper_[iRow] = rowUpper[new_row];
    // Cannot guarantee to create unique names, so name is blank
    if (have_names) lp.row_names_[iRow] = "";
  }
  return HighsStatus::kOk;
}

void appendToMatrix(HighsLp& lp, const HighsInt num_vec,
                    const HighsInt num_new_vec, const HighsInt num_new_nz,
                    const HighsInt* XAstart, const HighsInt* XAindex,
                    const double* XAvalue) {
  // Append packed vectors to a matrix
  // Determine the new number of vectors in the matrix and resize the
  // starts accordingly.
  HighsInt new_num_vec = num_vec + num_new_vec;
  lp.Astart_.resize(new_num_vec + 1);
  // If adding vectors to an empty LP then introduce the start for the
  // fictitious vector 0
  if (num_vec == 0) lp.Astart_[0] = 0;

  // Determine the current number of nonzeros and the new number of nonzeros
  HighsInt current_num_nz = lp.Astart_[num_vec];
  HighsInt new_num_nz = current_num_nz + num_new_nz;

  // Append the starts of the new vectors
  if (num_new_nz) {
    // Nontrivial number of nonzeros being added, so use XAstart
    assert(XAstart != NULL);
    for (HighsInt vec = 0; vec < num_new_vec; vec++)
      lp.Astart_[num_vec + vec] = current_num_nz + XAstart[vec];
  } else {
    // No nonzeros being added, so XAstart may be null, but entries of
    // zero are implied.
    for (HighsInt vec = 0; vec < num_new_vec; vec++)
      lp.Astart_[num_vec + vec] = current_num_nz;
  }
  lp.Astart_[num_vec + num_new_vec] = new_num_nz;

  // If no nonzeros are being added then there's nothing else to do
  if (num_new_nz <= 0) return;

  // Adding a non-trivial matrix: resize the matrix arrays accordingly
  lp.Aindex_.resize(new_num_nz);
  lp.Avalue_.resize(new_num_nz);
  // Copy in the new indices and values
  for (HighsInt el = 0; el < num_new_nz; el++) {
    lp.Aindex_[current_num_nz + el] = XAindex[el];
    lp.Avalue_[current_num_nz + el] = XAvalue[el];
  }
}

HighsStatus appendColsToLpMatrix(HighsLp& lp, const HighsInt num_new_col,
                                 const HighsInt num_new_nz,
                                 const HighsInt* XAstart,
                                 const HighsInt* XAindex,
                                 const double* XAvalue) {
  if (num_new_col < 0) return HighsStatus::kError;
  if (num_new_col == 0) return HighsStatus::kOk;
  // Check that nonzeros aren't being appended to a matrix with no rows
  if (num_new_nz > 0 && lp.numRow_ <= 0) return HighsStatus::kError;
  // Adding a positive number of columns to a matrix
  if (lp.orientation_ == MatrixOrientation::kNone) {
    // LP is currently empty, store the matrix column-wise
    assert(lp.numCol_ == 0 && lp.numRow_ == 0);
    lp.orientation_ = MatrixOrientation::kColwise;
  } else {
    // Ensure that the matrix is stored column-wise
    setOrientation(lp);
  }
  // Determine the new number of columns in the matrix and resize the
  // starts accordingly.
  HighsInt new_num_col = lp.numCol_ + num_new_col;
  lp.Astart_.resize(new_num_col + 1);
  // If adding columns to an empty LP then introduce the start for the
  // fictitious column 0
  if (lp.numCol_ == 0) lp.Astart_[0] = 0;

  // Determine the current number of nonzeros and the new number of nonzeros
  HighsInt current_num_nz = lp.Astart_[lp.numCol_];
  HighsInt new_num_nz = current_num_nz + num_new_nz;

  // Append the starts of the new columns
  if (num_new_nz) {
    // Nontrivial number of nonzeros being added, so use XAstart
    assert(XAstart != NULL);
    for (HighsInt col = 0; col < num_new_col; col++)
      lp.Astart_[lp.numCol_ + col] = current_num_nz + XAstart[col];
  } else {
    // No nonzeros being added, so XAstart may be null, but entries of
    // zero are implied.
    for (HighsInt col = 0; col < num_new_col; col++)
      lp.Astart_[lp.numCol_ + col] = current_num_nz;
  }
  lp.Astart_[lp.numCol_ + num_new_col] = new_num_nz;

  // If no nonzeros are being added then there's nothing else to do
  if (num_new_nz <= 0) return HighsStatus::kOk;

  // Adding a non-trivial matrix: resize the column-wise matrix arrays
  // accordingly
  lp.Aindex_.resize(new_num_nz);
  lp.Avalue_.resize(new_num_nz);
  // Copy in the new indices and values
  for (HighsInt el = 0; el < num_new_nz; el++) {
    lp.Aindex_[current_num_nz + el] = XAindex[el];
    lp.Avalue_[current_num_nz + el] = XAvalue[el];
  }
  return HighsStatus::kOk;
}

HighsStatus appendRowsToLpMatrix(HighsLp& lp, const HighsInt num_new_row,
                                 const HighsInt num_new_nz,
                                 const HighsInt* XARstart,
                                 const HighsInt* XARindex,
                                 const double* XARvalue) {
  if (num_new_row < 0) return HighsStatus::kError;
  if (num_new_row == 0) return HighsStatus::kOk;
  // Check that nonzeros aren't being appended to a matrix with no columns
  if (num_new_nz > 0 && lp.numCol_ <= 0) return HighsStatus::kError;
  // Adding a positive number of rows to a matrix
  HighsInt current_num_nz = 0;
  if (lp.orientation_ == MatrixOrientation::kNone) {
    // LP is currently empty, store the matrix row-wise
    assert(lp.numCol_ == 0 && lp.numRow_ == 0);
    lp.orientation_ = MatrixOrientation::kRowwise;
  } else if (lp.orientation_ == MatrixOrientation::kColwise) {
    assert(lp.numCol_ > 0);
    assert((HighsInt)lp.Astart_.size() >= lp.numCol_);
    current_num_nz = lp.Astart_[lp.numCol_];
    if (current_num_nz == 0) {
      // Matrix is currently empty and stored column-wise. It can be
      // converted trivially to row-wise storage so that rows can be
      // added easily.
      //
      // It's possible that the model could have columns and (empty)
      // rows - hence the assignment of zero starts for rows
      // 0...lp.numRow_.
      //
      // However, this allows efficient handling of the (common) case
      // where a modeller defines variables without constraints, and
      // then constraints one-by-one.
      lp.orientation_ = MatrixOrientation::kRowwise;
      lp.Astart_.assign(lp.numRow_ + 1, 0);
    }
  }
  if (lp.orientation_ == MatrixOrientation::kRowwise) {
    appendToMatrix(lp, lp.numRow_, num_new_row, num_new_nz, XARstart, XARindex,
                   XARvalue);
  } else {
    // Storing the matrix column-wise, so have to insert the new rows
    assert(lp.orientation_ == MatrixOrientation::kColwise);
    vector<HighsInt> Alength;
    Alength.assign(lp.numCol_, 0);
    for (HighsInt el = 0; el < num_new_nz; el++) Alength[XARindex[el]]++;
    // Determine the new number of nonzeros and resize the column-wise matrix
    // arrays
    HighsInt new_num_nz = current_num_nz + num_new_nz;
    lp.Aindex_.resize(new_num_nz);
    lp.Avalue_.resize(new_num_nz);
    // Append the new rows
    // Shift the existing columns to make space for the new entries
    HighsInt new_el = new_num_nz;
    for (HighsInt col = lp.numCol_ - 1; col >= 0; col--) {
      HighsInt start_col_plus_1 = new_el;
      new_el -= Alength[col];
      for (HighsInt el = lp.Astart_[col + 1] - 1; el >= lp.Astart_[col]; el--) {
        new_el--;
        lp.Aindex_[new_el] = lp.Aindex_[el];
        lp.Avalue_[new_el] = lp.Avalue_[el];
      }
      lp.Astart_[col + 1] = start_col_plus_1;
    }
    assert(new_el == 0);
    // Insert the new entries
    for (HighsInt row = 0; row < num_new_row; row++) {
      HighsInt first_el = XARstart[row];
      HighsInt last_el =
          (row < num_new_row - 1 ? XARstart[row + 1] : num_new_nz);
      for (HighsInt el = first_el; el < last_el; el++) {
        HighsInt col = XARindex[el];
        new_el = lp.Astart_[col + 1] - Alength[col];
        Alength[col]--;
        lp.Aindex_[new_el] = lp.numRow_ + row;
        lp.Avalue_[new_el] = XARvalue[el];
      }
    }
  }
  return HighsStatus::kOk;
}

HighsStatus deleteLpCols(const HighsLogOptions& log_options, HighsLp& lp,
                         const HighsIndexCollection& index_collection) {
  HighsInt new_num_col;
  HighsStatus call_status;
  call_status =
      deleteColsFromLpVectors(log_options, lp, new_num_col, index_collection);
  if (call_status != HighsStatus::kOk) return call_status;
  call_status = deleteColsFromLpMatrix(log_options, lp, index_collection);
  if (call_status != HighsStatus::kOk) return call_status;
  lp.numCol_ = new_num_col;
  return HighsStatus::kOk;
}

HighsStatus deleteColsFromLpVectors(
    const HighsLogOptions& log_options, HighsLp& lp, HighsInt& new_num_col,
    const HighsIndexCollection& index_collection) {
  HighsStatus return_status = HighsStatus::kOk;
  if (!assessIndexCollection(log_options, index_collection))
    return interpretCallStatus(HighsStatus::kError, return_status,
                               "assessIndexCollection");
  HighsInt from_k;
  HighsInt to_k;
  if (!limitsForIndexCollection(log_options, index_collection, from_k, to_k))
    return interpretCallStatus(HighsStatus::kError, return_status,
                               "limitsForIndexCollection");
  if (index_collection.is_set_) {
    // For deletion by set it must be increasing
    if (!increasingSetOk(index_collection.set_,
                         index_collection.set_num_entries_, 0, lp.numCol_ - 1,
                         true))
      return HighsStatus::kError;
  }
  // Initialise new_num_col in case none is removed due to from_k > to_k
  new_num_col = lp.numCol_;
  if (from_k > to_k) return HighsStatus::kOk;

  HighsInt delete_from_col;
  HighsInt delete_to_col;
  HighsInt keep_from_col;
  HighsInt keep_to_col = -1;
  HighsInt current_set_entry = 0;

  HighsInt col_dim = lp.numCol_;
  new_num_col = 0;
  bool have_names = lp.col_names_.size();
  for (HighsInt k = from_k; k <= to_k; k++) {
    updateIndexCollectionOutInIndex(index_collection, delete_from_col,
                                    delete_to_col, keep_from_col, keep_to_col,
                                    current_set_entry);
    // Account for the initial columns being kept
    if (k == from_k) new_num_col = delete_from_col;
    if (delete_to_col >= col_dim - 1) break;
    assert(delete_to_col < col_dim);
    for (HighsInt col = keep_from_col; col <= keep_to_col; col++) {
      lp.colCost_[new_num_col] = lp.colCost_[col];
      lp.colLower_[new_num_col] = lp.colLower_[col];
      lp.colUpper_[new_num_col] = lp.colUpper_[col];
      if (have_names) lp.col_names_[new_num_col] = lp.col_names_[col];
      new_num_col++;
    }
    if (keep_to_col >= col_dim - 1) break;
  }
  lp.colCost_.resize(new_num_col);
  lp.colLower_.resize(new_num_col);
  lp.colUpper_.resize(new_num_col);
  if (have_names) lp.col_names_.resize(new_num_col);
  return HighsStatus::kOk;
}

HighsStatus deleteColsFromLpMatrix(
    const HighsLogOptions& log_options, HighsLp& lp,
    const HighsIndexCollection& index_collection) {
  HighsStatus return_status = HighsStatus::kOk;
  if (!assessIndexCollection(log_options, index_collection))
    return interpretCallStatus(HighsStatus::kError, return_status,
                               "assessIndexCollection");
  HighsInt from_k;
  HighsInt to_k;
  if (!limitsForIndexCollection(log_options, index_collection, from_k, to_k))
    return interpretCallStatus(HighsStatus::kError, return_status,
                               "limitsForIndexCollection");
  if (index_collection.is_set_) {
    // For deletion by set it must be increasing
    if (!increasingSetOk(index_collection.set_,
                         index_collection.set_num_entries_, 0, lp.numCol_ - 1,
                         true))
      return HighsStatus::kError;
  }
  if (from_k > to_k) return HighsStatus::kOk;

  HighsInt delete_from_col;
  HighsInt delete_to_col;
  HighsInt keep_from_col;
  HighsInt keep_to_col = -1;
  HighsInt current_set_entry = 0;

  HighsInt col_dim = lp.numCol_;
  HighsInt new_num_col = 0;
  HighsInt new_num_nz = 0;
  for (HighsInt k = from_k; k <= to_k; k++) {
    updateIndexCollectionOutInIndex(index_collection, delete_from_col,
                                    delete_to_col, keep_from_col, keep_to_col,
                                    current_set_entry);
    if (k == from_k) {
      // Account for the initial columns being kept
      new_num_col = delete_from_col;
      new_num_nz = lp.Astart_[delete_from_col];
    }
    // Ensure that the starts of the deleted columns are zeroed to
    // avoid redundant start information for columns whose indices
    // are't used after the deletion takes place. In particular, if
    // all columns are deleted then something must be done to ensure
    // that the matrix isn't magially recreated by increasing the
    // number of columns from zero when there are no rows in the LP.
    for (HighsInt col = delete_from_col; col <= delete_to_col; col++)
      lp.Astart_[col] = 0;
    // Shift the starts - both in place and value - to account for the
    // columns and nonzeros removed
    const HighsInt keep_from_el = lp.Astart_[keep_from_col];
    for (HighsInt col = keep_from_col; col <= keep_to_col; col++) {
      lp.Astart_[new_num_col] = new_num_nz + lp.Astart_[col] - keep_from_el;
      new_num_col++;
    }
    for (HighsInt el = keep_from_el; el < lp.Astart_[keep_to_col + 1]; el++) {
      lp.Aindex_[new_num_nz] = lp.Aindex_[el];
      lp.Avalue_[new_num_nz] = lp.Avalue_[el];
      new_num_nz++;
    }
    if (keep_to_col >= col_dim - 1) break;
  }
  // Ensure that the start of the spurious last column is zeroed so
  // that it doesn't give a positive number of matrix entries if the
  // number of columns in the LP is increased when there are no rows
  // in the LP.
  lp.Astart_[lp.numCol_] = 0;
  lp.Astart_[new_num_col] = new_num_nz;
  lp.Astart_.resize(new_num_col + 1);
  lp.Aindex_.resize(new_num_nz);
  lp.Avalue_.resize(new_num_nz);
  return HighsStatus::kOk;
}

HighsStatus deleteLpRows(const HighsLogOptions& log_options, HighsLp& lp,
                         const HighsIndexCollection& index_collection) {
  HighsStatus return_status = HighsStatus::kOk;
  HighsStatus call_status;
  HighsInt new_num_row;
  call_status =
      deleteRowsFromLpVectors(log_options, lp, new_num_row, index_collection);
  return_status = interpretCallStatus(call_status, return_status,
                                      "deleteRowsFromLpVectors");
  if (return_status == HighsStatus::kError) return return_status;
  call_status = deleteRowsFromLpMatrix(log_options, lp, index_collection);
  return_status =
      interpretCallStatus(call_status, return_status, "deleteRowsFromLpMatrix");
  if (return_status == HighsStatus::kError) return return_status;
  lp.numRow_ = new_num_row;
  return HighsStatus::kOk;
}

HighsStatus deleteRowsFromLpVectors(
    const HighsLogOptions& log_options, HighsLp& lp, HighsInt& new_num_row,
    const HighsIndexCollection& index_collection) {
  HighsStatus return_status = HighsStatus::kOk;
  if (!assessIndexCollection(log_options, index_collection))
    return interpretCallStatus(HighsStatus::kError, return_status,
                               "assessIndexCollection");
  HighsInt from_k;
  HighsInt to_k;
  if (!limitsForIndexCollection(log_options, index_collection, from_k, to_k))
    return interpretCallStatus(HighsStatus::kError, return_status,
                               "limitsForIndexCollection");
  if (index_collection.is_set_) {
    // For deletion by set it must be increasing
    if (!increasingSetOk(index_collection.set_,
                         index_collection.set_num_entries_, 0, lp.numRow_ - 1,
                         true))
      return HighsStatus::kError;
  }
  // Initialise new_num_row in case none is removed due to from_k > to_k
  new_num_row = lp.numRow_;
  if (from_k > to_k) return HighsStatus::kOk;

  HighsInt delete_from_row;
  HighsInt delete_to_row;
  HighsInt keep_from_row;
  HighsInt keep_to_row = -1;
  HighsInt current_set_entry = 0;

  HighsInt row_dim = lp.numRow_;
  new_num_row = 0;
  bool have_names = (HighsInt)lp.row_names_.size() > 0;
  for (HighsInt k = from_k; k <= to_k; k++) {
    updateIndexCollectionOutInIndex(index_collection, delete_from_row,
                                    delete_to_row, keep_from_row, keep_to_row,
                                    current_set_entry);
    if (k == from_k) {
      // Account for the initial rows being kept
      new_num_row = delete_from_row;
    }
    if (delete_to_row >= row_dim - 1) break;
    assert(delete_to_row < row_dim);
    for (HighsInt row = keep_from_row; row <= keep_to_row; row++) {
      lp.rowLower_[new_num_row] = lp.rowLower_[row];
      lp.rowUpper_[new_num_row] = lp.rowUpper_[row];
      if (have_names) lp.row_names_[new_num_row] = lp.row_names_[row];
      new_num_row++;
    }
    if (keep_to_row >= row_dim - 1) break;
  }
  lp.rowLower_.resize(new_num_row);
  lp.rowUpper_.resize(new_num_row);
  if (have_names) lp.row_names_.resize(new_num_row);
  return HighsStatus::kOk;
}

HighsStatus deleteRowsFromLpMatrix(
    const HighsLogOptions& log_options, HighsLp& lp,
    const HighsIndexCollection& index_collection) {
  HighsStatus return_status = HighsStatus::kOk;
  if (!assessIndexCollection(log_options, index_collection))
    return interpretCallStatus(HighsStatus::kError, return_status,
                               "assessIndexCollection");
  HighsInt from_k;
  HighsInt to_k;
  if (!limitsForIndexCollection(log_options, index_collection, from_k, to_k))
    return interpretCallStatus(HighsStatus::kError, return_status,
                               "limitsForIndexCollection");
  if (index_collection.is_set_) {
    // For deletion by set it must be increasing
    if (!increasingSetOk(index_collection.set_,
                         index_collection.set_num_entries_, 0, lp.numRow_ - 1,
                         true))
      return HighsStatus::kError;
  }
  if (from_k > to_k) return HighsStatus::kOk;

  HighsInt delete_from_row;
  HighsInt delete_to_row;
  HighsInt keep_from_row;
  HighsInt row_dim = lp.numRow_;
  HighsInt keep_to_row = -1;
  HighsInt current_set_entry = 0;

  // Set up a row mask to indicate the new row index of kept rows and
  // -1 for deleted rows so that the kept entries in the column-wise
  // matrix can be identified and have their correct row index.
  vector<HighsInt> new_index;
  new_index.resize(lp.numRow_);
  HighsInt new_num_row = 0;
  bool mask = index_collection.is_mask_;
  const HighsInt* row_mask = index_collection.mask_;
  if (!mask) {
    keep_to_row = -1;
    current_set_entry = 0;
    for (HighsInt k = from_k; k <= to_k; k++) {
      updateIndexCollectionOutInIndex(index_collection, delete_from_row,
                                      delete_to_row, keep_from_row, keep_to_row,
                                      current_set_entry);
      if (k == from_k) {
        // Account for any initial rows being kept
        for (HighsInt row = 0; row < delete_from_row; row++) {
          new_index[row] = new_num_row;
          new_num_row++;
        }
      }
      for (HighsInt row = delete_from_row; row <= delete_to_row; row++) {
        new_index[row] = -1;
      }
      for (HighsInt row = keep_from_row; row <= keep_to_row; row++) {
        new_index[row] = new_num_row;
        new_num_row++;
      }
      if (keep_to_row >= row_dim - 1) break;
    }
  } else {
    for (HighsInt row = 0; row < lp.numRow_; row++) {
      if (row_mask[row]) {
        new_index[row] = -1;
      } else {
        new_index[row] = new_num_row;
        new_num_row++;
      }
    }
  }
  HighsInt new_num_nz = 0;
  for (HighsInt col = 0; col < lp.numCol_; col++) {
    HighsInt from_el = lp.Astart_[col];
    lp.Astart_[col] = new_num_nz;
    for (HighsInt el = from_el; el < lp.Astart_[col + 1]; el++) {
      HighsInt row = lp.Aindex_[el];
      HighsInt new_row = new_index[row];
      if (new_row >= 0) {
        lp.Aindex_[new_num_nz] = new_row;
        lp.Avalue_[new_num_nz] = lp.Avalue_[el];
        new_num_nz++;
      }
    }
  }
  lp.Astart_[lp.numCol_] = new_num_nz;
  lp.Astart_.resize(lp.numCol_ + 1);
  lp.Aindex_.resize(new_num_nz);
  lp.Avalue_.resize(new_num_nz);
  return HighsStatus::kOk;
}

HighsStatus changeLpMatrixCoefficient(HighsLp& lp, const HighsInt row,
                                      const HighsInt col,
                                      const double new_value) {
  if (row < 0 || row > lp.numRow_) return HighsStatus::kError;
  if (col < 0 || col > lp.numCol_) return HighsStatus::kError;
  HighsInt changeElement = -1;
  for (HighsInt el = lp.Astart_[col]; el < lp.Astart_[col + 1]; el++) {
    if (lp.Aindex_[el] == row) {
      changeElement = el;
      break;
    }
  }
  if (changeElement < 0) {
    changeElement = lp.Astart_[col + 1];
    HighsInt new_num_nz = lp.Astart_[lp.numCol_] + 1;
    lp.Aindex_.resize(new_num_nz);
    lp.Avalue_.resize(new_num_nz);
    for (HighsInt i = col + 1; i <= lp.numCol_; i++) lp.Astart_[i]++;
    for (HighsInt el = new_num_nz - 1; el > changeElement; el--) {
      lp.Aindex_[el] = lp.Aindex_[el - 1];
      lp.Avalue_[el] = lp.Avalue_[el - 1];
    }
  }
  lp.Aindex_[changeElement] = row;
  lp.Avalue_[changeElement] = new_value;

  return HighsStatus::kOk;
}

HighsStatus changeLpIntegrality(const HighsLogOptions& log_options, HighsLp& lp,
                                const HighsIndexCollection& index_collection,
                                const vector<HighsVarType>& new_integrality) {
  HighsStatus return_status = HighsStatus::kOk;
  // Check parameters for technique and, if OK set the loop limits
  if (!assessIndexCollection(log_options, index_collection))
    return interpretCallStatus(HighsStatus::kError, return_status,
                               "assessIndexCollection");
  HighsInt from_k;
  HighsInt to_k;
  if (!limitsForIndexCollection(log_options, index_collection, from_k, to_k))
    return interpretCallStatus(HighsStatus::kError, return_status,
                               "limitsForIndexCollection");
  if (from_k > to_k) return HighsStatus::kOk;

  const bool& interval = index_collection.is_interval_;
  const bool& mask = index_collection.is_mask_;
  const HighsInt* col_set = index_collection.set_;
  const HighsInt* col_mask = index_collection.mask_;

  // Change the integrality to the user-supplied integrality, according to the
  // technique
  HighsInt lp_col;
  HighsInt usr_col = -1;
  // May be adding integrality to a pure LP for which lp.integrality_
  // is of size 0.
  lp.integrality_.resize(lp.numCol_);
  for (HighsInt k = from_k; k < to_k + 1; k++) {
    if (interval || mask) {
      lp_col = k;
    } else {
      lp_col = col_set[k];
    }
    HighsInt col = lp_col;
    if (interval) {
      usr_col++;
    } else {
      usr_col = k;
    }
    if (mask && !col_mask[col]) continue;
    lp.integrality_[col] = new_integrality[usr_col];
  }
  return HighsStatus::kOk;
}

HighsStatus changeLpCosts(const HighsLogOptions& log_options, HighsLp& lp,
                          const HighsIndexCollection& index_collection,
                          const vector<double>& new_col_cost) {
  HighsStatus return_status = HighsStatus::kOk;
  // Check parameters for technique and, if OK set the loop limits
  if (!assessIndexCollection(log_options, index_collection))
    return interpretCallStatus(HighsStatus::kError, return_status,
                               "assessIndexCollection");
  HighsInt from_k;
  HighsInt to_k;
  if (!limitsForIndexCollection(log_options, index_collection, from_k, to_k))
    return interpretCallStatus(HighsStatus::kError, return_status,
                               "limitsForIndexCollection");
  if (from_k > to_k) return HighsStatus::kOk;

  const bool& interval = index_collection.is_interval_;
  const bool& mask = index_collection.is_mask_;
  const HighsInt* col_set = index_collection.set_;
  const HighsInt* col_mask = index_collection.mask_;

  // Change the costs to the user-supplied costs, according to the technique
  HighsInt lp_col;
  HighsInt usr_col = -1;
  for (HighsInt k = from_k; k < to_k + 1; k++) {
    if (interval || mask) {
      lp_col = k;
    } else {
      lp_col = col_set[k];
    }
    HighsInt col = lp_col;
    if (interval) {
      usr_col++;
    } else {
      usr_col = k;
    }
    if (mask && !col_mask[col]) continue;
    lp.colCost_[col] = new_col_cost[usr_col];
  }
  return HighsStatus::kOk;
}

HighsStatus changeLpColBounds(const HighsLogOptions& log_options, HighsLp& lp,
                              const HighsIndexCollection& index_collection,
                              const vector<double>& new_col_lower,
                              const vector<double>& new_col_upper) {
  return changeBounds(log_options, lp.colLower_, lp.colUpper_, index_collection,
                      new_col_lower, new_col_upper);
}

HighsStatus changeLpRowBounds(const HighsLogOptions& log_options, HighsLp& lp,
                              const HighsIndexCollection& index_collection,
                              const vector<double>& new_row_lower,
                              const vector<double>& new_row_upper) {
  return changeBounds(log_options, lp.rowLower_, lp.rowUpper_, index_collection,
                      new_row_lower, new_row_upper);
}

HighsStatus changeBounds(const HighsLogOptions& log_options,
                         vector<double>& lower, vector<double>& upper,
                         const HighsIndexCollection& index_collection,
                         const vector<double>& new_lower,
                         const vector<double>& new_upper) {
  HighsStatus return_status = HighsStatus::kOk;
  // Check parameters for technique and, if OK set the loop limits
  if (!assessIndexCollection(log_options, index_collection))
    return interpretCallStatus(HighsStatus::kError, return_status,
                               "assessIndexCollection");
  HighsInt from_k;
  HighsInt to_k;
  if (!limitsForIndexCollection(log_options, index_collection, from_k, to_k))
    return interpretCallStatus(HighsStatus::kError, return_status,
                               "limitsForIndexCollection");
  if (from_k > to_k) return HighsStatus::kOk;

  const bool& interval = index_collection.is_interval_;
  const bool& mask = index_collection.is_mask_;
  const HighsInt* ix_set = index_collection.set_;
  const HighsInt* ix_mask = index_collection.mask_;

  // Change the bounds to the user-supplied bounds, according to the technique
  HighsInt lp_ix;
  HighsInt usr_ix = -1;
  for (HighsInt k = from_k; k < to_k + 1; k++) {
    if (interval || mask) {
      lp_ix = k;
    } else {
      lp_ix = ix_set[k];
    }
    HighsInt ix = lp_ix;
    if (interval) {
      usr_ix++;
    } else {
      usr_ix = k;
    }
    if (mask && !ix_mask[ix]) continue;
    lower[ix] = new_lower[usr_ix];
    upper[ix] = new_upper[usr_ix];
  }
  return HighsStatus::kOk;
}

HighsInt getNumInt(const HighsLp& lp) {
  HighsInt num_int = 0;
  if (lp.integrality_.size()) {
    for (HighsInt iCol = 0; iCol < lp.numCol_; iCol++)
      if (lp.integrality_[iCol] == HighsVarType::kInteger) num_int++;
  }
  return num_int;
}

HighsStatus getLpCosts(const HighsLp& lp, const HighsInt from_col,
                       const HighsInt to_col, double* XcolCost) {
  if (from_col < 0 || to_col >= lp.numCol_) return HighsStatus::kError;
  if (from_col > to_col) return HighsStatus::kOk;
  for (HighsInt col = from_col; col < to_col + 1; col++)
    XcolCost[col - from_col] = lp.colCost_[col];
  return HighsStatus::kOk;
}

HighsStatus getLpColBounds(const HighsLp& lp, const HighsInt from_col,
                           const HighsInt to_col, double* XcolLower,
                           double* XcolUpper) {
  if (from_col < 0 || to_col >= lp.numCol_) return HighsStatus::kError;
  if (from_col > to_col) return HighsStatus::kOk;
  for (HighsInt col = from_col; col < to_col + 1; col++) {
    if (XcolLower != NULL) XcolLower[col - from_col] = lp.colLower_[col];
    if (XcolUpper != NULL) XcolUpper[col - from_col] = lp.colUpper_[col];
  }
  return HighsStatus::kOk;
}

HighsStatus getLpRowBounds(const HighsLp& lp, const HighsInt from_row,
                           const HighsInt to_row, double* XrowLower,
                           double* XrowUpper) {
  if (from_row < 0 || to_row >= lp.numRow_) return HighsStatus::kError;
  if (from_row > to_row) return HighsStatus::kOk;
  for (HighsInt row = from_row; row < to_row + 1; row++) {
    if (XrowLower != NULL) XrowLower[row - from_row] = lp.rowLower_[row];
    if (XrowUpper != NULL) XrowUpper[row - from_row] = lp.rowUpper_[row];
  }
  return HighsStatus::kOk;
}

// Get a single coefficient from the matrix
HighsStatus getLpMatrixCoefficient(const HighsLp& lp, const HighsInt Xrow,
                                   const HighsInt Xcol, double* val) {
  if (Xrow < 0 || Xrow >= lp.numRow_) return HighsStatus::kError;
  if (Xcol < 0 || Xcol >= lp.numCol_) return HighsStatus::kError;

  HighsInt get_el = -1;
  for (HighsInt el = lp.Astart_[Xcol]; el < lp.Astart_[Xcol + 1]; el++) {
    if (lp.Aindex_[el] == Xrow) {
      get_el = el;
      break;
    }
  }
  if (get_el < 0) {
    *val = 0;
  } else {
    *val = lp.Avalue_[get_el];
  }
  return HighsStatus::kOk;
}

// Methods for reporting an LP, including its row and column data and matrix
//
// Report the whole LP
void reportLp(const HighsLogOptions& log_options, const HighsLp& lp,
              const HighsLogType report_level) {
  reportLpBrief(log_options, lp);
  if ((HighsInt)report_level >= (HighsInt)HighsLogType::kDetailed) {
    reportLpColVectors(log_options, lp);
    reportLpRowVectors(log_options, lp);
    if ((HighsInt)report_level >= (HighsInt)HighsLogType::kVerbose)
      reportLpColMatrix(log_options, lp);
  }
}

// Report the LP briefly
void reportLpBrief(const HighsLogOptions& log_options, const HighsLp& lp) {
  reportLpDimensions(log_options, lp);
  reportLpObjSense(log_options, lp);
}

// Report the LP dimensions
void reportLpDimensions(const HighsLogOptions& log_options, const HighsLp& lp) {
  HighsInt lp_num_nz;
  if (lp.numCol_ == 0)
    lp_num_nz = 0;
  else
    lp_num_nz = lp.Astart_[lp.numCol_];
  highsLogDev(log_options, HighsLogType::kInfo,
              "LP has %" HIGHSINT_FORMAT " columns, %" HIGHSINT_FORMAT " rows",
              lp.numCol_, lp.numRow_);
  HighsInt num_int = getNumInt(lp);
  if (num_int) {
    highsLogDev(log_options, HighsLogType::kInfo,
                ", %" HIGHSINT_FORMAT " nonzeros and %" HIGHSINT_FORMAT
                " integer columns\n",
                lp_num_nz, num_int);
  } else {
    highsLogDev(log_options, HighsLogType::kInfo,
                " and %" HIGHSINT_FORMAT " nonzeros\n", lp_num_nz, num_int);
  }
}

// Report the LP objective sense
void reportLpObjSense(const HighsLogOptions& log_options, const HighsLp& lp) {
  if (lp.sense_ == ObjSense::kMinimize)
    highsLogDev(log_options, HighsLogType::kInfo,
                "Objective sense is minimize\n");
  else if (lp.sense_ == ObjSense::kMaximize)
    highsLogDev(log_options, HighsLogType::kInfo,
                "Objective sense is maximize\n");
  else
    highsLogDev(log_options, HighsLogType::kInfo,
                "Objective sense is ill-defined as %" HIGHSINT_FORMAT "\n",
                lp.sense_);
}

std::string getBoundType(const double lower, const double upper) {
  std::string type;
  if (highs_isInfinity(-lower)) {
    if (highs_isInfinity(upper)) {
      type = "FR";
    } else {
      type = "UB";
    }
  } else {
    if (highs_isInfinity(upper)) {
      type = "LB";
    } else {
      if (lower < upper) {
        type = "BX";
      } else {
        type = "FX";
      }
    }
  }
  return type;
}

// Report the vectors of LP column data
void reportLpColVectors(const HighsLogOptions& log_options, const HighsLp& lp) {
  if (lp.numCol_ <= 0) return;
  std::string type;
  HighsInt count;
  bool have_integer_columns = getNumInt(lp);
  bool have_col_names = lp.col_names_.size();

  highsLogDev(log_options, HighsLogType::kVerbose,
              "  Column        Lower        Upper         Cost       "
              "Type        Count");
  if (have_integer_columns)
    highsLogDev(log_options, HighsLogType::kVerbose, "  Discrete");
  if (have_col_names)
    highsLogDev(log_options, HighsLogType::kVerbose, "  Name");
  highsLogDev(log_options, HighsLogType::kVerbose, "\n");

  for (HighsInt iCol = 0; iCol < lp.numCol_; iCol++) {
    type = getBoundType(lp.colLower_[iCol], lp.colUpper_[iCol]);
    count = lp.Astart_[iCol + 1] - lp.Astart_[iCol];
    highsLogDev(log_options, HighsLogType::kVerbose,
                "%8" HIGHSINT_FORMAT
                " %12g %12g %12g         %2s %12" HIGHSINT_FORMAT "",
                iCol, lp.colLower_[iCol], lp.colUpper_[iCol], lp.colCost_[iCol],
                type.c_str(), count);
    if (have_integer_columns) {
      std::string integer_column = "";
      if (lp.integrality_[iCol] == HighsVarType::kInteger) {
        if (lp.colLower_[iCol] == 0 && lp.colUpper_[iCol] == 1) {
          integer_column = "Binary";
        } else {
          integer_column = "Integer";
        }
      }
      highsLogDev(log_options, HighsLogType::kVerbose, "  %-8s",
                  integer_column.c_str());
    }
    if (have_col_names)
      highsLogDev(log_options, HighsLogType::kVerbose, "  %-s",
                  lp.col_names_[iCol].c_str());
    highsLogDev(log_options, HighsLogType::kVerbose, "\n");
  }
}

// Report the vectors of LP row data
void reportLpRowVectors(const HighsLogOptions& log_options, const HighsLp& lp) {
  if (lp.numRow_ <= 0) return;
  std::string type;
  vector<HighsInt> count;
  bool have_row_names = lp.row_names_.size();

  count.resize(lp.numRow_, 0);
  if (lp.numCol_ > 0) {
    for (HighsInt el = 0; el < lp.Astart_[lp.numCol_]; el++)
      count[lp.Aindex_[el]]++;
  }

  highsLogDev(log_options, HighsLogType::kVerbose,
              "     Row        Lower        Upper       Type        Count");
  if (have_row_names)
    highsLogDev(log_options, HighsLogType::kVerbose, "  Name");
  highsLogDev(log_options, HighsLogType::kVerbose, "\n");

  for (HighsInt iRow = 0; iRow < lp.numRow_; iRow++) {
    type = getBoundType(lp.rowLower_[iRow], lp.rowUpper_[iRow]);
    std::string name = "";
    highsLogDev(log_options, HighsLogType::kVerbose,
                "%8" HIGHSINT_FORMAT
                " %12g %12g         %2s %12" HIGHSINT_FORMAT "",
                iRow, lp.rowLower_[iRow], lp.rowUpper_[iRow], type.c_str(),
                count[iRow]);
    if (have_row_names)
      highsLogDev(log_options, HighsLogType::kVerbose, "  %-s",
                  lp.row_names_[iRow].c_str());
    highsLogDev(log_options, HighsLogType::kVerbose, "\n");
  }
}

// Report the LP column-wise matrix
void reportLpColMatrix(const HighsLogOptions& log_options, const HighsLp& lp) {
  if (lp.numCol_ <= 0) return;
  if (lp.numRow_) {
    // With postitive number of rows, can assume that there are index and value
    // vectors to pass
    reportMatrix(log_options, "Column", lp.numCol_, lp.Astart_[lp.numCol_],
                 &lp.Astart_[0], &lp.Aindex_[0], &lp.Avalue_[0]);
  } else {
    // With no rows, can's assume that there are index and value vectors to pass
    reportMatrix(log_options, "Column", lp.numCol_, lp.Astart_[lp.numCol_],
                 &lp.Astart_[0], NULL, NULL);
  }
}

void reportMatrix(const HighsLogOptions& log_options, const std::string message,
                  const HighsInt num_col, const HighsInt num_nz,
                  const HighsInt* start, const HighsInt* index,
                  const double* value) {
  if (num_col <= 0) return;
  highsLogDev(log_options, HighsLogType::kVerbose,
              "%6s Index              Value\n", message.c_str());
  for (HighsInt col = 0; col < num_col; col++) {
    highsLogDev(log_options, HighsLogType::kVerbose,
                "    %8" HIGHSINT_FORMAT " Start   %10" HIGHSINT_FORMAT "\n",
                col, start[col]);
    HighsInt to_el = (col < num_col - 1 ? start[col + 1] : num_nz);
    for (HighsInt el = start[col]; el < to_el; el++)
      highsLogDev(log_options, HighsLogType::kVerbose,
                  "          %8" HIGHSINT_FORMAT " %12g\n", index[el],
                  value[el]);
  }
  highsLogDev(log_options, HighsLogType::kVerbose,
              "             Start   %10" HIGHSINT_FORMAT "\n", num_nz);
}

void analyseLp(const HighsLogOptions& log_options, const HighsLp& lp,
               const std::string message) {
  vector<double> min_colBound;
  vector<double> min_rowBound;
  vector<double> colRange;
  vector<double> rowRange;
  min_colBound.resize(lp.numCol_);
  min_rowBound.resize(lp.numRow_);
  colRange.resize(lp.numCol_);
  rowRange.resize(lp.numRow_);
  for (HighsInt col = 0; col < lp.numCol_; col++)
    min_colBound[col] = min(fabs(lp.colLower_[col]), fabs(lp.colUpper_[col]));
  for (HighsInt row = 0; row < lp.numRow_; row++)
    min_rowBound[row] = min(fabs(lp.rowLower_[row]), fabs(lp.rowUpper_[row]));
  for (HighsInt col = 0; col < lp.numCol_; col++)
    colRange[col] = lp.colUpper_[col] - lp.colLower_[col];
  for (HighsInt row = 0; row < lp.numRow_; row++)
    rowRange[row] = lp.rowUpper_[row] - lp.rowLower_[row];

  printf("\n%s model data: Analysis\n", message.c_str());
  analyseVectorValues(log_options, "Column costs", lp.numCol_, lp.colCost_);
  analyseVectorValues(log_options, "Column lower bounds", lp.numCol_,
                      lp.colLower_);
  analyseVectorValues(log_options, "Column upper bounds", lp.numCol_,
                      lp.colUpper_);
  analyseVectorValues(log_options, "Column min abs bound", lp.numCol_,
                      min_colBound);
  analyseVectorValues(log_options, "Column range", lp.numCol_, colRange);
  analyseVectorValues(log_options, "Row lower bounds", lp.numRow_,
                      lp.rowLower_);
  analyseVectorValues(log_options, "Row upper bounds", lp.numRow_,
                      lp.rowUpper_);
  analyseVectorValues(log_options, "Row min abs bound", lp.numRow_,
                      min_rowBound);
  analyseVectorValues(log_options, "Row range", lp.numRow_, rowRange);
  analyseVectorValues(log_options, "Matrix sparsity", lp.Astart_[lp.numCol_],
                      lp.Avalue_, true, lp.model_name_);
  analyseMatrixSparsity(log_options, "Constraint matrix", lp.numCol_,
                        lp.numRow_, lp.Astart_, lp.Aindex_);
  analyseModelBounds(log_options, "Column", lp.numCol_, lp.colLower_,
                     lp.colUpper_);
  analyseModelBounds(log_options, "Row", lp.numRow_, lp.rowLower_,
                     lp.rowUpper_);
}

void analyseScaledLp(const HighsLogOptions& log_options,
                     const HighsScale& scale, const HighsLp& scaled_lp) {
  if (!scale.is_scaled) return;
  analyseVectorValues(log_options, "Column scaling factors", scaled_lp.numCol_,
                      scale.col);
  analyseVectorValues(log_options, "Row    scaling factors", scaled_lp.numRow_,
                      scale.row);
  analyseLp(log_options, scaled_lp, "Scaled");
}

void writeSolutionToFile(FILE* file, const HighsLp& lp, const HighsBasis& basis,
                         const HighsSolution& solution, const bool pretty) {
  const bool have_value = solution.value_valid;
  const bool have_dual = solution.dual_valid;
  const bool have_basis = basis.valid;
  vector<double> use_col_value;
  vector<double> use_row_value;
  vector<double> use_col_dual;
  vector<double> use_row_dual;
  vector<HighsBasisStatus> use_col_status;
  vector<HighsBasisStatus> use_row_status;
  if (have_value) {
    use_col_value = solution.col_value;
    use_row_value = solution.row_value;
  }
  if (have_dual) {
    use_col_dual = solution.col_dual;
    use_row_dual = solution.row_dual;
  }
  if (have_basis) {
    use_col_status = basis.col_status;
    use_row_status = basis.row_status;
  }
  if (!have_value && !have_dual && !have_basis) return;
  if (pretty) {
    writeModelBoundSol(file, true, lp.numCol_, lp.colLower_, lp.colUpper_,
                       lp.col_names_, use_col_value, use_col_dual,
                       use_col_status);
    writeModelBoundSol(file, false, lp.numRow_, lp.rowLower_, lp.rowUpper_,
                       lp.row_names_, use_row_value, use_row_dual,
                       use_row_status);
  } else {
    fprintf(file,
            "%" HIGHSINT_FORMAT " %" HIGHSINT_FORMAT
            " : Number of columns and rows for primal or dual solution "
            "or basis\n",
            lp.numCol_, lp.numRow_);
    if (have_value) {
      fprintf(file, "T");
    } else {
      fprintf(file, "F");
    }
    fprintf(file, " Primal solution\n");
    if (have_dual) {
      fprintf(file, "T");
    } else {
      fprintf(file, "F");
    }
    fprintf(file, " Dual solution\n");
    if (have_basis) {
      fprintf(file, "T");
    } else {
      fprintf(file, "F");
    }
    fprintf(file, " Basis\n");
    fprintf(file, "Columns\n");
    for (HighsInt iCol = 0; iCol < lp.numCol_; iCol++) {
      if (have_value) fprintf(file, "%g", use_col_value[iCol]);
      if (have_dual) fprintf(file, "%g", use_col_dual[iCol]);
      if (have_basis)
        fprintf(file, " %" HIGHSINT_FORMAT "", (HighsInt)use_col_status[iCol]);
      fprintf(file, " \n");
    }
    fprintf(file, "Rows\n");
    for (HighsInt iRow = 0; iRow < lp.numRow_; iRow++) {
      if (have_value) fprintf(file, "%g", use_row_value[iRow]);
      if (have_dual) fprintf(file, "%g", use_row_dual[iRow]);
      if (have_basis)
        fprintf(file, " %" HIGHSINT_FORMAT "", (HighsInt)use_row_status[iRow]);
      fprintf(file, " \n");
    }
  }
}

HighsStatus writeBasisFile(const HighsLogOptions& log_options,
                           const HighsBasis& basis,
                           const std::string filename) {
  HighsStatus return_status = HighsStatus::kOk;
  if (basis.valid == false) {
    highsLogUser(log_options, HighsLogType::kError,
                 "writeBasisFile: Cannot write an invalid basis\n");
    return HighsStatus::kError;
  }
  std::ofstream outFile(filename);
  if (outFile.fail()) {
    highsLogUser(log_options, HighsLogType::kError,
                 "writeBasisFile: Cannot open writeable file \"%s\"\n",
                 filename.c_str());
    return HighsStatus::kError;
  }
  outFile << "HiGHS Version " << HIGHS_VERSION_MAJOR << std::endl;
  outFile << basis.col_status.size() << " " << basis.row_status.size()
          << std::endl;
  for (const auto& status : basis.col_status) {
    outFile << (HighsInt)status << " ";
  }
  outFile << std::endl;
  for (const auto& status : basis.row_status) {
    outFile << (HighsInt)status << " ";
  }
  outFile << std::endl;
  outFile << std::endl;
  outFile.close();
  return return_status;
}

HighsStatus readBasisFile(const HighsLogOptions& log_options, HighsBasis& basis,
                          const std::string filename) {
  // Reads a basis file, returning an error if what's read is
  // inconsistent with the sizes of the HighsBasis passed in
  HighsStatus return_status = HighsStatus::kOk;
  std::ifstream inFile(filename);
  if (inFile.fail()) {
    highsLogUser(log_options, HighsLogType::kError,
                 "readBasisFile: Cannot open readable file \"%s\"\n",
                 filename.c_str());
    return HighsStatus::kError;
  }
  std::string string_highs, string_version;
  HighsInt highs_version_number;
  inFile >> string_highs >> string_version >> highs_version_number;
  if (highs_version_number == 1) {
    HighsInt numCol, numRow;
    inFile >> numCol >> numRow;
    HighsInt basis_numCol = (HighsInt)basis.col_status.size();
    HighsInt basis_numRow = (HighsInt)basis.row_status.size();
    if (numCol != basis_numCol) {
      highsLogUser(log_options, HighsLogType::kError,
                   "readBasisFile: Basis file is for %" HIGHSINT_FORMAT
                   " columns, not %" HIGHSINT_FORMAT "\n",
                   numCol, basis_numCol);
      return HighsStatus::kError;
    }
    if (numRow != basis_numRow) {
      highsLogUser(log_options, HighsLogType::kError,
                   "readBasisFile: Basis file is for %" HIGHSINT_FORMAT
                   " rows, not %" HIGHSINT_FORMAT "\n",
                   numRow, basis_numRow);
      return HighsStatus::kError;
    }
    HighsInt int_status;
    for (HighsInt iCol = 0; iCol < numCol; iCol++) {
      inFile >> int_status;
      basis.col_status[iCol] = (HighsBasisStatus)int_status;
    }
    for (HighsInt iRow = 0; iRow < numRow; iRow++) {
      inFile >> int_status;
      basis.row_status[iRow] = (HighsBasisStatus)int_status;
    }
    if (inFile.eof()) {
      highsLogUser(
          log_options, HighsLogType::kError,
          "readBasisFile: Reached end of file before reading complete basis\n");
      return_status = HighsStatus::kError;
    }
  } else {
    highsLogUser(log_options, HighsLogType::kError,
                 "readBasisFile: Cannot read basis file for HiGHS version "
                 "%" HIGHSINT_FORMAT "\n",
                 highs_version_number);
    return_status = HighsStatus::kError;
  }
  inFile.close();
  return return_status;
}

HighsStatus calculateColDuals(const HighsLp& lp, HighsSolution& solution) {
  assert(solution.row_dual.size() > 0);
  if (!isSolutionRightSize(lp, solution)) return HighsStatus::kError;

  solution.col_dual.assign(lp.numCol_, 0);

  for (HighsInt col = 0; col < lp.numCol_; col++) {
    for (HighsInt i = lp.Astart_[col]; i < lp.Astart_[col + 1]; i++) {
      const HighsInt row = lp.Aindex_[i];
      assert(row >= 0);
      assert(row < lp.numRow_);
      // @FlipRowDual -= became +=
      solution.col_dual[col] += solution.row_dual[row] * lp.Avalue_[i];
    }
    solution.col_dual[col] += lp.colCost_[col];
  }

  return HighsStatus::kOk;
}

HighsStatus calculateRowValues(const HighsLp& lp, HighsSolution& solution) {
  // assert(solution.col_value.size() > 0);
  if (int(solution.col_value.size()) != lp.numCol_) return HighsStatus::kError;

  solution.row_value.clear();
  solution.row_value.assign(lp.numRow_, 0);

  for (HighsInt col = 0; col < lp.numCol_; col++) {
    for (HighsInt i = lp.Astart_[col]; i < lp.Astart_[col + 1]; i++) {
      const HighsInt row = lp.Aindex_[i];
      assert(row >= 0);
      assert(row < lp.numRow_);

      solution.row_value[row] += solution.col_value[col] * lp.Avalue_[i];
    }
  }

  return HighsStatus::kOk;
}

bool isBoundInfeasible(const HighsLogOptions& log_options, const HighsLp& lp) {
  HighsInt num_bound_infeasible = 0;
  for (HighsInt iCol = 0; iCol < lp.numCol_; iCol++)
    if (lp.colUpper_[iCol] < lp.colLower_[iCol]) num_bound_infeasible++;
  for (HighsInt iRow = 0; iRow < lp.numRow_; iRow++)
    if (lp.rowUpper_[iRow] < lp.rowLower_[iRow]) num_bound_infeasible++;
  if (num_bound_infeasible > 0)
    highsLogUser(log_options, HighsLogType::kInfo,
                 "Model infeasible due to %" HIGHSINT_FORMAT
                 " inconsistent bound(s)\n",
                 num_bound_infeasible);
  return num_bound_infeasible > 0;
}

bool isColDataNull(const HighsLogOptions& log_options,
                   const double* usr_col_cost, const double* usr_col_lower,
                   const double* usr_col_upper) {
  bool null_data = false;
  null_data =
      doubleUserDataNotNull(log_options, usr_col_cost, "column costs") ||
      null_data;
  null_data = doubleUserDataNotNull(log_options, usr_col_lower,
                                    "column lower bounds") ||
              null_data;
  null_data = doubleUserDataNotNull(log_options, usr_col_upper,
                                    "column upper bounds") ||
              null_data;
  return null_data;
}

bool isRowDataNull(const HighsLogOptions& log_options,
                   const double* usr_row_lower, const double* usr_row_upper) {
  bool null_data = false;
  null_data =
      doubleUserDataNotNull(log_options, usr_row_lower, "row lower bounds") ||
      null_data;
  null_data =
      doubleUserDataNotNull(log_options, usr_row_upper, "row upper bounds") ||
      null_data;
  return null_data;
}

bool isMatrixDataNull(const HighsLogOptions& log_options,
                      const HighsInt* usr_matrix_start,
                      const HighsInt* usr_matrix_index,
                      const double* usr_matrix_value) {
  bool null_data = false;
  null_data =
      intUserDataNotNull(log_options, usr_matrix_start, "matrix starts") ||
      null_data;
  null_data =
      intUserDataNotNull(log_options, usr_matrix_index, "matrix indices") ||
      null_data;
  null_data =
      doubleUserDataNotNull(log_options, usr_matrix_value, "matrix values") ||
      null_data;
  return null_data;
}

HighsStatus transformIntoEqualityProblem(const HighsLp& lp,
                                         HighsLp& equality_lp) {
  // Copy lp.
  equality_lp = lp;

  // Add slacks for each row with more than one bound.
  std::vector<double> rhs(lp.numRow_, 0);

  for (HighsInt row = 0; row < lp.numRow_; row++) {
    assert(equality_lp.Astart_[equality_lp.numCol_] ==
           (HighsInt)equality_lp.Avalue_.size());
    assert((HighsInt)equality_lp.Aindex_.size() ==
           (HighsInt)equality_lp.Avalue_.size());
    const HighsInt nnz = equality_lp.Astart_[equality_lp.numCol_];

    if (lp.rowLower_[row] <= -kHighsInf && lp.rowUpper_[row] >= kHighsInf) {
      // free row
      equality_lp.Astart_.push_back(nnz + 1);
      equality_lp.Aindex_.push_back(row);
      equality_lp.Avalue_.push_back(1.0);

      equality_lp.numCol_++;
      equality_lp.colLower_.push_back(-kHighsInf);
      equality_lp.colUpper_.push_back(kHighsInf);
      equality_lp.colCost_.push_back(0);
    } else if (lp.rowLower_[row] > -kHighsInf &&
               lp.rowUpper_[row] >= kHighsInf) {
      // only lower bound
      rhs[row] = lp.rowLower_[row];

      equality_lp.Astart_.push_back(nnz + 1);
      equality_lp.Aindex_.push_back(row);
      equality_lp.Avalue_.push_back(-1.0);

      equality_lp.numCol_++;
      equality_lp.colLower_.push_back(0);
      equality_lp.colUpper_.push_back(kHighsInf);
      equality_lp.colCost_.push_back(0);
    } else if (lp.rowLower_[row] <= -kHighsInf &&
               lp.rowUpper_[row] < kHighsInf) {
      // only upper bound
      rhs[row] = lp.rowUpper_[row];

      equality_lp.Astart_.push_back(nnz + 1);
      equality_lp.Aindex_.push_back(row);
      equality_lp.Avalue_.push_back(1.0);

      equality_lp.numCol_++;
      equality_lp.colLower_.push_back(0);
      equality_lp.colUpper_.push_back(kHighsInf);
      equality_lp.colCost_.push_back(0);
    } else if (lp.rowLower_[row] > -kHighsInf &&
               lp.rowUpper_[row] < kHighsInf &&
               lp.rowLower_[row] != lp.rowUpper_[row]) {
      // both lower and upper bound that are different
      double rhs_value, coefficient;
      double difference = lp.rowUpper_[row] - lp.rowLower_[row];
      if (fabs(lp.rowLower_[row]) < fabs(lp.rowUpper_[row])) {
        rhs_value = lp.rowLower_[row];
        coefficient = -1;
      } else {
        rhs_value = lp.rowUpper_[row];
        coefficient = 1;
      }
      rhs[row] = rhs_value;

      equality_lp.Astart_.push_back(nnz + 1);
      equality_lp.Aindex_.push_back(row);
      equality_lp.Avalue_.push_back(coefficient);

      equality_lp.numCol_++;
      equality_lp.colLower_.push_back(0);
      equality_lp.colUpper_.push_back(difference);
      equality_lp.colCost_.push_back(0);
    } else if (lp.rowLower_[row] == lp.rowUpper_[row]) {
      // equality row
      rhs[row] = lp.rowLower_[row];
    } else {
#ifdef HiGHSDEV
      printf(
          "transformIntoEqualityProblem: Unknown row type when adding slacks");
#endif
      return HighsStatus::kError;
    }
  }
  equality_lp.rowLower_ = rhs;
  equality_lp.rowUpper_ = rhs;
  equality_lp.integrality_.assign(equality_lp.numCol_,
                                  HighsVarType::kContinuous);
  return HighsStatus::kOk;
}

// Given (P) returns (D) for the pair
// (P)
//    min c'x st Ax=b
//     st l <= x <= u
// (D)
//    max b'y + l'zl - u'zu
//     st A'y + zl - zu = c
//        y free, zl >=0, zu >= 0
HighsStatus dualizeEqualityProblem(const HighsLp& lp, HighsLp& dual) {
  std::vector<double> colCost = lp.colCost_;
  if (lp.sense_ != ObjSense::kMinimize) {
    for (HighsInt col = 0; col < lp.numCol_; col++)
      colCost[col] = -colCost[col];
  }

  assert(lp.rowLower_ == lp.rowUpper_);

  const HighsInt ncols = lp.numRow_;
  const HighsInt nrows = lp.numCol_;

  dual.numRow_ = nrows;
  dual.rowLower_ = colCost;
  dual.rowUpper_ = colCost;

  // Add columns (y)
  dual.numCol_ = ncols;
  dual.colLower_.resize(ncols);
  dual.colUpper_.resize(ncols);
  dual.colCost_.resize(ncols);

  for (HighsInt col = 0; col < ncols; col++) {
    dual.colLower_[col] = -kHighsInf;
    dual.colUpper_[col] = kHighsInf;
    // cost b'y
    dual.colCost_[col] = lp.rowLower_[col];
  }

  // Get transpose of A
  HighsInt i, k;
  vector<HighsInt> iwork(lp.numRow_, 0);
  dual.Astart_.resize(lp.numRow_ + 1, 0);
  HighsInt AcountX = lp.Aindex_.size();
  dual.Aindex_.resize(AcountX);
  dual.Avalue_.resize(AcountX);
  for (HighsInt k = 0; k < AcountX; k++) iwork.at(lp.Aindex_.at(k))++;
  for (i = 1; i <= lp.numRow_; i++)
    dual.Astart_.at(i) = dual.Astart_.at(i - 1) + iwork.at(i - 1);
  for (i = 0; i < lp.numRow_; i++) iwork.at(i) = dual.Astart_.at(i);
  for (HighsInt iCol = 0; iCol < lp.numCol_; iCol++) {
    for (k = lp.Astart_.at(iCol); k < lp.Astart_.at(iCol + 1); k++) {
      HighsInt iRow = lp.Aindex_.at(k);
      HighsInt iPut = iwork.at(iRow)++;
      dual.Aindex_.at(iPut) = iCol;
      dual.Avalue_.at(iPut) = lp.Avalue_[k];
    }
  }

  // Add columns (zl)
  for (HighsInt col = 0; col < lp.numCol_; col++) {
    if (lp.colLower_[col] > -kHighsInf) {
      const HighsInt nnz = dual.Astart_[dual.numCol_];

      dual.colLower_.push_back(0);
      dual.colUpper_.push_back(kHighsInf);

      dual.colCost_.push_back(lp.colLower_[col]);

      // Add constaints
      dual.Astart_.push_back(nnz + 1);
      dual.Aindex_.push_back(col);
      dual.Avalue_.push_back(1.0);

      dual.numCol_++;
    }
  }

  // Add columns (zu)
  for (HighsInt col = 0; col < lp.numCol_; col++) {
    if (lp.colUpper_[col] < kHighsInf) {
      const HighsInt nnz = dual.Astart_[dual.numCol_];

      dual.colLower_.push_back(0);
      dual.colUpper_.push_back(kHighsInf);

      dual.colCost_.push_back(-lp.colUpper_[col]);

      // Add constaints
      dual.Astart_.push_back(nnz + 1);
      dual.Aindex_.push_back(col);
      dual.Avalue_.push_back(-1.0);

      dual.numCol_++;
    }
  }

  dual.sense_ = ObjSense::kMinimize;
  for (HighsInt col = 0; col < dual.numCol_; col++) {
    dual.colCost_[col] = -dual.colCost_[col];
  }

  dual.model_name_ = lp.model_name_ + "_dualized";

#ifdef HiGHSDEV
  printf("Dualized equality LP\n");
#endif

  return HighsStatus::kOk;
}

void reportPresolveReductions(const HighsLogOptions& log_options,
                              const HighsLp& lp, const HighsLp& presolve_lp) {
  HighsInt num_col_from = lp.numCol_;
  HighsInt num_row_from = lp.numRow_;
  HighsInt num_els_from = lp.Astart_[num_col_from];
  HighsInt num_col_to = presolve_lp.numCol_;
  HighsInt num_row_to = presolve_lp.numRow_;
  HighsInt num_els_to;
  if (num_col_to) {
    num_els_to = presolve_lp.Astart_[num_col_to];
  } else {
    num_els_to = 0;
  }
  char elemsignchar = '-';
  HighsInt elemdelta = num_els_from - num_els_to;
  if (num_els_from < num_els_to) {
    elemdelta = -elemdelta;
    elemsignchar = '+';
  }
  highsLogUser(
      log_options, HighsLogType::kInfo,
      "Presolve : Reductions: rows %" HIGHSINT_FORMAT "(-%" HIGHSINT_FORMAT
      "); columns %" HIGHSINT_FORMAT "(-%" HIGHSINT_FORMAT
      "); "
      "elements %" HIGHSINT_FORMAT "(%c%" HIGHSINT_FORMAT ")\n",
      num_row_to, (num_row_from - num_row_to), num_col_to,
      (num_col_from - num_col_to), num_els_to, elemsignchar, elemdelta);
}

void reportPresolveReductions(const HighsLogOptions& log_options,
                              const HighsLp& lp, const bool presolve_to_empty) {
  HighsInt num_col_from = lp.numCol_;
  HighsInt num_row_from = lp.numRow_;
  HighsInt num_els_from = lp.Astart_[num_col_from];
  HighsInt num_col_to;
  HighsInt num_row_to;
  HighsInt num_els_to;
  std::string message;
  if (presolve_to_empty) {
    num_col_to = 0;
    num_row_to = 0;
    num_els_to = 0;
    message = "- Reduced to empty";
  } else {
    num_col_to = num_col_from;
    num_row_to = num_row_from;
    num_els_to = num_els_from;
    message = "- Not reduced";
  }
  highsLogUser(log_options, HighsLogType::kInfo,
               "Presolve : Reductions: rows %" HIGHSINT_FORMAT
               "(-%" HIGHSINT_FORMAT "); columns %" HIGHSINT_FORMAT
               "(-%" HIGHSINT_FORMAT
               "); "
               "elements %" HIGHSINT_FORMAT "(-%" HIGHSINT_FORMAT ") %s\n",
               num_row_to, (num_row_from - num_row_to), num_col_to,
               (num_col_from - num_col_to), num_els_to,
               (num_els_from - num_els_to), message.c_str());
}

bool isLessInfeasibleDSECandidate(const HighsLogOptions& log_options,
                                  const HighsLp& lp) {
  HighsInt max_col_num_en = -1;
  const HighsInt max_allowed_col_num_en = 24;
  const HighsInt max_assess_col_num_en =
      std::max(HighsInt{9}, max_allowed_col_num_en);
  const HighsInt max_average_col_num_en = 6;
  vector<HighsInt> col_length_k;
  col_length_k.resize(1 + max_assess_col_num_en, 0);
  bool LiDSE_candidate = true;
  bool all_unit_nonzeros = true;
  for (HighsInt col = 0; col < lp.numCol_; col++) {
    // Check limit on number of entries in the column has not been breached
    HighsInt col_num_en = lp.Astart_[col + 1] - lp.Astart_[col];
    max_col_num_en = std::max(col_num_en, max_col_num_en);
    if (col_num_en > max_assess_col_num_en) {
#ifdef HiGHSDEV
      if (LiDSE_candidate)
        printf("Column %" HIGHSINT_FORMAT " has %" HIGHSINT_FORMAT
               " > %" HIGHSINT_FORMAT " entries so LP is not LiDSE candidate\n",
               col, col_num_en, max_allowed_col_num_en);
      LiDSE_candidate = false;
#else
      LiDSE_candidate = false;
      return LiDSE_candidate;
#endif
    } else {
      col_length_k[col_num_en]++;
    }
    for (HighsInt en = lp.Astart_[col]; en < lp.Astart_[col + 1]; en++) {
      double value = lp.Avalue_[en];
      // All nonzeros must be +1 or -1
      if (fabs(value) != 1) {
        all_unit_nonzeros = false;
#ifdef HiGHSDEV
        if (LiDSE_candidate)
          printf("Column %" HIGHSINT_FORMAT " has entry %" HIGHSINT_FORMAT
                 " with value %g so LP is not LiDSE "
                 "candidate\n",
                 col, en - lp.Astart_[col], value);
        LiDSE_candidate = false;
#else
        LiDSE_candidate = false;
        return LiDSE_candidate;
#endif
      }
    }
  }
#ifdef HiGHSDEV
  /*
  printf("LP has\n");
  HighsInt to_num_en = std::min(max_assess_col_num_en, max_col_num_en);
  for (HighsInt col_num_en = 0; col_num_en < to_num_en+1; col_num_en++)
    printf("%7" HIGHSINT_FORMAT " columns of count %1" HIGHSINT_FORMAT "\n",
  col_length_k[col_num_en], col_num_en);
  */
#endif
  double average_col_num_en = lp.Astart_[lp.numCol_];
  average_col_num_en = average_col_num_en / lp.numCol_;
  LiDSE_candidate =
      LiDSE_candidate && average_col_num_en <= max_average_col_num_en;
  std::string logic0 = "has";
  if (!all_unit_nonzeros) logic0 = "does not have";
  std::string logic1 = "is not";
  if (LiDSE_candidate) logic1 = "is";
  highsLogUser(log_options, HighsLogType::kInfo,
               "LP %s %s all |entries|=1; max column count = %" HIGHSINT_FORMAT
               " (limit %" HIGHSINT_FORMAT
               "); average "
               "column count = %0.2g (limit %" HIGHSINT_FORMAT
               "): So %s a candidate for LiDSE\n",
               lp.model_name_.c_str(), logic0.c_str(), max_col_num_en,
               max_allowed_col_num_en, average_col_num_en,
               max_average_col_num_en, logic1.c_str());
#ifdef HiGHSDEV
  HighsInt int_average_col_num_en = average_col_num_en;
  printf("grep_count_distrib,%s,%" HIGHSINT_FORMAT ",%" HIGHSINT_FORMAT
         ",%" HIGHSINT_FORMAT "\n",
         lp.model_name_.c_str(), max_col_num_en, int_average_col_num_en,
         LiDSE_candidate);
#endif
  return LiDSE_candidate;
}

HighsStatus setOrientation(HighsLp& lp,
                           const MatrixOrientation desired_orientation) {
  if (desired_orientation == MatrixOrientation::kNone)
    return HighsStatus::kError;
  if (lp.orientation_ == desired_orientation) return HighsStatus::kOk;
  if (lp.numCol_ == 0 && lp.numRow_ == 0) {
    // No rows or columns, so either orientation is possible and has
    // identical data: just requires the start of the fictitious
    // row/column 0
    lp.Astart_.assign(1, 0);
    lp.orientation_ = desired_orientation;
  } else {
    // Any LP with positive numbers of rows or columns must have an orientation
    assert(lp.orientation_ != MatrixOrientation::kNone);
    if (desired_orientation == MatrixOrientation::kColwise) {
      ensureColWise(lp);
    } else {
      ensureRowWise(lp);
    }
  }
  assert(lp.orientation_ == desired_orientation);
  return HighsStatus::kOk;
}

void ensureColWise(HighsLp& lp) {
  // Should only call this is orientation is ROWWISE
  assert(lp.orientation_ == MatrixOrientation::kRowwise);
  HighsInt num_nz;
  bool empty_matrix = lp.numCol_ == 0 || lp.numRow_ == 0;
  if (!empty_matrix) {
    // Matrix is probably non-empty
    assert((HighsInt)lp.Astart_.size() >= lp.numRow_ + 1);
    num_nz = lp.Astart_[lp.numRow_];
    assert(num_nz >= 0);
    assert((HighsInt)lp.Aindex_.size() >= num_nz);
    assert((HighsInt)lp.Avalue_.size() >= num_nz);
    empty_matrix = num_nz == 0;
    if (!empty_matrix) {
      // Matrix is non-empty, so transpose it
      vector<HighsInt>& ARstart = lp.Astart_;
      vector<HighsInt>& ARindex = lp.Aindex_;
      vector<double>& ARvalue = lp.Avalue_;
      vector<HighsInt> Astart;
      vector<HighsInt> Aindex;
      vector<double> Avalue;
      Astart.resize(lp.numCol_ + 1);
      Aindex.resize(num_nz);
      Avalue.resize(num_nz);
      vector<HighsInt> Alength;
      Alength.assign(lp.numCol_, 0);
      for (HighsInt iEl = ARstart[0]; iEl < num_nz; iEl++)
        Alength[ARindex[iEl]]++;
      Astart[0] = 0;
      for (HighsInt iCol = 0; iCol < lp.numCol_; iCol++)
        Astart[iCol + 1] = Astart[iCol] + Alength[iCol];
      for (HighsInt iRow = 0; iRow < lp.numRow_; iRow++) {
        for (HighsInt iEl = ARstart[iRow]; iEl < ARstart[iRow + 1]; iEl++) {
          HighsInt iCol = ARindex[iEl];
          HighsInt iCol_el = Astart[iCol];
          Aindex[iCol_el] = iRow;
          Avalue[iCol_el] = ARvalue[iEl];
          Astart[iCol]++;
        }
      }
      Astart[0] = 0;
      for (HighsInt iCol = 0; iCol < lp.numCol_; iCol++)
        Astart[iCol + 1] = Astart[iCol] + Alength[iCol];
      assert(Astart[lp.numCol_] == num_nz);
      // Now update the LP's matrix
      lp.Astart_ = Astart;
      lp.Aindex_ = Aindex;
      lp.Avalue_ = Avalue;
    }
  }
  if (empty_matrix) {
    // Matrix is empty, so set up empty column-wise structure
    lp.Astart_.assign(lp.numCol_ + 1, 0);
    lp.Aindex_.clear();
    lp.Avalue_.clear();
  }
  assert((HighsInt)lp.Astart_.size() >= lp.numCol_ + 1);
  num_nz = lp.Astart_[lp.numCol_];
  assert(num_nz >= 0);
  assert((HighsInt)lp.Aindex_.size() >= num_nz);
  assert((HighsInt)lp.Avalue_.size() >= num_nz);
  lp.orientation_ = MatrixOrientation::kColwise;
}

void ensureRowWise(HighsLp& lp) {
  // Should only call this is orientation is COLWISE
  assert(lp.orientation_ == MatrixOrientation::kColwise);
  HighsInt num_nz;
  bool empty_matrix = lp.numCol_ == 0 || lp.numRow_ == 0;
  if (!empty_matrix) {
    // Matrix is probably non-empty
    assert((HighsInt)lp.Astart_.size() >= lp.numCol_ + 1);
    num_nz = lp.Astart_[lp.numCol_];
    assert(num_nz >= 0);
    assert((HighsInt)lp.Aindex_.size() >= num_nz);
    assert((HighsInt)lp.Avalue_.size() >= num_nz);
    empty_matrix = num_nz == 0;
    if (!empty_matrix) {
      // Matrix is non-empty, so transpose it
      vector<HighsInt>& Astart = lp.Astart_;
      vector<HighsInt>& Aindex = lp.Aindex_;
      vector<double>& Avalue = lp.Avalue_;
      vector<HighsInt> ARstart;
      vector<HighsInt> ARindex;
      vector<double> ARvalue;
      ARstart.resize(lp.numRow_ + 1);
      ARindex.resize(num_nz);
      ARvalue.resize(num_nz);
      vector<HighsInt> ARlength;
      ARlength.assign(lp.numRow_, 0);
      for (HighsInt iEl = Astart[0]; iEl < num_nz; iEl++)
        ARlength[Aindex[iEl]]++;
      ARstart[0] = 0;
      for (HighsInt iRow = 0; iRow < lp.numRow_; iRow++)
        ARstart[iRow + 1] = ARstart[iRow] + ARlength[iRow];
      for (HighsInt iCol = 0; iCol < lp.numCol_; iCol++) {
        for (HighsInt iEl = Astart[iCol]; iEl < Astart[iCol + 1]; iEl++) {
          HighsInt iRow = Aindex[iEl];
          HighsInt iRow_el = ARstart[iRow];
          ARindex[iRow_el] = iCol;
          ARvalue[iRow_el] = Avalue[iEl];
          ARstart[iRow]++;
        }
      }
      ARstart[0] = 0;
      for (HighsInt iRow = 0; iRow < lp.numRow_; iRow++)
        ARstart[iRow + 1] = ARstart[iRow] + ARlength[iRow];
      assert(ARstart[lp.numRow_] == num_nz);
      // Now update the LP's matrix
      lp.Astart_ = ARstart;
      lp.Aindex_ = ARindex;
      lp.Avalue_ = ARvalue;
    }
  }
  if (empty_matrix) {
    // Matrix is empty, so set up empty row-wise structure
    lp.Astart_.assign(lp.numRow_ + 1, 0);
    lp.Aindex_.clear();
    lp.Avalue_.clear();
  }
  assert((HighsInt)lp.Astart_.size() >= lp.numRow_ + 1);
  num_nz = lp.Astart_[lp.numRow_];
  assert(num_nz >= 0);
  assert((HighsInt)lp.Aindex_.size() >= num_nz);
  assert((HighsInt)lp.Avalue_.size() >= num_nz);
  lp.orientation_ = MatrixOrientation::kRowwise;
}
