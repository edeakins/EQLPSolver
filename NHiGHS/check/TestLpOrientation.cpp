#include "Avgas.h"
#include "Highs.h"
#include "catch.hpp"
#include "lp_data/HighsLpUtils.h"
#include "util/HighsUtils.h"

const bool dev_run = false;

// No commas in test case name.
TEST_CASE("LP-orientation", "[lp_orientation]") {
  Avgas avgas;
  const HighsInt avgas_num_col = 8;
  const HighsInt avgas_num_row = 10;
  HighsInt num_row = 0;
  HighsInt num_row_nz = 0;
  vector<double> rowLower;
  vector<double> rowUpper;
  vector<HighsInt> ARstart;
  vector<HighsInt> ARindex;
  vector<double> ARvalue;

  for (HighsInt row = 0; row < avgas_num_row; row++)
    avgas.row(row, num_row, num_row_nz, rowLower, rowUpper, ARstart, ARindex,
              ARvalue);

  ARstart.resize(num_row + 1);
  ARstart[num_row] = num_row_nz;

  HighsInt num_col = 0;
  HighsInt num_col_nz = 0;
  vector<double> colCost;
  vector<double> colLower;
  vector<double> colUpper;
  vector<HighsInt> Astart;
  vector<HighsInt> Aindex;
  vector<double> Avalue;
  for (HighsInt col = 0; col < avgas_num_col; col++)
    avgas.col(col, num_col, num_col_nz, colCost, colLower, colUpper, Astart,
              Aindex, Avalue);
  Astart.resize(num_col + 1);
  Astart[num_col] = num_col_nz;
  assert(num_col_nz == num_row_nz);

  double optimal_objective_function_value = -7.75;
  Highs highs;
  if (!dev_run) {
    highs.setOptionValue("output_flag", false);
  } else {
    highs.setOptionValue("log_dev_level", kHighsLogDevLevelVerbose);
  }
  const HighsLp& highs_lp = highs.getLp();
  const HighsInfo& info = highs.getInfo();

  REQUIRE(highs_lp.orientation_ == MatrixOrientation::kNone);

  // Set up the LP externally
  HighsLp lp;
  lp.numCol_ = num_col;
  lp.numRow_ = num_row;
  lp.colCost_ = colCost;
  lp.colLower_ = colLower;
  lp.colUpper_ = colUpper;
  lp.rowLower_ = rowLower;
  lp.rowUpper_ = rowUpper;
  lp.Astart_ = Astart;
  lp.Aindex_ = Aindex;
  lp.Avalue_ = Avalue;
  lp.orientation_ = MatrixOrientation::kColwise;
  highs.passModel(lp);
  highs.run();
  REQUIRE(info.objective_function_value == optimal_objective_function_value);

  // Make the external LP row-wise then pass and solve it
  REQUIRE(setOrientation(lp, MatrixOrientation::kRowwise) == HighsStatus::kOk);
  highs.passModel(lp);
  highs.run();
  REQUIRE(info.objective_function_value == optimal_objective_function_value);

  // Make the external LP col-wise then pass and solve it
  REQUIRE(setOrientation(lp) == HighsStatus::kOk);
  highs.passModel(lp);
  highs.run();
  REQUIRE(info.objective_function_value == optimal_objective_function_value);

  // Clear the internal LP
  highs.clearModel();
  REQUIRE(highs.addCols(num_col, &colCost[0], &colLower[0], &colUpper[0], 0,
                        NULL, NULL, NULL) == HighsStatus::kOk);
  REQUIRE(highs.addRows(num_row, &rowLower[0], &rowUpper[0], num_row_nz,
                        &ARstart[0], &ARindex[0],
                        &ARvalue[0]) == HighsStatus::kOk);
  highs.run();
  REQUIRE(info.objective_function_value == optimal_objective_function_value);

  // Clear the internal LP
  highs.clearModel();
  highs.addCols(num_col, &colCost[0], &colLower[0], &colUpper[0], 0, NULL, NULL,
                NULL);
  vector<double> one_row_Lower;
  vector<double> one_row_Upper;
  vector<HighsInt> one_row_start;
  vector<HighsInt> one_row_index;
  vector<double> one_row_value;
  for (HighsInt row = 0; row < avgas_num_row; row++) {
    HighsInt one_row_numnz = 0;
    HighsInt one_row_numrow = 0;
    avgas.row(row, one_row_numrow, one_row_numnz, one_row_Lower, one_row_Upper,
              one_row_start, one_row_index, one_row_value);
    REQUIRE(highs.addRows(1, &one_row_Lower[0], &one_row_Upper[0],
                          one_row_numnz, &one_row_start[0], &one_row_index[0],
                          &one_row_value[0]) == HighsStatus::kOk);
  }
  highs.run();
  REQUIRE(info.objective_function_value == optimal_objective_function_value);

  //  return_status = highs.writeModel("");
}
