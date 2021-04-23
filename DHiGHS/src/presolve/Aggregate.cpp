#include "Aggregate.h"
using namespace std;

HighsAggregate::HighsAggregate(HighsLp& lp, const struct eq_part& ep, HighsSolution& solution, HighsBasis& basis){
	// From the original lp
  iter = 0;
	numRow = lp.numRow_;
	numCol = lp.numCol_;
	numTot = numRow + numCol;
  nnz = lp.nnz_;
	rowLower.assign(lp.rowLower_.begin(), lp.rowLower_.end());
	rowUpper.assign(lp.rowUpper_.begin(), lp.rowUpper_.end());
  // for (int i = 0; i < lp.colLower_.size(); ++i){
  //   lp.colUpper_[i] = +HIGHS_CONST_INF;
  // }
	colUpper.assign(lp.colUpper_.begin(), lp.colUpper_.end());
	colLower.assign(lp.colLower_.begin(), lp.colLower_.end());
	colCost.assign(lp.colCost_.begin(), lp.colCost_.end());
	Avalue.assign(lp.Avalue_.begin(), lp.Avalue_.end());
	Aindex.assign(lp.Aindex_.begin(), lp.Aindex_.end());
	Astart.assign(lp.Astart_.begin(), lp.Astart_.end());
	// //Equitable partition info
	partition = ep;
  cell.resize(numTot);
  cellFront.resize(numTot);
  cellSize.resize(numTot);
  labels.resize(numTot);
  previousCell.resize(numTot);
  previousCellSize.resize(numTot);
  previousCellFront.resize(numTot);
  previousLabels.resize(numTot);
  colsToReps.assign(numCol, -1);
  prevColsToReps.assign(numCol, -1);
  repsToCols.assign(numCol, -1);
  prevRepsToCols.assign(numCol, -1);
  rowsToReps.assign(numRow, -1);
  prevRowsToReps.assign(numRow, -1);
  repsToRows.assign(numRow, -1);
  prevRepsToRows.assign(numRow, -1);
  lastSolveRow.resize(numRow);
  lastSolveCol.resize(numCol);
  row.resize(numRow);
  prevCol.resize(numCol);
  col.resize(numCol);
  cellMarked.assign(numCol, false);
  linked.resize(numCol);
	// Previous solution
	col_value = (solution.col_value);
	row_value = (solution.row_value);
	// Previous basis
	col_status = (basis.col_status);
	row_status = (basis.row_status);
  nonBasicCol.resize(numCol);
  nonBasicRow.resize(numRow);
  // New Lp info
  maxLinkCols = ep.nsplits;
  maxLinkSpace = 3 * maxLinkCols;
  alp = (HighsLp *)calloc(1, sizeof(HighsLp));
  alpBasis = (HighsBasis *)calloc(1, sizeof(HighsBasis));
  alp->colCost_.resize(numCol + maxLinkCols);
  alp->Avalue_.resize(nnz + maxLinkSpace);
  alp->Aindex_.resize(nnz + maxLinkSpace);
  alp->Astart_.resize(numCol + maxLinkCols + 1);
  alp->colUpper_.resize(numCol + maxLinkCols);
  alp->colLower_.resize(numCol + maxLinkCols);
  alp->rowLower_.resize(numRow + maxLinkCols);
  alp->rowUpper_.resize(numRow + maxLinkCols);
  alp->linkers.resize(maxLinkCols);
  alp->linkLower_.resize(maxLinkCols);
  alp->linkUpper_.resize(maxLinkCols);
  rowRepsToFix.assign(numRow, false);
  rowRepsValue.resize(numRow);
  colRepsToFix.assign(numRow, false);
  colRepsValue.resize(numRow);
  rowRepsScale.resize(numRow);
  alpBasis->col_status.resize(numCol + maxLinkCols);
  alpBasis->row_status.resize(numRow + maxLinkCols);
  parent.resize(maxLinkCols);
  child.resize(maxLinkCols);
  coeff.assign(numRow, 0);
  linkARstart.resize(maxLinkCols + 1);
  linkARindex.resize(maxLinkSpace);
  linkARvalue.resize(maxLinkSpace);
  linkAlength.assign(numCol + maxLinkCols, 0);
  // linkLB.assign(maxLinkCols, 0);
  // linkUB.assign(maxLinkCols, 0);
  // coeff.assign(numTot);
  AindexPacked_.resize(nnz);
  // Translate fronts array to colors for vertices
  translateFrontsToColors();
  packVectors();
  foldObj();
  foldMatrix();
  fixMatrix();
  foldRhsInit();
  foldBndsInit();
}

bool HighsAggregate::update(const struct eq_part& ep, const HighsSolution& solution, const HighsBasis& basis){
  // equtiable partition update
  // previousPartition = partition;
  partition = ep;
  // Previous solution
	col_value = (solution.col_value);
	row_value = (solution.row_value);
	// Previous basis
	col_status = (basis.col_status);
	row_status = (basis.row_status);
  reset();
  saveRowsAndColsFromLastSolve();
  translateFrontsToColors();
  identifyLinks();
  if (!numLinkers_) return false;
  packVectors();
  foldObj();
  foldMatrix();
  fixMatrix();
  if (!iter) {findRowRepsToFix(); findColRepsToFix();}
  foldRhs();
  foldBnds();
  setRowBasis();
  setColBasis();
  createLinkRows();
  addCols();
  addRows();
  ++iter;
  return true;
}

void HighsAggregate::translateFrontsToColors(){
  int i = 0, c = 0, r = 0, rep, cf, colCounter = 0, rowCounter = 0;
  map<int, int> fCell;
  // map fronts to colors and count numCol_ and numRow_
  previousNumCol_ = numCol_;
  numCol_ = numRow_ = numTot_ = 0;
  prevCol = col;
  int numColCell = 0;
  for (i = 0; i < numTot; ++i){
    if (fCell.insert(pair<int, int>(partition.fronts[i], c)).second){
      ++c;
      partition.fronts[i] < numCol ? ++numCol_ : ++numRow_;
      ++numTot_;
    }
  }
  for (i = 0; i < numTot; ++i){
    labels[i] = partition.labels[i];
    cell[i] = fCell.find(partition.fronts[i])->second;
    // if (cell[i] == 0)
    //   std::cout << "color of node i: " << i << " = " << cell[i] << std::endl; 
    cellFront[fCell.find(partition.fronts[i])->second] = partition.fronts[i];
    ++cellSize[fCell.find(partition.fronts[i])->second];
  }
  alp->numCol_ = numCol_;
  alp->numRow_ = numRow_;
  for (i = 0; i < numCol_; ++i){
    cf = cellFront[i];
    rep = labels[cf];
    repsToCols[rep] = colCounter;
    colsToReps[colCounter++] = rep;
  }
  for (i = 0; i < numRow_; ++i){
    cf = cellFront[i + numCol_];
    rep = labels[cf];
    repsToRows[rep - numCol] = rowCounter;
    rowsToReps[rowCounter++] = rep - numCol;
  }
  for (i = 0; i < numCol; ++i){
    cf = cellFront[cell[i]];
    rep = labels[cf];
    c = repsToCols[rep];
    col[i] = c;
  }
  for (i = numCol; i < numTot; ++i){
    cf = cellFront[cell[i]];
    rep = labels[cf];
    r = repsToRows[rep - numCol];
    row[i - numCol] = r;
  }
}

void HighsAggregate::packVectors(){
  int i, r;
  for (i = 0; i < nnz; ++i){
    r = row[Aindex[i]];
    AindexPacked_[i] = r;
  }
}

// Find the rows to fix from the original aggregate
void HighsAggregate::findRowRepsToFix(){
  int oldNCol = col_value.size();
  vector<bool> fixRow(numRow, false);
  for (int i = 0; i < row_value.size(); ++i){
    int rep = previousLabels[previousCellFront[i + oldNCol]];
    if (std::fabs(row_value[i] - rowLower[rep - numCol] * previousCellSize[i + oldNCol]) < 1e-6 ||
        std::fabs(row_value[i] - rowUpper[rep - numCol] * previousCellSize[i + oldNCol]) < 1e-6){
          fixRow[i] = true;
    }
  }
  for (int i = numCol; i < numTot; ++i){
    if (fixRow[previousCell[i] - oldNCol]){
      rowRepsToFix[i - numCol] = true;
      rowRepsValue[i - numCol] = row_value[previousCell[i] - oldNCol];
      rowRepsScale[i - numCol] = previousCellSize[previousCell[i]];
    }
  }
}

// Find the col reps to be fixed from original aggregation
void HighsAggregate::findColRepsToFix(){
  vector<bool> fixCol(numCol, false);
  for (int i = 0; i < col_value.size(); ++i){
    int rep = previousLabels[previousCellFront[i]];
    if (std::fabs(col_value[i] - colLower[rep]) < 1e-6 ||
        std::fabs(col_value[i] - colUpper[rep]) < 1e-6){
          fixCol[i] = true;
    }
  }
  for (int i = 0; i < numCol; ++i){
    if (fixCol[previousCell[i]]){
      colRepsToFix[i] = true;
      colRepsValue[i] = col_value[previousCell[i]];
    }
  }
}

// Fold the objective
void HighsAggregate::foldObj(){
  int i, rep;
  for (i = 0; i < numCol_; ++i){
    rep = partition.labels[cellFront[i]];
    alp->colCost_[i] = colCost[rep] * cellSize[i];
  }
  alp->sense_ = 1;
}

// Fold the matrix down based on current ep
void HighsAggregate::foldMatrix(){
  int i, j, r, c, rep, idx = 0;
  alp->Astart_[0] = 0;
  for (i = 0; i < numCol_; ++i){
    rep = colsToReps[i];
    for (j = Astart[rep]; j < Astart[rep + 1]; ++j){
      r = AindexPacked_[j];
      coeff[r] += Avalue[j];
    }
    for (j = 0; j < numRow_; ++j){
      if (coeff[j]){
        c = cell[rep];
        alp->Avalue_[alp->nnz_] = coeff[j] * cellSize[c];
        alp->Aindex_[alp->nnz_++] = j;
        coeff[j] = 0;
      }
    }
    alp->Astart_[i + 1] = alp->nnz_;
  }
}

void HighsAggregate::fixMatrix(){
  int i;
  for (i = alp->numCol_; i < numCol + maxLinkCols; ++i)
    alp->Astart_[i + 1] = alp->Astart_[i];
}

void HighsAggregate::saveRowsAndColsFromLastSolve(){
  lastSolveCol = col;
  lastSolveRow = row;
  prevRepsToCols = repsToCols;
  prevColsToReps = colsToReps;
  prevRepsToRows = repsToRows;
  prevRowsToReps = rowsToReps;
}

void HighsAggregate::reset(){
  int i;
  alp->nnz_ = 0;
  alp->pivotTime = 0;
  alp->invertTime = 0;
  alp->unfoldIter = 0;
  previousNumCol_ = alp->numCol_ - numLinkers_;
  for (i = 0; i < numTot; ++i){
    previousCell[i] = cell[i];
    previousCellSize[i] = cellSize[i];
    previousCellFront[i] = cellFront[i];
    previousLabels[i] = labels[i];
    cellSize[i] = 0;
  }
  for (i = 0; i < numCol; ++i)
    linked[i] = 0;
  for (i = 0; i < alp->nnz_; ++i)
    alp->Avalue_[i] = 0;
  for (i = 0; i < numRow; ++i)  
    nonBasicRow[i] = false;
  for (i = 0; i < numCol; ++i)
    nonBasicCol[i] = false;
  for (i = 0; i < maxLinkSpace; ++i){
    linkARindex[i] = 0;
    linkARvalue[i] = 0;
  }
  for (i = 0; i < maxLinkCols + 1; ++i)
    linkARstart[i] = 0;
  for (i = 0; i < maxLinkCols; ++i){
    parent[i] = -1;
    child[i] = -1;
  }
  numLinkers_ = 0;
}

// Fold the row bounds based on current ep
void HighsAggregate::foldRhsInit(){
  int i, rep;
  for (i = 0; i < numRow_; ++i){
    rep = partition.labels[cellFront[i + numCol_]] - numCol;
    alp->rowLower_[i] = std::max(rowLower[rep] * cellSize[i + numCol_], -HIGHS_CONST_INF);
    alp->rowUpper_[i] = std::min(rowUpper[rep] * cellSize[i + numCol_], HIGHS_CONST_INF);
  }
}

void HighsAggregate::foldRhs(){
  int i, rep, rowIdx, cellIdx;
  for (i = 0; i < numRow_; ++i){
    rep = partition.labels[cellFront[i + numCol_]];
    rowIdx = previousCell[rep] - previousNumCol_;
    bool fixed = rowRepsToFix[rep - numCol];
    cellIdx = previousCell[rep];
    if (rowRepsToFix[rep - numCol]){
      double rhs = (double)rowRepsValue[rep - numCol]/rowRepsScale[rep - numCol];
      alp->rowLower_[i] = alp->rowUpper_[i] = rhs * cellSize[i + numCol_];
    }
    // if (fabs(row_value[rowIdx] - (rowLower[rep - numCol] * previousCellSize[cellIdx])) < 1e-6){
    //   alp->rowLower_[i] = rowLower[rep - numCol] * cellSize[i + numCol_];
    //   alp->rowUpper_[i] = rowLower[rep - numCol] * cellSize[i + numCol_];  
    // }
    // else if (fabs(row_value[rowIdx] - (rowUpper[rep - numCol] * previousCellSize[cellIdx])) < 1e-6){
    //   alp->rowLower_[i] = rowUpper[rep - numCol] * cellSize[i + numCol_];
    //   alp->rowUpper_[i] = rowUpper[rep - numCol] * cellSize[i + numCol_];  
    // }
    else{
      alp->rowLower_[i] = std::max(rowLower[rep - numCol] * cellSize[i + numCol_], -HIGHS_CONST_INF);
      alp->rowUpper_[i] = std::min(rowUpper[rep - numCol] * cellSize[i + numCol_], HIGHS_CONST_INF);  
    }
  }
}

// Set row basis
void HighsAggregate::setRowBasis(){
  int i, rep, rowIdx;
  alpBasis->numRow_ = 0;
  for (i = 0; i < numRow_; ++i){
    rep = partition.labels[cellFront[i + numCol_]];
    rowIdx = previousCell[rep] - previousNumCol_;
    if ((row_status[rowIdx] == HighsBasisStatus::UPPER ||
        row_status[rowIdx] == HighsBasisStatus::LOWER ||
        row_status[rowIdx] == HighsBasisStatus::NONBASIC) &&
        !nonBasicRow[rowIdx]){
      alpBasis->row_status[i] = row_status[rowIdx];
      nonBasicRow[rowIdx] = true;
      ++alpBasis->numRow_;
    }
    else if ((row_status[rowIdx] == HighsBasisStatus::UPPER ||
            row_status[rowIdx] == HighsBasisStatus::LOWER ||
            row_status[rowIdx] == HighsBasisStatus::NONBASIC) &&
            nonBasicRow[rowIdx]){
      alpBasis->row_status[i] = HighsBasisStatus::BASIC;
      ++alpBasis->numRow_;
    }
    else{
      alpBasis->row_status[i] = row_status[rowIdx];
      ++alpBasis->numRow_;
    }
  }
}

// Fold the row bounds based on current ep
void HighsAggregate::foldBndsInit(){
  int i, rep;
  for (i = 0; i < numCol_; ++i){
    rep = partition.labels[cellFront[i]];
    alp->colLower_[i] = colLower[rep];
    alp->colUpper_[i] = colUpper[rep];
  }
}

void HighsAggregate::foldBnds(){
  int i, rep, colIdx;
  for (i = 0; i < numCol_; ++i){
    rep = partition.labels[cellFront[i]];
    colIdx = previousCell[rep];
    bool fixed = colRepsToFix[rep];
    if (colRepsToFix[rep]){
      double bnd = (double)colRepsValue[rep];
      alp->colLower_[i] = alp->colUpper_[i] = bnd;
    }
    // if (fabs(col_value[colIdx] - (colLower[rep])) < 1e-6){
    //   alp->colLower_[i] = colLower[rep];
    //   alp->colUpper_[i] = colLower[rep];
    // }
    // else if (fabs(col_value[colIdx] - (colUpper[rep])) < 1e-6){
    //   alp->colLower_[i] = colUpper[rep];
    //   alp->colUpper_[i] = colUpper[rep];
    // }
    else{
      alp->colLower_[i] = colLower[rep];
      alp->colUpper_[i] = colUpper[rep];
    }
  }
}

// Set col basis
void HighsAggregate::setColBasis(){
  int i, rep, colIdx;
  alpBasis->numCol_ = 0;
  for (i = 0; i < numCol_; ++i){
    rep = partition.labels[cellFront[i]];
    colIdx = previousCell[rep];
    if ((col_status[colIdx] == HighsBasisStatus::UPPER ||
        col_status[colIdx] == HighsBasisStatus::LOWER ||
        col_status[colIdx] == HighsBasisStatus::NONBASIC) &&
        !nonBasicCol[colIdx]){
      alpBasis->col_status[i] = col_status[colIdx];
      nonBasicCol[colIdx] = true;
      ++alpBasis->numCol_;
    }
    else if ((col_status[colIdx] == HighsBasisStatus::UPPER ||
            col_status[colIdx] == HighsBasisStatus::LOWER ||
            col_status[colIdx] == HighsBasisStatus::NONBASIC) &&
            nonBasicCol[colIdx]){
      alpBasis->col_status[i] = HighsBasisStatus::BASIC;
      ++alpBasis->numCol_;
    }
    else{
      alpBasis->col_status[i] = col_status[colIdx];
      ++alpBasis->numCol_;
    }
  }
}

void HighsAggregate::addRows(){
  int i, stop = numLinkers_;
  for (i = 0; i < maxLinkCols; ++i){
    if (i < stop){
      alp->rowLower_[alp->numRow_] = 0;
      alp->rowUpper_[alp->numRow_] = 0;
      alpBasis->row_status[alp->numRow_++] = HighsBasisStatus::LOWER;
      ++alpBasis->numRow_;
    }
  }
  alp->numRow_ -= stop;
  appendRowsToMatrix();
}

void HighsAggregate::appendRowsToMatrix(){
  int num_new_row = numLinkers_;
  int num_new_nz = numLinkers_ * 3;
  int new_num_nz = alp->nnz_ + num_new_nz;
  for (int el = 0; el < num_new_nz; el++) linkAlength[linkARindex[el]]++;
  // Append the new rows
  // Shift the existing columns to make space for the new entries
  int new_el = new_num_nz;
  for (int col = alp->numCol_ - 1; col >= 0; col--) {
    int start_col_plus_1 = new_el;
    new_el -= linkAlength[col];
    for (int el = alp->Astart_[col + 1] - 1; el >= alp->Astart_[col]; el--) {
      new_el--;
      alp->Aindex_[new_el] = alp->Aindex_[el];
      alp->Avalue_[new_el] = alp->Avalue_[el];
    }
    alp->Astart_[col + 1] = start_col_plus_1;
  }
  for (int i = alp->numCol_; i < numCol + maxLinkCols; ++i)
    alp->Astart_[i + 1] = alp->Astart_[i];
  assert(new_el == 0);
  // Insert the new entries
  for (int row = 0; row < num_new_row; row++) {
    int first_el = linkARstart[row];
    int last_el = (row < num_new_row - 1 ? linkARstart[row + 1] : num_new_nz);
    for (int el = first_el; el < last_el; el++) {
      int col = linkARindex[el];
      new_el = alp->Astart_[col + 1] - linkAlength[col];
      linkAlength[col]--;
      alp->Aindex_[new_el] = alp->numRow_ + row;
      alp->Avalue_[new_el] = linkARvalue[el];
    }
  }
  alp->numRow_ += num_new_row;
}

void HighsAggregate::addCols(){
  int i, stop = numLinkers_;
  for (i = 0; i < maxLinkCols; ++i){
    if (i < stop){
      alpBasis->col_status[alp->numCol_] = HighsBasisStatus::LOWER;
      alp->colLower_[alp->numCol_] = 0;
      alp->colUpper_[alp->numCol_] = 0;
      alp->linkers[i] = alp->numCol_;
      ++alp->numCol_;
      ++alpBasis->numCol_;
    }
  }
  alp->numLinkers = numLinkers_;
}

void HighsAggregate::identifyLinks(){
  int i, idx = 0;
  for (i = 0; i < numCol; ++i){
    if (partition.parents[i] >= 0){
      int x1 = cell[partition.labels[partition.parents[i]]];
      int x2 = cell[partition.labels[i]];
      parent[numLinkers_] = x1; child[numLinkers_] = x2;
      int rep1 = partition.labels[cellFront[x1]];
      int rep2 = partition.labels[cellFront[x2]];
      ++numLinkers_;
      // alp->linkLower_[numLinkers_] = colLower[rep1] - colUpper[rep2];
      // alp->linkUpper_[numLinkers_++] = colUpper[rep1] - colLower[rep2]; 
    }
  }
  // int lCnt = 0;
  // for (i = 0; i < numCol_; ++i){
  //   if (previousCellSize[i] != cellSize[i] && !linked[cell[i]]
  //       && cell[i] != previousCell[i]){
  //     // parent[numLinkers_] = previousCell[i];
  //     // child[numLinkers_++] = cell[i];
  //     lCnt++;
  //     linked[cell[i]] = true;
  //   }
  // }
  // std::cout << "lCnt: " << lCnt << std::endl;
  // std::cin.get();
}

void HighsAggregate::createLinkRows(){
  int i, offset = 0, idx, stop = numLinkers_;
  linkARstart[0] = 0;
  for (i = 0; i < maxLinkCols; ++i){
    if (i < stop){
      linkARindex[offset] = parent[i]; linkARvalue[offset] = 1;
      linkARindex[offset + 1] = child[i]; linkARvalue[offset + 1] = -1;
      linkARindex[offset + 2] = alp->numCol_ + i; linkARvalue[offset + 2] = -1;
      linkARstart[i + 1] = linkARstart[i] + 3;
      offset += 3;
    }
  }
  for (i = stop; i < maxLinkCols; ++i)
    linkARstart[i + 1] = linkARstart[i];
}

HighsLp* HighsAggregate::getAlp(){
  return alp;
}

HighsBasis* HighsAggregate::getBasis(){
  return alpBasis;
}

// HighsBasis& HighsAggregate::getAlpBasis(){
//   alpBasis.col_status = col_status_;
//   alpBasis.row_status = row_status_;
// 	return alpBasis;
// }

// HighsLp& HighsAggregate::getAlp(){
// 	alp.numRow_ = (numRow_);
// 	alp.numCol_ = (numCol_);
// 	alp.numInt_ = 0;
// 	alp.nnz_ = Avalue_.size();
// 	alp.linkers.assign(linkers.begin(), linkers.end());
// 	//alp.artificialVariables.assign(artificialVariables.begin(), artificialVariables.end());
// 	alp.Astart_.assign(Astart_.begin(), Astart_.end());
// 	alp.Aindex_.assign(Aindex_.begin(), Aindex_.end());
// 	alp.Avalue_.assign(Avalue_.begin(), Avalue_.end());
// 	alp.colUpper_.assign(colUpper_.begin(), colUpper_.end());
// 	alp.colLower_.assign(colLower_.begin(), colLower_.end());
// 	alp.rowUpper_.assign(rowUpper_.begin(), rowUpper_.end());
// 	alp.rowLower_.assign(rowLower_.begin(), rowLower_.end());
// 	alp.colCost_.assign(colCost_.begin(), colCost_.end());
// 	alp.model_name_ = (model_name_);
// 	alp.lp_name_ = (lp_name_);
// 	alp.sense_ = 1;
// 	alp.offset_ = 0;
// 	return alp;
// }

// void HighsAggregate::aggregate(){
//   aggregateColBounds();
//   aggregateRowBounds();
//   aggregateCVector();
//   aggregateAMatrix();
//   setColBasis();
//   setRowBasis();
//   appendLinkersToLp();
// }

// void HighsAggregate::aggregateColBounds(){
//   colUpper_.assign(numCol_, 0);
//   colLower_.assign(numCol_, 0);
//   if (col_value.size()){
//     for (int i = 0; i < numCol_; ++i){
//       int rep = C[i].front();
//       int pCol = previousColumnColoring[rep];
//       if (fabs(col_value[pCol] - colUpper[rep]) < 1e-6){
//         colUpper_[i] = colUpper[rep];
//         colLower_[i] = colUpper[rep];
//       }
//       else if (fabs(col_value[pCol] - colLower[rep]) < 1e-6){
//         colUpper_[i] = colLower[rep];
//         colLower_[i] = colLower[rep];
//       }
//       else{
//         colLower_[i] = colLower[rep];
//         colUpper_[i] = colUpper[rep];
//       }
// 	  }
//   }
//   else{
//     for (int i = 0; i < numCol_; ++i){
//       int rep = C[i].front();
//       colLower_[i] = colLower[rep];
//       colUpper_[i] = colUpper[rep];
//     }
//   }
// }

// void HighsAggregate::aggregateRowBounds(){
//   rowUpper_.assign(numRow_, 0);
//   rowLower_.assign(numRow_, 0);
//   if (row_value.size()){
//     for (int i = 0; i < numRow_; ++i){
//       int rep = C[i + numCol].front() - numCol;
//       int pCol = previousRowColoring[rep] - numCol;
//       if (fabs(row_value[pCol] - (rowLower[rep] * prevC[pCol + numCol].size())) < 1e-6){
//         rowLower_[i] = rowLower[rep] * C[i + numCol].size();
//         rowUpper_[i] = rowLower[rep] * C[i + numCol].size();
//       }
//       else if (fabs(row_value[pCol] - (rowUpper[rep] * prevC[pCol + numCol].size())) < 1e-6){
//         rowLower_[i] = rowUpper[rep] * C[i + numCol].size();
//         rowUpper_[i] = rowUpper[rep] * C[i + numCol].size();
//       }
//       else{
//         rowLower_[i] = rowLower[rep] * C[i + numCol].size();
//         rowUpper_[i] = rowUpper[rep] * C[i + numCol].size();
//       }
//     }
//   }
//   else{
//     for (int i = 0; i < numRow_; ++i){
//       int rep = C[i + numCol].front() - numCol;
//       rowLower_[i] = rowLower[rep] * C[i + numCol].size();
//       rowUpper_[i] = rowUpper[rep] * C[i + numCol].size();
//     }
//   }
// }

// void HighsAggregate::aggregateCVector(){
//   colCost_.assign(numCol_, 0);
// 	for (int i = 0; i < numCol_; ++i){
//     int rep = C[i].front();
//     colCost_[i] = C[i].size() * colCost[rep];
//   }
// }

// void HighsAggregate::aggregateAMatrix(){
// 	Astart_.push_back(0);
// 	for (int i = 0; i < numCol_; ++i){
//     vector<double> coeff(numRow_, 0);
// 		int rep = C[i].front();
// 		for (int j = Astart[rep]; j < Astart[rep + 1]; ++j){
//       int rowIdx = AindexSub[j] - numCol;
//       coeff[rowIdx] += Avalue[j];
//     }
//     for (int j = 0; j < coeff.size(); ++j){
//       if (coeff[j]){
//         Avalue_.push_back(coeff[j] * C[i].size());
//         Aindex_.push_back(j);
//       }
//     }
//     Astart_.push_back(Avalue_.size());
// 	}
// }

// void HighsAggregate::setColBasis(){
//   if (col_status.size()){
//     col_status_.assign(numCol_, HighsBasisStatus::LOWER);
//     vector<bool> nonBasics(numCol, false);
//     for (int i = 0; i < numCol_; ++i){
//       int rep = C[i].front();
//       int pCol = previousColumnColoring[rep];
//       if ((col_status[pCol] == HighsBasisStatus::LOWER
//       || col_status[pCol] == HighsBasisStatus::UPPER
//       || col_status[pCol] == HighsBasisStatus::NONBASIC) 
//       && !nonBasics[pCol]){
//         col_status_[i] = col_status[pCol];
//         nonBasics[pCol] = true;
//       }
//       else if ((col_status[pCol] == HighsBasisStatus::LOWER
//       || col_status[pCol] == HighsBasisStatus::UPPER
//       || col_status[pCol] == HighsBasisStatus::NONBASIC) 
//       && nonBasics[pCol])
//         col_status_[i] = HighsBasisStatus::BASIC;
//       else
//         col_status_[i] = col_status[pCol];
//     }
//   }
//   else return;
// }

// void HighsAggregate::setRowBasis(){
//   if (row_status.size()){
//     row_status_.assign(numRow_, HighsBasisStatus::BASIC);
//     vector<bool> nonBasics(numRow, false);
//     for (int i = 0; i < numRow_; ++i){
//       int rep = C[i + numCol].front() - numCol;
// 			int pCol = previousRowColoring[rep] - numCol;
// 			if ((row_status[pCol] == HighsBasisStatus::LOWER
// 			|| row_status[pCol] == HighsBasisStatus::UPPER
//       || row_status[pCol] == HighsBasisStatus::NONBASIC)
// 			&& !nonBasics[pCol]){
// 				row_status_[i] = row_status[pCol];
// 				nonBasics[pCol] = true;
// 			}
// 			else if ((row_status[pCol] == HighsBasisStatus::LOWER
// 			|| row_status[pCol] == HighsBasisStatus::UPPER
//       || row_status[pCol] == HighsBasisStatus::NONBASIC)
// 			&& nonBasics[pCol])
// 				row_status_[i] = HighsBasisStatus::BASIC;
// 			else
// 				row_status_[i] = row_status[pCol];
//     }
//   }
//   else return;
// }

// void HighsAggregate::appendLinkersToLp(){
//   numPairs = 0;
//   for (int i = 0; i < parentPartition.size(); ++i)
//     if (parentPartition[i] > -1) numPairs++;
//   if (!numPairs) return;
//   int num_new_col = numPairs;
// 	int num_new_row = numPairs;
// 	int num_new_nz = numPairs * 3;  
// 	vector<double> linkColCost, linkColBounds, linkRowBounds, linkAvalue, linkARvalue;
// 	vector<int> linkAstart, linkAindex, linkARstart, linkARindex;
// 	linkARstart.push_back(0);
//   int nnz = Astart_[numCol_];
//   numPairs = 0;
// 	for (int i = 0; i < parentPartition.size(); ++i){
//     if (parentPartition[i] < 0) continue;
// 		int x0 = parentPartition[i], x1 = i, r = numCol_ + numPairs;
//     linkers.push_back(r);
// 		linkARvalue.push_back(1); linkARvalue.push_back(-1); linkARvalue.push_back(-1);
// 		linkARindex.push_back(x0); linkARindex.push_back(x1); linkARindex.push_back(r);
// 		linkARstart.push_back(linkARvalue.size());
// 		linkRowBounds.push_back(0);
// 		linkColBounds.push_back(0);
// 		linkColCost.push_back(0);
//     Astart_.push_back(nnz);
//     numPairs++;
// 	}
// 	// transpose(linkAstart, linkAindex, linkAvalue, linkARstart, linkARindex, linkARvalue);
// 	appendColsToLpVectors(num_new_col, linkColCost, linkColBounds, linkColBounds);
//   numCol_ += numPairs;
// 	// appendColsToLpMatrix(num_new_col, num_new_nz, linkAstart, linkAindex, linkAvalue);
// 	appendRowsToLpVectors(num_new_row, linkRowBounds, linkRowBounds);
//   appendRowsToLpMatrix(num_new_row, num_new_nz, linkARstart, linkARindex, linkARvalue);
//   numRow_ += numPairs;
//   numTot_ = numRow_ + numCol_;
// }

// void HighsAggregate::appendColsToLpVectors(const int num_new_col,
//                                   vector<double>& XcolCost,
//                                   vector<double>& XcolLower,
//                                   vector<double>& XcolUpper) {
//   int new_num_col = numCol_ + num_new_col;
//   colCost_.resize(new_num_col);
//   colLower_.resize(new_num_col);
//   colUpper_.resize(new_num_col);
//   col_status_.resize(new_num_col);
//   for (int new_col = 0; new_col < num_new_col; new_col++) {
//     int iCol = numCol_ + new_col;
//     colCost_[iCol] = XcolCost[new_col];
//     colLower_[iCol] = XcolLower[new_col];
//     colUpper_[iCol] = XcolUpper[new_col];
//     col_status_[iCol] = HighsBasisStatus::LOWER;
//     // Cannot guarantee to create unique names, so name is blank
//   }
// }

// // void HighsAggregate::appendColsToLpMatrix(const int num_new_col,
// //                                  const int num_new_nz, vector<int>& XAstart,
// //                                  vector<int>& XAindex, vector<double>& XAvalue) {
// //   // Determine the new number of columns in the matrix and resize the
// //   // starts accordingly. 
// //   int new_num_col = numCol_ + num_new_col;
// //   Astart_.resize(new_num_col + 1);
// //   // If adding columns to an empty LP then introduce the start for the
// //   // fictitious column 0
// //   if (numCol_ == 0) Astart_[0] = 0;

// //   // Determine the current number of nonzeros and the new number of nonzeros
// //   int current_num_nz = Astart_[numCol_];
// //   int new_num_nz = current_num_nz + num_new_nz;

// //   // Append the starts of the new columns
// //   if (num_new_nz) {
// //   // Nontrivial number of nonzeros being added, so use XAstart
// //     for (int col = 0; col < num_new_col; col++)
// //       Astart_[numCol_ + col] = current_num_nz + XAstart[col];
// //   } else {
// //     // No nonzeros being added, so XAstart may be null, but entries of
// //     // zero are implied.
// //     for (int col = 0; col < num_new_col; col++)
// //       Astart_[numCol_ + col] = current_num_nz;    
// //   }
// //   Astart_[numCol_ + num_new_col] = new_num_nz;

// //   // Adding a non-trivial matrix: resize the column-wise matrix arrays
// //   // accordingly
// //   Aindex_.resize(new_num_nz);
// //   Avalue_.resize(new_num_nz);
// //   // Copy in the new indices and values
// //   for (int el = 0; el < num_new_nz; el++) {
// //     Aindex_[current_num_nz + el] = XAindex[el];
// //     Avalue_[current_num_nz + el] = XAvalue[el];
// //   }
// // }

// void HighsAggregate::appendRowsToLpMatrix(const int num_new_row,
//                                  const int num_new_nz, vector<int>& XARstart,
//                                  vector<int>& XARindex, vector<double>& XARvalue) {
 
//   int current_num_nz = Astart_[numCol_];
//   vector<int> Alength;
//   Alength.assign(numCol_ + numPairs, 0);
//   for (int el = 0; el < num_new_nz; el++) Alength[XARindex[el]]++;
//   // Determine the new number of nonzeros and resize the column-wise matrix
//   // arrays
//   int new_num_nz = current_num_nz + num_new_nz;
//   Aindex_.resize(new_num_nz);
//   Avalue_.resize(new_num_nz);

//   // Append the new rows
//   // Shift the existing columns to make space for the new entries
//   int new_el = new_num_nz;
//   for (int col = numCol_ - 1; col >= 0; col--) {
//     int start_col_plus_1 = new_el;
//     new_el -= Alength[col];
//     for (int el = Astart_[col + 1] - 1; el >= Astart_[col]; el--) {
//       new_el--;
//       Aindex_[new_el] = Aindex_[el];
//       Avalue_[new_el] = Avalue_[el];
//     }
//     Astart_[col + 1] = start_col_plus_1;
//   }
//   assert(new_el == 0);

//   // Insert the new entries
//   for (int row = 0; row < num_new_row; row++) {
//     int first_el = XARstart[row];
//     int last_el = (row < num_new_row - 1 ? XARstart[row + 1] : num_new_nz);
//     for (int el = first_el; el < last_el; el++) {
//       int col = XARindex[el];
//       new_el = Astart_[col + 1] - Alength[col];
//       Alength[col]--;
//       Aindex_[new_el] = numRow_ + row;
//       Avalue_[new_el] = XARvalue[el];
//     }
//   }
// }

// void HighsAggregate::appendRowsToLpVectors(const int num_new_row,
//                                   vector<double>& XrowLower,
//                                   vector<double>& XrowUpper) {
//   int new_num_row = numRow_ + num_new_row;
//   rowLower_.resize(new_num_row);
//   rowUpper_.resize(new_num_row);
//   row_status_.resize(new_num_row);
//   for (int new_row = 0; new_row < num_new_row; new_row++) {
//     int iRow = numRow_ + new_row;
//     rowLower_[iRow] = XrowLower[new_row];
//     rowUpper_[iRow] = XrowUpper[new_row];
//     row_status_[iRow] = HighsBasisStatus::LOWER;
//     // Cannot guarantee to create unique names, so name is blank
//   }
// }

// void HighsAggregate::transpose(vector<int>& xAstart, vector<int>& xAindex, vector<double>& xAvalue,
// 					vector<int>& xARstart, vector<int>& xARindex, vector<double> &xARvalue){
// 	vector<int> iwork(numCol_ + numPairs, 0);
// 	xAstart.resize(numPairs + 1, 0);
// 	int AcountX = numPairs * 3;
// 	xAindex.resize(AcountX);
// 	xAvalue.resize(AcountX);
// 	for (int k = 0; k < AcountX; k++) iwork.at(xARindex.at(k))++;
// 	for (int i = 1; i <= numPairs; i++)
// 		xAstart.at(i) = xAstart.at(i - 1) + iwork.at(i - 1);
// 	for (int i = 0; i < numPairs; i++) iwork.at(i) = xAstart.at(i);
// 	for (int iRow = 0; iRow < numPairs; iRow++) {
// 		for (int k = xARstart.at(iRow); k < xARstart.at(iRow + 1); k++) {
// 		int iCol = xARindex.at(k);
// 		int iPut = iwork.at(iCol)++;
// 		xAindex.at(iPut) = iRow;
// 		xAvalue.at(iPut) = xARvalue[k];
// 		}
// 	}					
// }

// // void HighsAggregate::initialAggregateAMatrix(){
// // 	HighsTimer timer;
// // 	bool run_highs_clock_already_running = timer.runningRunHighsClock();
// //   	if (!run_highs_clock_already_running) timer.startRunHighsClock();
// //   	double initial_time = timer.readRunHighsClock();
// // 	int i, j;
// // 	numCol_ = 0;
// // 	numRow_ = 0;
// // 	for (i = 0; i < C.size(); ++i){
// // 		if (C[i].size() && i < numCol){
// // 			numCol_++;
// // 			realNumCol_++;
// // 		}
// // 		else if(C[i].size()){
// // 			numRow_++;
// // 			realNumRow_++;
// // 		}
// // 	}
// // 	numTot_ = numCol_ + numRow_;
// // 	//std::cout << "counter numRow and numCol time: " <<  timer.readRunHighsClock() - initial_time << std::endl;
// // 	//std::cin.get();
// // 	initial_time = timer.readRunHighsClock();
// // 	Astart_.push_back(0);
// // 	coeff.clear();
// // 	coeffIdx.clear();
// // 	for (i = 0; i < numCol_; ++i){
// // 		rowCoeff(i);
// // 		//std::cout << i << std::endl;
// // 		for (j = 0; j < coeff.size(); ++j){
// // 			Avalue_.push_back(coeff[j]);
// // 			Aindex_.push_back(coeffIdx[j]);
// // 		}
// // 		Astart_.push_back(Avalue_.size());
// // 	}
// // 	//std::cout << "built Amatrix only time: " << timer.readRunHighsClock() - initial_time << std::endl;
// // 	initial_time = timer.readRunHighsClock();
// // 	colUpper_.assign(numCol_, HIGHS_CONST_INF); 
// // 	colLower_.assign(numCol_, -HIGHS_CONST_INF);
// // 	rowUpper_.assign(numRow_, HIGHS_CONST_INF); 
// // 	rowLower_.assign(numRow_, -HIGHS_CONST_INF);
// // 	colCost_.assign(numCol_, 0);
// // 	//std::cout << "built A matrix and all row/column vectors time: " << timer.readRunHighsClock() - initial_time << std::endl;
// // 	//std::cin.get();
// // }

// // void HighsAggregate::examinePartition(){
// // 	HighsTimer timer;
// // 	bool run_highs_clock_already_running = timer.runningRunHighsClock();
// //   	if (!run_highs_clock_already_running) timer.startRunHighsClock();
// //   	double initial_time = timer.readRunHighsClock();
// // 	int i, j;
// // 	for (i = 0; i < C.size(); ++i){
// // 		if (C[i].size() && i < numCol){
// // 			numCol_++;
// // 			realNumCol_++;
// // 		}
// // 		else if(C[i].size()){
// // 			numRow_++;
// // 			realNumRow_++;
// // 		}
// // 	}
// // 	numCol_ += linkingPairs.size();
// // 	numRow_ += linkingPairs.size();
// // 	Astart_.assign(numCol_ + 1, 0);
// // 	colUpper_.assign(numCol_, HIGHS_CONST_INF); 
// // 	colLower_.assign(numCol_, -HIGHS_CONST_INF);
// // 	rowUpper_.assign(numRow_, HIGHS_CONST_INF); 
// // 	rowLower_.assign(numRow_, -HIGHS_CONST_INF);
// // 	colCost_.assign(numCol_, 0);
// // 	//std::cout << "examine part and row/col vectors set up time: " << timer.readRunHighsClock() - initial_time << std::endl;
// // 	//std::cin.get();
// // }

// // void HighsAggregate::collectColumns(){
// // 	HighsTimer timer;
// // 	bool run_highs_clock_already_running = timer.runningRunHighsClock();
// //   	if (!run_highs_clock_already_running) timer.startRunHighsClock();
// //   	double initial_time = timer.readRunHighsClock();
// // 	int i, j, x1, x2, rIdx, aIdx, rowIdx;
// // 	for (i = 0; i < numCol_; ++i){
// // 		if (i < realNumCol_){
// // 			rowCoeff(i);
// // 			for (j = 0; j < coeff.size(); ++j){
// // 				Avalue_.push_back(coeff[j]);
// // 				Aindex_.push_back(coeffIdx[j]);
// // 			}
// // 		}
// // 		Astart_[i + 1] = Avalue_.size();
// // 	}
// // 	setRhsAndBounds();
// // 	//std::cout << "built Amatrix and set bounds times: " << timer.readRunHighsClock() - initial_time << std::endl;
// // 	//std::cin.get();
// // 	initial_time = timer.readRunHighsClock();
// // 	rIdx = realNumCol_;
// // 	rowIdx = realNumRow_;
// // 	coeff.clear();
// // 	coeffIdx.clear();
// // 	for (i = 0; i < linkingPairs.size(); ++i){
// // 		x1 = linkingPairs[i].first;
// // 		x2 = linkingPairs[i].second;
// // 		coeff.push_back(1), coeff.push_back(-1), coeff.push_back(-1);
// // 		coeffIdx.push_back(x1), coeffIdx.push_back(x2); coeffIdx.push_back(rIdx);
// // 		linkers.push_back(rIdx);
// // 		appendLinkersToAMatrix(coeff, coeffIdx, rowIdx, rIdx);
// // 		appendLinkersToColBounds(rIdx);
// // 		appendLinkersToRowRhs(rowIdx);
// // 		rIdx++;
// // 		rowIdx++;
// // 		coeff.clear();
// // 		coeffIdx.clear();
// // 	}
// // 	//std::cout << "linkers added time: " << timer.readRunHighsClock() - initial_time << std::endl;
// // 	//std::cin.get();
// // }

// // void HighsAggregate::findLpBasis(){
// // 	alpBasis.col_status.resize(numCol_);
// // 	alpBasis.row_status.resize(numRow_);
// // 	vector<bool> nonBasics(numTot, false);
// // 	for (int i = 0; i < numCol_; ++i){
// // 		if (i < numCol_ - linkingPairs.size()){
// // 			int rep = C[i].front();
// // 			int pCol = previousColumnColoring[rep];
// // 			if ((col_status[pCol] == HighsBasisStatus::LOWER
// // 			|| col_status[pCol] == HighsBasisStatus::UPPER) 
// // 			&& !nonBasics[pCol]){
// // 				alpBasis.col_status[i] = col_status[pCol];
// // 				nonBasics[pCol] = true;
// // 			}
// // 			else if ((col_status[pCol] == HighsBasisStatus::LOWER
// // 			|| col_status[pCol] == HighsBasisStatus::UPPER) 
// // 			&& nonBasics[pCol])
// // 				alpBasis.col_status[i] = HighsBasisStatus::BASIC;
// // 			else
// // 				alpBasis.col_status[i] = col_status[pCol];
// // 		}
// // 		else
// // 			alpBasis.col_status[i] = HighsBasisStatus::NONBASIC;
// // 	}
// // 	for (int i = 0; i < numRow_; ++i){
// // 		if (i < numRow_ - linkingPairs.size()){
// // 			int rep = C[i + numCol].front() - numCol;
// // 			int pCol = previousRowColoring[rep];
// // 			if ((row_status[pCol - numCol] == HighsBasisStatus::LOWER
// // 			|| row_status[pCol - numCol] == HighsBasisStatus::UPPER)
// // 			&& !nonBasics[pCol]){
// // 				alpBasis.row_status[i] = row_status[pCol - numCol];
// // 				nonBasics[pCol] = true;
// // 			}
// // 			else if ((row_status[pCol - numCol] == HighsBasisStatus::LOWER
// // 			|| row_status[pCol - numCol] == HighsBasisStatus::UPPER)
// // 			&& nonBasics[pCol])
// // 				alpBasis.row_status[i] = HighsBasisStatus::BASIC;
// // 			else
// // 				alpBasis.row_status[i] = row_status[pCol - numCol];
			
// // 		}
// // 		else
// // 			alpBasis.row_status[i] = HighsBasisStatus::NONBASIC;
// // 	}
// // }

// // void HighsAggregate::aggregateCVector(){
// // 	colCost_.assign(numCol_, 0);
// // 	for (int i = 0; i < numCol_; ++i){
// // 		if (i < numCol_ - linkingPairs.size()){
// // 			for (int j = 0; j < C[i].size(); ++j)
// // 				colCost_[i] += colCost[C[i][j]];
// // 		}
// // 	}
// // }

// // void HighsAggregate::appendLinkersToAMatrix(vector<double>& row, vector<int> idx, int rowIdx, int rIdx){
// // 	for (int i = 0; i < row.size(); ++i){
// // 		int colStartIdx = Astart_[idx[i] + 1];
// // 		for (int j = idx[i] + 1; j < Astart_.size(); ++j)
// // 			Astart_[j]++;
// // 		Aindex_.insert(Aindex_.begin() + colStartIdx, rowIdx);
// // 		Avalue_.insert(Avalue_.begin() + colStartIdx, row[i]);
// // 	}
// // }

// // void HighsAggregate::appendLinkersToRowRhs(int rowIdx){
// // 	rowLower_[rowIdx] = 0;
// // 	rowUpper_[rowIdx] = 0;
// // }

// // void HighsAggregate::appendLinkersToColBounds(int rIdx){
// // 	colLower_[rIdx] = 0;
// // 	colUpper_[rIdx] = 0;
// // }

// // void HighsAggregate::setInitialRhsAndBounds(){
// // 	for (int i = 0; i < numRow_; ++i){
// // 		int rep = C[i + numCol].front();
// // 		rowLower_[i] = rowLower[rep - numCol];
// // 		rowUpper_[i] = rowUpper[rep - numCol];
// // 	}
// // 	for (int i = 0; i < numCol_; ++i){
// // 		int rep = C[i].front();
// // 		colLower_[i] = colLower[rep];
// // 		colUpper_[i] = colUpper[rep];
// // 	}
// // }

// // void HighsAggregate::setRhsAndBounds(){
// // 	for (int i = 0; i < realNumRow_; ++i){
// // 		int rep = C[i + numCol].front() - numCol;
// // 		int pCol = previousRowColoring[rep] - numCol;
// // 		if (row_value[pCol] == rowLower[rep]){
// // 			rowLower_[i] = rowLower[rep];
// // 			rowUpper_[i] = rowLower[rep];
// // 		}
// // 		else if (row_value[pCol] == rowUpper[rep]){
// // 			rowLower_[i] = rowUpper[rep];
// // 			rowUpper_[i] = rowUpper[rep];
// // 		}
// // 		else{
// // 			rowLower_[i] = rowLower[rep];
// // 			rowUpper_[i] = rowUpper[rep];
// // 		}
// // 	}
// // 	for (int i = 0; i < realNumCol_; ++i){
// // 		int rep = C[i].front();
// // 		int pCol = previousColumnColoring[rep];
// // 		if (col_value[pCol] == colUpper[rep]){
// // 			colUpper_[i] = colUpper[rep];
// // 			colLower_[i] = colUpper[rep];
// // 		}
// // 		else if (col_value[pCol] == colLower[rep]){
// // 			colUpper_[i] = colLower[rep];
// // 			colLower_[i] = colLower[rep];
// // 		}
// // 		else{
// // 			colLower_[i] = colLower[rep];
// // 			colUpper_[i] = colUpper[rep];
// // 		}
// // 	}
// // }

// // void HighsAggregate::rowCoeff(int column){
// // 	coeff.clear();
// // 	coeffIdx.clear();
// // 	int i, j, var, rep, colWeight = 0;
// // 	for (i = 0; i < realNumRow_; ++i){
// // 		rep = C[i + numCol].front();
// // 		for (j = 0; j < adjListLab[rep].size(); ++j){
// // 			var = adjListLab[rep][j];
// // 			if (color[var] == column){
// // 				colWeight += adjListWeight[rep][j]; 
// // 			}
// // 		}
// // 		if (colWeight){
// // 			coeff.push_back(colWeight);
// // 			coeffIdx.push_back(i);
// // 		}
// // 		colWeight = 0;
// // 	}
// // }

