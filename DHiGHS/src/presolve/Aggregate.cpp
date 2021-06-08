#include "Aggregate.h"
using namespace std;

int cCounter = 0;
int rCounter = 0;
int numBasic = 0;

HighsAggregate::HighsAggregate(HighsLp& lp, const struct eq_part* ep, HighsSolution& solution, HighsBasis& basis, 
int numRefinements){
	// To solve and was the LP solved
  numRef = numRefinements;
  solve = true;
  solved = true;
  // From the original lp
  iter = 0;
	numRow = lp.numRow_;
	numCol = lp.numCol_;
	numTot = numRow + numCol;
  nnz = lp.nnz_;
  prevBasis = basis;
  prevSol = solution;
	rowLower.assign(lp.rowLower_.begin(), lp.rowLower_.end());
	rowUpper.assign(lp.rowUpper_.begin(), lp.rowUpper_.end());
  fixed.assign(numCol, false);
  // for (int i = 0; i < lp.colLower_.size(); ++i){
  //   lp.colUpper_[i] = +HIGHS_CONST_INF;
  // }
  // for (int i = 0; i < lp.colLower_.size(); ++i){
  //   lp.colLower_[i] = -HIGHS_CONST_INF;
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
  isRep.assign(numTot, false);
  prevRow.resize(numRow);
  row.resize(numRow);
  prevCol.resize(numCol);
  col.resize(numCol);
  // cellMarked.assign(numCol, false);
  newCellForRep.assign(numTot, -1);
  cellContainsOldRep.assign(numTot, false);
  cellToCol.assign(numCol, -1);
  cellToRow.assign(numRow, -1);
  linked.resize(numCol);
	// Previous solution
	col_value = prevSol.col_value;
	row_value = prevSol.row_value;
	// Previous basis
	col_status = prevBasis.col_status;
	row_status = prevBasis.row_status;
  col_status_.assign(numCol, HighsBasisStatus::BASIC);
  row_status_.assign(numRow, HighsBasisStatus::BASIC);
  nonBasicCol.resize(numCol);
  nonBasicRow.resize(numRow);
  // New Lp info
  // numLinkers_ = 0;
  parents.resize(numCol);
  countNumLinkers();
  skip.assign(numLinkers_, false);
  maxLinkSpace = numLinkers_ * 3;
  // countNumRefinements();
  alp = (HighsLp *)calloc(1, sizeof(HighsLp));
  alpBasis = (HighsBasis *)calloc(1, sizeof(HighsBasis));
  alp->colCost_.resize(numCol + numLinkers_);
  alp->Avalue_.resize(nnz + maxLinkSpace);
  alp->Aindex_.resize(nnz + maxLinkSpace);
  alp->Astart_.resize(numCol + numLinkers_ + 1);
  alp->colUpper_.resize(numCol + numLinkers_);
  alp->colLower_.resize(numCol + numLinkers_);
  alp->rowLower_.resize(numRow + numLinkers_);
  alp->rowUpper_.resize(numRow + numLinkers_);
  alp->linkers.resize(numLinkers_);
  alp->linkLower_.resize(numLinkers_);
  alp->linkUpper_.resize(numLinkers_);
  rowRepsToFix.assign(numRow, false);
  rowRepsValue.resize(numRow);
  colRepsToFix.assign(numCol, false);
  colRepsValue.resize(numCol);
  rowRepsScale.resize(numRow);
  alpBasis->col_status.resize(numCol + numLinkers_);
  alpBasis->row_status.resize(numRow + numLinkers_);
  parent.resize(numLinkers_);
  child.resize(numLinkers_);
  coeff.assign(numRow, 0);
  linkARstart.resize(numLinkers_ + 1);
  linkARindex.resize(maxLinkSpace);
  linkARvalue.resize(maxLinkSpace);
  linkAlength.assign(numCol + numLinkers_, 0);
  // linkLB.assign(numLinkers_, 0);
  // linkUB.assign(numLinkers_, 0);
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

int HighsAggregate::update(HighsSolution& solution, HighsBasis& basis){
  savePartition();
  if (solved){
    prevBasis = basis;
    prevSol = solution;
    saveRowsAndColsFromLastSolve();
    clearLp();
    clearLinks();
    findNonbasicCols();
    findNonbasicRows();
    if (!iter) {findRowRepsToFix(); findColRepsToFix();}
  }
  ++iter;
  translateFrontsToColors();
  identifyLinks();
  if (!numLinkers_) return 0;
  (numLinkers_ == numLinkers_) ? solve = true : solve = false;
  solved = solve;
  if (solve){
    packVectors();
    foldObj();
    foldMatrix();
    fixMatrix();
    foldRhs();
    foldBnds();
    // findNonbasicRows();
    // findNonbasicCols();
    setRowBasis();
    setColBasis();
    createLinkRows();
    addCols();
    addRows();
    return 2;
  }
  else{
    return 1;
  }
  // reset();
  // saveRowsAndColsFromLastSolve();
  // translateFrontsToColors();
  // identifyLinks();
  // if (!numLinkers_) return false;
  // packVectors();
  // foldObj();
  // foldMatrix();
  // fixMatrix();
  
  // foldRhs();
  // foldBnds();
  // findNonbasicRows();
  // findNonbasicCols();
  // setRowBasis();
  // setColBasis();
  // createLinkRows();
  // addCols();
  // addRows();
  // ++iter;
  // return true;
}

void HighsAggregate::lift(HighsSolution &solution, HighsBasis& basis){
  col_value = solution.col_value;
  row_value = solution.row_value;
  col_status = basis.col_status;
  row_status = basis.row_status;
  numLinkers_ = 0; 
  liftColBasis();
  liftRowBasis();
  liftBnd();
  liftRhs();
  countNumLinkers();
  liftObjective();
  liftAMatrix();
  fixAstart();
  makeLinks();
  createLinkRows();
  addCols();
  addRows();
}

void HighsAggregate::countNumRefinements(){
  int i;
  for (i = 0; i < numRef; ++i){
    numLinkers_ += partition[i].nsplits;
  }
  maxLinkSpace = numLinkers_ * 3;
}

void HighsAggregate::translateFrontsToColors(){
  int i = 0, c = 0, r = 0, rep, cRep, pRep, v, pf, cf, colCounter = 0, rowCounter = 0;
  map<int, int> fCell;
  vector<bool> cellMarked(numTot, false);
  // map fronts to colors and count numCol_ and numRow_
  previousNumCol_ = numCol_;
  previousNumRow_ = numRow_;
  numTot_ = 0;
  
  int numColCell = 0;
  for (i = 0; i < numTot; ++i){
    if (fCell.insert(pair<int, int>(partition[iter].fronts[i], c)).second){
      ++c;
      ++numTot_;
    }
  }
  for (i = 0; i < numTot; ++i){
    labels[i] = partition[iter].labels[i];
    if (i < numCol) parents[i] = partition[iter].parents[i];
    cell[i] = fCell.find(partition[iter].fronts[i])->second;
    cellFront[fCell.find(partition[iter].fronts[i])->second] = partition[iter].fronts[i];
    ++cellSize[fCell.find(partition[iter].fronts[i])->second];
  }
  /* New mapping code: Go through col by col in ascending order,
  when a cell is counted it is marked as true, go to next col that
  that is in a new cell and mark that cell */
  for (i = 0; i < numTot; ++i){
    c = cell[i];
    if (cellMarked[c]){ 
      i < numCol ? col[i] = cellToCol[c] : row[i - numCol] = cellToRow[c - numCol_];
      continue;
    }
    if (i < numCol){
      repsToCols[i] = numCol_;
      colsToReps[numCol_] = i;
      cellToCol[c] = numCol_;
      col[i] = numCol_++;
      cellMarked[c] = true;
    }
    else{
      repsToRows[i - numCol] = numRow_;
      rowsToReps[numRow_] = i - numCol;
      cellToRow[c - numCol_] = numRow_;
      row[i - numCol] = numRow_++;
      cellMarked[c] = true;
    }
  }
  alp->numCol_ = numCol_;
  alp->numRow_ = numRow_;
  prevCol = col;
  prevRow = row;
}

// New lifting funcs
void HighsAggregate::liftObjective(){
  int i;
  for (i = 0; i < numCol; ++i)
    alp->colCost_[i] = colCost[i];
}

void HighsAggregate::liftBnd(){
  int i;
  HighsBasisStatus status, basic = HighsBasisStatus::BASIC;
  for (i = 0; i < numCol; ++i){
    int pCol = col[i];
    int pRep = colsToReps[pCol];
    double pVal = col_value[pCol];
    status = col_status[pCol];
    double ub = colUpper[pRep];
    double lb = colLower[pRep];
    if (fabs(pVal - ub) < 1e-6 ||
        fabs(pVal - lb) < 1e-6){
      alp->colUpper_[i] = alp->colLower_[i] = pVal;
      // fixed[i] = true;
      if (status != basic) fixed[i] = true;
    }
    else{
      alp->colUpper_[i] = colUpper[i];
      alp->colLower_[i] = colLower[i];
    }
  }
  // for (i = 0; i < numCol; ++i){
  //   alp->colLower_[i] = colLower[i];
  //   alp->colUpper_[i] = colUpper[i];
  // }
  // for (i = 0; i < col_value.size(); ++i){
  //   int pCol = i;
  //   int pRep = colsToReps[pCol];
  //   double pVal = col_value[pCol];
  //   status = col_status[pCol];
  //   double ub = colUpper[pRep];
  //   double lb = colLower[pRep];
  //   if (fabs(pVal - ub) < 1e-6 ||
  //       fabs(pVal - lb) < 1e-6){
  //     alp->colUpper_[i] = alp->colLower_[i] = pVal;
  //     // fixed[i] = true;
  //     if (status != basic) fixed[i] = true;
  //   }
  // }
}

void HighsAggregate::liftRhs(){
  int i;
  for (i = 0; i < numRow; ++i){
    int pRow = row[i];
    int pRep = rowsToReps[pRow];
    double pVal = row_value[pRow];
    double ub = rowUpper[pRep];
    double lb = rowLower[pRep];
    int c = cell[i + numCol];
    if (fabs(pVal - ub * cellSize[c]) < 1e-6 ||
        fabs(pVal - lb * cellSize[c]) < 1e-6){
      alp->rowUpper_[i] = alp->rowLower_[i] = (double)pVal/cellSize[c];
    }
    else{
      alp->rowUpper_[i] = rowUpper[i];
      alp->rowLower_[i] = rowLower[i];
    }
  }
  // for (i = 0; i < numRow; ++i){
  //   alp->rowLower_[i] = rowLower[i];
  //   alp->rowUpper_[i] = rowUpper[i];
  // }
  // for (i = 0; i < row_value.size(); ++i){
  //   int pRow = i;
  //   int pRep = rowsToReps[pRow];
  //   double pVal = row_value[pRow];
  //   double ub = rowUpper[pRep];
  //   double lb = rowLower[pRep];
  //   int c = cell[i + numCol];
  //   if (fabs(pVal - ub * cellSize[c]) < 1e-6 ||
  //       fabs(pVal - lb * cellSize[c]) < 1e-6){
  //     alp->rowUpper_[i] = alp->rowLower_[i] = (double)pVal/cellSize[c];
  //   }
  // }
}

void HighsAggregate::liftColBasis(){
  int i; 
  HighsBasisStatus status, basic = HighsBasisStatus::BASIC;
  alpBasis->numCol_ = numCol;
  // for (i = 0; i < numCol; ++i)
  //   alpBasis->col_status[i] = basic;
  // numBasic = numCol;
  // for (i = 0; i < numCol_; ++i){
  //   status = col_status[i];
  //   if (status != basic){
  //     numBasic--;
  //     int pRep = colsToReps[i];
  //     alpBasis->col_status[pRep] = status;
  //   }
  // }
  for (i = 0; i < numCol; ++i){
    int rep = i;
    int pCol = prevCol[i];
    status = col_status[pCol];
    if (status == basic) numBasic++;
    alpBasis->col_status[i] = status;
  }
}

void HighsAggregate::liftRowBasis(){
  int i;
  HighsBasisStatus status, basic = HighsBasisStatus::BASIC;
  alpBasis->numRow_ = numRow;
  for (i = 0; i < numRow; ++i)
    alpBasis->row_status[i] = basic;
  numBasic += numRow;
  for (i = 0; i < numRow_; ++i){
    status = row_status[i];
    if (status != basic){
      numBasic--;
      int pRep = rowsToReps[i];
      alpBasis->row_status[pRep] = status;
    }
  }
}

void HighsAggregate::liftAMatrix(){
  int i;
  numCol_ = numCol;
  numRow_ = numRow;
  alp->numCol_ = numCol_;
  alp->numRow_ = numRow_;
  alp->Astart_[0] = Astart[0];
  for (i = 1; i <= numCol_; ++i)
    alp->Astart_[i] = Astart[i];
  for (i = 0; i < nnz; ++i){
    alp->Avalue_[i] = Avalue[i];
    alp->Aindex_[i] = Aindex[i];
  }
  alp->nnz_ = nnz;
}

void HighsAggregate::fixAstart(){
  int i;
  for (i = alp->numCol_; i < numCol + numLinkers_; ++i)
    alp->Astart_[i + 1] = alp->Astart_[i];
}

void HighsAggregate::countNumLinkers(){
  int i;
  for (i = 0; i < numCol; ++i){
    int p = partition[0].parents[i];
    if (p > -1 && !fixed[p]) numLinkers_++;
    parents[i] = partition[0].parents[i];
  }
}

void HighsAggregate::makeLinks(){
  int i, p, c, n = 0, cnt = 0;
  for (i = 0; i < numCol; ++i){
    if (parents[i] > -1){
      p = parents[i];
      bool fix = fixed[p];
      if (fix) continue;
      c = i;
      parent[n] = p;
      child[n++] = i;
      // bool fixedP = fixed[p];
      // bool fixedC = fixed[c];
      // parent[n] = p;
      // child[n] = c;
      // if (fixedP && fixedC){
      //   cnt++;
      //   skip[n++] = true;
      // }
      // else{
      //   skip[n++] = false;
      // }
    }
  }
  alp->skip = skip;
}

// New Lifting funcs

void HighsAggregate::packVectors(){
  int i, r, rep, c, cf, offset;
  for (i = 0; i < nnz; ++i){
    offset = Aindex[i] + numCol;
    c = cell[offset];
    r = cellToRow[c - numCol_];
    // if (r == -1) 
    //   std::cin.get();
    AindexPacked_[i] = r;
  }
}

// Find the rows to fix from the original aggregate
void HighsAggregate::findRowRepsToFix(){
  int oldNCol = col_value.size();
  vector<bool> fixRow(numRow, false);
  for (int i = 0; i < row_value.size(); ++i){
    int rep = prevRowsToReps[i];
    int pCell = previousCell[rep + numCol];
    if (std::fabs(row_value[i] - rowLower[rep] * previousCellSize[pCell]) < 1e-6 ||
        std::fabs(row_value[i] - rowUpper[rep] * previousCellSize[pCell]) < 1e-6){
          fixRow[i] = true;
    }
  }
  for (int i = numCol; i < numTot; ++i){
    int pRow = prevRow[i - numCol];
    int pCell = previousCell[i];
    if (fixRow[pRow]){
      rowRepsToFix[i - numCol] = true;
      rowRepsValue[i - numCol] = row_value[pRow];
      rowRepsScale[i - numCol] = previousCellSize[pCell];
    }
  }
}

// Find which new rows should be basic based off representatives
void HighsAggregate::findNonbasicRows(){
  int i;
  for (i = 0; i < previousNumSolveRow_; ++i){
    if (row_status[i] != HighsBasisStatus::BASIC){
      // int con = i;
      // int pRep = prevRowsToReps[con];
      // int c = cell[pRep + numCol];
      // int cCon = cellToRow[c - numCol_];
      // nonBasicRow[cRow] = true;
      row_status_[i] = row_status[i];
    }
  }
}

void HighsAggregate::findNonbasicCols(){
  int i;
  for (i = 0; i < previousNumSolveCol_; ++i){
    if (col_status[i] != HighsBasisStatus::BASIC){
      // int var = i;
      // int pRep = prevColsToReps[var];
      // int c = cell[pRep];
      // int cVar = cellToCol[c];
      col_status_[i] = col_status[i];
    }
  }
}

// Find the col reps to be fixed from original aggregation
void HighsAggregate::findColRepsToFix(){
  vector<bool> fixCol(numCol, false);
  for (int i = 0; i < col_value.size(); ++i){
    int rep = prevColsToReps[i];
    if (std::fabs(col_value[i] - colLower[rep]) < 1e-6 ||
        std::fabs(col_value[i] - colUpper[rep]) < 1e-6){
          fixCol[i] = true;
    }
  }
  for (int i = 0; i < numCol; ++i){
    int pCol = prevCol[i];
    if (fixCol[pCol]){
      colRepsToFix[i] = true;
      colRepsValue[i] = col_value[pCol];
    }
  }
}

// Fold the objective
void HighsAggregate::foldObj(){
  int i, c, rep;
  for (i = 0; i < numCol_; ++i){
    rep = colsToReps[i];
    c = cell[rep];
    alp->colCost_[i] = colCost[rep] * cellSize[c];
  }
  alp->sense_ = 1;
}

// Fold the matrix down based on current ep
void HighsAggregate::foldMatrix(){
  int i, j, r, c, cf, l, rep;
  alp->Astart_[0] = 0;
  for (i = 0; i < numCol_; ++i){
    rep = colsToReps[i];
    // cf = cellFront[i];
    // rep = labels[cf];
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
  for (i = alp->numCol_; i < numCol + numLinkers_; ++i)
    alp->Astart_[i + 1] = alp->Astart_[i];
}

void HighsAggregate::saveRowsAndColsFromLastSolve(){
  lastSolveCol = col;
  lastSolveRow = row;
  prevRepsToCols = repsToCols;
  prevColsToReps = colsToReps;
  prevRepsToRows = repsToRows;
  prevRowsToReps = rowsToReps;
  col_status = prevBasis.col_status;
	row_status = prevBasis.row_status;
  col_value = prevSol.col_value;
	row_value = prevSol.row_value;
  previousNumSolveCol_ = alp->numCol_ - numLinkers_;
  previousNumSolveRow_ = alp->numRow_ - numLinkers_;
}

void HighsAggregate::clearLp(){
  int i;
  alp->nnz_ = 0;
  alp->pivotTime = 0;
  alp->invertTime = 0;
  alp->unfoldIter = 0;
  // previousNumCol_ = alp->numCol_ - numLinkers_;
  alpBasis->numCol_ = 0;
  alpBasis->numRow_ = 0;
  for (i = 0; i < alp->nnz_; ++i)
    alp->Avalue_[i] = 0;
  for (i = 0; i < numRow; ++i){  
    nonBasicRow[i] = false;
    row_status_[i] = HighsBasisStatus::BASIC;
    cellToRow[i] = -1;
  }
  for (i = 0; i < numCol; ++i){
    nonBasicCol[i] = false;
    col_status_[i] = HighsBasisStatus::BASIC;
    cellToCol[i] = -1;
  }
}

void HighsAggregate::savePartition(){
  int i;
  for (i = 0; i < numTot; ++i){
    previousCell[i] = cell[i];
    previousCellSize[i] = cellSize[i];
    previousCellFront[i] = cellFront[i];
    previousLabels[i] = labels[i];
    cellSize[i] = 0;
    cellContainsOldRep[i] = false;
  }
  prevCol = col;
  prevRow = row;
}

void HighsAggregate::clearLinks(){
  int i;
  for (i = 0; i < numCol; ++i)
    linked[i] = 0;
  for (i = 0; i < maxLinkSpace; ++i){
    linkARindex[i] = 0;
    linkARvalue[i] = 0;
  }
  for (i = 0; i < numLinkers_ + 1; ++i)
    linkARstart[i] = 0;
  for (i = 0; i < numLinkers_; ++i){
    parent[i] = -1;
    child[i] = -1;
  }
  numLinkers_ = 0;
}

void HighsAggregate::reset(){}

// Fold the row bounds based on current ep
void HighsAggregate::foldRhsInit(){
  int i, rep, cellIdx;
  for (i = 0; i < numRow_; ++i){
    rep = rowsToReps[i];
    cellIdx = cell[rep + numCol];
    alp->rowLower_[i] = std::max(rowLower[rep] * cellSize[cellIdx], -HIGHS_CONST_INF);
    alp->rowUpper_[i] = std::min(rowUpper[rep] * cellSize[cellIdx], HIGHS_CONST_INF);
  }
}

void HighsAggregate::foldRhs(){
  int i, rep, rowIdx, cellIdx;
  for (i = 0; i < numRow_; ++i){
    rep = rowsToReps[i];
    rowIdx = prevRow[rep];
    bool fixed = rowRepsToFix[rep];
    cellIdx = cell[rep + numCol];
    if (rowRepsToFix[rep]){
      double rhs = (double)rowRepsValue[rep]/rowRepsScale[rep];
      alp->rowLower_[i] = alp->rowUpper_[i] = rhs * cellSize[cellIdx];
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
      alp->rowLower_[i] = std::max(rowLower[rep] * cellSize[cellIdx], -HIGHS_CONST_INF);
      alp->rowUpper_[i] = std::min(rowUpper[rep] * cellSize[cellIdx], HIGHS_CONST_INF);  
    }
  }
}

// Set row basis
void HighsAggregate::setRowBasis(){
  int i, rep, rowIdx;
  // alpBasis->numRow_ = 0;
  // for (i = 0; i < numRow_; ++i){
  //   if (i == 2619){
  //     std::cout << "problem row: " << i << std::endl;
  //   }
  //   rep = rowsToReps[i];
  //   rowIdx = prevRow[rep];
  //   if ((row_status[rowIdx] == HighsBasisStatus::UPPER ||
  //       row_status[rowIdx] == HighsBasisStatus::LOWER ||
  //       row_status[rowIdx] == HighsBasisStatus::NONBASIC) &&
  //       !nonBasicRow[rowIdx]){
  //     alpBasis->row_status[i] = row_status[rowIdx];
  //     nonBasicRow[rowIdx] = true;
  //     ++alpBasis->numRow_;
  //   }
  //   else if ((row_status[rowIdx] == HighsBasisStatus::UPPER ||
  //           row_status[rowIdx] == HighsBasisStatus::LOWER ||
  //           row_status[rowIdx] == HighsBasisStatus::NONBASIC) &&
  //           nonBasicRow[rowIdx]){
  //     alpBasis->row_status[i] = HighsBasisStatus::BASIC;
  //     ++alpBasis->numRow_;
  //   }
  //   else{
  //     alpBasis->row_status[i] = row_status[rowIdx];
  //     ++alpBasis->numRow_;
  //   }
  // }
  for (i = 0; i < numRow_; ++i){
    alpBasis->row_status[i] = row_status_[i];
    ++alpBasis->numRow_;
  }
}

// Fold the row bounds based on current ep
void HighsAggregate::foldBndsInit(){
  int i, rep;
  for (i = 0; i < numCol_; ++i){
    rep = colsToReps[i];
    alp->colLower_[i] = colLower[rep];
    alp->colUpper_[i] = colUpper[rep];
  }
}

void HighsAggregate::foldBnds(){
  int i, rep, colIdx;
  for (i = 0; i < numCol_; ++i){
    rep = colsToReps[i];
    colIdx = prevCol[rep];
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
  // alpBasis->numCol_ = 0;
  // for (i = 0; i < numCol_; ++i){
  //   rep = colsToReps[i];
  //   colIdx = prevCol[rep];
  //   if ((col_status[colIdx] == HighsBasisStatus::UPPER ||
  //       col_status[colIdx] == HighsBasisStatus::LOWER ||
  //       col_status[colIdx] == HighsBasisStatus::NONBASIC) &&
  //       !nonBasicCol[colIdx]){
  //     alpBasis->col_status[i] = col_status[colIdx];
  //     nonBasicCol[colIdx] = true;
  //     ++alpBasis->numCol_;
  //   }
  //   else if ((col_status[colIdx] == HighsBasisStatus::UPPER ||
  //           col_status[colIdx] == HighsBasisStatus::LOWER ||
  //           col_status[colIdx] == HighsBasisStatus::NONBASIC) &&
  //           nonBasicCol[colIdx]){
  //     alpBasis->col_status[i] = HighsBasisStatus::BASIC;
  //     ++alpBasis->numCol_;
  //   }
  //   else{
  //     alpBasis->col_status[i] = col_status[colIdx];
  //     ++alpBasis->numCol_;
  //   }
  // }
  for (i = 0; i < numCol_; ++i){
    alpBasis->col_status[i] = col_status_[i];
    ++alpBasis->numCol_;
  }
}

void HighsAggregate::addRows(){
  int i, stop = numLinkers_;
  for (i = 0; i < numLinkers_; ++i){
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
  for (int i = alp->numCol_; i < numCol + numLinkers_; ++i)
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
  alp->numLinkers_ = num_new_row;
}

void HighsAggregate::addCols(){
  int i, stop = numLinkers_;
  for (i = 0; i < numLinkers_; ++i){
    if (i < stop){
      alpBasis->col_status[alp->numCol_] = HighsBasisStatus::ZERO;
      alp->colLower_[alp->numCol_] = -HIGHS_CONST_INF;
      alp->colUpper_[alp->numCol_] = +HIGHS_CONST_INF;
      alp->linkers[i] = alp->numCol_;
      ++alp->numCol_;
      ++alpBasis->numCol_;
    }
  }
  alp->numLinkers = numLinkers_;
}

void HighsAggregate::identifyLinks(){
  int i, idx = 0;
  for (i = previousNumCol_; i < numCol_; ++i){
    int rep1 = colsToReps[i];
    int rep2 = colsToReps[prevCol[rep1]];
    int pCol = prevCol[rep1];
    parent[numLinkers_] = pCol;
    child[numLinkers_++] = i;
  }
}

void HighsAggregate::createLinkRows(){
  int i, offset = 0, idx, stop = numLinkers_;
  linkARstart[0] = 0;
  for (i = 0; i < numLinkers_; ++i){
    if (i < stop){
      linkARindex[offset] = parent[i]; linkARvalue[offset] = 1;
      linkARindex[offset + 1] = child[i]; linkARvalue[offset + 1] = -1;
      linkARindex[offset + 2] = alp->numCol_ + i; linkARvalue[offset + 2] = -1;
      linkARstart[i + 1] = linkARstart[i] + 3;
      offset += 3;
    }
  }
  for (i = stop; i < numLinkers_; ++i)
    linkARstart[i + 1] = linkARstart[i];
}

HighsLp* HighsAggregate::getAlp(){
  return alp;
}

HighsBasis* HighsAggregate::getBasis(){
  return alpBasis;
}



