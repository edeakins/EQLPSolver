#include "Aggregate.h"
using namespace std;

int cCounter = 0;
int rCounter = 0;
int numBasic = 0;

HighsAggregate::HighsAggregate(HighsLp& lp, const struct lpPartition* ep){
	/* Make elp and partition point to lp and ep */
  elp = &lp;
  partition = ep;
}

void HighsAggregate::allocateAlp(){
  /* Allocate Space for alp/alpBasis and its data structs */
  alp = (HighsLp *)calloc(1, sizeof(HighsLp));
  alp->numCol_ = alpNumCol_ = partition->numCol_;
  alp->numRow_ = alpNumRow_ = partition->numRow_;
  alpNumTot_ = alpNumCol_ + alpNumRow_;
  alp->colCost_.resize(alpNumCol_);
  alp->colLower_.resize(alpNumCol_);
  alp->colUpper_.resize(alpNumCol_);
  alp->rowLower_.resize(alpNumRow_);
  alp->rowUpper_.resize(alpNumRow_);
  alp->Astart_.resize(alpNumCol_ + 1);
  alp->Aindex_.reserve(elpNnz_);
  alp->Avalue_.reserve(elpNnz_);
  coeff.assign(numRow, 0);
  /* Allocate space for partition info */
  cell = partition->cell;
  cellFront = partition->cellFront;
  cellSize = partition->cellSize;
  labels = partition->labels;
  parents = partition->parents;
  prevRow.resize(numRow);
  row = partition->row;
  prevCol.resize(numCol);
  col = partition->col;
  cellToCol = partition->cellToCol;
  cellToRow = partition->cellToRow;
  colsToReps = partition->colsToReps;
  repsToCols = partition->repsToCols;
  rowsToReps = partition->rowsToReps;
  repsToRows = partition->repsToRows;
  AindexPacked_.resize(elpNnz_);
  /* Allocate space for linker storage info */
  fixed.assign(numCol, false);
  parents.resize(numCol);
  parent.resize(numLinkers_);
  child.resize(numLinkers_);
  linkARstart.resize(numLinkers_ + 1);
  linkARindex.resize(maxLinkSpace);
  linkARvalue.resize(maxLinkSpace);
  linkAlength.assign(numCol + numLinkers_, 0);
}

void HighsAggregate::allocateElp(){
   /* Resize original LP to contain full extended LP
  info at level N partition */
  elpNnz_ = elp->nnz_;
  elpNumRow_ = elp->numRow_;
  elpNumCol_ = elp->numCol_;
  elpNumTot_ = elpNumRow_ + elpNumCol_;
  elpNumResCols_ = elpNumResRows_ = partition->numResCol_;
  elpNumResNnz_ = partition->numResCol_ * 3;
  elp->colCost_.resize(elpNumCol_ + elpNumResCols_);
  elp->colLower_.resize(elpNumCol_ + elpNumResCols_);
  elp->colUpper_.resize(elpNumCol_ + elpNumResCols_);
  elp->rowLower_.resize(elpNumRow_ + elpNumResRows_);
  elp->rowUpper_.resize(elpNumRow_ + elpNumResRows_);
  elp->Astart_.resize(elpNumCol_ + elpNumResCols_ + 1);
  elp->Aindex_.resize(elpNnz_ + elpNumResNnz_);
  elp->Avalue_.resize(elpNnz_ + elpNumResNnz_);
  elpBasis = (HighsBasis *)calloc(1, sizeof(HighsBasis));
  elpBasis->col_status.resize(elpNumCol_ + elpNumResCols_);
  elpBasis->row_status.resize(elpNumRow_ + elpNumResCols_);
}

void HighsAggregate::fold(){
  packVectors();
  foldObj();
  foldMatrix();
  foldRhs();
  foldBnd();
}

void HighsAggregate::lift(HighsSolution &solution, HighsBasis& basis){
  col_value = solution.col_value;
  row_value = solution.row_value;
  col_status = basis.col_status;
  row_status = basis.row_status;
  liftColBasis();
  liftRowBasis();
  liftBnd();
  liftRhs();
  countNumLinkers();
  liftObjective();
  // liftAMatrix();
  fixAstart();
  makeLinks();
  createLinkRows();
  addCols();
  addRows();
}

// New lifting funcs
void HighsAggregate::liftObjective(){
  int i;
  for (i = 0; i < elpNumCol_; ++i)
    alp->colCost_[i] = elp->colCost_[i];
}

void HighsAggregate::liftBnd(){
  int i;
  HighsBasisStatus status, basic = HighsBasisStatus::BASIC;
  for (i = 0; i < elpNumCol_; ++i){
    int pCol = col[i];
    int pRep = colsToReps[pCol];
    double pVal = col_value[pCol];
    status = col_status[pCol];
    double ub = elp->colUpper_[pRep];
    double lb = elp->colLower_[pRep];
    if (fabs(pVal - ub) < 1e-6 ||
        fabs(pVal - lb) < 1e-6){
      elp->colUpper_[i] = elp->colLower_[i] = pVal;
      // fixed[i] = true;
      if (status != basic) fixed[i] = true;
    }
  }
}

void HighsAggregate::liftRhs(){
  int i;
  for (i = 0; i < elpNumRow_; ++i){
    int pRow = row[i];
    int pRep = rowsToReps[pRow];
    double pVal = row_value[pRow];
    double ub = elp->rowUpper_[pRep];
    double lb = elp->rowLower_[pRep];
    int c = cell[i + elpNumCol_];
    if (fabs(pVal - ub * cellSize[c]) < 1e-6 ||
        fabs(pVal - lb * cellSize[c]) < 1e-6){
      elp->rowUpper_[i] = elp->rowLower_[i] = (double)pVal/cellSize[c];
    }
  }
}

void HighsAggregate::liftColBasis(){
  int i; 
  HighsBasisStatus status, basic = HighsBasisStatus::BASIC;
  elpBasis->numCol_ = elpNumCol_;
  for (i = 0; i < elpNumCol_; ++i){
    int rep = i;
    int pCol = prevCol[i];
    status = col_status[pCol];
    elpBasis->col_status[i] = status;
  }
}

void HighsAggregate::liftRowBasis(){
  int i;
  HighsBasisStatus status, basic = HighsBasisStatus::BASIC;
  elpBasis->numRow_ = elpNumRow_;
  for (i = 0; i < elpNumRow_; ++i)
    elpBasis->row_status[i] = basic;
  for (i = 0; i < alpNumRow_; ++i){
    status = row_status[i];
    if (status != basic){
      int pRep = rowsToReps[i];
      elpBasis->row_status[pRep] = status;
    }
  }
}

void HighsAggregate::fixAstart(){
  int i;
  for (i = elpNumCol_; i < elpNumCol_ + elpNumResCols_; ++i)
    elp->Astart_[i + 1] = elp->Astart_[i];
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
    }
  }
  alp->skip = skip;
}

void HighsAggregate::packVectors(){
  int i, r, rep, c, cf, offset;
  for (i = 0; i < elpNnz_; ++i){
    offset = elp->Aindex_[i] + elpNumCol_;
    c = cell[offset];
    r = cellToRow[c - alpNumCol_];
    AindexPacked_[i] = r;
  }
}

// Fold the objective
void HighsAggregate::foldObj(){
  int i, c, rep;
  for (i = 0; i < alpNumCol_; ++i){
    rep = colsToReps[i];
    c = cell[rep];
    alp->colCost_[i] = elp->colCost_[rep] * cellSize[c];
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
    for (j = elp->Astart_[rep]; j < elp->Astart_[rep + 1]; ++j){
      r = AindexPacked_[j];
      coeff[r] += elp->Avalue_[j];
    }
    for (j = 0; j < numRow_; ++j){
      if (coeff[j]){
        c = cell[rep];
        alp->Avalue_.push_back(coeff[j] * cellSize[c]);
        alp->Aindex_.push_back(j);
        alp->nnz_++;
        coeff[j] = 0;
      }
    }
    alp->Astart_[i + 1] = alp->nnz_;
  }
}

// Fold the row bounds based on current ep
void HighsAggregate::foldRhs(){
  int i, rep, cellIdx;
  for (i = 0; i < alpNumRow_; ++i){
    rep = rowsToReps[i];
    cellIdx = cell[rep + elpNumCol_];
    alp->rowLower_[i] = std::max(elp->rowLower_[rep] * cellSize[cellIdx], -HIGHS_CONST_INF);
    alp->rowUpper_[i] = std::min(elp->rowUpper_[rep] * cellSize[cellIdx], HIGHS_CONST_INF);
  }
}

// Fold the row bounds based on current ep
void HighsAggregate::foldBnd(){
  int i, rep;
  for (i = 0; i < alpNumCol_; ++i){
    rep = colsToReps[i];
    alp->colLower_[i] = elp->colLower_[rep];
    alp->colUpper_[i] = elp->colUpper_[rep];
  }
}

void HighsAggregate::addRows(){
  int i, stop = numLinkers_;
  for (i = 0; i < elpNumResCols_; ++i){
    if (i < stop){
      elp->rowLower_[elp->numRow_] = 0;
      elp->rowUpper_[elp->numRow_] = 0;
      elpBasis->row_status[elp->numRow_++] = HighsBasisStatus::LOWER;
      ++elpBasis->numRow_;
    }
  }
  elp->numRow_ -= stop;
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



