#include "Aggregate.h"
using namespace std;

HighsAggregate::HighsAggregate(HighsLp& lp, const HighsEquitable& ep, HighsSolution& solution, HighsBasis& basis, HighsTableau& tableau, bool flag){
	// From the original lp
	numRow = lp.numRow_;
	numCol = lp.numCol_;
	numTot = numRow + numCol;
	masterIter = lp.masterIter;
	rowLower.assign(lp.rowLower_.begin(), lp.rowLower_.end());
	rowUpper.assign(lp.rowUpper_.begin(), lp.rowUpper_.end());
	colUpper.assign(lp.colUpper_.begin(), lp.colUpper_.end());
	colLower.assign(lp.colLower_.begin(), lp.colLower_.end());
	colCost.assign(lp.colCost_.begin(), lp.colCost_.end());
	Avalue.assign(lp.Avalue_.begin(), lp.Avalue_.end());
	Aindex.assign(lp.Aindex_.begin(), lp.Aindex_.end());
	Astart.assign(lp.Astart_.begin(), lp.Astart_.end());
	realNumRow = lp.realNumRow_;
	realNumCol = lp.realNumCol_;
	//Equitable partition info
	previousRowColoring.assign(ep.previousRowColoring.begin(), ep.previousRowColoring.end());
	previousColumnColoring.assign(ep.previousColumnColoring.begin(), ep.previousColumnColoring.end());
	numRow_ = ep.cCol - numCol;
	numCol_ = ep.vCol
	numTot_ = numRow_ + numCol_;
	C.assign(ep.C.begin(), ep.C.end());
	prevC.assign(ep.prevC.begin(), ep.prevC.end());
	color.assign(ep.color.begin(), ep.color.end());
	linkingPairs.assign(ep.linkingPairs.begin(), ep.linkingPairs.end());
	Xindex.assign(ep.Xindex_.begin(), ep.Xindex_.end());
	// Previous solution
	col_value = (solution.col_value);
	row_value = (solution.row_value);
	// Previous basis
	col_status = (basis.col_status);
	row_status = (basis.row_status);
	// flag status to find linkers or not
	flag_ = flag;
	aggregateAMatrix();
}

// HighsBasis& HighsAggregate::getAlpBasis(){
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
// 	alp.numRealRows = numRow_ - numLinkers_;
// 	alp.numRealCols = numCol_ - numLinkers_;
// 	return alp;
// }

void HighsAggregate::aggregate(){
	if (!flag_){
		aggregateAMatrix();
		// initialAggregateAMatrix();
		// setInitialRhsAndBounds();
		// aggregateCVector();
		return;
	}
	else{
		examinePartition();
		collectColumns();
		findLpBasis();
		aggregateCVector();
	}
}

void HighsAggregate::aggregateColBounds(){

}

void HighsAggregate::aggregateRowBounds(){

}

void HighsAggregate::aggregateCVector(){
	
}

void HighsAggregate::aggregateAMatrix(){
	Astart_.push_back(0);
	for (int i = 0; i < numCol_; ++i){
		int rep = C[i]->front();
		int idx = Xindex[Astart[rep]]; 
		double weight = Avalue[Astart[rep]];
		for (int j = Astart[rep] + 1; j < Astart[rep + 1]; ++j){
			if (idx = Xindex[j])
				weight += Avalue[j];
			else{
				Avalue_.push_back(weight);
				Aindex_.push_back(idx - numCol);
				idx = Xindex[j];
				weight = Avalue[j];
			}
		}
		Aindex_.push_back(Avalue_.size());
	}
	int num_new_col = linkingPairs.size();
	int num_new_row = linkingPairs.size();
	int num_new_nz = linkingPairs.size() * 3;  
	vector<double> linkColCost, linkColBounds, linkRowBounds, linkAvalue, linkARvalue;
	vector<int> linkAstart, linkAindex, linkARstart, linkARindex;
	linkARstart.push_back(0);
	for (int i = 0; i < linkingPairs.size(); ++i){
		int x0 = linkingPairs[i].first, x1 = linkingPairs[i].second, r = numCol_ + i;
		linkARvalue.push_back(1); linkARvalue.push_back(-1); linkARvalue.push_back(-1);
		linkARindex.push_back(x0); linkARindex.push_back(x1); linkARindex.push_back(r);
		linkARstart.push_back(linkARvalue.size());
		linkRowBounds.push_back(0);
		linkColBounds.push_back(0);
		linkColCost.push_back(0);
	}
	transpose(linkAstart, linkAindex, linkAvalue, linkARstart, linkARindex, linkARvalue);
	appendColsToLpVectors(num_new_col, linkColCost, linkColBounds, linkColBounds);
	appendColsToLpMatrix(num_new_col, num_new_nz, linkAstart, linkAindex, linkAvalue);
	appendRowsToLpVectors(num_new_row, linkRowBounds, linkRowBounds);
}

void HighsAggregate::appendColsToLpVectors(const int num_new_col,
                                  const double* XcolCost,
                                  const double* XcolLower,
                                  const double* XcolUpper) {
  int new_num_col = numCol_ + num_new_col;
  colCost_.resize(new_num_col);
  colLower_.resize(new_num_col);
  colUpper_.resize(new_num_col);
  bool have_names = col_names_.size();
  for (int new_col = 0; new_col < num_new_col; new_col++) {
    int iCol = numCol_ + new_col;
    colCost_[iCol] = XcolCost[new_col];
    colLower_[iCol] = XcolLower[new_col];
    colUpper_[iCol] = XcolUpper[new_col];
    // Cannot guarantee to create unique names, so name is blank
  }
}

void HighsAggregate::appendColsToLpMatrix(const int num_new_col,
                                 const int num_new_nz, const int* XAstart,
                                 const int* XAindex, const double* XAvalue) {
  if (num_new_col < 0) return HighsStatus::Error;
  if (num_new_col == 0) return HighsStatus::OK;
  // Check that nonzeros aren't being appended to a matrix with no rows
  if (num_new_nz > 0 && numRow_ <= 0) return HighsStatus::Error;
  // Determine the new number of columns in the matrix and resize the
  // starts accordingly. 
  int new_num_col = numCol_ + num_new_col;
  Astart_.resize(new_num_col + 1);
  // If adding columns to an empty LP then introduce the start for the
  // fictitious column 0
  if (numCol_ == 0) Astart_[0] = 0;

  // Determine the current number of nonzeros and the new number of nonzeros
  int current_num_nz = Astart_[numCol_];
  int new_num_nz = current_num_nz + num_new_nz;

  // Append the starts of the new columns
  if (num_new_nz) {
  // Nontrivial number of nonzeros being added, so use XAstart
    assert(XAstart != NULL);
    for (int col = 0; col < num_new_col; col++)
      Astart_[numCol_ + col] = current_num_nz + XAstart[col];
  } else {
    // No nonzeros being added, so XAstart may be null, but entries of
    // zero are implied.
    for (int col = 0; col < num_new_col; col++)
      Astart_[numCol_ + col] = current_num_nz;    
  }
  Astart_[numCol_ + num_new_col] = new_num_nz;

  // If no nonzeros are being added then there's nothing else to do
  if (num_new_nz <= 0) return HighsStatus::OK;

  // Adding a non-trivial matrix: resize the column-wise matrix arrays
  // accordingly
  Aindex_.resize(new_num_nz);
  Avalue_.resize(new_num_nz);
  // Copy in the new indices and values
  for (int el = 0; el < num_new_nz; el++) {
    Aindex_[current_num_nz + el] = XAindex[el];
    Avalue_[current_num_nz + el] = XAvalue[el];
  }
}

void HighsAggregate::appendRowsToLpMatrix(const int num_new_row,
                                 const int num_new_nz, const int* XARstart,
                                 const int* XARindex, const double* XARvalue) {
 
  int current_num_nz = Astart_[numCol_];
  vector<int> Alength;
  Alength.assign(numCol_, 0);
  for (int el = 0; el < num_new_nz; el++) Alength[XARindex[el]]++;
  // Determine the new number of nonzeros and resize the column-wise matrix
  // arrays
  int new_num_nz = current_num_nz + num_new_nz;
  Aindex_.resize(new_num_nz);
  Avalue_.resize(new_num_nz);

  // Append the new rows
  // Shift the existing columns to make space for the new entries
  int new_el = new_num_nz;
  for (int col = numCol_ - 1; col >= 0; col--) {
    int start_col_plus_1 = new_el;
    new_el -= Alength[col];
    for (int el = Astart_[col + 1] - 1; el >= Astart_[col]; el--) {
      new_el--;
      Aindex_[new_el] = Aindex_[el];
      Avalue_[new_el] = Avalue_[el];
    }
    Astart_[col + 1] = start_col_plus_1;
  }
  assert(new_el == 0);

  // Insert the new entries
  for (int row = 0; row < num_new_row; row++) {
    int first_el = XARstart[row];
    int last_el = (row < num_new_row - 1 ? XARstart[row + 1] : num_new_nz);
    for (int el = first_el; el < last_el; el++) {
      int col = XARindex[el];
      new_el = Astart_[col + 1] - Alength[col];
      Alength[col]--;
      Aindex_[new_el] = numRow_ + row;
      Avalue_[new_el] = XARvalue[el];
    }
  }
}

void HighsAggregate::appendRowsToLpVectors(const int num_new_row,
                                  const double* XrowLower,
                                  const double* XrowUpper) {
  int new_num_row = numRow_ + num_new_row;
  rowLower_.resize(new_num_row);
  rowUpper_.resize(new_num_row);
  bool have_names = row_names_.size();

  for (int new_row = 0; new_row < num_new_row; new_row++) {
    int iRow = numRow_ + new_row;
    rowLower_[iRow] = XrowLower[new_row];
    rowUpper_[iRow] = XrowUpper[new_row];
    // Cannot guarantee to create unique names, so name is blank
  }
}

void HighsAggregate::transpose(vector<int>& xAstart, vector<int>& xAindex, vector<double>& xAvalue,
					vector<int>& xARstart, vector<int>& xARindex, vector<double> &xARvalue){
	vector<int> iwork(numCol_, 0);
	xAstart.resize(linkingPairs.size() + 1, 0);
	int AcountX = linkingPairs.size() * 3;
	xAindex.resize(AcountX);
	xAvalue.resize(AcountX);
	for (int k = 0; k < AcountX; k++) iwork.at(xARindex_.at(k))++;
	for (int i = 1; i <= linkingPairs.size(); i++)
		xAstart.at(i) = xAstart.at(i - 1) + iwork.at(i - 1);
	for (int i = 0; i < linkingPairs.size(); i++) iwork.at(i) = xAstart.at(i);
	for (int iRow = 0; iRow < linkingPairs.size(); iRow++) {
		for (int k = xARstart.at(iRow); k < xARstart.at(iRow + 1); k++) {
		int iCol = xARindex.at(k);
		int iPut = iwork.at(iCol)++;
		xAindex.at(iPut) = iRow;
		xAvalue.at(iPut) = xARvalue[k];
		}
	}					
}

// void HighsAggregate::initialAggregateAMatrix(){
// 	HighsTimer timer;
// 	bool run_highs_clock_already_running = timer.runningRunHighsClock();
//   	if (!run_highs_clock_already_running) timer.startRunHighsClock();
//   	double initial_time = timer.readRunHighsClock();
// 	int i, j;
// 	numCol_ = 0;
// 	numRow_ = 0;
// 	for (i = 0; i < C.size(); ++i){
// 		if (C[i].size() && i < numCol){
// 			numCol_++;
// 			realNumCol_++;
// 		}
// 		else if(C[i].size()){
// 			numRow_++;
// 			realNumRow_++;
// 		}
// 	}
// 	numTot_ = numCol_ + numRow_;
// 	//std::cout << "counter numRow and numCol time: " <<  timer.readRunHighsClock() - initial_time << std::endl;
// 	//std::cin.get();
// 	initial_time = timer.readRunHighsClock();
// 	Astart_.push_back(0);
// 	coeff.clear();
// 	coeffIdx.clear();
// 	for (i = 0; i < numCol_; ++i){
// 		rowCoeff(i);
// 		//std::cout << i << std::endl;
// 		for (j = 0; j < coeff.size(); ++j){
// 			Avalue_.push_back(coeff[j]);
// 			Aindex_.push_back(coeffIdx[j]);
// 		}
// 		Astart_.push_back(Avalue_.size());
// 	}
// 	//std::cout << "built Amatrix only time: " << timer.readRunHighsClock() - initial_time << std::endl;
// 	initial_time = timer.readRunHighsClock();
// 	colUpper_.assign(numCol_, HIGHS_CONST_INF); 
// 	colLower_.assign(numCol_, -HIGHS_CONST_INF);
// 	rowUpper_.assign(numRow_, HIGHS_CONST_INF); 
// 	rowLower_.assign(numRow_, -HIGHS_CONST_INF);
// 	colCost_.assign(numCol_, 0);
// 	//std::cout << "built A matrix and all row/column vectors time: " << timer.readRunHighsClock() - initial_time << std::endl;
// 	//std::cin.get();
// }

// void HighsAggregate::examinePartition(){
// 	HighsTimer timer;
// 	bool run_highs_clock_already_running = timer.runningRunHighsClock();
//   	if (!run_highs_clock_already_running) timer.startRunHighsClock();
//   	double initial_time = timer.readRunHighsClock();
// 	int i, j;
// 	for (i = 0; i < C.size(); ++i){
// 		if (C[i].size() && i < numCol){
// 			numCol_++;
// 			realNumCol_++;
// 		}
// 		else if(C[i].size()){
// 			numRow_++;
// 			realNumRow_++;
// 		}
// 	}
// 	numCol_ += linkingPairs.size();
// 	numRow_ += linkingPairs.size();
// 	Astart_.assign(numCol_ + 1, 0);
// 	colUpper_.assign(numCol_, HIGHS_CONST_INF); 
// 	colLower_.assign(numCol_, -HIGHS_CONST_INF);
// 	rowUpper_.assign(numRow_, HIGHS_CONST_INF); 
// 	rowLower_.assign(numRow_, -HIGHS_CONST_INF);
// 	colCost_.assign(numCol_, 0);
// 	//std::cout << "examine part and row/col vectors set up time: " << timer.readRunHighsClock() - initial_time << std::endl;
// 	//std::cin.get();
// }

// void HighsAggregate::collectColumns(){
// 	HighsTimer timer;
// 	bool run_highs_clock_already_running = timer.runningRunHighsClock();
//   	if (!run_highs_clock_already_running) timer.startRunHighsClock();
//   	double initial_time = timer.readRunHighsClock();
// 	int i, j, x1, x2, rIdx, aIdx, rowIdx;
// 	for (i = 0; i < numCol_; ++i){
// 		if (i < realNumCol_){
// 			rowCoeff(i);
// 			for (j = 0; j < coeff.size(); ++j){
// 				Avalue_.push_back(coeff[j]);
// 				Aindex_.push_back(coeffIdx[j]);
// 			}
// 		}
// 		Astart_[i + 1] = Avalue_.size();
// 	}
// 	setRhsAndBounds();
// 	//std::cout << "built Amatrix and set bounds times: " << timer.readRunHighsClock() - initial_time << std::endl;
// 	//std::cin.get();
// 	initial_time = timer.readRunHighsClock();
// 	rIdx = realNumCol_;
// 	rowIdx = realNumRow_;
// 	coeff.clear();
// 	coeffIdx.clear();
// 	for (i = 0; i < linkingPairs.size(); ++i){
// 		x1 = linkingPairs[i].first;
// 		x2 = linkingPairs[i].second;
// 		coeff.push_back(1), coeff.push_back(-1), coeff.push_back(-1);
// 		coeffIdx.push_back(x1), coeffIdx.push_back(x2); coeffIdx.push_back(rIdx);
// 		linkers.push_back(rIdx);
// 		appendLinkersToAMatrix(coeff, coeffIdx, rowIdx, rIdx);
// 		appendLinkersToColBounds(rIdx);
// 		appendLinkersToRowRhs(rowIdx);
// 		rIdx++;
// 		rowIdx++;
// 		coeff.clear();
// 		coeffIdx.clear();
// 	}
// 	//std::cout << "linkers added time: " << timer.readRunHighsClock() - initial_time << std::endl;
// 	//std::cin.get();
// }

// void HighsAggregate::findLpBasis(){
// 	alpBasis.col_status.resize(numCol_);
// 	alpBasis.row_status.resize(numRow_);
// 	vector<bool> nonBasics(numTot, false);
// 	for (int i = 0; i < numCol_; ++i){
// 		if (i < numCol_ - linkingPairs.size()){
// 			int rep = C[i].front();
// 			int pCol = previousColumnColoring[rep];
// 			if ((col_status[pCol] == HighsBasisStatus::LOWER
// 			|| col_status[pCol] == HighsBasisStatus::UPPER) 
// 			&& !nonBasics[pCol]){
// 				alpBasis.col_status[i] = col_status[pCol];
// 				nonBasics[pCol] = true;
// 			}
// 			else if ((col_status[pCol] == HighsBasisStatus::LOWER
// 			|| col_status[pCol] == HighsBasisStatus::UPPER) 
// 			&& nonBasics[pCol])
// 				alpBasis.col_status[i] = HighsBasisStatus::BASIC;
// 			else
// 				alpBasis.col_status[i] = col_status[pCol];
// 		}
// 		else
// 			alpBasis.col_status[i] = HighsBasisStatus::NONBASIC;
// 	}
// 	for (int i = 0; i < numRow_; ++i){
// 		if (i < numRow_ - linkingPairs.size()){
// 			int rep = C[i + numCol].front() - numCol;
// 			int pCol = previousRowColoring[rep];
// 			if ((row_status[pCol - numCol] == HighsBasisStatus::LOWER
// 			|| row_status[pCol - numCol] == HighsBasisStatus::UPPER)
// 			&& !nonBasics[pCol]){
// 				alpBasis.row_status[i] = row_status[pCol - numCol];
// 				nonBasics[pCol] = true;
// 			}
// 			else if ((row_status[pCol - numCol] == HighsBasisStatus::LOWER
// 			|| row_status[pCol - numCol] == HighsBasisStatus::UPPER)
// 			&& nonBasics[pCol])
// 				alpBasis.row_status[i] = HighsBasisStatus::BASIC;
// 			else
// 				alpBasis.row_status[i] = row_status[pCol - numCol];
			
// 		}
// 		else
// 			alpBasis.row_status[i] = HighsBasisStatus::NONBASIC;
// 	}
// }

// void HighsAggregate::aggregateCVector(){
// 	colCost_.assign(numCol_, 0);
// 	for (int i = 0; i < numCol_; ++i){
// 		if (i < numCol_ - linkingPairs.size()){
// 			for (int j = 0; j < C[i].size(); ++j)
// 				colCost_[i] += colCost[C[i][j]];
// 		}
// 	}
// }

// void HighsAggregate::appendLinkersToAMatrix(vector<double>& row, vector<int> idx, int rowIdx, int rIdx){
// 	for (int i = 0; i < row.size(); ++i){
// 		int colStartIdx = Astart_[idx[i] + 1];
// 		for (int j = idx[i] + 1; j < Astart_.size(); ++j)
// 			Astart_[j]++;
// 		Aindex_.insert(Aindex_.begin() + colStartIdx, rowIdx);
// 		Avalue_.insert(Avalue_.begin() + colStartIdx, row[i]);
// 	}
// }

// void HighsAggregate::appendLinkersToRowRhs(int rowIdx){
// 	rowLower_[rowIdx] = 0;
// 	rowUpper_[rowIdx] = 0;
// }

// void HighsAggregate::appendLinkersToColBounds(int rIdx){
// 	colLower_[rIdx] = 0;
// 	colUpper_[rIdx] = 0;
// }

// void HighsAggregate::setInitialRhsAndBounds(){
// 	for (int i = 0; i < numRow_; ++i){
// 		int rep = C[i + numCol].front();
// 		rowLower_[i] = rowLower[rep - numCol];
// 		rowUpper_[i] = rowUpper[rep - numCol];
// 	}
// 	for (int i = 0; i < numCol_; ++i){
// 		int rep = C[i].front();
// 		colLower_[i] = colLower[rep];
// 		colUpper_[i] = colUpper[rep];
// 	}
// }

// void HighsAggregate::setRhsAndBounds(){
// 	for (int i = 0; i < realNumRow_; ++i){
// 		int rep = C[i + numCol].front() - numCol;
// 		int pCol = previousRowColoring[rep] - numCol;
// 		if (row_value[pCol] == rowLower[rep]){
// 			rowLower_[i] = rowLower[rep];
// 			rowUpper_[i] = rowLower[rep];
// 		}
// 		else if (row_value[pCol] == rowUpper[rep]){
// 			rowLower_[i] = rowUpper[rep];
// 			rowUpper_[i] = rowUpper[rep];
// 		}
// 		else{
// 			rowLower_[i] = rowLower[rep];
// 			rowUpper_[i] = rowUpper[rep];
// 		}
// 	}
// 	for (int i = 0; i < realNumCol_; ++i){
// 		int rep = C[i].front();
// 		int pCol = previousColumnColoring[rep];
// 		if (col_value[pCol] == colUpper[rep]){
// 			colUpper_[i] = colUpper[rep];
// 			colLower_[i] = colUpper[rep];
// 		}
// 		else if (col_value[pCol] == colLower[rep]){
// 			colUpper_[i] = colLower[rep];
// 			colLower_[i] = colLower[rep];
// 		}
// 		else{
// 			colLower_[i] = colLower[rep];
// 			colUpper_[i] = colUpper[rep];
// 		}
// 	}
// }

// void HighsAggregate::rowCoeff(int column){
// 	coeff.clear();
// 	coeffIdx.clear();
// 	int i, j, var, rep, colWeight = 0;
// 	for (i = 0; i < realNumRow_; ++i){
// 		rep = C[i + numCol].front();
// 		for (j = 0; j < adjListLab[rep].size(); ++j){
// 			var = adjListLab[rep][j];
// 			if (color[var] == column){
// 				colWeight += adjListWeight[rep][j]; 
// 			}
// 		}
// 		if (colWeight){
// 			coeff.push_back(colWeight);
// 			coeffIdx.push_back(i);
// 		}
// 		colWeight = 0;
// 	}
// }

