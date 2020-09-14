#include "Aggregate.h"
using namespace std;

HighsAggregate::HighsAggregate(HighsLp& lp, const HighsEquitable& ep, HighsSolution& solution, HighsBasis& basis, HighsTableau& tableau, bool flag){
	// From the original lp
	numRow = lp.numRow_;
	numCol = lp.numCol_;
	numTot = numRow + numCol;
	impliedNumRow = lp.addNumRow_;
	masterIter = lp.masterIter;
	rowLower.assign(lp.rowLower_.begin(), lp.rowLower_.end());
	rowUpper.assign(lp.rowUpper_.begin(), lp.rowUpper_.end());
	colUpper.assign(lp.colUpper_.begin(), lp.colUpper_.end());
	colLower.assign(lp.colLower_.begin(), lp.colLower_.end());
	colCost.assign(lp.colCost_.begin(), lp.colCost_.end());
	Avalue.assign(lp.Avalue_.begin(), lp.Avalue_.end());
	Aindex.assign(lp.Aindex_.begin(), lp.Aindex_.end());
	Astart.assign(lp.Astart_.begin(), lp.Astart_.end());
	impliedARstart.assign(lp.addARstart_.begin(), lp.addARstart_.end());
	impliedARvalue.assign(lp.addARvalue_.begin(), lp.addARvalue_.end());
	impliedARindex.assign(lp.addARindex_.begin(), lp.addARindex_.end());
	ARtableauValue.assign(tableau.ARtableauValue.begin(), tableau.ARtableauValue.end());
	ARtableauIndex.assign(tableau.ARtableauIndex.begin(), tableau.ARtableauIndex.end());
	ARtableauStart.assign(tableau.ARtableauStart.begin(), tableau.ARtableauStart.end());
	ARreducedRHS.assign(tableau.ARreducedRHS.begin(), tableau.ARreducedRHS.end());
	rowColor.assign(lp.rowColor.begin(), lp.rowColor.end());
	activeColorHistory.assign(lp.activeColorHistory.begin(), lp.activeColorHistory.end());
	realNumRow = lp.realNumRow_;
	realNumCol = lp.realNumCol_;
	//Equitable partition info
	previousRowColoring.assign(ep.previousRowColoring.begin(), ep.previousRowColoring.end());
	previousColumnColoring.assign(ep.previousColumnColoring.begin(), ep.previousColumnColoring.end());
	C.assign(ep.C.begin(), ep.C.end());
	prevC.assign(ep.prevC.begin(), ep.prevC.end());
	color.assign(ep.color.begin(), ep.color.end());
	adjListLab.assign(ep.adjListLab.begin(), ep.adjListLab.end());
	adjListWeight.assign(ep.adjListWeight.begin(), ep.adjListWeight.end());
	linkingPairs.assign(ep.linkingPairs.begin(), ep.linkingPairs.end());
	partSize.assign(ep.partSize.begin(), ep.partSize.end());
	previousPartSize.assign(ep.previousPartSize.begin(), ep.previousPartSize.end());
	commonLinkers = ep.commonLinkers;
	// Used for setting active set
	activeBounds_.assign(numCol, false);
	activeConstraints_.assign(numRow, false);
	// Preivous solution
	col_value = (solution.col_value);
	row_value = (solution.row_value);
	// Previous basis
	col_status = (basis.col_status);
	row_status = (basis.row_status);
	// flag status to find linkers or not
	flag_ = flag;
	numSplits.assign(numRow, 0);
	// if (ARtableauStart.size()){
	// 	for (int i = 0; i < ARtableauStart.size() - 1; ++i){
	// 		cout << "prev row: " << i << endl;
	// 		for (int j = ARtableauStart[i]; j < ARtableauStart[i + 1]; ++j){
	// 			cout << ARtableauValue[j] << "x_" << ARtableauIndex[j] << " ";
	// 		}
	// 		cout << endl;
	// 	}
	// 	cin.get();
	// }
	aggregateAMatrix();
}

HighsBasis& HighsAggregate::getAlpBasis(){
	return alpBasis;
}

HighsLp& HighsAggregate::getAlp(){
	alp.numRow_ = (numRow_);
	alp.numCol_ = (numCol_);
	alp.numInt_ = 0;
	alp.nnz_ = Avalue_.size();
	alp.linkers.assign(linkers.begin(), linkers.end());
	//alp.artificialVariables.assign(artificialVariables.begin(), artificialVariables.end());
	alp.Astart_.assign(Astart_.begin(), Astart_.end());
	alp.Aindex_.assign(Aindex_.begin(), Aindex_.end());
	alp.Avalue_.assign(Avalue_.begin(), Avalue_.end());
	alp.colUpper_.assign(colUpper_.begin(), colUpper_.end());
	alp.colLower_.assign(colLower_.begin(), colLower_.end());
	alp.rowUpper_.assign(rowUpper_.begin(), rowUpper_.end());
	alp.rowLower_.assign(rowLower_.begin(), rowLower_.end());
	alp.colCost_.assign(colCost_.begin(), colCost_.end());
	alp.model_name_ = (model_name_);
	alp.lp_name_ = (lp_name_);
	alp.sense_ = 1;
	alp.offset_ = 0;
	alp.numRealRows = numRow_ - numLinkers_;
	alp.numRealCols = numCol_ - numLinkers_;
	return alp;
}

void HighsAggregate::aggregateAMatrix(){
	if (!flag_){
		initialAggregateAMatrix();
		setInitialRhsAndBounds();
		aggregateCVector();
		return;
	}
	else{
		examinePartition();
		collectColumns();
		findLpBasis();
		// aggregateCVector();
	}
}

void HighsAggregate::initialAggregateAMatrix(){
	int i, j;
	numCol_ = 0;
	numRow_ = 0;
	for (i = 0; i < C.size(); ++i){
		if (C[i].size() && i < numCol){
			numCol_++;
			realNumCol_++;
		}
		else if(C[i].size()){
			numRow_++;
			realNumRow_++;
		}
	}
	numTot_ = numCol_ + numRow_;
	Astart_.push_back(0);
	for (i = 0; i < numCol_; ++i){
		vector<double> coeff = rowCoeff(i);
		for (j = 0; j < coeff.size(); ++j){
			if (coeff[j]){
				Avalue_.push_back(coeff[j]);
				Aindex_.push_back(j);
			}
		}
		Astart_.push_back(Avalue_.size());
	}
	colUpper_.assign(numCol_, HIGHS_CONST_INF); 
	colLower_.assign(numCol_, -HIGHS_CONST_INF);
	rowUpper_.assign(numRow_, HIGHS_CONST_INF); 
	rowLower_.assign(numRow_, -HIGHS_CONST_INF);
	colCost_.assign(numCol_, 0);
}

void HighsAggregate::examinePartition(){
	int i, j;
	for (i = 0; i < C.size(); ++i){
		if (C[i].size() && i < numCol){
			realNumCol_++;
			numCol_++;
		}
		else if (C[i].size()){
			realNumRow_++;
			numRow_++;
		}
	}
	numCol_ += linkingPairs.size();
	numRow_ += linkingPairs.size();
	Astart_.assign(numCol_ + 1, 0);
	colUpper_.assign(numCol_, HIGHS_CONST_INF); 
	colLower_.assign(numCol_, -HIGHS_CONST_INF);
	rowUpper_.assign(numRow_, HIGHS_CONST_INF); 
	rowLower_.assign(numRow_, -HIGHS_CONST_INF);
	colCost_.assign(numCol_, 0);
}

void HighsAggregate::collectColumns(){
	int i, j, x1, x2, rIdx, aIdx, rowIdx;
	for (i = 0; i < numCol_; ++i){
		if (i < realNumCol_){
			vector<double> coeff = rowCoeff(i);
			for (j = 0; j < coeff.size(); ++j){
				if (coeff[j]){
					Avalue_.push_back(coeff[j]);
					Aindex_.push_back(j);
				}
			}
		}
		Astart_[i + 1] = Avalue_.size();
	}
	setRhsAndBounds();
	rIdx = realNumCol_;
	rowIdx = realNumRow_;
	for (i = 0; i < linkingPairs.size(); ++i){
		vector<double> coeff(numCol_);
		x1 = linkingPairs[i].first;
		x2 = linkingPairs[i].second;
		coeff[x1] = 1, coeff[x2] = -1, coeff[rIdx] = -1;
		linkers.push_back(rIdx);
		appendLinkersToAMatrix(coeff, rowIdx, rIdx);
		appendLinkersToColBounds(rIdx);
		appendLinkersToRowRhs(rowIdx);
		rIdx++;
		rowIdx++;
	}
}

void HighsAggregate::findLpBasis(){
	alpBasis.col_status.resize(numCol_);
	alpBasis.row_status.resize(numRow_);
	vector<bool> nonBasics(numTot, false);
	for (int i = 0; i < numCol_; ++i){
		if (i < numCol_ - linkingPairs.size()){
			int rep = C[i].front();
			int pCol = previousColumnColoring[rep];
			if ((col_status[pCol] == HighsBasisStatus::LOWER
			|| col_status[pCol] == HighsBasisStatus::UPPER) 
			&& !nonBasics[pCol]){
				alpBasis.col_status[i] = col_status[pCol];
				nonBasics[pCol] = true;
			}
			else if ((col_status[pCol] == HighsBasisStatus::LOWER
			|| col_status[pCol] == HighsBasisStatus::UPPER) 
			&& nonBasics[pCol])
				alpBasis.col_status[i] = HighsBasisStatus::BASIC;
			else
				alpBasis.col_status[i] = col_status[pCol];
		}
		else
			alpBasis.col_status[i] = HighsBasisStatus::NONBASIC;
	}
	for (int i = 0; i < numRow_; ++i){
		if (i < numRow_ - linkingPairs.size()){
			int rep = C[i + numCol].front() - numCol;
			int pCol = previousRowColoring[rep];
			if ((row_status[pCol - numCol] == HighsBasisStatus::LOWER
			|| row_status[pCol - numCol] == HighsBasisStatus::UPPER)
			&& !nonBasics[pCol]){
				alpBasis.row_status[i] = row_status[pCol - numCol];
				nonBasics[pCol] = true;
			}
			else if ((row_status[pCol - numCol] == HighsBasisStatus::LOWER
			|| row_status[pCol - numCol] == HighsBasisStatus::UPPER)
			&& nonBasics[pCol])
				alpBasis.row_status[i] = HighsBasisStatus::BASIC;
			else
				alpBasis.row_status[i] = row_status[pCol - numCol];
			
		}
		else
			alpBasis.row_status[i] = HighsBasisStatus::NONBASIC;
	}
}

void HighsAggregate::aggregateCVector(){
	colCost_.assign(numCol_, 0);
	for (int i = 0; i < numCol_; ++i){
		for (int j = 0; j < C[i].size(); ++j)
			colCost_[i] += colCost[C[i][j]];
	}
}

void HighsAggregate::appendLinkersToAMatrix(vector<double>& row, int rowIdx, int rIdx){
	for (int i = 0; i < row.size(); ++i){
		if (row[i]){
			int colStartIdx = Astart_[i + 1];
			for (int j = i + 1; j < Astart_.size(); ++j)
				Astart_[j]++;
			Aindex_.insert(Aindex_.begin() + colStartIdx, rowIdx);
			Avalue_.insert(Avalue_.begin() + colStartIdx, row[i]);
		}
	}
}

void HighsAggregate::appendLinkersToRowRhs(int rowIdx){
	rowLower_[rowIdx] = 0;
	rowUpper_[rowIdx] = 0;
}

void HighsAggregate::appendLinkersToColBounds(int rIdx){
	colLower_[rIdx] = 0;
	colUpper_[rIdx] = 0;
}

void HighsAggregate::setInitialRhsAndBounds(){
	for (int i = 0; i < numRow_; ++i){
		int rep = C[i + numCol].front();
		rowLower_[i] = rowLower[rep - numCol];
		rowUpper_[i] = rowUpper[rep - numCol];
	}
	for (int i = 0; i < numCol_; ++i){
		int rep = C[i].front();
		colLower_[i] = colLower[rep];
		colUpper_[i] = colUpper[rep];
	}
}

void HighsAggregate::setRhsAndBounds(){
	for (int i = 0; i < realNumRow_; ++i){
		int rep = C[i + numCol].front() - numCol;
		int pCol = previousRowColoring[rep] - numCol;
		if (row_value[pCol] == rowLower[rep]){
			rowLower_[i] = rowLower[rep];
			rowUpper_[i] = rowLower[rep];
		}
		else if (row_value[pCol] == rowUpper[rep]){
			rowLower_[i] = rowUpper[rep];
			rowUpper_[i] = rowUpper[rep];
		}
		else{
			rowLower_[i] = rowLower[rep];
			rowUpper_[i] = rowUpper[rep];
		}
	}
	for (int i = 0; i < realNumCol_; ++i){
		int rep = C[i].front();
		int pCol = previousColumnColoring[rep];
		if (col_value[pCol] == colUpper[rep]){
			colUpper_[i] = colUpper[rep];
			colLower_[i] = colUpper[rep];
		}
		else if (col_value[pCol] == colLower[rep]){
			colUpper_[i] = colLower[rep];
			colLower_[i] = colLower[rep];
		}
		else{
			colLower_[i] = colLower[rep];
			colUpper_[i] = colUpper[rep];
		}
	}
}

vector<double> HighsAggregate::rowCoeff(int column){
	vector<double> coeffs(realNumRow_, 0);
	int i, j, var, rep, vWeight = 0, colWeight = 0;
	for (i = 0; i < realNumRow_; ++i){
		rep = C[i + numCol].front();
		for (j = 0; j < adjListLab[rep].size(); ++j){
			var = adjListLab[rep][j];
			if (color[var] == column){
				colWeight += adjListWeight[rep][j]; 
			}
		}
		coeffs[i] = colWeight;
		colWeight = 0;
	}
	return coeffs;
}

