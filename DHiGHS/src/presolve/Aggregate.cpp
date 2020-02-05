#include "Aggregate.h"
using namespace std;

HighsAggregate::HighsAggregate(const HighsLp& lp, const HighsEquitable& ep, HighsSolution* solution, HighsBasis* basis, bool flag){
	// From the original lp 
	numRow = lp.numRow_;
	numCol = lp.numCol_;
	rowLower.assign(lp.rowLower_.begin(), lp.rowLower_.end());
	rowUpper.assign(lp.rowUpper_.begin(), lp.rowUpper_.end());
	colUpper.assign(lp.colUpper_.begin(), lp.colUpper_.end());
	colLower.assign(lp.colLower_.begin(), lp.colUpper_.end());
	colCost.assign(lp.colCost_.begin(), lp.colCost_.end());
	Avalue.assign(lp.Avalue_.begin(), lp.Avalue_.end());
	Aindex.assign(lp.Aindex_.begin(), lp.Aindex_.end());
	Astart.assign(lp.Astart_.begin(), lp.Astart_.end());
	//Equitable partition info
	previousRowColoring.assign(ep.previousRowColoring.begin(), ep.previousRowColoring.end());
	previousColumnColoring.assign(ep.previousColumnColoring.begin(), ep.previousColumnColoring.end());
	C.assign(ep.C.begin(), ep.C.end());
	prevC.assign(ep.prevC.begin(), ep.prevC.end());
	color.assign(ep.color.begin(), ep.color.end());
	if (!flag_){
		adjListLab.assign(ep.adjListLab.begin(), ep.adjListLab.end());
		adjListWeight.assign(ep.adjListWeight.begin(), ep.adjListWeight.end());
	}
	linkingPairs.assign(ep.linkingPairs.begin(), ep.linkingPairs.end());
	// Used for setting active set
	activeBounds_.assign(numCol, false);
	activeConstraints_.assign(numRow, false);
	// Preivous solution
	col_value = (solution->col_value);
	row_value = (solution->row_value);
	// Previous basis
	col_status = (basis->col_status);
	row_status = (basis->row_status);
	// flag status to find linkers or not 
	flag_ = flag;
	aggregateAMatrix();
}

// void HighsAggregate::setup(HighsLp& lp){
// 	// From the original lp 
// 	numRow = lp.numRow_;
// 	numCol = lp.numCol_;
// 	rowLower.assign(lp.rowLower_.begin(), lp.rowLower_.end());
// 	rowUpper.assign(lp.rowUpper_.begin(), lp.rowUpper_.end());
// 	colUpper.assign(lp.colUpper_.begin(), lp.colUpper_.end());
// 	colLower.assign(lp.colLower_.begin(), lp.colUpper_.end());
// 	colCost.assign(lp.colCost_.begin(), lp.colCost_.end());
// 	Avalue.assign(lp.Avalue_.begin(), lp.Avalue_.end());
// 	Aindex.assign(lp.Aindex_.begin(), lp.Aindex_.end());
// 	Astart.assign(lp.Astart_.begin(), lp.Astart_.end());
// }

HighsLp& HighsAggregate::getAlp(){
	alp.numRow_ = (numRow_);
	alp.numCol_ = (numCol_);
	alp.numInt_ = 0;
	alp.nnz_ = Avalue_.size();
	alp.Astart_.assign(Astart_.begin(), Astart_.end());
	alp.Aindex_.assign(Aindex_.begin(), Aindex_.end());
	alp.Avalue_.assign(Avalue_.begin(), Avalue_.end());
	alp.colUpper_.assign(colUpper_.begin(), colUpper_.end());
	alp.colLower_.assign(colLower_.begin(), colLower_.end());
	alp.rowUpper_.assign(rowUpper_.begin(), rowUpper_.end());
	alp.rowLower_.assign(rowLower_.begin(), rowLower_.end());
	alp.colCost_ = (colCost_);
	alp.model_name_ = (model_name_);
	alp.lp_name_ = (lp_name_);
	alp.sense_ = 1;
	alp.offset_ = 0;
	return alp;
}

void HighsAggregate::aggregateAMatrix(){
	int i, j;
	numCol_ = 0;
	numRow_ = 0;
	for (i = 0; i < C.size(); ++i){
		if (C[i].size() && i < numCol)
			numCol_++;
		else if (C[i].size())
			numRow_++;
	}
	numTot_ = numCol_ + numRow_;
	colCost_.assign(numCol_, 1);
	Astart_.push_back(0);
	for (i = 0; i < numCol_; ++i){
		doubleVec coeff = rowCoeff(i);
		for (j = 0; j < coeff.size(); ++j){
			if (coeff[j]){
				Avalue_.push_back(coeff[j]);
				Aindex_.push_back(j);
			}
		}
		Astart_.push_back(Avalue_.size());
	}
	if (!flag_){
		findPreviousBasisForRows();
		findPreviousBasisForColumns();
		setAggregateRealRowsRhs();
		setAggregateRealColsBounds();
		createRowWiseAMatrix();
		// createAlp();
		return;
	}
	// else{
	// 	findPreviousBasisForRows();
	// 	findPreviousBasisForColumns();
	// 	setAggregateRealRowsRhs();
	// 	setAggregateRealColsBounds();
	// 	findMissingBasicColumns();
	// 	//appendLinkerSubMatrix();
	// }
}

void HighsAggregate::createRowWiseAMatrix(){
    int AcountX = Astart_[numCol_];
    ARindex_.resize(AcountX);
    ARvalue_.resize(AcountX);
    // Build row copy - pointers
    ARstart_.resize(numRow_ + 1);
    AR_Nend_.assign(numRow_, 0);
    for (int k = 0; k < AcountX; ++k)
        AR_Nend_[Aindex_[k]]++;
    ARstart_[0] = 0;
    for (int i = 1; i <= numRow_; ++i)
        ARstart_[i] = ARstart_[i - 1] + AR_Nend_[i - 1];
    for (int i = 0; i < numRow; ++i)
        AR_Nend_[i] = ARstart_[i];
    // Build row copy - elements
    for (int iCol = 0; iCol < numCol_; ++iCol) {
        for (int k = Astart_[iCol]; k < Astart_[iCol + 1]; ++k) {
            int iRow = Aindex_[k];
            int iPut = AR_Nend_[iRow]++;
            ARindex_[iPut] = iCol;
            ARvalue_[iPut] = Avalue_[k];
		}
	}
}

// void HighsAggregate::findMissingBasicColumns(){
// 	createRowWiseAMatrix();
// 	int i, j;
// 	int rowIdx = 0;
// 	int numRowsToTest = linkingPairs.size() + numActiveBounds_ + numActiveRows_;
// 	cout << "numRowsToTest: " << numRowsToTest << endl;
// 	vector<vector<double> > AM(numRowsToTest, vector<double>(numCol_, 0.0));
// 	cout << "matrix initilization" << endl;
// 	// for (i = 0; i < AM.size(); ++i){
// 	// 	cout << "row: " << i << " {";
// 	// 	for (j = 0; j < AM[i].size(); ++j){
// 	// 		cout << AM[i][j] << " ";
// 	// 	}
// 	// 	cout << "}" << endl;
// 	// }
// 	for (i = 0; i < numRow_; ++i){
// 		if (activeConstraints_[i]){
// 			for (j = ARstart_[i]; j < ARstart_[i + 1]; ++j)
// 				AM[rowIdx][ARindex_[j]] = ARvalue_[j];
// 			rowIdx++;
// 		}
// 	}
// 	for (i = 0; i < numCol_; ++i){
// 		if (activeBounds_[i]){
// 			AM[rowIdx][i] = 1;
// 			rowIdx++;
// 		}
// 	}
// 	for (i = 0; i < linkingPairs.size(); ++i){
// 		AM[rowIdx][linkingPairs[i].first] = 1;
// 		AM[rowIdx][linkingPairs[i].second] = -1;
// 		rowIdx++;
// 	}
// 	cout << "matrix before gramSchmidt" << endl;
// 	// for (i = 0; i < AM.size(); ++i){
// 	// 	cout << "row: " << i << " {";
// 	// 	for (j = 0; j < AM[i].size(); ++j){
// 	// 		cout << AM[i][j] << " ";
// 	// 	}
// 	// 	cout << "}" << endl;
// 	// }
// 	QR.gramSchmidt(AM);
// 	cout << "matrix after gramSchmidt" << endl;
// 	// for (i = 0; i < AM.size(); ++i){
// 	// 	cout << "row: " << i << " {";
// 	// 	for (j = 0; j < AM[i].size(); ++j){
// 	// 		cout << AM[i][j] << " ";
// 	// 	}
// 	// 	cout << "}" << endl;
// 	// }
// 	for (i = numActiveRows_ + numActiveBounds_; i < numRowsToTest; ++i)
// 		numLinkers_ += !dependanceCheck(AM[i]);
// 	for (i = numActiveRows_; i < numActiveRows_ + numActiveBounds_; ++i){
// 		if (dependanceCheck(AM[i]))
// 			startingBasicColumns_.push_back(potentialBasicColumns_[i - numActiveRows_]);
// 	}
// 	for (i = 0; i < numActiveRows_; ++i){
// 		if (dependanceCheck(AM[i]))
// 			startingBasicRows_.push_back(potentialBasicRows_[i]);
// 	}
// 	numRow_ += numLinkers_;
// 	numCol_ += numLinkers_;
// 	numTot_ = numRow_ + numCol_;
// }

// HighsBasis& HighsAggregate::getAlpBasis(){
// 	alpBasis.col_status.resize(numCol_);
// 	alpBasis.row_status.resize(numRow_);
// 	int previousColumnColor = -1;
// 	int previousRowColor = -1;
// 	int rep = -1;
// 	for (int i = 0; i < numCol_; ++i){
// 		if (i < numCol_ - numLinkers_){
// 			rep = C[i].front();
// 			previousColumnColor = previousColumnColoring[rep];
// 			alpBasis.col_status[i] = previousColumnInfo[previousColumnColor];
// 		}
// 		else{
// 			alpBasis.col_status[i] = HighsBasisStatus::LOWER;
// 			colLower_.push_back(0);
// 			colUpper_.push_back(0);
// 		}
// 		// cout << "colLower: " << colLower_[i] << endl;
// 		// cout << "colUpper: " << colUpper_[i] << endl;
// 		// cin.get();
// 	}
// 	for (int i = 0; i < numRow_; ++i){
// 		if (i < numRow_ - numLinkers_){
// 			rep = C[i + numCol].front() - numCol;
// 			previousRowColor = previousRowColoring[rep];
// 			alpBasis.row_status[i] = previousRowInfo[previousRowColor];
// 		}
// 		else{
// 			alpBasis.row_status[i] = HighsBasisStatus::LOWER;
// 			rowLower_.push_back(0);
// 			rowUpper_.push_back(0);
// 		}
// 		cout << "rowLower: " << rowLower_[i] << endl;
// 		cout << "rowUpper: " << rowUpper_[i] << endl;
// 		cin.get();

// 	}
// 	for (int i = 0; i < startingBasicColumns_.size(); ++i)
// 		alpBasis.col_status[startingBasicColumns_[i]] = HighsBasisStatus::BASIC;
// 	for (int i = 0; i < startingBasicRows_.size(); ++i)
// 		alpBasis.row_status[startingBasicRows_[i]] = HighsBasisStatus::BASIC;
// }

// bool HighsAggregate::dependanceCheck(vector<double> &v){
// 	for (int i = 0; i < v.size(); ++i){
// 		if (fabs(v[i]) > 1e-5)
// 			return false;
// 	}
// 	return true;
// }

void HighsAggregate::setAggregateRealRowsRhs(){
	int previousRowColor;
	int rep;
	int i, j;
	for (i = 0; i < numRow_; ++i){
		rep = C[i + numCol].front() - numCol;
		previousRowColor = previousRowColoring[rep];
		if (previousRowColor == -1){
			rowLower_.push_back(rowLower[rep]);
			rowUpper_.push_back(rowUpper[rep]);
		}
		else if (previousRowInfo[previousRowColor] == HighsBasisStatus::LOWER){
			rowLower_.push_back(previousRowValue[previousRowColor]);
			rowUpper_.push_back(previousRowValue[previousRowColor]);
			activeConstraints_[i] = true;
			potentialBasicRows_.push_back(i);
			numActiveRows_++;
		}
		else{
			rowLower_.push_back(rowLower[rep]);
			rowUpper_.push_back(rowUpper[rep]);
		}
	}
}

void HighsAggregate::setAggregateRealColsBounds(){
	int i, j;
	int rep, previousColumnColor;
	for (i = 0; i < numCol_; ++i){
		rep = C[i].front();
 		previousColumnColor = previousColumnColoring[rep];
 		if (previousColumnColor == -1){
			colLower_.push_back(colLower[rep]);
			colUpper_.push_back(colUpper[rep]);
		}
		else if (previousColumnInfo[previousColumnColor] == HighsBasisStatus::LOWER ||
				 previousColumnInfo[previousColumnColor] == HighsBasisStatus::UPPER){
			colUpper_.push_back(previousColumnValue[previousColumnColor]);
			colLower_.push_back(previousColumnValue[previousColumnColor]);
			activeBounds_[i] = true;
			potentialBasicColumns_.push_back(i);
			numActiveBounds_++;
		}
		else{
			colUpper_.push_back(colUpper[rep]);
			colLower_.push_back(colLower[rep]);
		}
	}
}

void HighsAggregate::findPreviousBasisForRows(){
	int i, j, rep, previousRowColor;
	double rhs;
	for (i = 0; i < row_value.size(); ++i){
		rep = prevC[i + numCol].front();
		rhs = min(fabs(rowUpper[rep]), fabs(rowLower[rep]));
		if (row_status[i] == HighsBasisStatus::LOWER){
			previousRowInfo.insert(pair<int, HighsBasisStatus>(i + numCol, HighsBasisStatus::LOWER));
			previousRowValue.insert(pair<int, double>(i + numCol, row_value[i]));
		}
		else if (rhs == row_value[i]){
			previousRowInfo.insert(pair<int, HighsBasisStatus>(i + numCol, HighsBasisStatus::LOWER));
			previousRowValue.insert(pair<int, double>(i + numCol, row_value[i]));
		}
		else if (row_status[i] == HighsBasisStatus::UPPER){
			previousRowInfo.insert(pair<int, HighsBasisStatus>(i + numCol, HighsBasisStatus::UPPER));
			previousRowValue.insert(pair<int, double>(i + numCol, row_value[i]));
		}
		else if (row_status[i] == HighsBasisStatus::BASIC){
			previousRowInfo.insert(pair<int, HighsBasisStatus>(i + numCol, HighsBasisStatus::BASIC));
			previousRowValue.insert(pair<int, double>(i + numCol, row_value[i]));
		}
		else if (row_status[i] == HighsBasisStatus::NONBASIC){
			previousRowInfo.insert(pair<int, HighsBasisStatus>(i + numCol, HighsBasisStatus::NONBASIC));
			previousRowValue.insert(pair<int, double>(i + numCol, row_value[i]));
		}
		else if (row_status[i] == HighsBasisStatus::SUPER){
			previousRowInfo.insert(pair<int, HighsBasisStatus>(i + numCol, HighsBasisStatus::SUPER));
			previousRowValue.insert(pair<int, double>(i + numCol, row_value[i]));
		}
		else if (row_status[i] == HighsBasisStatus::ZERO){
			previousRowInfo.insert(pair<int, HighsBasisStatus>(i + numCol, HighsBasisStatus::ZERO));
			previousRowValue.insert(pair<int, double>(i + numCol, row_value[i]));
		}
	}
}

void HighsAggregate::findPreviousBasisForColumns(){
	int i, j, rep;
	for (i = 0; i < col_value.size(); ++i){
		if (col_status[i] == HighsBasisStatus::UPPER){
			previousColumnInfo.insert(pair<int, HighsBasisStatus>(i, HighsBasisStatus::UPPER));
			previousColumnValue.insert(pair<int, double>(i, col_value[i]));
		}
		else if (col_status[i] == HighsBasisStatus::LOWER){
			previousColumnInfo.insert(pair<int, HighsBasisStatus>(i, HighsBasisStatus::LOWER));
			previousColumnValue.insert(pair<int, double>(i, col_value[i]));
		}
		else if (col_status[i] == HighsBasisStatus::BASIC){
			previousColumnInfo.insert(pair<int, HighsBasisStatus>(i, HighsBasisStatus::BASIC));
			previousColumnValue.insert(pair<int, double>(i, col_value[i]));
		}
		else if (col_status[i] == HighsBasisStatus::SUPER){
			previousColumnInfo.insert(pair<int, HighsBasisStatus>(i, HighsBasisStatus::SUPER));
			previousColumnValue.insert(pair<int, double>(i, col_value[i]));
		}
		else if (col_status[i] == HighsBasisStatus::ZERO){
			previousColumnInfo.insert(pair<int, HighsBasisStatus>(i, HighsBasisStatus::ZERO));
			previousColumnValue.insert(pair<int, double>(i, col_value[i]));
		}
		else if (col_status[i] == HighsBasisStatus::BASIC){
			previousColumnInfo.insert(pair<int, HighsBasisStatus>(i, HighsBasisStatus::NONBASIC));
			previousColumnValue.insert(pair<int, double>(i, col_value[i]));
		}
	} 
}

doubleVec HighsAggregate::rowCoeff(int column){
	doubleVec coeffs(numRow_, 0);
	int i, j, var, rep, vWeight = 0, colWeight = 0;
	for (i = 0; i < numRow_; ++i){
		rep = C[i + numCol].front();
		for (j = 0; j < adjListLab[rep].size(); ++j){
			var = adjListLab[rep][j];
			if (color[var] == column){
				vWeight = adjListWeight[rep][j];
				colWeight += vWeight;
			}
		}
		coeffs[i] = colWeight;
		colWeight = 0;
	}
	return coeffs;	
}