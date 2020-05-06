#include "Aggregate.h"
using namespace std;

HighsAggregate::HighsAggregate(HighsLp& lp, const HighsEquitable& ep, HighsSolution& solution, HighsBasis& basis, bool flag){
	// From the original lp
	numRow = lp.numRow_;
	numCol = lp.numCol_;
	impliedNumRow = lp.addNumRow_;
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
	//Equitable partition info
	previousRowColoring.assign(ep.previousRowColoring.begin(), ep.previousRowColoring.end());
	previousColumnColoring.assign(ep.previousColumnColoring.begin(), ep.previousColumnColoring.end());
	C.assign(ep.C.begin(), ep.C.end());
	prevC.assign(ep.prevC.begin(), ep.prevC.end());
	color.assign(ep.color.begin(), ep.color.end());
	adjListLab.assign(ep.adjListLab.begin(), ep.adjListLab.end());
	adjListWeight.assign(ep.adjListWeight.begin(), ep.adjListWeight.end());
	linkingPairs.assign(ep.linkingPairs.begin(), ep.linkingPairs.end());
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
	aggregateAMatrix();
	createImpliedRows(lp);
}


void HighsAggregate::createImpliedRows(HighsLp& lp){
	int i,j, dom, slav;
	if (equalColors.size() && !lp.addARstart_.size())
		lp.addARstart_.push_back(0);
	for (i = 0; i < equalColors.size(); ++i){
		dom = equalColors[i].first;
		slav = equalColors[i].second;
		lp.addNumRow_++;
		for (j = 0; j < C[dom].size(); ++j){
			lp.addARvalue_.push_back((double)1/C[dom].size());
			lp.addARindex_.push_back(C[dom][j]);
		}
		for (j = 0; j < C[slav].size(); ++j){
			lp.addARvalue_.push_back((double)-1/C[slav].size());
			lp.addARindex_.push_back(C[slav][j]);
		}
		lp.addARstart_.push_back(lp.addARvalue_.size());
	}
}

HighsLp& HighsAggregate::getAlp(){
	alp.numRow_ = (numRow_);
	alp.numCol_ = (numCol_);
	alp.numInt_ = 0;
	alp.nnz_ = Avalue_.size();
	alp.linkers.assign(linkers.begin(), linkers.end());
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
	return alp;
}

void HighsAggregate::aggregateAMatrix(){
	int i, j, rep, idx;
	numCol_ = 0;
	numRow_ = 0;
	for (i = 0; i < C.size(); ++i){
		if (C[i].size() && i < numCol){
			numCol_++;
			realNumCol_++;
		}
		else if (C[i].size()){
			numRow_++;
			realNumRow_++;
			rep = C[i].front();
			idx = previousRowColoring[rep - numCol] - numCol;
			if (idx < 0) continue;
			numSplits[previousRowColoring[rep - numCol] - numCol]++;
		}
	}
	numTot_ = numCol_ + numRow_;
	realNumTot_ = numTot_;
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
	if (!flag_){
		findPreviousBasisForRows();
		findPreviousBasisForColumns();
		setAggregateRealRowsRhs();
		setAggregateRealColsBounds();
		aggregateCVector();
		return;
	}
	else{
		findPreviousBasisForRows();
		findPreviousBasisForColumns();
		setAggregateRealRowsRhs();
		setAggregateRealColsBounds();
		eraseLinkersIfNotNeeded();
		//getAggImpliedRows();
		findMissingBasicColumns();
	}
}

void HighsAggregate::aggregateCVector(){
	colCost_.assign(numCol_, 0);
	for (int i = 0; i < numCol_; ++i){
		int rep = C[i].front();
		colCost_[i] = colCost[rep];
	}
}

void HighsAggregate::appendLinkersToAMatrix(vector<double>& row){
	linkers.push_back(numCol_);
	for (int i = 0; i < row.size(); ++i){
		if (row[i]){
			int colStartIdx = Astart_[i + 1];
			for (int j = i + 1; j <= numCol_; j++)
				Astart_[j]++;
			Aindex_.insert(Aindex_.begin() + colStartIdx, numRow_);
			Avalue_.insert(Avalue_.begin() + colStartIdx, row[i]);
		}
	}
	Astart_.push_back(Astart_[numCol_] + 1);
	Avalue_.push_back(-1);
	Aindex_.push_back(numRow_);
	numCol_++;
	numRow_++;
}

void HighsAggregate::appendLinkersToRowRhs(){
	rowLower_.push_back(0);
	rowUpper_.push_back(0);
}

void HighsAggregate::appendLinkersToColBounds(){
	colLower_.push_back(0);
	colUpper_.push_back(0);
}

void HighsAggregate::createRowWiseAMatrix(){
    int AcountX = Astart_[numCol_];
    ARindex_.resize(AcountX);
    ARvalue_.resize(AcountX);
    // Build row copy - pointers
    ARstart_.assign(numRow_ + 1, 0);
    AR_Nend_.assign(numRow_, 0);
    for (int k = 0; k < AcountX; ++k)
        AR_Nend_[Aindex_[k]]++;
    for (int i = 1; i <= numRow_; ++i)
        ARstart_[i] = ARstart_[i - 1] + AR_Nend_[i - 1];
    for (int i = 0; i < numRow_; ++i)
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
	ARstartSub_ = ARstart_;
	ARindexSub_ = ARindex_;
	ARvalueSub_ = ARvalue_;
}

void HighsAggregate::findMissingBasicColumns(){
	createRowWiseAMatrix();
	int i, j, parent, child, currSize, prevSize;
	currSize = linkingPairs.size();
	prevSize = 0;
	while (currSize != prevSize && linkingPairs.size()){
		prevSize = currSize;
		for (i = 0; i < partsForGS.size(); ++i){
			//if (numSplits[partsForGS[i] - numCol] <= 1) continue;
			doGramSchmidt(partsForGS[i]);
			if (!linkingPairs.size())
				break;
		}
		currSize = linkingPairs.size();
	}
	for (i = 0; i < linkingPairs.size(); ++i){
		parent = linkingPairs[i].first;
		child = linkingPairs[i].second;
		vector<double> linkRow(realNumCol_,  0);
		linkRow[parent] = 1;
		linkRow[child] = -1;
		numLinkers_++;
		appendLinkersToAMatrix(linkRow);
		appendLinkersToColBounds();
		appendLinkersToRowRhs();
	}
	// for (i = numActiveRows_; i < numActiveRows_ + numActiveBounds_; ++i){
	// 	if (dependanceCheck(QRmat[i]))
	// 		startingBasicColumns_.push_back(potentialBasicColumns_[i - numActiveRows_]);
	// }
	numTot_ = numRow_ + numCol_;
	colCost_.assign(numCol_, 0);
}

void HighsAggregate::doGramSchmidt(int oldPart){
	int i, j, k, x, rep, prevColor, domLink, slavLink;
	int rowIdx = 0;
	int numNonLinkRows = numSplits[oldPart - numCol];
	int impliedRowsIdx = aggImpliedRows.size() + numNonLinkRows;
	int numRowsToTest = linkingPairs.size() + impliedRowsIdx;
	vector<int> currentRows;
	vector<vector<double> > AM(numRowsToTest, vector<double>(numCol_, 0.0));
	for (i = 0; i < numRow_; ++i){
		rep = C[i + numCol].front();
		prevColor = previousRowColoring[rep - numCol];
		if ((activeConstraints_[i] && prevColor == oldPart)){
			currentRows.push_back(i);
			for (j = ARstartSub_[i]; j < ARstartSub_[i + 1]; ++j){
				AM[rowIdx][ARindexSub_[j]] = activeBounds_[ARindexSub_[j]] ? 0 : ARvalueSub_[j];
			}
			rowIdx++;
		}
	}
	for (i = 0; i < aggImpliedRows.size(); ++i){
		for (j = 0; j < aggImpliedRows[i].size(); ++j){
			AM[rowIdx][j] = aggImpliedRows[i][j];
		}
		rowIdx++;
	}
	for (i = 0; i < linkingPairs.size(); ++i){
		AM[rowIdx][linkingPairs[i].first] = 1;
		AM[rowIdx][linkingPairs[i].second] = -1;
		rowIdx++;
	}
	vector<vector<double> > QRmat = AM;
	for (i = 0; i < QRmat.size(); ++i){
		cout << "[ ";
		for (j = 0; j < QRmat[i].size(); ++j){
			cout << QRmat[i][j] << " ";
		}
		cout << "]" << endl;
	}
	QR.gramSchmidt(QRmat);
	cout << endl;
	for (i = 0; i < QRmat.size(); ++i){
		cout << "[ ";
		for (j = 0; j < QRmat[i].size(); ++j){
			cout << QRmat[i][j] << " ";
		}
		cout << "]" << endl;
	}
	cin.get();
	for (i = 0; i < numNonLinkRows; ++i)
		if (dependanceCheck(QRmat[i]))
			startingBasicRows_.push_back(currentRows[i]);
	for (i = impliedRowsIdx; i < numRowsToTest; ++i)
		if (dependanceCheck(QRmat[i]))
			linkIsNeeded[i - impliedRowsIdx] = false;
	i = 0;
	while (i < linkingPairs.size()){
		if (!linkIsNeeded[i]){
			domLink = linkingPairs[i].first;
			slavLink = linkingPairs[i].second;
			editRowWiseMatrix(domLink, slavLink, i);
			equalColors.push_back(pair<int, int>(domLink, slavLink));
			linkingPairs.erase(linkingPairs.begin() + i);
			linkIsNeeded.erase(linkIsNeeded.begin() + i);
		}
		else
			++i;
	}
}

void HighsAggregate::editRowWiseMatrix(int domLink, int slavLink, int idx){
	int i, j, k;
	bool isDom, isSlav;
	vector<int> ARstartTemp;
	vector<int> ARindexTemp;
	vector<double> ARvalueTemp;
	vector<bool> need(numRow_, false);
	vector<double> sub;
	vector<double> coeff(realNumCol_, 0);
	ARstartTemp.push_back(0);
	for (i = 0; i < numRow_; ++i){
		isDom = false, isSlav = false;
		for (j = ARstartSub_[i]; j < ARstartSub_[i + 1]; ++j){
			if (ARindexSub_[j] == domLink) isDom = true;
			if (ARindexSub_[j] == slavLink) isSlav = true;
		}
		if (isDom && isSlav)
			need[i] = true;
	}
	for (i = 0; i < numRow_; ++i){
		for (j = ARstartSub_[i]; j < ARstartSub_[i + 1]; ++j){
			if (ARindexSub_[j] != slavLink){
				ARindexTemp.push_back(ARindexSub_[j]);
				ARvalueTemp.push_back(ARvalueSub_[j]);
			}
			else{
				if (need[i])
					sub.push_back(ARvalueSub_[j]);
				else{
					ARindexTemp.push_back(ARindexSub_[j]);
					ARvalueTemp.push_back(ARvalueSub_[j]);
				}
			}
		}
		ARstartTemp.push_back(ARvalueTemp.size());
	}
	k= 0;	
	for (i = 0; i < numRow_; ++i){
		if (need[i]){
			for (j = ARstartTemp[i]; j < ARstartTemp[i + 1]; ++j){
				if (ARindexTemp[j] == domLink){
					ARvalueTemp[j] += sub[k];
					++k;
					break;
				}
			}
		}
	}
	for (i = 0; i < idx; ++i){
		coeff[linkingPairs[i].first] = 1;
		coeff[linkingPairs[i].second] = -1;
	}
	coeff[domLink] = 1;
	coeff[slavLink] = -1;
	aggImpliedRows.push_back(coeff);
	ARstartSub_ = ARstartTemp;
	ARvalueSub_ = ARvalueTemp;
	ARindexSub_ = ARindexTemp;
}

HighsBasis& HighsAggregate::getAlpBasis(){
	alpBasis.col_status.resize(numCol_);
	alpBasis.row_status.resize(numRow_);
	int previousColumnColor = -1;
	int previousRowColor = -1;
	int rep = -1;
	int numBasicVars = 0;
	for (int i = 0; i < numCol_; ++i){
		if (i < numCol_ - numLinkers_){
			rep = C[i].front();
			previousColumnColor = previousColumnColoring[rep];
			alpBasis.col_status[i] = previousColumnInfo[previousColumnColor];
			if (previousColumnInfo[previousColumnColor] == HighsBasisStatus::BASIC)
				numBasicVars++;
		}
		else
			alpBasis.col_status[i] = HighsBasisStatus::LOWER;
	}
	for (int i = 0; i < numRow_; ++i){
		if (i < numRow_ - numLinkers_){
			rep = C[i + numCol].front() - numCol;
			previousRowColor = previousRowColoring[rep];
			alpBasis.row_status[i] = previousRowInfo[previousRowColor];
			if (previousRowInfo[previousRowColor] == HighsBasisStatus::BASIC)
				numBasicVars++;
		}
		else
			alpBasis.row_status[i] = HighsBasisStatus::LOWER;
	}
	for (int i = 0; i < startingBasicColumns_.size(); ++i){
		if (numBasicVars < numRow_){
			alpBasis.col_status[startingBasicColumns_[i]] = HighsBasisStatus::BASIC;
			numBasicVars++;
		}
		else
			break;
	}
	for (int i = 0; i < startingBasicRows_.size(); ++i){
		if (numBasicVars < numRow_){
			alpBasis.row_status[startingBasicRows_[i]] = HighsBasisStatus::BASIC;
			numBasicVars++;
		}
		else 
			break;
	}
	return alpBasis;
}

bool HighsAggregate::dependanceCheck(vector<double> &v){
	bool cond = true;
	for (int i = 0; i < v.size(); ++i){
		if (fabs(v[i]) > 1e-6)
			cond = false;
		else
			v[i] = 0;
	}
	return cond;
}

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
	int rep, previousColumnColor, oldRep;
	double lb, ub;
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
		else if (previousColumnInfo[previousColumnColor] == HighsBasisStatus::BASIC &&
			     previousColumnValue[previousColumnColor] == colLower[prevC[previousColumnColor].front()]){
			colUpper_.push_back(previousColumnValue[previousColumnColor]);
			colLower_.push_back(previousColumnValue[previousColumnColor]);
			activeBounds_[i] = true;
			potentialBasicColumns_.push_back(i);
			numActiveBounds_++;
		}
		else if (previousColumnInfo[previousColumnColor] == HighsBasisStatus::BASIC &&
			     previousColumnValue[previousColumnColor] == colUpper[prevC[previousColumnColor].front()]){
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

bool HighsAggregate::varIsBounded(pair<int, int> link){
	int rep = C[link.first].front();
	int previousColumnColor = previousColumnColoring[rep];
	if (previousColumnInfo[previousColumnColor] == HighsBasisStatus::LOWER ||
		previousColumnInfo[previousColumnColor] == HighsBasisStatus::UPPER)
		return true;
	else if (previousColumnInfo[previousColumnColor] == HighsBasisStatus::BASIC &&
			 previousColumnValue[previousColumnColor] == colLower[prevC[previousColumnColor].front()])	
		return true;
	else if (previousColumnInfo[previousColumnColor] == HighsBasisStatus::BASIC &&
			 previousColumnValue[previousColumnColor] == colUpper[prevC[previousColumnColor].front()])
		return true;
	return false;
}

void HighsAggregate::eraseLinkersIfNotNeeded(){
	int i = 0;
	while (i < linkingPairs.size()){
		if (varIsBounded(linkingPairs[i]))
			linkingPairs.erase(linkingPairs.begin() + i);
		else
			++i;
	}
	linkIsNeeded.assign(linkingPairs.size(), true);
}

void HighsAggregate::findPreviousBasisForRows(){
	int i, rep, previousRowColor;
	double rhs;
	for (i = 0; i < row_value.size(); ++i){
		if (prevC[i + numCol].size()){
			rep = prevC[i + numCol].front();
			rhs = min(fabs(rowUpper[rep - numCol]), fabs(rowLower[rep - numCol]));
			if (row_status[i] == HighsBasisStatus::LOWER){
				previousRowInfo.insert(pair<int, HighsBasisStatus>(i + numCol, HighsBasisStatus::LOWER));
				previousRowValue.insert(pair<int, double>(i + numCol, row_value[i]));
				if (!(prevC[i + numCol].size() == 1))
					partsForGS.push_back(i + numCol);
			}
			else if (rhs == row_value[i]){
				previousRowInfo.insert(pair<int, HighsBasisStatus>(i + numCol, col_status[i]));
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
}

void HighsAggregate::findPreviousBasisForColumns(){
	int i, rep;
	double lb, ub;
	for (i = 0; i < col_value.size(); ++i){
		if (prevC[i].size()){
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
}

vector<double> HighsAggregate::rowCoeff(int column){
	vector<double> coeffs(numRow_, 0);
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

vector<double> HighsAggregate::aggregateImpliedRow(int impliedRow){
	int i, j, weight = 0;
	vector<double> coeff(realNumCol_, 0);
	for (i = impliedARstart[impliedRow]; i < impliedARstart[impliedRow + 1]; ++i){
		for (j = 0; j < realNumCol_; ++j){
			if (color[impliedARindex[i]] == j){
				cout << "index: " << impliedARindex[i] << " color: " << color[impliedARindex[i]] << endl;
				coeff[j] += impliedARvalue[i];
			}
		}
	}
	cout << "coeff" << endl;
	return coeff;
}

void HighsAggregate::getAggImpliedRows(){
	int i;
	for (i = 0; i < impliedNumRow; ++i){
		aggImpliedRows.push_back(aggregateImpliedRow(i));
	}
}

