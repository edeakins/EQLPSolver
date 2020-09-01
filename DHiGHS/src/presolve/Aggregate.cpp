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
	trackRowColors(lp);
	createImpliedRows(lp);
}

void HighsAggregate::createImpliedRows(HighsLp& lp){
	int i,j;
	if (!aggImpliedRows.size())
		return;
	if (aggImpliedRows.size() && !lp.addARstart_.size())
		lp.addARstart_.push_back(0);
	for (i = 0; i < aggImpliedRows.size(); ++i){
		for (j = 0; j < aggImpliedRows[i].size(); ++j){
			if (aggImpliedRows[i][j]){
				lp.addARindex_.push_back(j);
				lp.addARvalue_.push_back(aggImpliedRows[i][j]);
			}
		}
		lp.addARstart_.push_back(lp.addARvalue_.size());
	}
}

void HighsAggregate::trackRowColors(HighsLp& lp){
	lp.rowColor = rowColor_;
	lp.realNumRow_ = numLiftedRow_;
	lp.realNumCol_ = realNumCol_;
	lp.masterIter++;
	lp.activeColorHistory = activeColorHistory_;
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
	alp.numRealRows = numRow_ - numLinkers_;
	alp.numRealCols = numCol_ - numLinkers_;
	return alp;
}

void HighsAggregate::aggregateAMatrix(){
	if (!flag_){
		initialAggregateAMatrix();
		setAggregateRealRowsRhs();
		setAggregateRealColsBounds();
		aggregateCVector();
		return;
	}
	else{
		examinePartition();
		collectRowsForGS();
		collectPartsForGS();
		createRowWiseAMatrix();
		findPreviousBasisForRows();
		findPreviousBasisForColumns();
		setAggregateRealRowsRhs();
		setAggregateRealColsBounds();
		liftTableau();
		eraseLinkersIfNotNeeded();
		initGSMatricesAndGraphs();
		findMissingBasicColumns();
		transposeMatrix();
		// aggregateCVector();
	}
}

void HighsAggregate::examinePartition(){
	int i, j, k, idx, rep, previousColor, oldColor;
	for (i = 0; i < C.size(); ++i){
		if (C[i].size() && i < numCol){
			numCol_++;
			realNumCol_++;
			rep = C[i].front();
			oldColor = previousColumnColoring[rep];
			scale.push_back((double)partSize[i]/previousPartSize[oldColor]);
		}
		else if (C[i].size()){
			numRowGS_++;
			realNumRow_++;
			rep = C[i].front();
			idx = previousRowColoring[rep - numCol] - numCol;
			if (idx < 0) continue;
			numSplits[previousRowColoring[rep - numCol] - numCol]++;
		}
	}
}

void HighsAggregate::collectRowsForGS(){
	int i, j;
	AstartGS_.push_back(0);
	for (i = 0; i < realNumCol_; ++i){
		vector<double> coeff = rowCoeff(i);
		for (j = 0; j < coeff.size(); ++j){
			if (coeff[j]){
				AvalueGS_.push_back(coeff[j]);
				AindexGS_.push_back(j);
			}
		}
		AstartGS_.push_back(AvalueGS_.size());
	}
}

void HighsAggregate::createRowWiseAMatrix(){
    int AcountXSub = AstartGS_[numCol_];
    ARindexGS_.resize(AcountXSub);
    ARvalueGS_.resize(AcountXSub);
    // Build row copy - pointers
    ARstartGS_.assign(realNumRow_ + 1, 0);
    AR_NendGS_.assign(realNumRow_, 0);
    for (int k = 0; k < AcountXSub; ++k)
        AR_NendGS_[AindexGS_[k]]++;
    for (int i = 1; i <= realNumRow_; ++i)
        ARstartGS_[i] = ARstartGS_[i - 1] + AR_NendGS_[i - 1];
    for (int i = 0; i < realNumRow_; ++i)
        AR_NendGS_[i] = ARstartGS_[i];
    // Build row copy - elements
    for (int iCol = 0; iCol < realNumCol_; ++iCol) {
        for (int k = AstartGS_[iCol]; k < AstartGS_[iCol + 1]; ++k) {
            int iRow = AindexGS_[k];
            int iPut = AR_NendGS_[iRow]++;
            ARindexGS_[iPut] = iCol;
            ARvalueGS_[iPut] = AvalueGS_[k];
		}
	}
	for (int i = 0; i < realNumRow_; ++i){
		cout << "R" << i << ": ";
		for (int j = ARstartGS_[i]; j < ARstartGS_[i + 1]; ++j){
			cout << ARvalueGS_[j] << "x_" << ARindexGS_[j] << " ";
		}
		cout << endl;
	}
	cin.get();
}

void HighsAggregate::initialAggregateAMatrix(){
	int i, j;
	numCol_ = 0;
	numRow_ = 0;
	for (i = 0; i < C.size(); ++i){
		if (C[i].size() && i < numCol){
			realNumCol_++;
			numCol_++;
		}
		else if(C[i].size()){
			numLiftedRow_++;
			realNumRow_++;
			numRow_++;
			rowColor_.push_back(i);
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
}

void HighsAggregate::liftTableau(){
	int i, j, k, previousRowColor, previousColColor, rep, varRep, rowCol, start = 0;
	vector<bool> liftedColors(prevC.size(), false);
	for (i = 0; i < prevC.size(); ++i)
		if (prevC[i].size() && i < numCol)
			start++;
	ARstart_.push_back(0);
	for (i = 0; i < realNumRow_; ++i){
		// cout << "row: " << i << endl;
		// cout << "ARstartGS i: " << ARstartGS_[i] << " ARstartGS_ i + 1: " << ARstartGS_[i + 1] << endl;
		// cin.get();

		rep = C[i + numCol].front() - numCol;
		previousRowColor = previousRowColoring[rep];
		// cout << "i: " << i << endl;
		// cout << "rep: " << rep << endl;
		// cout << "prevColor: " << previousRowColor << endl;
		// cout << "rowColor: " << rowCol << endl;
		// cin.get();
		if (!previousRowInfo.count(previousRowColor))
			continue;
		if (previousRowInfo[previousRowColor] == HighsBasisStatus::BASIC){
			// cout << "previous Row Color: " << previousRowColor << endl;
			for (j = ARstartGS_[i]; j < ARstartGS_[i + 1]; ++j){
				ARvalue_.push_back(ARvalueGS_[j]);
				ARindex_.push_back(ARindexGS_[j]);
			}
			ARstart_.push_back(ARvalue_.size());
			numLiftedRow_++;
			numRow_++;
			rowColor_.push_back(color[rep + numCol]);
		}
		else if (liftedColors[previousRowColor]){ 
			// cout << "previousRowColor: " << previousRowColor << " skipped" << endl;
			continue;
		}
		else if (previousRowInfo[previousRowColor] != HighsBasisStatus::BASIC
			&& (ARtableauStart[i] == ARtableauStart[i + 1])){
			// cout << "previous Row Color: " << previousRowColor << endl;
			// cout << "next time this color comes it should be skipped" << endl;
			for (j = ARstartGS_[i]; j < ARstartGS_[i +1]; ++j){
				ARvalue_.push_back(ARvalueGS_[j]);
				ARindex_.push_back(ARindexGS_[j]);
			}
			ARstart_.push_back(ARvalue_.size());
			liftedColors[previousRowColor] = true;
			numLiftedRow_++;
			numRow_++;
			rowColor_.push_back(color[rep + numCol]);
		}
		else{
			for (j = ARtableauStart[i]; j < ARtableauStart[i + 1]; ++j){
				ARvalue_.push_back(ARtableauValue[j] * scale[ARtableauIndex[j]]);
				ARindex_.push_back(ARtableauIndex[j]);
				for (k = start; k < realNumCol_; ++k){
					varRep = C[k].front();
					previousColColor = previousColumnColoring[varRep];
					if (previousColColor == ARtableauIndex[j]){
						ARvalue_.push_back(ARtableauValue[j] * scale[k]);
						ARindex_.push_back(k);
					}
				}
			}
			ARstart_.push_back(ARvalue_.size());
			liftedColors[previousRowColor] = true;
			numLiftedRow_++;
			numRow_++;
			rowColor_.push_back(color[rep + numCol]);
		}
	}
	numRowAfterImp_ = numLiftedRow_;
	if (impliedARstart.size()){
		for (i = 0; i < impliedARstart.size() - 1; ++i){
			for (j = impliedARstart[i]; j < impliedARstart[i + 1]; ++j){
				ARvalue_.push_back(impliedARvalue[j]);
				ARindex_.push_back(impliedARindex[j]);
			}
			ARstart_.push_back(ARvalue_.size());
			numRow_++;
			numRowAfterImp_++;
			rowLower_.push_back(0);
			rowUpper_.push_back(0);
		}
	}	
	// cout << "numRow_: " << numRow_ << endl;
	// cin.get();
}

void HighsAggregate::transposeMatrix(){
	cout << "numRow_: " << numRow_ << endl;
	cout << "final row wise matrix" << endl;
	for (int i = 0; i < numRow_; ++i){
		cout << "row: " << i << endl;
		cout << rowLower_[i] << " <= ";
		for (int j = ARstart_[i]; j < ARstart_[i + 1]; ++j){
			cout << ARvalue_[j] << "x_" << ARindex_[j] << " ";
		}
		cout << " <= " << rowUpper_[i];
		cout << "\n" << endl;
	}
	int i, j, parent, child;
	int _numRow_ = numRow_;
	int _numCol_ = numCol_;
    int AcountX = ARstart_[_numRow_];
    Aindex_.resize(AcountX);
    Avalue_.resize(AcountX);
    // Build row copy - pointers
    Astart_.assign(_numCol_ + 1, 0);
    A_Nend_.assign(_numCol_, 0);
    for (int k = 0; k < AcountX; ++k)
        A_Nend_[ARindex_[k]]++;
    for (int i = 1; i <= _numCol_; ++i)
        Astart_[i] = Astart_[i - 1] + A_Nend_[i - 1];
    for (int i = 0; i < _numCol_; ++i)
        A_Nend_[i] = Astart_[i];
    // Build row copy - elements
    for (int iRow = 0; iRow < _numRow_; ++iRow) {
        for (int k = ARstart_[iRow]; k < ARstart_[iRow + 1]; ++k) {
            int iCol = ARindex_[k];
            int iPut = A_Nend_[iCol]++;
            Aindex_[iPut] = iRow;
            Avalue_[iPut] = ARvalue_[k];
		}
	}
	for (i = 0; i < linkingPairs.size(); ++i){
		if (linkIsNeeded[i]){
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
	}
	numTot_ = numRow_ + numCol_;
	colCost_.assign(numCol_, 0);
}

void HighsAggregate::aggregateCVector(){
	colCost_.assign(numCol_, 0);
	for (int i = 0; i < numCol_; ++i){
		for (int j = 0; j < C[i].size(); ++j)
			colCost_[i] += colCost[C[i][j]];
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

void HighsAggregate::initGSMatricesAndGraphs(){
	int i, j, k, rep, prevCol, nodeCnt = 0;
	/* Start constructing matrices for gram schmidt process.  We go through 
	and separate rows based on their previous colors and then add them to 
	the respective matrix for that previous color class. */
	for (i = 0; i < partsForGS.size(); ++i){
		linkGraph tempGraph;
		vector<vector<double> > tempMat;
		impliedIdx[partsForGS[i]] = 0;
		cLGraphs[partsForGS[i]] = tempGraph; 
		QRMatrices[partsForGS[i]] = tempMat;
		for (j = 0; j < linkingPairs.size(); ++j){
			int x1 = linkingPairs[j].first;
			int x2 = linkingPairs[j].second;
			varToLink[x1].push_back(j + numTot);
			varToLink[x2].push_back(j + numTot);
			cLGraphs[partsForGS[i]].addNode(j + numTot);
		}
		for (j = 0; j < numRowGS_; ++j){
			rep = C[j + numCol].front();
			prevCol = previousRowColoring[rep - numCol];
			if (prevCol = partsForGS[i]){
				impliedIdx[partsForGS[i]]++;
				nodeCnt++;
			}
		}
		cLGraphs[partsForGS[i]].initializeMat(nodeCnt + linkingPairs.size());
		nodeCnt = 0;
	}
	for (i = 0; i < numRowGS_; ++i){
		rep = C[i + numCol].front();
		prevCol = previousRowColoring[rep - numCol];
		if (colorsForGS[prevCol - numCol]){
			vector<double> temp(numCol_ + originalNumLinkers);
			cLGraphs[prevCol].addNode(i + numCol);
			for (j = ARstartGS_[i]; j < ARstartGS_[i + 1]; ++j){
				if (varToLink.count(ARindexGS_[j])){
					for (k = 0; k < varToLink[ARindexGS_[j]].size(); ++k)
						cLGraphs[prevCol].addEdge(i + numCol, varToLink[ARindexGS_[j]][k]);
				}
				temp[ARindexGS_[j]] = ARvalueGS_[j];
			}
			QRMatrices[prevCol].push_back(temp);
		}
	}
	/* Loop over the broken up matrices for gram schmidt and add
	the rows corresponing to linking constraints and linking 
	variables.  The connected component for a particular "old" 
	color class tells us which linking rows should actually be added
	to the particular gram schmidt matrix */
	for (i = 0; i < partsForGS.size(); ++i){
		cLGraphs[partsForGS[i]].connectedComponents();
		vector<int> comp = cLGraphs[partsForGS[i]].components[0];
		for (j = 0; j < comp.size(); ++j){
			if (comp[j] >= numTot){
				vector<double> temp(numCol_ + originalNumLinkers);
				int x1 = linkingPairs[comp[j] - numTot].first;
				int x2 = linkingPairs[comp[j] - numTot].second;
				int r = comp[j] - numTot + numCol_;
				temp[x1] = 1; temp[x2] = -1; temp[r] = -1;
				QRMatrices[partsForGS[i]].push_back(temp);
				linksInMatrix[partsForGS[i]].push_back(comp[j] - numTot);
			}
		}
	}
}

void HighsAggregate::findMissingBasicColumns(){
	int i, j, idx = 0, rep, prevCol, previous = 0, current = numLinkers;
	while (current != previous && numLinkers){
		previous = current;
		if (partsForGS.size() == 1){
			cout << "only one part" << endl;
			for (i = 0; i < partsForGS.size(); ++i){
				doGramSchmidt(partsForGS[i], i);
			}
			break;
		}
		for (i = 0; i < partsForGS.size(); ++i){
			doGramSchmidt(partsForGS[i], i);
			if (!linkingPairs.size())
				break;
		}
		current = numLinkers;
	}
}

void HighsAggregate::doGramSchmidt(int oldPart, int idx){
	cout << "oldPart: " << oldPart << endl;
	cin.get();
	int i, j;
	int testRowStart = impliedIdx[oldPart];
	int linkColIdx = numCol_;
	vector<double> temp;
	vector<int> remainingLinks = linksInMatrix[oldPart];
	vector<vector<double> > AM = QRMatrices[oldPart];
	vector<vector<double> >& AMref = QRMatrices[oldPart];
	cout << "Before GS" << endl;
	cout << " A = [ " << endl; 
	for (i = 0; i < AM.size(); ++i){
		for (j = 0; j < realNumCol_ + originalNumLinkers; ++j){
			cout << AM[i][j] << " ";
		}
		cout << ";" << endl;
	}
	cout << "]" << endl;
	QR.gramSchmidt(AM, linkColIdx);
	cout << endl;
	cout << "After GS" << endl;
	cout << " A = [ " << endl; 
	for (i = 0; i < AM.size(); ++i){
		for (j = 0; j < realNumCol_ + originalNumLinkers; ++j){
			cout << AM[i][j] << " ";
		}
		cout << ";" << endl;
	}
	cout << "]" << endl;
	cin.get();
	for (i = testRowStart; i < AM.size(); ++i){
		int linkIdx = remainingLinks[i - testRowStart] + realNumCol_;
		if (dependanceCheck(AM[i], testRowStart)){ // linkIdx)){
			numLinkers--;
			linkIsNeeded[linkIdx - realNumCol_] = false;
			temp = createImpliedLinkRows(AMref[i], linkIdx);
			AMref.erase // Working Here 
		}
	}
	//QRIndexUpdate[idx] = aggImpliedRows.size() + numNonLinkRows;
	// for (i = 0; i < currentRows.size(); ++i){
	// 	if (i == currentRows.size() - 1){
	// 		cout << "R" << currentRows[i] << " zero out " << endl;
	// 		break;
	// 	}
	// 	cout << "R" << currentRows[i] << ", ";
	// }
	for (i = 0; i < linkingPairs.size(); ++i){
		if (!linkIsNeeded[i] and !linkIsErased[i]){
			// int dom = linkingPairs[i].first;
			// int slav = linkingPairs[i].second;
			// editRowWiseMatrix(dom, slav);
			linkIsErased[i] = true;
		}
	}
}

void HighsAggregate::editRowWiseMatrix(int domLink, int slavLink){
	// int i, j, k;
	// bool isDom, isSlav;
	// vector<int> ARstartTemp;
	// vector<int> ARindexTemp;
	// vector<double> ARvalueTemp;
	// vector<bool> need(numRow_, false);
	// vector<double> sub;
	// ARstartTemp.push_back(0);
	// for (i = 0; i < numRow_; ++i){
	// 	isDom = false, isSlav = false;
	// 	for (j = ARstartGS_[i]; j < ARstartGS_[i + 1]; ++j){
	// 		if (ARindexGS_[j] == domLink) isDom = true;
	// 		if (ARindexGS_[j] == slavLink) isSlav = true;
	// 	}
	// 	if (isDom && isSlav)
	// 		need[i] = true;
	// }
	// for (i = 0; i < numRow_; ++i){
	// 	for (j = ARstartGS_[i]; j < ARstartGS_[i + 1]; ++j){
	// 		if (ARindexGS_[j] != slavLink){
	// 			ARindexTemp.push_back(ARindexGS_[j]);
	// 			ARvalueTemp.push_back(ARvalueGS_[j]);
	// 		}
	// 		else{
	// 			if (need[i])
	// 				sub.push_back(ARvalueGS_[j]);
	// 			else{
	// 				ARindexTemp.push_back(ARindexGS_[j]);
	// 				ARvalueTemp.push_back(ARvalueGS_[j]);
	// 			}
	// 		}
	// 	}
	// 	ARstartTemp.push_back(ARvalueTemp.size());
	// }
	// k= 0;	
	// for (i = 0; i < numRow_; ++i){
	// 	if (need[i]){
	// 		for (j = ARstartTemp[i]; j < ARstartTemp[i + 1]; ++j){
	// 			if (ARindexTemp[j] == domLink){
	// 				ARvalueTemp[j] += sub[k];
	// 				++k;
	// 				break;
	// 			}
	// 		}
	// 	}
	// }	
	// ARstartGS_ = ARstartTemp;
	// ARindexGS_ = ARindexTemp;
	// ARvalueGS_ = ARvalueTemp;
}

vector<double> HighsAggregate::createImpliedLinkRows(vector<double>& linkRow, int linkIdx){
	int i, v1, v2;
	vector<double> temp(linkRow.size(), 0);
	for (i = realNumCol_; i < linkRow.size(); ++i){
		if (linkRow[i]){
			v1 = linkingPairs[i - realNumCol_].first;
			v2 = linkingPairs[i - realNumCol_].second;
			temp[v1] += linkRow[i] * 1;
			temp[v2] += linkRow[i] * -1;
		}
	}
	for (i = 0; i < temp.size(); ++i){
		if (temp[i]){
			ARvalue_.push_back(temp[i]);
			ARindex_.push_back(i);
		}
	}
	ARstart_.push_back(ARvalue_.size());
	rowLower_.push_back(0);
	rowUpper_.push_back(0);
	numRow_++;
	aggImpliedRows.push_back(temp);
	return temp;
	// cout << "implied rows being added" << endl;
	
}

HighsBasis& HighsAggregate::getAlpBasis(){
	alpBasis.col_status.resize(numCol_);
	alpBasis.row_status.resize(numRow_);
	int previousColumnColor, previousRowColor, rep, idx = 0, numBasicVars = 0;
	vector<bool> liftedColors(prevC.size(), false);
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
	for (int i = 0; i < realNumRow_; ++i){
		rep = C[i + numCol].front() - numCol;
		previousRowColor = previousRowColoring[rep];
		// cout << "previousRowColor: " << previousRowColor << endl;
		if (!previousRowInfo.count(previousRowColor))
			continue;
		else if (previousRowInfo[previousRowColor] == HighsBasisStatus::BASIC){ 
			// cout << "row: " << i << " basic" << endl; 
			// cout << "row: " << i << " status: " << (int)previousRowInfo[previousRowColor] << endl; 
			numBasicVars ++;
			alpBasis.row_status[idx] = previousRowInfo[previousRowColor];
			++idx;
		}
		else if (liftedColors[previousRowColor]){
			// cout << "row: " << i << " skipped" << endl; 
			continue;
		}
		else{
			// cout << "row: " << i << " lower" << endl; 
			// cout << "row: " << i << " status: " << (int)previousRowInfo[previousRowColor] << endl; 
			liftedColors[previousRowColor] = true;
			alpBasis.row_status[idx] = previousRowInfo[previousRowColor];
			++idx;
		}
	}
	for (int i = numLiftedRow_; i < numRowAfterImp_; ++i){
		alpBasis.row_status[i] = HighsBasisStatus::LOWER;
	}
	for (int i = numRowAfterImp_; i < numRow_; ++i){
		alpBasis.row_status[i] = HighsBasisStatus::LOWER;
	}
	// for (int i = 0; i < startingBasicColumns_.size(); ++i){
	// 	if (numBasicVars < numRow_){
	// 		alpBasis.col_status[startingBasicColumns_[i]] = HighsBasisStatus::BASIC;
	// 		numBasicVars++;
	// 	}
	// 	else
	// 		break;
	// }
	// for (int i = 0; i < startingBasicRows_.size(); ++i){
	// 	if (numBasicVars < numRow_){
	// 		alpBasis.row_status[startingBasicRows_[i]] = HighsBasisStatus::BASIC;
	// 		numBasicVars++;
	// 	}
	// 	else 
	// 		break;
	// }
	return alpBasis;
}

bool HighsAggregate::dependanceCheck(vector<double> &v){ // int impliedRowsIdx){ //int linkIdx){
	// cout << "linkIdx: " << linkIdx - realNumCol_ << endl;
	bool cond = true;
	for (int i = 0; i < realNumCol_; ++i){
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
	int i, j, rhs;
	vector<bool> liftedColors(prevC.size(), false);
	for (i = 0; i < realNumRow_; ++i){
		rep = C[i + numCol].front() - numCol;
		previousRowColor = previousRowColoring[rep];
		cout << "previousRowColor: " << previousRowColor << endl;
		if (previousRowColor == -1){
			rowLower_.push_back(rowLower[rep]);
			rowUpper_.push_back(rowUpper[rep]);
		}
		else if (!previousRowInfo.count(previousRowColor)){
			// cout << "previous Row Color not there: " << previousRowColor << endl;
			continue;
		}
		else if (previousRowInfo[previousRowColor] == HighsBasisStatus::BASIC){
			rowLower_.push_back(rowLower[rep]);
			rowUpper_.push_back(rowUpper[rep]);
		}
		else if (liftedColors[previousRowColor]){
			//activeConstraints_[i] = true;
			continue;
		}
		else if (previousRowInfo[previousRowColor] != HighsBasisStatus::BASIC
			     && (ARtableauStart[i] == ARtableauStart[i + 1])
			     && i < ARtableauStart.size() - 1){
			rhs = min(fabs(rowLower[rep]), fabs(rowUpper[rep]));
			rowLower_.push_back(rhs);
			rowUpper_.push_back(rhs);
			potentialBasicRows_.push_back(i);
			//activeConstraints_[i] = true;
			numActiveRows_++;
			liftedColors[previousRowColor] = true;
		}
		else{
			rowLower_.push_back(ARreducedRHS[previousRowColor - numCol]);
			rowUpper_.push_back(ARreducedRHS[previousRowColor - numCol]);
			//activeConstraints_[i] = true;
			potentialBasicRows_.push_back(i);
			numActiveRows_++;
			liftedColors[previousRowColor] = true;
		}
	}
	// for (int i = 0; i < rowUpper_.size(); ++i){
	// 	cout << "row: " << i << " upper: " << rowUpper_[i] << endl;
	// 	cout << "row: " << i << " lower: " << rowLower_[i] << endl;
	// }
	// cin.get();
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
			     fabs(previousColumnValue[previousColumnColor] - colLower[prevC[previousColumnColor].front()])
			     <= 1e-5){
			colUpper_.push_back(previousColumnValue[previousColumnColor]);
			colLower_.push_back(previousColumnValue[previousColumnColor]);
			activeBounds_[i] = true;
			potentialBasicColumns_.push_back(i);
			numActiveBounds_++;
		}
		else if (previousColumnInfo[previousColumnColor] == HighsBasisStatus::BASIC &&
			     fabs(previousColumnValue[previousColumnColor] - colUpper[prevC[previousColumnColor].front()])
			      <= 1e-5){
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
			 fabs(previousColumnValue[previousColumnColor] - colLower[prevC[previousColumnColor].front()])
			 <= 1e-5)	
		return true;
	else if (previousColumnInfo[previousColumnColor] == HighsBasisStatus::BASIC &&
			 fabs(previousColumnValue[previousColumnColor] - colUpper[prevC[previousColumnColor].front()])
			 <= 1e-5)
		return true;
	return false;
}

void HighsAggregate::eraseLinkersIfNotNeeded(){
	int i = 0, rep, prevColor;
	while (i < linkingPairs.size()){
		if (varIsBounded(linkingPairs[i])){
			rep = C[linkingPairs[i].first].front();
			prevColor = previousColumnColoring[rep];
				if (previousColumnInfo[prevColor] == HighsBasisStatus::BASIC){
				ARvalue_.push_back(1);
				ARindex_.push_back(linkingPairs[i].first);
				ARvalue_.push_back(-1);
				ARindex_.push_back(linkingPairs[i].second);
				ARstart_.push_back(ARvalue_.size());
				rowLower_.push_back(0);
				rowUpper_.push_back(0);
				numRow_++;
				vector<double> temp(realNumCol_, 0);
				temp[linkingPairs[i].first] = 1; temp[linkingPairs[i].second] = -1;
				aggImpliedRows.push_back(temp);
			}
			linkingPairs.erase(linkingPairs.begin() + i);
		}
		else
			++i;
	}
	linkIsNeeded.assign(linkingPairs.size(), true);
	linkIsErased.assign(linkingPairs.size(), false);
	originalNumLinkers = linkingPairs.size();
	numLinkers = linkingPairs.size();
}

void HighsAggregate::findPreviousBasisForRows(){
	int i, rep, previousRowColor, realRowColor;
	double rhs;
	for (i = 0; i < realNumRow; ++i){
		realRowColor = rowColor[i];
		// cout << "realRowColor: " << realRowColor << endl;
		if (prevC[realRowColor].size()){
			rep = prevC[realRowColor].front();
			// cout << "rep: " << rep << endl;
			previousRowColor = previousRowColoring[rep - numCol];
			// cout << "previousRowColor: " << previousRowColor << endl;
			rhs = min(fabs(rowUpper[rep - numCol]), fabs(rowLower[rep - numCol]));
			if (row_status[i] == HighsBasisStatus::LOWER){
				previousRowInfo.insert(pair<int, HighsBasisStatus>(previousRowColor, HighsBasisStatus::LOWER));
				previousRowValue.insert(pair<int, double>(previousRowColor, row_value[i]));
				// if (!(prevC[i + numCol].size() == 1))
				// 	partsForGS.push_back(previousRowColor);
			}
			else if (rhs == row_value[i]){
				previousRowInfo.insert(pair<int, HighsBasisStatus>(previousRowColor, row_status[i]));
				previousRowValue.insert(pair<int, double>(previousRowColor, row_value[i]));
			}
			else if (row_status[i] == HighsBasisStatus::UPPER){
				previousRowInfo.insert(pair<int, HighsBasisStatus>(previousRowColor, HighsBasisStatus::UPPER));
				previousRowValue.insert(pair<int, double>(previousRowColor, row_value[i]));
			}
			else if (row_status[i] == HighsBasisStatus::BASIC){
				previousRowInfo.insert(pair<int, HighsBasisStatus>(previousRowColor, HighsBasisStatus::BASIC));
				previousRowValue.insert(pair<int, double>(previousRowColor, row_value[i]));
			}
			else if (row_status[i] == HighsBasisStatus::NONBASIC){
				previousRowInfo.insert(pair<int, HighsBasisStatus>(previousRowColor, HighsBasisStatus::NONBASIC));
				previousRowValue.insert(pair<int, double>(previousRowColor, row_value[i]));
			}
			else if (row_status[i] == HighsBasisStatus::SUPER){
				previousRowInfo.insert(pair<int, HighsBasisStatus>(previousRowColor, HighsBasisStatus::SUPER));
				previousRowValue.insert(pair<int, double>(previousRowColor, row_value[i]));
			}
			else if (row_status[i] == HighsBasisStatus::ZERO){
				previousRowInfo.insert(pair<int, HighsBasisStatus>(previousRowColor, HighsBasisStatus::ZERO));
				previousRowValue.insert(pair<int, double>(previousRowColor, row_value[i]));
			}
		}
	}
}

void HighsAggregate::collectPartsForGS(){
	activeColorHistory_.assign(numRow, false);
	colorsForGS.assign(numRow, false);
	vector<bool> colorAdded(numRow, false);
	int i, j, rep, rhs, rCol, prevColor;
	if (masterIter < 2){
		for (i = 0; i < realNumRow; ++i){
			if (row_status[i] != HighsBasisStatus::BASIC){
				if (prevC[i].size() > 1 && (prevC[i].size() != C[i].size())){
					partsForGS.push_back(i + numCol);
					colorsForGS[i] = true;
				}
				activeColorHistory_[i] = true;
				activeConstraints_[i] = true;
			}
		}
		for (i = 0; i < realNumRow_; ++i){
			rep = C[i + numCol].front();
			prevColor = previousRowColoring[rep - numCol];
			if (activeColorHistory_[prevColor - numCol]){
				activeColorHistory_[i] = true;
				activeConstraints_[i] = true;
			}
		}
	}
	else{
		for (i = 0; i < realNumRow; ++i){
			rCol = rowColor[i];
			cout << "rCol: " << rCol << endl;
			if (row_status[i] != HighsBasisStatus::BASIC){
				if (prevC[rCol].size() > 1 && (prevC[rCol].size() != C[rCol].size())){
					partsForGS.push_back(rCol);
					colorAdded[rCol - numCol] = true;
					colorsForGS[rCol - numCol] = true;
				}
				activeColorHistory_[rCol - numCol] = true;
				activeConstraints_[rCol - numCol] = true;
			} 
		}
		for (i = 0; i < activeColorHistory.size(); ++i){
			if (!activeColorHistory[i])
				continue;
			else{
				if (!colorAdded[i]){
					if (prevC[i + numCol].size() > 1 && (prevC[i + numCol].size() != C[i + numCol].size())){
						partsForGS.push_back(i + numCol);
						colorsForGS[i] = true;
					}
					activeColorHistory_[i] = true;
					activeConstraints_[i] = true;
				}
			}
		}
		for (i = 0; i < realNumRow_; ++i){
			rep = C[i + numCol].front();
			//cout << "rep: " << rep << endl;
			prevColor = previousRowColoring[rep - numCol];
			//cout << "prevColor: " << prevColor << endl;
			if (activeColorHistory_[prevColor - numCol]){
				//cout << "i: " << i << endl;
				activeColorHistory_[i] = true;
				activeConstraints_[i] = true;
			}
		}
	}
}

void HighsAggregate::findPreviousBasisForColumns(){
	int i, rep;
	double lb, ub;
	for (i = 0; i < realNumCol; ++i){
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

vector<double> HighsAggregate::aggregateImpliedRow(int impliedRow){
	int i, j, weight = 0;
	vector<double> coeff(realNumCol_, 0);
	for (i = impliedARstart[impliedRow]; i < impliedARstart[impliedRow + 1]; ++i){
		for (j = 0; j < realNumCol_; ++j){
			if (color[impliedARindex[i]] == j){
				coeff[j] += impliedARvalue[i];
			}
		}
	}
	return coeff;
}

void HighsAggregate::getAggImpliedRows(){
	int i;
	for (i = 0; i < impliedNumRow; ++i){
		aggImpliedRows.push_back(aggregateImpliedRow(i));
	}
}

