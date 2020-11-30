#include "Aggregate.hpp"

AggregateLp::AggregateLp(EquitablePartition& ep){
    // Original Lp info
    nRows = ep.nRows;
	nCols = ep.nCols;
	numTot = nRows + nCols;
	colCost.assign(ep.colCost.begin(), ep.colCost.end());
    colLower.assign(ep.colLower.begin(), ep.colLower.end());
    colUpper.assign(ep.colUpper.begin(), ep.colUpper.end());
    rowLower.assign(ep.rowLower.begin(), ep.rowLower.end());
    rowUpper.assign(ep.rowUpper.begin(), ep.rowUpper.end());
    Avalue.assign(ep.Avalue.begin(), ep.Avalue.end());
    Aindex.assign(ep.Aindex.begin(), ep.Aindex.end());
    AindexP.assign(ep.AindexP.begin(), ep.AindexP.end());
    Astart.assign(ep.Astart.begin(), ep.Astart.end());

    // EP info
    C.assign(ep.C.begin(), ep.C.end());
    Csize.assign(ep.Csize.begin(), ep.Csize.end());
    color.assign(ep.color.begin(), ep.color.end());
}

void AggregateLp::updateMasterLpAndEp(EquitablePartition& ep, int _nC, int _nR,
                            int _nnz, vector<int>& As, vector<int>& Ai,
                            vector<double>& Av, vector<double>& rL, vector<double>& rU){
   // update original Lp info
    _nRows = _nR;
    _nCols = _nC;
	_numTot = _nRows + _nCols;
    rowLower.assign(rL.begin(), rL.end());
    rowUpper.assign(rU.begin(), rU.end());
    Avalue.assign(Av.begin(), Av.end());
    Aindex.assign(Ai.begin(), Ai.end());
    Astart.assign(As.begin(), As.end());

    // update EP info
    C.assign(ep.C.begin(), ep.C.end());
    Csize.assign(ep.Csize.begin(), ep.Csize.end());
    color.assign(ep.color.begin(), ep.color.end());
}

void AggregateLp::clear(){
    Avalue_.clear();
    Astart_.clear();
    Aindex_.clear();
}

void AggregateLp::findDimensions(){
    for (int i = tempNCols_; i < nCols; ++i)
        if (Csize[i])
            nCols_++;
    for (int i = tempNRows_; i < nRows; ++i)
        if (Csize[i])
            nRows_++;
    colLower_.resize(nCols_);
    colUpper_.resize(nCols_);
    colCost_.resize(nCols_);
    rowLower_.resize(nRows_);
    rowUpper_.resize(nRows_); 
}

void AggregateLp::scanForCuts(){
    cut.clear();
    cutIdx.clear();
    for (int i = 0; i < AindexP.size(); ++i){
        if (AindexP[i] - nCols < nRows_)
            continue;
        // Working here for cuts need to do set insert and make mapping of cut to its aggregated index
    }
}

void AggregateLp::aggregateColBounds(){
    for (int i = tempNCols_; i < nCols_; ++i){
        int rep = C[i].front();
        colLower_[i] = colLower[rep];
        colUpper_[i] = colUpper[rep];
    }
}

void AggregateLp::aggregateRowBounds(){
     for (int i = tempNRows_; i < nRows_; ++i){
        int rep = C[i + nCols].front() - nCols;
        rowLower_[i] = rowLower[rep];
        rowUpper_[i] = rowUpper[rep];
    }
}

void AggregateLp::aggregateAMatrix(){
	Astart_.push_back(0);
	for (int i = 0; i < nCols_; ++i){
        vector<double> coeff(nRows_ + _nRows - nRows, 0);
		int rep = C[i].front();
		for (int j = Astart[rep]; j < Astart[rep + 1]; ++j){
            int rowIdx = AindexP[j] - nCols;
            if (rowIdx < nRows_)
                coeff[rowIdx] += Avalue[j];
            else{
                coeff[cutIdx[AindexP[j]]] += Avalue[j];
            }
        }
        for (int j = 0; j < coeff.size(); ++j){
            if (coeff[j]){
                Avalue_.push_back(coeff[j] * C[i].size());
                Aindex_.push_back(j);
            }
        }
    Astart_.push_back(Avalue_.size());
	}
}

void AggregateLp::aggregate(){
    clear();
    findDimensions();
    aggregateColBounds();
    aggregateRowBounds();
    aggregateAMatrix();
    aggregateCostVector();
    tempNCols_ = nCols_;
    tempNRows_ = nRows_;
}

void AggregateLp::aggregateCostVector(){
    colCost_.assign(nCols_, 0);
	for (int i = 0; i < nCols_; ++i){
        int rep = C[i].front();
        colCost_[i] = C[i].size() * colCost[rep];
    }
}

vector<double>& AggregateLp::getColUpper(){
    return colUpper_;
}

vector<double>& AggregateLp::getColLower(){
    return colLower_;
}

vector<double>& AggregateLp::getColCost(){
    return colCost_;
}

vector<double>& AggregateLp::getRowUpper(){
    return rowUpper_;
}

vector<double>& AggregateLp::getRowLower(){
    return rowLower_;
}
vector<double>& AggregateLp::getAvalue(){
    return Avalue_;
}

vector<int>& AggregateLp::getAindex(){
    return Aindex_;
}   

vector<int>& AggregateLp::getAstart(){
    return Astart_;
}
