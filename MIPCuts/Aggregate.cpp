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
    _rowLower.assign(rL.begin(), rL.end());
    _rowUpper.assign(rU.begin(), rU.end());
    _Avalue.assign(Av.begin(), Av.end());
    _Aindex.assign(Ai.begin(), Ai.end());
    _Astart.assign(As.begin(), As.end());

    // update EP info
    C.assign(ep.C.begin(), ep.C.end());
    Csize.assign(ep.Csize.begin(), ep.Csize.end());
    color.assign(ep.color.begin(), ep.color.end());
    AindexP.assign(ep.AindexP.begin(), ep.AindexP.end());
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
}

void AggregateLp::scanForCuts(){
    cut.clear();
    cutIdx.clear();
    nCuts = nRows_;
    for (int i = 0; i < _Aindex.size(); ++i){
        if (_Aindex[i] >= nRows){
            pair<set<int>::iterator, bool> ret = cut.insert(_Aindex[i]);
            if (ret.second){
                cutIdx[_Aindex[i]] = nCuts;
                nCuts++;
            }
        }
    }
    rowLower_.resize(nCuts);
    rowUpper_.resize(nCuts);
}

void AggregateLp::aggregateColBounds(){
    for (int i = tempNCols_; i < nCols_; ++i){
        int rep = C[i].front();
        colLower_[i] = colLower[rep];
        colUpper_[i] = colUpper[rep];
    }
}

void AggregateLp::aggregateRowBounds(){
    for (int i = 0; i < nRows_; ++i){
        int rep = C[i + nCols].front() - nCols;
        int orbitSize = Csize[color[rep + nCols]];
        rowLower_[i] = rowLower[rep] * orbitSize;
        rowUpper_[i] = rowUpper[rep] * orbitSize;
    }
    for (int i = nRows_; i < nCuts; ++i){
        rowLower_[i] = _rowLower[i - nRows_ + nRows];
        rowUpper_[i] = _rowUpper[i - nRows_ + nRows];
    }
}

void AggregateLp::aggregateAMatrix(){
	Astart_.push_back(0);
	for (int i = 0; i < nCols_; ++i){
        vector<double> coeff(nRows_);
		int rep = C[i].front();
		for (int j = Astart[rep]; j < Astart[rep + 1]; ++j){
            int rowIdx = AindexP[j] - nCols;
            coeff[rowIdx] += Avalue[j];
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

void AggregateLp::addCutsToAggregate(){
    if (!_Astart.size())
        return;
    for (int i = 0; i < nCols_; ++i){
        int rep = C[i].front();
        // vector<double> coeff(nRows_ + _nRows - nRows);
        for (int j = _Astart[rep]; j < _Astart[rep + 1]; ++j){
            if (_Aindex[j] >= nRows){
                int aggColEnd = Astart_[i + 1];
                int cIdx = cutIdx[_Aindex[j]];
                double value = _Avalue[j];
                Aindex_.insert(Aindex_.begin() + aggColEnd, cIdx);
                Avalue_.insert(Avalue_.begin() + aggColEnd, value);
                for (int k = i + 1; k <= nCols_; ++k){
                    Astart_[k]++;
                }
            }
        }
    }
}

void AggregateLp::aggregate(){
    clear();
    findDimensions();
    scanForCuts();
    aggregateColBounds();
    aggregateRowBounds();
    aggregateAMatrix();
    addCutsToAggregate();
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

int AggregateLp::getNumCol(){
    return nCols_;
}

int AggregateLp::getNumRow(){
    if (!_Astart.size())
        return nRows_;
    return nRows_ + _nRows - nRows;
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
