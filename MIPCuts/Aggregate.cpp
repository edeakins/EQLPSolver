#include "Aggregate.hpp"

AggregateLp::AggregateLp(EquitablePartition& ep){
    // Original Lp info
    nRows = ep.nRows;
	nCols = ep.nCols;
	numTot = nRows + nCols;
	rowLower = ep.rowLower;
	rowUpper = ep.rowUpper;
	colUpper = ep.colUpper;
	colLower = ep.colLower;
	colCost = ep.colCost;
	Avalue = ep.Avalue;
	Aindex = ep.Aindex;
    AindexP = ep.AindexP;
	Astart = ep.Astart;

    // EP info
    C = ep.C;
    Csize = ep.Csize;
    color = ep.color;
}

void AggregateLp::updatePacking(vector<int>& packed){
    AindexP = packed;
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

void AggregateLp::aggregateColBounds(){
    for (int i = tempNCols_; i < nCols_; ++i){
        int rep = C[i].front();
        colLower_[i] = colLower[rep];
        colUpper_[i] = colUpper[rep];
    }
}

void AggregateLp::aggregateRowBounds(){
     for (int i = tempNRows_; i < nRows_; ++i){
        int rep = C[i].front() - nCols;
        rowLower_[i] = rowLower[rep];
        rowUpper_[i] = rowUpper[rep];
    }
}

void AggregateLp::aggregateAMatrix(){
	Astart_.push_back(0);
	for (int i = 0; i < nCols_; ++i){
        vector<double> coeff(nRows_, 0);
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
