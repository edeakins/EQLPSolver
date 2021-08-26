#include "OCAggregate.h"

HighsLp* HighsOCAggregate::allocate(HighsLp* lp, OCpartition* partition, HighsBasis* b, HighsSolution* s){
    olp = lp;
    ep = partition;
    basis = b;
    solution = s;
    numCol = lp->numCol_;
    numRow = lp->numRow_;
    numTot = numCol + numRow;
    nnz = lp->Avalue_.size();
    numResiduals = numCol - partition->ncsplits;
    // Allocate new lp container
    elp = (HighsLp*)calloc(1, sizeof(HighsLp));
    elp->colCost_.resize(numCol + numResiduals);
    elp->colUpper_.resize(numCol + numResiduals);
    elp->colLower_.resize(numCol + numResiduals);
    elp->rowUpper_.resize(numRow + numResiduals);
    elp->rowLower_.resize(numRow + numResiduals);
    elp->Astart_.resize(numCol + numResiduals);
    elp->Aindex_.resize(nnz + numResiduals * 3);
    elp->Avalue_.resize(nnz + numResiduals * 3);
    // Allocate col and row pointers
    col.resize(numCol);
    pcol.resize(numCol);
    row.resize(numRow);
    prow.resize(numRow);
    // Allocate temp column storage
    columnI.resize(nnz);
    columnX.resize(numRow, 0.0);
    columnF.resize(numRow);
    buildColPointers();
    buildRowPointers();
    buildLp();
}

HighsLp* HighsOCAggregate::buildLp(){
    buildObj();
    buildAmatrix();
    buildRhs();
    buildBnds();
    return elp;
}

void HighsOCAggregate::buildAmatrix(){
    int i, j, ci, f, flen, xf, xi,
    nnz = 0, start = 0;
    double xlen, xv;
    for (i = 0; i < numCol; i += ep->len[i] + 1){
        xi = col[ep->label[i]];
        xf = ep->front[ep->label[i]];
        xlen = ep->len[xf] + 1;
        for (j = olp->Astart_[ep->label[i]]; j < olp->Astart_[ep->label[i] + 1]; ++j){
            ci = row[olp->Aindex_[j]];
            f = ep->front[olp->Aindex_[j] + numCol];
            flen = ep->len[f];
            if (!columnF[ci]++){ 
                columnI[nnz++] = ci;
            }
            columnX[ci] += flen > 0 ? olp->Avalue_[j] : (olp->Avalue_[j] * xlen);
        }
        for (j = start; j < nnz; ++j){
            elp->Aindex_[j] = columnI[j];
            elp->Avalue_[j] = columnX[columnI[j]];
            columnF[columnI[j]] = 0;
            columnX[columnI[j]] = 0;
            columnI[j] = 0;
        }
        elp->Astart_[xi + 1] = nnz;
        start = nnz;
    }
    elp->nnz_ = nnz;
    elp->numCol_ = ep->ncsplits;
    elp->numRow_ = ep->nrsplits;
}

void HighsOCAggregate::buildObj(){
    int i, x;
    double c, xlen;
    for (i = 0; i < numCol; i += ep->len[i] + 1){
        xlen = ep->len[i] + 1;
        c = olp->colCost_[ep->label[i]];
        x = col[ep->label[i]];
        elp->colCost_[x] = c * xlen;
    }
    elp->sense_ = olp->sense_;
}

void HighsOCAggregate::buildRhs(){
    if (solution->row_value.size())
        buildRhsFromSolution();
    else
        buildRhsFromScratch();
}

void HighsOCAggregate::buildRhsFromScratch(){
    int i;
    for (i = numCol; i < numTot; i += ep->len[i] + 1){
        elp->rowLower_[row[ep->label[i] - numCol]] = 
            olp->rowLower_[ep->label[i] - numCol];
        elp->rowUpper_[row[ep->label[i] - numCol]] = 
            olp->rowUpper_[ep->label[i] - numCol];
    }
}

void HighsOCAggregate::buildRhsFromSolution(){

}

void HighsOCAggregate::buildBnds(){
    if (solution->col_value.size())
        buildBndsFromSolution();
    else
        buildBndsFromScratch();
}

void HighsOCAggregate::buildBndsFromScratch(){
    int i;
    for (i = 0; i < numCol; i += ep->len[i] + 1){
        elp->colLower_[col[ep->label[i]]] = 
            olp->colLower_[ep->label[i]];
        elp->colUpper_[col[ep->label[i]]] = 
            olp->colUpper_[ep->label[i]];
    }
}

void HighsOCAggregate::buildBndsFromSolution(){

}

void HighsOCAggregate::buildColPointers(){
    int i, j, cnt = 0;
    for (i = 0; i < numCol; i += ep->len[i] + 1){
        for (j = i; j <= i + ep->len[i]; ++j){
            col[ep->label[j]] = cnt;
        }
        ++cnt;
    }
}

void HighsOCAggregate::buildRowPointers(){
    int i, j, cnt = 0;
    for (i = numCol; i < numTot; i += ep->len[i] + 1){
        for (j = i; j <= i + ep->len[i]; ++j){
            row[ep->label[j] - numCol] = cnt;
        }
        ++cnt;
    }
}

