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
    numTotResiduals = numCol - partition->ncsplits;
    // Allocate new lp container
    elp = (HighsLp*)calloc(1, sizeof(HighsLp));
    elp->colCost_.resize(numCol + numTotResiduals);
    elp->colUpper_.resize(numCol + numTotResiduals);
    elp->colLower_.resize(numCol + numTotResiduals);
    elp->rowUpper_.resize(numRow + numTotResiduals);
    elp->rowLower_.resize(numRow + numTotResiduals);
    elp->Astart_.resize(numCol + numTotResiduals + 1);
    elp->Aindex_.resize(nnz + numTotResiduals * 3);
    elp->Avalue_.resize(nnz + numTotResiduals * 3);
    // Allocate col and row pointers
    col.resize(numCol);
    row.resize(numRow);
    // Allocate temp column storage
    columnI.resize(nnz);
    columnX.resize(numRow);
    columnF.resize(numRow);
    // Allocate temp row storage
    rowF.resize(numCol + numTotResiduals);
    // Allocate link storage
    parent.resize(numCol);
    child.resize(numCol);
    buildLp();
    return elp;
}

HighsLp* HighsOCAggregate::buildLp(){
    buildColPointers();
    buildRowPointers();
    buildObj();
    buildAmatrix();
    buildRhs();
    buildBnds();
    buildResiduals();
    pcol = col;
    prow = row;
    return elp;
}

HighsLp* HighsOCAggregate::buildLp(OCpartition* partition, HighsBasis* b, HighsSolution* s){
    ep = partition;
    basis = b;
    solution = s;
    buildColPointers();
    buildRowPointers();
    buildObj();
    buildAmatrix();
    buildRhs();
    buildBnds();
    buildResiduals();
    pcol = col;
    prow = row;
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
            // columnX[ci] += flen > 0 ? olp->Avalue_[j] : (olp->Avalue_[j] * xlen);
            columnX[ci] += olp->Avalue_[j] * xlen;
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
    for (i = elp->numCol_; i < olp->numCol_ + numTotResiduals; ++i)
        elp->Astart_[i + 1] = elp->Astart_[i];
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
        buildRhsFromScratch(); //buildRhsFromSolution();
    else
        buildRhsFromScratch();
}

void HighsOCAggregate::buildRhsFromScratch(){
    int i;
    double clen;
    for (i = numCol; i < numTot; i += ep->len[i] + 1){
        clen = ep->len[ep->front[ep->label[i]]] + 1;
        elp->rowLower_[row[ep->label[i] - numCol]] = 
            olp->rowLower_[ep->label[i] - numCol] * clen;
        elp->rowUpper_[row[ep->label[i] - numCol]] = 
            olp->rowUpper_[ep->label[i] - numCol] * clen;
    }
}

void HighsOCAggregate::buildRhsFromSolution(){

}

void HighsOCAggregate::buildBnds(){
    if (solution->col_value.size())
        buildBndsFromScratch(); //buildBndsFromSolution();
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

void HighsOCAggregate::buildResiduals(){
    if (!ep->level) return;
    buildResidualLinks();
    buildResidualCols();
    buildResidualRows();
    buildResidualSubMatrix();
}

void HighsOCAggregate::buildResidualLinks(){
    int i, x1, x2;
    numResiduals = 0;
    for (i = 0; i < numCol; i += ep->len[i] + 1){
        x1 = pcol[ep->label[i]];
        x2 = col[ep->label[i]];
        if (x1 == x2) continue;
        parent[numResiduals] = x1;
        child[numResiduals++] = x2;
    }
}

void HighsOCAggregate::buildResidualCols(){
    int i, p, c, idx = elp->numCol_;
    for (i = 0; i < numResiduals; ++i){
        p = parent[i];
        c = child[i];
        elp->colLower_[idx] = -HIGHS_CONST_INF;
        elp->colUpper_[idx++] = HIGHS_CONST_INF;
    }
}

void HighsOCAggregate::buildResidualRows(){
    int i, idx = elp->numRow_;
    for (i = 0; i < numResiduals; ++i)
        elp->rowLower_[idx] = elp->rowUpper_[idx++] = 0;
}

void HighsOCAggregate::buildResidualSubMatrix(){
    int i, j, start;
    int numNewRow = numResiduals;
    int numNewNz = numResiduals * 3;
    int newNumNz = elp->nnz_ + numNewNz;
    for (i = 0; i < numResiduals; ++i){
        rowF[parent[i]]++;
        rowF[child[i]]++;
        rowF[i + elp->numCol_]++;
    }
    int newEl = newNumNz;
    for (i = elp->numCol_ + numResiduals - 1; i >= 0; --i){
        start = newEl;
        newEl -= rowF[i];
        for (j = elp->Astart_[i + 1] - 1; j >= elp->Astart_[i]; --j){
            --newEl;
            elp->Aindex_[newEl] = elp->Aindex_[j];
            elp->Avalue_[newEl] = elp->Avalue_[j];
        }
        elp->Astart_[i + 1] = start;
    }
    for (i = numResiduals - 1; i >= 0; --i){
        elp->Avalue_[elp->Astart_[parent[i] + 1] - rowF[parent[i]]] = 1;
        elp->Avalue_[elp->Astart_[child[i] + 1] - rowF[child[i]]] = -1;
        elp->Avalue_[elp->Astart_[elp->numCol_ + 1] - rowF[elp->numCol_]] = -1;
        elp->Aindex_[elp->Astart_[parent[i] + 1] - rowF[parent[i]]--] = elp->numRow_;
        elp->Aindex_[elp->Astart_[child[i] + 1] - rowF[child[i]]--] = elp->numRow_;
        elp->Aindex_[elp->Astart_[elp->numCol_ + 1] - rowF[elp->numCol_++]--] = elp->numRow_++;
    }
    for (i = elp->numCol_; i < olp->numCol_ + numTotResiduals; ++i)
        elp->Astart_[i + 1] = elp->Astart_[i];
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

