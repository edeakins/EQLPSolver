#include "OCAggregate.h"

void HighsOCAggregate::allocate(HighsLp* lp, OCpartition* partition){
    olp = lp;
    ep = partition;
    numCol = lp->numCol_;
    numRow = lp->numRow_;
    numTot = numCol + numRow;
    nnz = lp->Avalue_.size();
    numTotResiduals = numCol - partition->ncsplits;
    // Allocate new lp container, new lp solution, new lp basis
    elp = (HighsLp*)calloc(1, sizeof(HighsLp));
    elpSolution = (HighsSolution*)calloc(1, sizeof(HighsSolution));
    elpBasis = (HighsBasis*)calloc(1, sizeof(HighsBasis));
    elp->colCost_.resize(numCol + numTotResiduals);
    elp->colUpper_.resize(numCol + numTotResiduals);
    elp->colLower_.resize(numCol + numTotResiduals);
    elp->rowUpper_.resize(numRow + numTotResiduals);
    elp->rowLower_.resize(numRow + numTotResiduals);
    elp->Astart_.resize(numCol + numTotResiduals + 1);
    elp->Aindex_.resize(nnz + numTotResiduals * 3);
    elp->Avalue_.resize(nnz + numTotResiduals * 3);
    elpSolution->col_value.resize(numCol + numTotResiduals);
    elpSolution->col_dual.resize(numCol + numTotResiduals);
    elpSolution->row_value.resize(numRow + numTotResiduals);
    elpSolution->row_dual.resize(numRow + numTotResiduals);
    elpBasis->col_status.resize(numCol + numTotResiduals, HighsBasisStatus::NONBASIC);
    elpBasis->row_status.resize(numRow + numTotResiduals, HighsBasisStatus::BASIC);
    // Allocate col and row pointers
    col.resize(numCol);
    colrep.resize(numCol);
    row.resize(numRow);
    rowrep.resize(numRow);
    // Allocate temp column storage
    columnI.resize(nnz);
    columnX.resize(numRow);
    columnF.resize(numRow);
    // Allocate temp row storage
    rowF.resize(numCol + numTotResiduals);
    rowNonbasic.resize(numRow);
    // Allocate link storage
    parent.resize(numCol);
    child.resize(numCol);
    buildLp();
}

void HighsOCAggregate::buildLp(){
    buildColPointers();
    buildRowPointers();
    buildObj();
    buildAmatrix();
    buildRhs();
    buildBnds();
    buildResiduals();
}

void HighsOCAggregate::buildLp(OCpartition* partition, HighsBasis* b,
                               HighsSolution* s, bool finish, bool extended){
    buildFinalLp = finish;
    ep = partition;
    basis = b;
    solution = s;
    pcol = col;
    pcolrep = colrep;
    prow = row;
    prowrep = rowrep;
    if (buildFinalLp && extended){
        buildFinalPointers();
        buildFinalObj();
        buildFinalAmatrix();
        buildFinalRhs();
        buildFinalBnds();
        buildFinalResiduals();
    }
    else if (buildFinalLp){
        buildFinalPointers();
        buildFinalObj();
        buildFinalAmatrix();
        buildFinalRhs();
        buildFinalBnds();
    }
    else{
        buildColPointers();
        buildRowPointers();
        buildObj();
        buildAmatrix();
        buildRhs();
        buildBnds();
        buildResiduals();
    }
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
    for (i = elp->numCol_; i < olp->numCol_ + numTotResiduals; ++i)
        elp->Astart_[i + 1] = elp->Astart_[i];
}

void HighsOCAggregate::buildFinalAmatrix(){
    int i, j;
    for (i = 0; i < numCol; ++i){
        for (j = olp->Astart_[i]; j < olp->Astart_[i + 1]; ++j){
            elp->Avalue_[j] = olp->Avalue_[j];
            elp->Aindex_[j] = olp->Aindex_[j];
        }
        elp->Astart_[i + 1] = olp->Astart_[i + 1];
    }
    elp->nnz_ = elp->Astart_[numCol];
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

void HighsOCAggregate::buildFinalObj(){
    int i;
    for (i = 0; i < numCol; ++i)
        elp->colCost_[i] = olp->colCost_[i];
    elp->sense_ = olp->sense_;
}

void HighsOCAggregate::buildRhs(){
    if (ep->level)
        buildRhsFromScratch(); //buildRhsFromSolution();
    else
        buildRhsFromScratch();
}

void HighsOCAggregate::buildFinalRhs(){
    // buildFinalRhsFromScratch(); 
    buildFinalRhsFromSolution();
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

void HighsOCAggregate::buildFinalRhsFromScratch(){
    int i;
    for (i = 0; i < numRow; ++i){
        elp->rowLower_[i] = olp->rowLower_[i];
        elp->rowUpper_[i] = olp->rowUpper_[i];
    }
}

void HighsOCAggregate::buildFinalRhsFromSolution(){
    int i, r, pr, pclen;
    double pv, lb, ub, lbDiff, ubDiff, tol = 1e-6;
    for (i = numCol; i < numTot; ++i){
        r = i - numCol;
        pr = prow[i - numCol];
        pclen = ep->len[ep->front[i]] + 1;
        pv = (double)solution->row_value[pr]/pclen;
        lb = olp->rowLower_[r];
        ub = olp->rowUpper_[r];
        lbDiff = std::fabs(pv - lb);
        ubDiff = std::fabs(pv - ub);
        if (ubDiff < tol){
            elp->rowLower_[r] = ub;
            elp->rowUpper_[r] = ub;
        }
        else if (lbDiff < tol){
            elp->rowLower_[r] = lb;
            elp->rowUpper_[r] = lb;
        }
        else{
            elp->rowLower_[r] = lb;
            elp->rowUpper_[r] = ub;
        }
    }
}

void HighsOCAggregate::buildBnds(){
    if (ep->level)
        buildBndsFromScratch(); //buildBndsFromSolution();
    else
        buildBndsFromScratch();
}

void HighsOCAggregate::buildFinalBnds(){
    // buildFinalBndsFromScratch(); 
    buildFinalBndsFromSolution();
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

void HighsOCAggregate::buildFinalBndsFromScratch(){
    int i;
    for (i = 0; i < numCol; ++i){
        elp->colLower_[i] = olp->colLower_[i];
        elp->colUpper_[i] = olp->colUpper_[i];
    }
}

void HighsOCAggregate::buildFinalBndsFromSolution(){
    int i, pc, c;
    double pv, ubDiff, lbDiff, lb, ub, tol = 1e-6;
    for (i = 0; i < numCol; ++i){
        c = i;
        pc = pcol[i];
        pv = solution->col_value[pc];
        ub = olp->colUpper_[c];
        lb = olp->colLower_[c];
        ubDiff = std::fabs(pv - ub);
        lbDiff = std::fabs(pv - lb);
        if (ubDiff < tol){
            elp->colLower_[c] = ub;
            elp->colUpper_[c] = ub;
        }
        else if (lbDiff < tol){
            elp->colLower_[c] = lb;
            elp->colUpper_[c] = lb;
        }
        else{
            elp->colLower_[c] = lb;
            elp->colUpper_[c] = ub;
        }
    }
}   

void HighsOCAggregate::buildResiduals(){
    if (!ep->level) return;
    buildResidualLinks();
    buildResidualCols();
    buildResidualRows();
    buildResidualSubMatrix();
}

void HighsOCAggregate::buildFinalResiduals(){
    buildFinalResidualLinks();
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

void HighsOCAggregate::buildFinalResidualLinks(){
    int i, x1, x2;
    numResiduals = 0;
    for (i = 0; i < numCol; ++i){
        x1 = pcol[i];
        x2 = i;
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
    for (i = 0; i < numResiduals; ++i){
        elp->rowLower_[idx] = elp->rowUpper_[idx] = 0;
        idx++;
    }
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
    int reverse = numResiduals - 1;
    for (i = numResiduals - 1; i >= 0; --i){
        elp->Avalue_[elp->Astart_[parent[reverse - i] + 1] - rowF[parent[reverse - i]]] = 1;
        elp->Avalue_[elp->Astart_[child[reverse - i] + 1] - rowF[child[reverse - i]]] = -1;
        elp->Avalue_[elp->Astart_[elp->numCol_ + 1] - rowF[elp->numCol_]] = -1;
        elp->Aindex_[elp->Astart_[parent[reverse - i] + 1] - rowF[parent[reverse - i]]--] = elp->numRow_;
        elp->Aindex_[elp->Astart_[child[reverse - i] + 1] - rowF[child[reverse - i]]--] = elp->numRow_;
        elp->Aindex_[elp->Astart_[elp->numCol_ + 1] - rowF[elp->numCol_++]--] = elp->numRow_++;
    }
    for (i = elp->numCol_; i < olp->numCol_ + numTotResiduals; ++i)
        elp->Astart_[i + 1] = elp->Astart_[i];
    elp->numResiduals_ = numResiduals;
}

void HighsOCAggregate::buildSolution(bool finish, bool extended){
    if (!finish && extended){
        buildColSolution();
        buildRowSolution();
        buildResidualColSolution();
        buildResidualRowSolution();
        elpSolution->numCol = elp->numCol_;
        elpSolution->numRow = elp->numRow_;
    }
    else if (finish && extended){
        buildFinalColSolution();
        buildFinalRowSolution();
        buildResidualColSolution();
        buildResidualRowSolution();
        elpSolution->numCol = elp->numCol_;
        elpSolution->numRow = elp->numRow_;
    }
    else{
        buildFinalColSolution();
        buildFinalRowSolution();
        elpSolution->numCol = elp->numCol_;
        elpSolution->numRow = elp->numRow_;
    }
}

void HighsOCAggregate::buildColSolution(){
    int i, c, pc;
    double ppv, pdv;
    for (i = 0; i < numCol; i += ep->len[i] + 1){
        c = col[ep->label[i]];
        pc = pcol[ep->label[i]];
        ppv = solution->col_value[pc];
        pdv = solution->col_dual[pc];
        elpSolution->col_value[c] = ppv;
        elpSolution->col_dual[c] = pdv;
    }
}

void HighsOCAggregate::buildResidualColSolution(){
    int i;
    for (i = elp->numCol_ - numResiduals; i < elp->numCol_; ++i){
        elpSolution->col_value[i] = 0;
        elpSolution->col_dual[i] = 0;
    }
}

void HighsOCAggregate::buildRowSolution(){
    int i, r, pr;
    double ppv, pdv;
    for (i = numCol; i < numTot; i += ep->len[i] + 1){
        r = row[ep->label[i] - numCol];
        pr = prow[ep->label[i] - numCol];
        ppv = solution->row_value[pr];
        pdv = solution->row_dual[pr];
        elpSolution->row_value[r] = ppv;
        elpSolution->row_dual[r] = pdv;
    }
}

void HighsOCAggregate::buildResidualRowSolution(){
    int i;
    for (i = elp->numRow_ - numResiduals; i < elp->numRow_; ++i){
        elpSolution->row_value[i] = 0;
        elpSolution->row_dual[i] = 0;
    }
}

void HighsOCAggregate::buildFinalColSolution(){
    int i, c, pc;
    double ppv, pdv;
    for (i = 0; i < numCol; ++i){
        c = i;
        pc = pcol[i];
        ppv = solution->col_value[pc];
        pdv = solution->col_dual[pc];
        elpSolution->col_value[c] = ppv;
        elpSolution->col_dual[c] = pdv;
    }
}

void HighsOCAggregate::buildFinalRowSolution(){
    int i, r, pr;
    double ppv, pdv, pclen;
    for (i = numCol; i < numTot; ++i){
        r = i - numCol;
        pr = prow[i - numCol];
        pclen = ep->len[ep->front[i]] + 1;
        ppv = solution->row_value[pr];
        pdv = solution->row_dual[pr];
        elpSolution->row_value[r] = (double)ppv/pclen;
        elpSolution->row_dual[r] = pdv;// (double)pdv/pclen;
    }
}

void HighsOCAggregate::buildBasis(bool finish){
    if (finish){
        buildFinalColBasis();
        buildFinalRowBasis();
        buildResidualColBasis();
        buildResidualRowBasis();
    }
    else{
        buildColBasis();
        buildRowBasis();
        buildResidualColBasis();
        buildResidualRowBasis();
    }
}

void HighsOCAggregate::buildColBasis(){
    
}

void HighsOCAggregate::buildRowBasis(){

}

void HighsOCAggregate::buildResidualColBasis(){
    int i;
    for (i = elp->numCol_ - elp->numResiduals_; i < elp->numCol_; ++i)
        elpBasis->col_status[i] = HighsBasisStatus::LOWER;
    elpBasis->numCol_ = elp->numCol_;
}

void HighsOCAggregate::buildResidualRowBasis(){
    int i;
    for (i = elp->numRow_ - elp->numResiduals_; i < elp->numRow_; ++i)
        elpBasis->row_status[i] = HighsBasisStatus::LOWER;
    elpBasis->numRow_ = elp->numRow_;
}

void HighsOCAggregate::buildFinalColBasis(){
    int i, c, pc;
    HighsBasisStatus status;;
    for (i = 0; i < numCol; ++i){
        c = i;
        pc = pcol[i];
        status = basis->col_status[pc];
        elpBasis->col_status[c] = status;
    }
}

void HighsOCAggregate::buildFinalRowBasis(){
    int i, r, pr, rep;
    HighsBasisStatus basic = HighsBasisStatus::BASIC;
    for (i = 0; i < basis->row_status.size(); ++i){
        if (basis->row_status[i] != basic){
            rep = prowrep[i];
            elpBasis->row_status[rep] = basis->row_status[i];
        }
    }
}

void HighsOCAggregate::buildColPointers(){
    int i, j, cnt = 0;
    for (i = 0; i < numCol; i += ep->len[i] + 1){
        colrep[cnt] = ep->label[i];
        for (j = i; j <= i + ep->len[i]; ++j){
            col[ep->label[j]] = cnt;
        }
        ++cnt;
    }
    elp->numCol_ = ep->ncsplits;
}

void HighsOCAggregate::buildRowPointers(){
    int i, j, cnt = 0;
    for (i = numCol; i < numTot; i += ep->len[i] + 1){
        rowrep[cnt] = ep->label[i] - numCol;
        for (j = i; j <= i + ep->len[i]; ++j){
            row[ep->label[j] - numCol] = cnt;
        }
        ++cnt;
    }
    elp->numRow_ = ep->nrsplits;
}

void HighsOCAggregate::buildFinalPointers(){
    int i;
    for (i = 0; i < numCol; ++i)
        col[i] = colrep[i] = i;
    for (i = 0; i < numRow; ++i)
        row[i] = rowrep[i] = i;
    elp->numCol_ = numCol;
    elp->numRow_ = numRow;
}

HighsLp* HighsOCAggregate::getLp(){
    return elp;
}

HighsSolution* HighsOCAggregate::getSolution(){
    return elpSolution;
}

HighsBasis* HighsOCAggregate::getBasis(){
    return elpBasis;
}

