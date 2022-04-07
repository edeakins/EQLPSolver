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
    agglp = (HighsLp*)calloc(1, sizeof(HighsLp));
    agglp->colCost_.resize(numCol);
    agglp->colUpper_.resize(numCol);
    agglp->colLower_.resize(numCol);
    agglp->rowUpper_.resize(numRow);
    agglp->rowLower_.resize(numRow);
    agglp->Astart_.resize(numCol + 1);
    agglp->Aindex_.resize(nnz);
    agglp->Avalue_.resize(nnz);
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
    col.assign(numCol, -1);
    colrep.assign(numCol, -1);
    row.assign(numRow, -1);
    rowrep.assign(numRow, -1);
    finalRowRep.resize(numRow, false);
    degenSlack.resize(numRow + numTotResiduals, 0);
    degenRow.resize(numRow, true);
    // Allocate temp column storage
    columnI.resize(nnz + numTotResiduals * 3);
    columnX.resize(numRow);
    columnF.resize(numRow);
    // Allocate temp row storage
    rowF.resize(numCol + numTotResiduals);
    rowNonbasic.resize(numRow);
    // Allocate link storage
    parent.resize(numCol);
    isParent.resize(numCol);
    parentStart.resize(numCol + 1);
    parentFreq.resize(numCol);
    parentRow.resize(numCol);
    child.resize(numCol);
    isChild.resize(numCol);
    childRow.resize(numCol);
    residualCol.resize(numCol);
    residualRow.resize(numCol);
    // Reserve mem for basicIndex and nonbasicFlag
    basicIndex.reserve(numRow + numTotResiduals);
    nonbasicFlag.reserve(numCol + numTotResiduals + numRow);
    rowPerm.assign(numRow, -1);
    rowUnperm.assign(numRow, -1);
    colPerm.assign(numCol + numTotResiduals, -1);
    colUnperm.assign(numCol + numTotResiduals, -1);
    frontCol.assign(numCol, -1);
    colFront.assign(numCol, -1);
    frontRow.assign(numRow + numCol, -1);
    rowFront.assign(numRow + numCol, -1);
    buildLp();
}

void HighsOCAggregate::buildLp(){
    buildColPointers();
    buildRowPointers();
    buildObj();
    buildAmatrix();
    buildRhs();
    buildBnds();
    // buildResiduals();
    epMinusOne = *ep;
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
    pFrontCol = frontCol;
    pcolCnt = colCnt;
    pFrontRow = frontRow;
    prowCnt = rowCnt;
    if (buildFinalLp && extended){
        buildFinalPointers();
        buildFinalObj();
        buildFinalAmatrix();
        buildFinalRhs();
        buildFinalBnds();
        buildFinalResiduals();
        buildBasis(true, true);
    }
    else if (buildFinalLp){
        buildFinalPointers();
        buildFinalObj();
        buildFinalAmatrix();
        buildFinalRhs();
        buildFinalBnds();
        buildBasis(true, false);
    }
    else{
        buildColPointers();
        buildRowPointers();
        buildResidualLinks();
        buildBasis(false, false);
        buildObj();
        buildObjExtended();
        buildAmatrix();
        buildAmatrixExtended();
        // checkAMatrix();
        // buildStandardMatrix();
        buildRhs();
        buildRhsExtended();
        buildBnds();
        buildBndsExtended();
        // printAMatrixToMatlabFormat();
        // setUpPreBasicIndex();
        // setUpPreNonbasicFlag();
        // setUpColAq();
        // setUpPreLU();
        // setUpPreMatrix();
        // factor();
        // buildSubLp();
        // buildSubLpBasis();
        // buildResiduals();
        // buildBasis(false, false);
        // buildSparseMat();
        // buildLUFactor();
    }
    epMinusOne = *ep;
}

// Build A matrix without linkers if they exist
void HighsOCAggregate::buildAmatrix(){
    int iCol, cf, crep, clen, c;
    int mIdx, ro, rf, r, rlen;
    int nnzStan = 0, start = 0, startStan = 0;
    nnz = 0;
    double xlen, xv;
    std::set<int>::iterator it;
    for (iCol = 0; iCol < colCnt; ++iCol){
        crep = colrep[iCol];
        cf = colFront[iCol];
        clen = ep->len[cf] + 1;
        c = iCol;
        for (mIdx = olp->Astart_[crep]; mIdx < olp->Astart_[crep + 1]; ++mIdx){
            ro = olp->Aindex_[mIdx];
            rf = ep->front[ro + numCol];
            r = frontRow[rf];
            rlen = ep->len[rf] + 1;
            if (!columnF[r]++){ 
                columnI[nnz++] = r;
                // nnzStan++;
            }
            columnX[r] += olp->Avalue_[mIdx] * clen;
        }
        for (mIdx = start; mIdx < nnz; ++mIdx){
            agglp->Aindex_[mIdx] = columnI[mIdx];
            // AdegenRIndex[j] = columnI[j];
            agglp->Avalue_[mIdx] = columnX[columnI[mIdx]];
            // AdegenRValue[j] = columnX[columnI[j]];
            columnF[columnI[mIdx]] = 0;
            columnX[columnI[mIdx]] = 0;
            columnI[mIdx] = 0;
        }
        agglp->Astart_[c + 1] = nnz;
        // AdegenRStart[xi + 1] = nnzStan;
        start = nnz;
    }
    agglp->numCol_ = colCnt;
    agglp->numRow_ = rowCnt;
    agglp->numResiduals_ = 0;
    agglp->nnz_ = nnz;
}

// Build A matrix with linkers if they exist
void HighsOCAggregate::buildAmatrixExtended(){
    int iCol, cf, crep, clen, c;
    int mIdx, ro, rf, r, rlen;
    int nnzStan = 0, start = 0, startStan = 0;
    nnz = 0;
    double xlen, xv;
    std::set<int>::iterator it;
    for (iCol = 0; iCol < colCnt; ++iCol){
        crep = colrep[iCol];
        cf = colFront[iCol];
        clen = ep->len[cf] + 1;
        c = iCol;
        for (mIdx = olp->Astart_[crep]; mIdx < olp->Astart_[crep + 1]; ++mIdx){
            ro = olp->Aindex_[mIdx];
            rf = ep->front[ro + numCol];
            r = frontRow[rf];
            rlen = ep->len[rf] + 1;
            if (!columnF[r]++){ 
                columnI[nnz++] = r;
                // nnzStan++;
            }
            columnX[r] += olp->Avalue_[mIdx] * clen;
        }
        for (mIdx = start; mIdx < nnz; ++mIdx){
            elp->Aindex_[mIdx] = columnI[mIdx];
            // AdegenRIndex[j] = columnI[j];
            elp->Avalue_[mIdx] = columnX[columnI[mIdx]];
            // AdegenRValue[j] = columnX[columnI[j]];
            columnF[columnI[mIdx]] = 0;
            columnX[columnI[mIdx]] = 0;
            columnI[mIdx] = 0;
        }
        // Add linking row entry for parent columns
        if (isParent[c]){
            for (mIdx = parentStart[c]; mIdx < parentStart[c + 1]; ++mIdx){
                elp->Aindex_[nnz] = mIdx;
                // AdegenRIndex[nnzStan] = j;
                elp->Avalue_[nnz++] = 1;
                // AdegenRValue[nnzStan++] = 1;
            }
        } 
        // Add linking row entry for child columns
        if (isChild[c]){
            elp->Aindex_[nnz] = childRow[c];
            // AdegenRIndex[nnzStan] = childRow[xi];
            elp->Avalue_[nnz++] = -1;
            // AdegenRValue[nnzStan++] = -1;
        }
        elp->Astart_[c + 1] = nnz;
        // AdegenRStart[xi + 1] = nnzStan;
        start = nnz;
    }
    // Add linking row entry for residual columns
    for (iCol = 0; iCol < numResiduals; ++iCol){
        elp->Aindex_[nnz] = residualRow[iCol];
        // AdegenRIndex[nnzStan] = residualRow[i];
        elp->Avalue_[nnz++] = -1;
        // AdegenRValue[nnzStan++] = -1;
        elp->Astart_[elp->numCol_ + iCol + 1] = nnz;
        // AdegenRStart[elp->numCol_ + i + 1] = nnzStan;
    }
    elp->numCol_ += numResiduals;
    elp->numRow_ += numResiduals;
    elp->numResiduals_ = numResiduals;
    elp->residuals_ = residualCol;
    elp->nnz_ = nnz;
    elpBasis->numCol_ = elp->numCol_;
    elpBasis->numRow_ = elp->numRow_;
}

void HighsOCAggregate::checkAMatrix(){
    int iCol, crep, iRow, rrep, aIndex;
    std::vector<int> compareElp;
    std::vector<int> compareOlp(numRow);
    for (iCol = 0; iCol < elp->numX_; ++iCol){
        int numColEntries = 0;
        compareElp.clear();
        crep = colrep[iCol];
        for (aIndex = elp->Astart_[iCol]; aIndex < elp->Astart_[iCol + 1]; ++aIndex){
            iRow = elp->Aindex_[aIndex];
            if (iRow >= elp->numS_) continue;
            rrep = rowrep[iRow] - numCol;
            compareElp.push_back(rrep);
            numColEntries++;
        }
        for (aIndex = olp->Astart_[crep]; aIndex < olp->Astart_[crep + 1]; ++aIndex){
            iRow = olp->Aindex_[aIndex];
            compareOlp[iRow] = 1;
        }
        for (int i = 0; i < numColEntries; ++i){
            int elpRep = compareElp[i];
            if (compareOlp[elpRep] != 1){
                std::cout << "bad col: " << iCol << std::endl;
                std::cout << "bad col rep: " << crep << std::endl;
                std::cout << "bad row: " << iRow << std::endl; 
                std::cin.get();
            }
        }
        for (int i = 0; i < numRow; ++i){
            compareOlp[i] = 0;
        }
    }
}

void HighsOCAggregate::printAMatrixToMatlabFormat(){
    // A Matrix
    std::vector<std::vector<double> > printMat;
    for (int iCol = 0; iCol < elp->numCol_ + elp->numRow_; ++iCol){
        std::vector<double> column(elp->numRow_, 0);
        if (iCol < elp->numCol_){
            for (int idx = elp->Astart_[iCol]; idx < elp->Astart_[iCol + 1]; ++idx){
                int iRow = elp->Aindex_[idx];
                int iRowValue = elp->Avalue_[idx];
                column[iRow] = iRowValue;
            }
        }
        else{
            int iRow = iCol - elp->numCol_;
            column[iRow] = -1;
        }
        printMat.push_back(column);
    }
    std::cout << "Standard Form A Matrix" << std::endl;
    std::cout << "[ ";
    for (int iCol = 0; iCol < elp->numCol_ + elp->numRow_; ++iCol){
        for (int iRow = 0; iRow < elp->numRow_; ++iRow){
            std::cout << printMat[iCol][iRow] << " ";
        }
        std::cout << ";" << std::endl;
    }
    std::cout << "]" << std::endl;
    // b vector
    std::cout << "Rhs" << std::endl;
    std::cout << "[ "; 
    for (int iRow = 0; iRow < elp->numRow_; ++iRow)
        std::cout << std::min(std::fabs(elp->rowLower_[iRow]), std::fabs(elp->rowUpper_[iRow])) << " ";
    std::cout << "]" << std::endl;
    // basis vectors
    std::cout << "Basic Indices" << std::endl;
    HighsBasisStatus basic = HighsBasisStatus::BASIC;
    std::vector<int> baseIndex;
    for (int iCol = 0; iCol < elp->numCol_; ++iCol)
        if (elpBasis->col_status[iCol] == basic) baseIndex.push_back(iCol);
    for (int iRow = 0; iRow < elp->numRow_; ++iRow)
        if (elpBasis->row_status[iRow] == basic) baseIndex.push_back(iRow + elp->numCol_);
    std::cout << "[ ";
    for (int i = 0; i < baseIndex.size(); ++i)
        std::cout << baseIndex[i] + 1 << " ";
    std::cout << "]" << std::endl;
    std::cout << "Objective Coefficients" << std::endl;
    std::cout << "[ ";
    for (int iCol = 0; iCol < elp->numCol_; ++iCol)
        std::cout << elp->colCost_[iCol] << " ";
    std::cout << "]" << std::endl;
}

// void HighsOCAggregate::buildStandardMatrix(){
//     int iCol, j, iRow, nnzStan = 0, nRow = 0;
//     double iVal;
//     std::vector<int> rowPerm(elp->numRow_, -1);
//     // Permute degenerate row indices
//     for (iRow = 0; iRow < elp->numRow_; ++iRow)
//         if (degenSlack[iRow]) rowPerm[iRow] = nRow++;
//     // Fill in degenerate and linking row sub matrix
//     for (iCol = 0; iCol < elp->numCol_; ++iCol){
//         for (j = elp->Astart_[iCol]; j < elp->Astart_[iCol + 1]; ++j){
//             iRow = rowPerm[elp->Aindex_[j]];
//             if (iRow < 0) continue;
//             iVal = elp->Avalue_[j];
//             AdegenRIndex[nnzStan] = iRow;
//             AdegenRValue[nnzStan++] = iVal;
//         }
//         AdegenRStart[iCol + 1] = nnzStan;
//     }
//     // // Add slacks to A matrix copy for LU factor to remove some r columns
//     // for (int i = 0; i < elp->numS_; ++i){
//     //     AdegenRIndex[nnzStan] = i;
//     //     AdegenRValue[nnzStan++] = 1;
//     //     AdegenRStart[elp->numCol_ + i + 1] = nnzStan;
//     // }
// }

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
    int iCol, rep, cf;
    double c, clen;
    agglp->colCost_.resize(colCnt);
    for (iCol = 0; iCol < colCnt; ++iCol){
        rep = colrep[iCol];
        cf = colFront[iCol];
        clen = ep->len[cf] + 1;
        c = olp->colCost_[rep];
        agglp->colCost_[iCol] = clen * c;
    }
    agglp->sense_ = olp->sense_;
}

void HighsOCAggregate::buildObjExtended(){
    int iCol, rep, cf;
    double c, clen;
    for (iCol = 0; iCol < colCnt; ++iCol){
        rep = colrep[iCol];
        cf = colFront[iCol];
        clen = ep->len[cf] + 1;
        c = olp->colCost_[rep];
        elp->colCost_[iCol] = clen * c;
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
        buildRhsFromSolution();
    else
        buildRhsFromScratch();
}

void HighsOCAggregate::buildRhsExtended(){
    if (ep->level)
        buildRhsFromSolutionExtended();
    else
        buildRhsFromScratchExtended();
}

void HighsOCAggregate::buildFinalRhs(){
    // buildFinalRhsFromScratch(); 
    buildFinalRhsFromSolution();
}

void HighsOCAggregate::buildRhsFromScratch(){
    int iRow, rf, rrep;
    double rlen;
    agglp->rowUpper_.resize(rowCnt);
    agglp->rowLower_.resize(rowCnt);
    for (iRow = 0; iRow < rowCnt; ++iRow){
        rrep = rowrep[iRow];
        rf = rowFront[iRow];
        rlen = ep->len[rf] + 1;
        agglp->rowLower_[iRow] = 
            olp->rowLower_[rrep - numCol] * rlen;
        agglp->rowUpper_[iRow] = 
            olp->rowUpper_[rrep - numCol] * rlen;
    }
}

void HighsOCAggregate::buildRhsFromScratchExtended(){
    int iRow, rf, rrep;
    double rlen;
    for (iRow = 0; iRow < rowCnt; ++iRow){
        rrep = rowrep[iRow];
        rf = rowFront[iRow];
        rlen = ep->len[rf] + 1;
        elp->rowLower_[iRow] = 
            olp->rowLower_[rrep - numCol] * rlen;
        elp->rowUpper_[iRow] = 
            olp->rowUpper_[rrep - numCol] * rlen;
    }
}

void HighsOCAggregate::buildRhsFromSolution(){
    int iRow, pr, rf, rpf, rlen, prlen, rep;
    double pv, lb, ub, lbDiff, ubDiff, tol = 1e-6;
    agglp->rowUpper_.resize(rowCnt);
    agglp->rowLower_.resize(rowCnt);
    for (iRow = 0; iRow < rowCnt; ++iRow){
        rep = rowrep[iRow];
        rf = ep->front[rep];
        rpf = epMinusOne.front[rep];
        rlen = ep->len[rf] + 1;
        pr = pFrontRow[rpf];
        prlen = epMinusOne.len[rpf] + 1;
        pv = (double)solution->row_value[pr]/prlen;
        lb = olp->rowLower_[rep - numCol];
        ub = olp->rowUpper_[rep - numCol];
        lbDiff = std::fabs(pv - lb);
        ubDiff = std::fabs(pv - ub);
        if (ubDiff < tol){
            agglp->rowLower_[iRow] = ub * rlen;
            agglp->rowUpper_[iRow] = ub * rlen;
        }
        else if (lbDiff < tol){
            agglp->rowLower_[iRow] = lb * rlen;
            agglp->rowUpper_[iRow] = lb * rlen;
        }
        else{
            agglp->rowLower_[iRow] = lb * rlen;
            agglp->rowUpper_[iRow] = ub * rlen;
        }
    }
}

void HighsOCAggregate::buildRhsFromSolutionExtended(){
    int iRow, pr, rf, rpf, rlen, prlen, rep;
    double pv, lb, ub, lbDiff, ubDiff, tol = 1e-6;
    for (iRow = 0; iRow < rowCnt; ++iRow){
        rep = rowrep[iRow];
        rf = ep->front[rep];
        rpf = epMinusOne.front[rep];
        rlen = ep->len[rf] + 1;
        pr = pFrontRow[rpf];
        prlen = epMinusOne.len[rpf] + 1;
        pv = (double)solution->row_value[pr]/prlen;
        lb = olp->rowLower_[rep - numCol];
        ub = olp->rowUpper_[rep - numCol];
        lbDiff = std::fabs(pv - lb);
        ubDiff = std::fabs(pv - ub);
        if (ubDiff < tol){
            elp->rowLower_[iRow] = ub * rlen;
            elp->rowUpper_[iRow] = ub * rlen;
        }
        else if (lbDiff < tol){
            elp->rowLower_[iRow] = lb * rlen;
            elp->rowUpper_[iRow] = lb * rlen;
        }
        else{
            elp->rowLower_[iRow] = lb * rlen;
            elp->rowUpper_[iRow] = ub * rlen;
        }
    }
    for (iRow = rowCnt; iRow < elp->numRow_; ++iRow){
        elp->rowLower_[iRow] = 0;
        elp->rowUpper_[iRow] = 0;
    }
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
        pr = prow[r];
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
        buildBndsFromSolution();
    else
        buildBndsFromScratch();
}

void HighsOCAggregate::buildBndsExtended(){
    if (ep->level)
        buildBndsFromSolutionExtended();
    else
        buildBndsFromScratchExtended();
}

void HighsOCAggregate::buildFinalBnds(){
    // buildFinalBndsFromScratch(); 
    buildFinalBndsFromSolution();
}

void HighsOCAggregate::buildBndsFromScratch(){
    int iCol, crep;
    agglp->colLower_.resize(colCnt);
    agglp->colUpper_.resize(colCnt);
    for (iCol = 0; iCol < colCnt; ++iCol){
        crep = colrep[iCol];
        agglp->colLower_[iCol] = 
            olp->colLower_[crep];
        agglp->colUpper_[iCol] = 
            olp->colUpper_[crep];
    }
}

void HighsOCAggregate::buildBndsFromScratchExtended(){
    int iCol, crep;
    for (iCol = 0; iCol < colCnt; ++iCol){
        crep = colrep[iCol];
        elp->colLower_[iCol] = 
            olp->colLower_[crep];
        elp->colUpper_[iCol] = 
            olp->colUpper_[crep];
    }
}

void HighsOCAggregate::buildBndsFromSolution(){
    int iCol, pc, cf, pcf, rep;
    double pv, ubDiff, lbDiff, lb, ub, tol = 1e-6;
    agglp->colLower_.resize(colCnt);
    agglp->colUpper_.resize(colCnt);
    for (iCol = 0; iCol < colCnt; ++iCol){
        rep = colrep[iCol];
        cf = ep->front[rep];
        pcf = epMinusOne.front[rep];
        pc = pFrontCol[pcf];
        pv = solution->col_value[pc];
        ub = olp->colUpper_[rep];
        lb = olp->colLower_[rep];
        ubDiff = std::fabs(pv - ub);
        lbDiff = std::fabs(pv - lb);
        if (ubDiff < tol){
            agglp->colLower_[iCol] = ub;
            agglp->colUpper_[iCol] = ub;
        }
        else if (lbDiff < tol){
            agglp->colLower_[iCol] = lb;
            agglp->colUpper_[iCol] = lb;
        }
        else{
            agglp->colLower_[iCol] = lb;
            agglp->colUpper_[iCol] = ub;
        }
    }
}

void HighsOCAggregate::buildBndsFromSolutionExtended(){
    int iCol, pc, cf, pcf, rep;
    double pv, ubDiff, lbDiff, lb, ub, tol = 1e-6;
    for (iCol = 0; iCol < colCnt; ++iCol){
        rep = colrep[iCol];
        cf = ep->front[rep];
        pcf = epMinusOne.front[rep];
        pc = pFrontCol[pcf];
        pv = solution->col_value[pc];
        ub = olp->colUpper_[rep];
        lb = olp->colLower_[rep];
        ubDiff = std::fabs(pv - ub);
        lbDiff = std::fabs(pv - lb);
        if (ubDiff < tol){
            elp->colLower_[iCol] = ub;
            elp->colUpper_[iCol] = ub;
        }
        else if (lbDiff < tol){
            elp->colLower_[iCol] = lb;
            elp->colUpper_[iCol] = lb;
        }
        else{
            elp->colLower_[iCol] = lb;
            elp->colUpper_[iCol] = ub;
        }
    }
    for (iCol = colCnt; iCol < elp->numCol_; ++iCol){
        elp->colLower_[iCol] = 0;
        elp->colUpper_[iCol] = 0;
    }
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
    int i, x1, x2, x0, nf, of, newNumRow = elp->numRow_;
    int rep, xOld, xNew, pcf, minParent = numCol, oldRep, newRep;
    HighsBasisStatus basic = HighsBasisStatus::BASIC;
    std::pair<std::set<std::pair<int, int> >::iterator, bool> ret;
    numResiduals = 0;
    std::fill(parent.begin(), parent.end(), 0);
    std::fill(isParent.begin(), isParent.end(), 0);
    std::fill(parentStart.begin(), parentStart.end(), 0);
    std::fill(parentFreq.begin(), parentFreq.end(), 0);
    std::fill(child.begin(), child.end(), 0);
    std::fill(isChild.begin(), isChild.end(), 0);
    std::fill(childRow.begin(), childRow.end(), 0);
    std::fill(residualRow.begin(), residualRow.end(), 0);
    int elpNumRow = elp->numRow_;
    int elpNumCol = elp->numCol_;
    for (int iCol = pcolCnt; iCol < colCnt; ++iCol){
        xNew = iCol;
        rep = colrep[iCol];
        newRep = colrep[iCol];
        pcf = epMinusOne.front[rep];
        xOld = pFrontCol[pcf];
        oldRep = colrep[xOld];
        if (basis->col_status[xOld] != basic) continue;
        if (xOld < minParent) minParent = xOld;
        isParent[xOld] = 1;
        isChild[xNew] = 1;
        // parent[]
        // parentFreq[xOld]++;
        // childRow[xNew] = elpNumRow;
        // residualCol[numResiduals] = elp->numCol_ + numResiduals;
        // residualRow[numResiduals++] = elpNumRow++;
        splitCells[xOld].push_back(xNew);
    }
    // parentStart[minParent] = elp->numRow_;
    // parentStart[minParent + 1] = parentStart[minParent] + parentFreq[minParent];
    // int lastStart = minParent + 1;
    // for (int iCol = minParent + 1; iCol < pcolCnt; ++iCol){
    //     if (!isParent[iCol]) continue;
    //     parentStart[iCol] = parentStart[lastStart];
    //     parentStart[iCol + 1] = parentStart[iCol] + parentFreq[iCol];
    //     lastStart = iCol + 1;
    // }
    for (const auto split : splitCells){
        parentStart[split.first] = elpNumRow;
        parentStart[split.first + 1] = elpNumRow + split.second.size();
        for (int idx = 0; idx < split.second.size(); ++idx){
            int iCol = split.second.at(idx);
            childRow[iCol] = elpNumRow;
            residualCol[numResiduals] = elpNumCol++;
            residualRow[numResiduals++] = elpNumRow++;
        }
    }
    splitCells.clear();
}

void HighsOCAggregate::buildFinalResidualLinks(){
    int i, oldCol, newCol, oldRep, newRep, of, nf;
    numResiduals = 0;
    for (i = 0; i < numCol; ++i){
        nf = i;
        of = ep->front[i];
        oldRep = ep->label[of];
        newRep = i;
        oldCol = pcol[of];
        newCol = col[nf];
        if (oldRep == newRep) continue;
        parent[numResiduals] = oldRep;
        child[numResiduals] = newRep;
        residualCol[numResiduals] = elp->numCol_ + numResiduals;
        numResiduals++;
    }
}

void HighsOCAggregate::buildResidualCols(){
    int i, p, c, idx = elp->numCol_;
    for (i = 0; i < numResiduals; ++i){
        p = parent[i];
        c = child[i];
        elp->colLower_[idx] = 0;
        elp->colUpper_[idx++] = 0;
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
    int i, j, start, nnz;
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
    elp->residuals_ = residualCol;
    elp->nnz_ = newNumNz;
    for (int i = elp->numX_; i < elp->numCol_; ++i){
        for (int j = elp->Astart_[i]; j < elp->Astart_[i + 1]; ++j){
            AdegenRIndex.push_back(elp->Aindex_[j]);
            AdegenRValue.push_back(elp->Avalue_[j]);
        }
        AdegenRStart.push_back(elp->Astart_[i + 1]);
    }
    // Fill in slacks
    nnz = elp->Astart_[elp->numCol_];
    for (int i = 0; i < elp->numS_; ++i){
        AdegenRIndex.push_back(i + elp->numCol_);
        AdegenRValue.push_back(1);
        AdegenRStart.push_back(nnz++);
    }
}

// Set up the basicIndex
void HighsOCAggregate::setUpPreBasicIndex(){
    int iCol, iRow;
    HighsBasisStatus basic = HighsBasisStatus::BASIC;
    for (iCol = 0; iCol < elp->numCol_; ++iCol)
        if (elpBasis->col_status[iCol] == basic)
            basicIndex.push_back(iCol);
    for (iRow = 0; iRow < elp->numRow_; ++iRow)
        if (elpBasis->row_status[iRow] == basic)
            basicIndex.push_back(iRow + elp->numCol_);
}

// Set up nonbasicFlag
void HighsOCAggregate::setUpPreNonbasicFlag(){
    int iCol, iRow;
    HighsBasisStatus basic = HighsBasisStatus::BASIC;
    for (iCol = 0; iCol < elp->numCol_; ++iCol){
        if (elpBasis->col_status[iCol] == basic)
            nonbasicFlag.push_back(0);
        else nonbasicFlag.push_back(1);
    }
    for (iRow = 0; iRow < elp->numRow_; ++iRow){
        if (elpBasis->row_status[iRow] == basic)
            nonbasicFlag.push_back(0);
        else nonbasicFlag.push_back(1);
    }
}

// Set up colAq for reduce r columns computation
void HighsOCAggregate::setUpColAq(){
    colAq.setup(elp->numRow_);
}

// Set up the pre LU factor to build B^-1 * A_r sub matrix for pivoting out basic, degenerate slacks
void HighsOCAggregate::setUpPreLU(){
    LU.setup(elp->numCol_, elp->numRow_, &elp->Astart_[0],
                 &elp->Aindex_[0], &elp->Avalue_[0],
                 &basicIndex[0], NULL);
}

void HighsOCAggregate::setUpPreMatrix(){
    matrix.setup(elp->numCol_, elp->numRow_, &elp->Astart_[0],
                 &elp->Aindex_[0], &elp->Avalue_[0],
                 &nonbasicFlag[0]);
}

void HighsOCAggregate::factor(){
    int rankDeficiency;
    rankDeficiency = LU.build();
}

void HighsOCAggregate::reduceColumn(int iCol){
    colAq.clear();
    colAq.packFlag = true;
    matrix.collect_aj(colAq, iCol, 1);
    LU.ftran(colAq, 0);
}

void HighsOCAggregate::buildSubLp(){
    int i, j, iCol, iRow, rIdx = 0, cIdx = 0;
    // Set numCol, numRow
    sublp->numCol_ = numResiduals;
    sublp->numRow_ = degenCnt;
    // Set up colCost
    sublp->colCost_.resize(numResiduals, 0);
    // Set up colLower, colUpper
    sublp->colLower_.assign(numResiduals, 0);
    sublp->colUpper_.assign(numResiduals, 0);
    // Set up rowLower, rowUpper
    sublp->rowLower_.assign(degenCnt, 0);
    sublp->rowUpper_.assign(degenCnt, 0);
    // Build Astart, Aindex, Avalue
    sublp->Astart_.push_back(0);
    for (i = 0; i < numResiduals; ++i){
        iCol = residualCol[i];
        colPerm[iCol] = cIdx;
        colUnperm[cIdx++] = iCol;
        reduceColumn(iCol);
        for (j = 0; j < colAq.count; ++j){
            iRow = colAq.index[j];
            if (!degenSlack[iRow]) continue;
            if (rowPerm[iRow] > -1){
                sublp->Aindex_.push_back(rowPerm[iRow]);
                sublp->Avalue_.push_back(colAq.array[iRow]);
            }
            else{
                rowPerm[iRow] = rIdx;
                rowUnperm[rIdx++] = iRow;
                sublp->Aindex_.push_back(rowPerm[iRow]);
                sublp->Avalue_.push_back(colAq.array[iRow]);
            }
        }
        sublp->Astart_.push_back(sublp->Aindex_.size());
    }
}

void HighsOCAggregate::buildSubLpBasis(){
    sublpBasis->col_status.assign(sublp->numCol_, HighsBasisStatus::NONBASIC);
    sublpBasis->row_status.assign(sublp->numRow_, HighsBasisStatus::BASIC);
    sublpBasis->numCol_ = sublp->numCol_;
    sublpBasis->numRow_ = sublp->numRow_;
}

void HighsOCAggregate::buildSparseMat(){
    // A.LoadFromArrays(elp->numRow_, elp->numCol_ + elp->numS_, &AdegenRStart[0], 
    //                 &AdegenRStart[1], &AdegenRIndex[0], &AdegenRValue[0]);
    int iRow, nRow = 0, dim;
    for (iRow = 0; iRow < elp->numRow_; ++iRow)
        if (degenSlack[iRow]) nRow++;
    // dim = nRow < elp->numCol_ ? elp->numCol_ : nRow;
    A.LoadFromArrays(nRow, elp->numCol_, &AdegenRStart[0], 
                    &AdegenRStart[1], &AdegenRIndex[0], &AdegenRValue[0]);
    AT = A.Transpose(A);
    // std::vector<int> copyRows;
    // for (int i = 0; i < elp->numRow_; ++i){
    //     if (i >= elp->numRow_ - elp->numResiduals_)
    //         copyRows.push_back(i);
    //     else if (elp->degenSlack_[i])
    //         copyRows.push_back(i);
    // }
    // ATsub = AT.CopyColumns(AT, copyRows);
    // Asub = ATsub.Transpose(ATsub);
}

void HighsOCAggregate::buildLUFactor(){
    int rankDeficiency, rank;
    int iRow, iCol, nRow = 0, dim, dimRow, dimPseudo, cnt = 0;
    basicIndex.resize(0);
    indepRows.resize(0);
    // Get fake dimension of B
    dim = AT.rows();
    // Set initial basic index for A'
    for (iCol = 0; iCol < A.cols(); ++iCol)
        basicIndex.push_back(iCol);
    
    // Set up pseudo LU factor for A'
    LU.setup(A.cols(), A.cols(), &A.colptr()[0],
                 &A.rowidx()[0], &A.values()[0],
                 &basicIndex[0], NULL);
    // Build pseudo LU and grab rank/rank deficiency
    rankDeficiency = LU.build();
    rank = dim - rankDeficiency;
    dimPseudo = dim - degenCnt;
    dimRow = degenCnt;
    // Grab independent columns
    for (iRow = 0; iRow < dim; ++iRow)
        if (basicIndex[iRow] < dimRow) ++cnt;
    std::cout << "Real Rows: " << cnt << std::endl;
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

void HighsOCAggregate::buildBasis(bool finish, bool extended){
    if (finish && extended){
        buildFinalColBasis();
        buildFinalRowBasis();
        buildResidualColBasis();
        buildResidualRowBasis();
    }
    else if(finish){
        buildFinalColBasis();
        buildFinalRowBasis();
    }
    else{
        buildColBasis();
        buildRowBasis();
        // buildResidualColBasis();
        // buildResidualRowBasis();
    }
}

void HighsOCAggregate::buildColBasis(){
    int iCol, pCol, pf, pc, crep, pcrep;
    HighsBasisStatus basic = HighsBasisStatus::BASIC, status;
    std::set<int> cols;
    int numOldBasic = 0;
    int numOldBasicToSplit = 0;
    int numNewBasic = 0;
    int numNonSingle = 0;
    std::fill(elpBasis->col_status.begin(), elpBasis->col_status.end(), basic);
    for (iCol = 0; iCol < colCnt; ++iCol){
        crep = colrep[iCol];
        pf = epMinusOne.front[crep];
        pCol = pFrontCol[pf];
        status = basis->col_status[pCol];
        elpBasis->col_status[iCol] = status;
    }
    for (iCol = colCnt; iCol < colCnt + numResiduals; ++iCol)
        elpBasis->col_status[iCol] = HighsBasisStatus::LOWER;
    // for (iCol = 0; iCol < pcolCnt; ++iCol){
    //     if (basis->col_status[iCol] != HighsBasisStatus::BASIC &&
    //         basis->col_status[iCol] != HighsBasisStatus::LOWER &&
    //         basis->col_status[iCol] != HighsBasisStatus::UPPER)
    //         std::cout << "iCol: " << iCol << std::endl;
    // }
    for (iCol = pcolCnt; iCol < colCnt; ++iCol){
        if (elpBasis->col_status[iCol] == HighsBasisStatus::BASIC)
            numNewBasic++;
    }
}

void HighsOCAggregate::buildRowBasis(){
    int iRow, r, pr, pf, rrep, rlen;
    int of, nf;
    HighsBasisStatus basic = HighsBasisStatus::BASIC, status;
    std::fill(elpBasis->row_status.begin(), elpBasis->row_status.end(), basic);
    int numNewBasic = 0;
    int numBasicToSplit = 0;
    int numBasicSplits = 0;
    int numNonBasicToSplit = 0;
    int numNonBasicSplits = 0;
    int numNonBasic = 0;
    for (iRow = 0; iRow < prowCnt; ++iRow){
        elpBasis->row_status[iRow] = basis->row_status[iRow];
    }
    for (iRow = rowCnt; iRow < rowCnt + numResiduals; ++iRow){
        elpBasis->row_status[iRow] = HighsBasisStatus::LOWER;
    }
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
    int i, c, oldRep, newRep, nf, of, oldCol, newCol;
    HighsBasisStatus status;;
    for (i = 0; i < numCol; ++i){
        nf = i;
        of = ep->front[i];
        oldRep = ep->label[of];
        newRep = i;
        oldCol = pcol[of];
        newCol = i;
        status = basis->col_status[oldCol];
        elpBasis->col_status[newCol] = status;
    }
    elpBasis->numCol_ = numCol;
}

void HighsOCAggregate::buildFinalRowBasis(){
    int i, r, pr, rep, of, nf, oldRep, newRep, oldRow, newRow, cnt = 0;
    HighsBasisStatus basic = HighsBasisStatus::BASIC, status;
    for (i = 0; i < basis->row_status.size(); ++i){
        if (basis->row_status[i] != basic){
            oldRep = prowrep[i];
            elpBasis->row_status[oldRep] = basis->row_status[i];
            finalRowRep[oldRep] = true;
        }
    }
    for (i = 0; i < numRow; ++i){
        if (finalRowRep[i]) continue;
        of = ep->front[i + numCol];
        nf = i + numCol;
        pr = prow[i];
        if (basis->row_status[pr] != basic){
            degenSlack[i] = 1;
            cnt++;
        }
    }
    // for (i = 0; i < numRow; ++i){
    //     r = i;
    //     pr = prow[i];
    //     status = basis->row_status[pr];
    //     elpBasis->row_status[r] = status;
    // }
    elpBasis->numRow_ = numRow;
    elp->degenSlack_ = degenSlack;
    elp->numDegenSlack_ = cnt;
}

void HighsOCAggregate::buildColPointers(){
    int i, j;
    // front.clear();
    std::set<int> fronts;
    std::set<int> testCols;
    std::vector<int> colFreq(numCol, 0);
    bool newFront;
    bool duplicatedRep;
    std::pair<std::set<int>::iterator,bool> newRep;
    newColReps.clear();
    for (i = 0; i < numCol; ++i){
        if (col[i] > -1){
            // colFreq[col[i]] ++;
            fronts.insert(ep->front[i]);
            frontCol[ep->front[i]] = col[i];
            colFront[col[i]] = ep->front[i];
            // colrep[colCnt] = i;
            continue;
        }
        newFront = fronts.insert(ep->front[i]).second;
        if (newFront){
            frontCol[ep->front[i]] = colCnt;
            colFront[colCnt] = ep->front[i];
            colrep[colCnt] = i;
            col[i] = colCnt++;
            // colFreq[col[i]]++;
        }
    }
    elp->numCol_ = ep->ncsplits;
    elp->numX_ = ep->ncsplits;
    // for (i = 0; i < colCnt; ++i){
    //     if (testCols.insert(colrep[i]).second)
    //         colFreq[colrep[i]]++;
    // }
    // for (i = 0; i < numCol; ++i){
    //     if (colFreq[i] > 1)
    //         std::cin.get();
    // }
}

void HighsOCAggregate::buildRowPointers(){
    int i, j;
    std::set<int> fronts;
    std::set<int> testRows;
    std::vector<int> rowFreq(numRow, 0);
    bool newFront;
    std::pair<std::set<int>::iterator,bool> newRep;
    newRowReps.clear();
    for (i = numCol; i < numTot; ++i){
        if (row[i - numCol] > -1){
            // rowFreq[row[i - numCol]]++;
            fronts.insert(ep->front[i]);
            frontRow[ep->front[i]] = row[i - numCol];
            rowFront[row[i - numCol]] = ep->front[i];
            // rowrep[rowCnt] = i;
            continue;
        }
        newFront = fronts.insert(ep->front[i]).second;
        if (newFront){
            frontRow[ep->front[i]] = rowCnt;
            rowFront[rowCnt] = ep->front[i];
            rowrep[rowCnt] = i;
            row[i - numCol] = rowCnt++;
            // rowFreq[row[i - numCol]]++;
        }
    }
    elp->numRow_ = ep->nrsplits;
    elp->numS_ = ep->nrsplits;
    // for (i = 0; i < rowCnt; ++i){
    //     if (testRows.insert(rowrep[i]).second)
    //         rowFreq[colrep[i]]++;
    // }
    // for (i = 0; i < numRow; ++i){
    //     if (rowFreq[i] > 1)
    //         std::cin.get();
    // }
}

void HighsOCAggregate::buildFinalPointers(){
    int i;
    for (i = 0; i < numCol; ++i)
        col[i] = colrep[i] = i;
    for (i = 0; i < numRow; ++i)
        row[i] = rowrep[i] = i;
    elp->numCol_ = numCol;
    elp->numX_ = numCol;
    elp->numRow_ = numRow;
    elp->numS_ = numRow;
}

HighsLp* HighsOCAggregate::getLp(){
    return elp;
}

HighsLp* HighsOCAggregate::getAggLp(){
    return agglp;
}

HighsSolution* HighsOCAggregate::getSolution(){
    return elpSolution;
}

HighsBasis* HighsOCAggregate::getBasis(){
    return elpBasis;
}

HighsLp* HighsOCAggregate::getSubLp(){
    return sublp;
}

HighsBasis* HighsOCAggregate::getSubBasis(){
    return sublpBasis;
}

std::vector<int>& HighsOCAggregate::getUnPerm(){
    return rowUnperm;
}
