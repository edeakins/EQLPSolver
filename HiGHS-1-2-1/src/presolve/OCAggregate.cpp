#include "OCAggregate.h"

void HighsOCAggregate::allocate(HighsLp* lp, OCPartition* partition){
    olp = lp;
    ep = partition;
    epMinusOne = (OCPartition*)calloc(1, sizeof(OCPartition));
    numCol = lp->num_col_;
    numRow = lp->num_row_;
    numTot = numCol + numRow;
    nnz = lp->a_matrix_.value_.size();
    numTotResiduals = numCol - partition->ncsplits;
    // Allocate new lp container, new lp solution, new lp basis
    agglp = (HighsLp*)calloc(1, sizeof(HighsLp));
    agglp->col_cost_.resize(numCol);
    agglp->col_upper_.resize(numCol);
    agglp->col_lower_.resize(numCol);
    agglp->row_upper_.resize(numRow);
    agglp->row_lower_.resize(numRow);
    agglp->a_matrix_.start_.resize(numCol + 1);
    agglp->a_matrix_.index_.resize(nnz);
    agglp->a_matrix_.value_.resize(nnz);
    elp = (HighsLp*)calloc(1, sizeof(HighsLp));
    elpBasis = (HighsBasis*)calloc(1, sizeof(HighsBasis));
    elp->col_cost_.resize(numCol + numTotResiduals);
    elp->col_upper_.resize(numCol + numTotResiduals);
    elp->col_lower_.resize(numCol + numTotResiduals);
    elp->row_upper_.resize(numRow + numTotResiduals);
    elp->row_lower_.resize(numRow + numTotResiduals);
    elp->a_matrix_.start_.resize(numCol + numTotResiduals + 1);
    elp->a_matrix_.index_.resize(nnz + numTotResiduals * 3);
    elp->a_matrix_.value_.resize(nnz + numTotResiduals * 3);
    elpBasis->col_status.resize(numCol + numTotResiduals, HighsBasisStatus::kNonbasic);
    elpBasis->row_status.resize(numRow + numTotResiduals, HighsBasisStatus::kBasic);
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
    copyPartition();
}

void HighsOCAggregate::buildLp(OCPartition* partition, HighsBasis* b,
                               HighsSolution* s, bool finish, bool extended){
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
    copyPartition();
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
        for (mIdx = olp->a_matrix_.start_[crep]; mIdx < olp->a_matrix_.start_[crep + 1]; ++mIdx){
            ro = olp->a_matrix_.index_[mIdx];
            rf = ep->front[ro + numCol];
            r = frontRow[rf];
            rlen = ep->len[rf] + 1;
            if (!columnF[r]++){ 
                columnI[nnz++] = r;
                // nnzStan++;
            }
            columnX[r] += olp->a_matrix_.value_[mIdx] * clen;
        }
        for (mIdx = start; mIdx < nnz; ++mIdx){
            agglp->a_matrix_.index_[mIdx] = columnI[mIdx];
            // AdegenRIndex[j] = columnI[j];
            agglp->a_matrix_.value_[mIdx] = columnX[columnI[mIdx]];
            // AdegenRValue[j] = columnX[columnI[j]];
            columnF[columnI[mIdx]] = 0;
            columnX[columnI[mIdx]] = 0;
            columnI[mIdx] = 0;
        }
        agglp->a_matrix_.start_[c + 1] = nnz;
        // AdegenRStart[xi + 1] = nnzStan;
        start = nnz;
    }
    agglp->num_col_ = colCnt;
    agglp->num_row_ = rowCnt;
    agglp->num_residual_cols_ = 0;
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
        for (mIdx = olp->a_matrix_.start_[crep]; mIdx < olp->a_matrix_.start_[crep + 1]; ++mIdx){
            ro = olp->a_matrix_.index_[mIdx];
            rf = ep->front[ro + numCol];
            r = frontRow[rf];
            rlen = ep->len[rf] + 1;
            if (!columnF[r]++){ 
                columnI[nnz++] = r;
                // nnzStan++;
            }
            columnX[r] += olp->a_matrix_.value_[mIdx] * clen;
        }
        for (mIdx = start; mIdx < nnz; ++mIdx){
            elp->a_matrix_.index_[mIdx] = columnI[mIdx];
            // AdegenRIndex[j] = columnI[j];
            elp->a_matrix_.value_[mIdx] = columnX[columnI[mIdx]];
            // AdegenRValue[j] = columnX[columnI[j]];
            columnF[columnI[mIdx]] = 0;
            columnX[columnI[mIdx]] = 0;
            columnI[mIdx] = 0;
        }
        // Add linking row entry for parent columns
        if (isParent[c]){
            for (mIdx = parentStart[c]; mIdx < parentStart[c + 1]; ++mIdx){
                elp->a_matrix_.index_[nnz] = mIdx;
                // AdegenRIndex[nnzStan] = j;
                elp->a_matrix_.value_[nnz++] = 1;
                // AdegenRValue[nnzStan++] = 1;
            }
        } 
        // Add linking row entry for child columns
        if (isChild[c]){
            elp->a_matrix_.index_[nnz] = childRow[c];
            // AdegenRIndex[nnzStan] = childRow[xi];
            elp->a_matrix_.value_[nnz++] = -1;
            // AdegenRValue[nnzStan++] = -1;
        }
        elp->a_matrix_.start_[c + 1] = nnz;
        // AdegenRStart[xi + 1] = nnzStan;
        start = nnz;
    }
    // Add linking row entry for residual columns
    for (iCol = 0; iCol < numResiduals; ++iCol){
        elp->a_matrix_.index_[nnz] = residualRow[iCol];
        // AdegenRIndex[nnzStan] = residualRow[i];
        elp->a_matrix_.value_[nnz++] = -1;
        // AdegenRValue[nnzStan++] = -1;
        elp->a_matrix_.start_[elp->num_col_ + iCol + 1] = nnz;
        // AdegenRStart[elp->num_col_ + i + 1] = nnzStan;
    }
    elp->num_col_ += numResiduals;
    elp->num_row_ += numResiduals;
    elp->num_residual_cols_ = numResiduals;
    elp->residual_cols_ = residualCol;
    elpBasis->num_col_ = elp->num_col_;
    elpBasis->num_row_ = elp->num_row_;
}

void HighsOCAggregate::checkAMatrix(){
    int iCol, crep, iRow, rrep, aIndex;
    std::vector<int> compareElp;
    std::vector<int> compareOlp(numRow);
    for (iCol = 0; iCol < elp->num_aggregate_cols_; ++iCol){
        int numColEntries = 0;
        compareElp.clear();
        crep = colrep[iCol];
        for (aIndex = elp->a_matrix_.start_[iCol]; aIndex < elp->a_matrix_.start_[iCol + 1]; ++aIndex){
            iRow = elp->a_matrix_.index_[aIndex];
            if (iRow >= elp->num_aggregate_rows_) continue;
            rrep = rowrep[iRow] - numCol;
            compareElp.push_back(rrep);
            numColEntries++;
        }
        for (aIndex = olp->a_matrix_.start_[crep]; aIndex < olp->a_matrix_.start_[crep + 1]; ++aIndex){
            iRow = olp->a_matrix_.index_[aIndex];
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
    for (int iCol = 0; iCol < elp->num_col_ + elp->num_row_; ++iCol){
        std::vector<double> column(elp->num_row_, 0);
        if (iCol < elp->num_col_){
            for (int idx = elp->a_matrix_.start_[iCol]; idx < elp->a_matrix_.start_[iCol + 1]; ++idx){
                int iRow = elp->a_matrix_.index_[idx];
                int iRowValue = elp->a_matrix_.value_[idx];
                column[iRow] = iRowValue;
            }
        }
        else{
            int iRow = iCol - elp->num_col_;
            column[iRow] = -1;
        }
        printMat.push_back(column);
    }
    std::cout << "Standard Form A Matrix" << std::endl;
    std::cout << "[ ";
    for (int iCol = 0; iCol < elp->num_col_ + elp->num_row_; ++iCol){
        for (int iRow = 0; iRow < elp->num_row_; ++iRow){
            std::cout << printMat[iCol][iRow] << " ";
        }
        std::cout << ";" << std::endl;
    }
    std::cout << "]" << std::endl;
    // b vector
    std::cout << "Rhs" << std::endl;
    std::cout << "[ "; 
    for (int iRow = 0; iRow < elp->num_row_; ++iRow)
        std::cout << std::min(std::fabs(elp->row_lower_[iRow]), std::fabs(elp->row_upper_[iRow])) << " ";
    std::cout << "]" << std::endl;
    // basis vectors
    std::cout << "Basic Indices" << std::endl;
    HighsBasisStatus basic = HighsBasisStatus::kBasic;
    std::vector<int> baseIndex;
    for (int iCol = 0; iCol < elp->num_col_; ++iCol)
        if (elpBasis->col_status[iCol] == basic) baseIndex.push_back(iCol);
    for (int iRow = 0; iRow < elp->num_row_; ++iRow)
        if (elpBasis->row_status[iRow] == basic) baseIndex.push_back(iRow + elp->num_col_);
    std::cout << "[ ";
    for (int i = 0; i < baseIndex.size(); ++i)
        std::cout << baseIndex[i] + 1 << " ";
    std::cout << "]" << std::endl;
    std::cout << "Objective Coefficients" << std::endl;
    std::cout << "[ ";
    for (int iCol = 0; iCol < elp->num_col_; ++iCol)
        std::cout << elp->col_cost_[iCol] << " ";
    std::cout << "]" << std::endl;
}

void HighsOCAggregate::buildObj(){
    int iCol, rep, cf;
    double c, clen;
    agglp->col_cost_.resize(colCnt);
    for (iCol = 0; iCol < colCnt; ++iCol){
        rep = colrep[iCol];
        cf = colFront[iCol];
        clen = ep->len[cf] + 1;
        c = olp->col_cost_[rep];
        agglp->col_cost_[iCol] = clen * c;
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
        c = olp->col_cost_[rep];
        elp->col_cost_[iCol] = clen * c;
    }
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

void HighsOCAggregate::buildRhsFromScratch(){
    int iRow, rf, rrep;
    double rlen;
    agglp->row_upper_.resize(rowCnt);
    agglp->row_lower_.resize(rowCnt);
    for (iRow = 0; iRow < rowCnt; ++iRow){
        rrep = rowrep[iRow];
        rf = rowFront[iRow];
        rlen = ep->len[rf] + 1;
        agglp->row_lower_[iRow] = 
            olp->row_lower_[rrep - numCol] * rlen;
        agglp->row_upper_[iRow] = 
            olp->row_upper_[rrep - numCol] * rlen;
    }
}

void HighsOCAggregate::buildRhsFromScratchExtended(){
    int iRow, rf, rrep;
    double rlen;
    for (iRow = 0; iRow < rowCnt; ++iRow){
        rrep = rowrep[iRow];
        rf = rowFront[iRow];
        rlen = ep->len[rf] + 1;
        elp->row_lower_[iRow] = 
            olp->row_lower_[rrep - numCol] * rlen;
        elp->row_upper_[iRow] = 
            olp->row_upper_[rrep - numCol] * rlen;
    }
}

void HighsOCAggregate::buildRhsFromSolution(){
    int iRow, pr, rf, rpf, rlen, prlen, rep;
    double pv, lb, ub, lbDiff, ubDiff, tol = 1e-6;
    agglp->row_upper_.resize(rowCnt);
    agglp->row_lower_.resize(rowCnt);
    for (iRow = 0; iRow < rowCnt; ++iRow){
        rep = rowrep[iRow];
        rf = ep->front[rep];
        rpf = epMinusOne->front[rep];
        rlen = ep->len[rf] + 1;
        pr = pFrontRow[rpf];
        prlen = epMinusOne->len[rpf] + 1;
        pv = (double)solution->row_value[pr]/prlen;
        lb = olp->row_lower_[rep - numCol];
        ub = olp->row_upper_[rep - numCol];
        lbDiff = std::fabs(pv - lb);
        ubDiff = std::fabs(pv - ub);
        if (ubDiff < tol){
            agglp->row_lower_[iRow] = ub * rlen;
            agglp->row_upper_[iRow] = ub * rlen;
        }
        else if (lbDiff < tol){
            agglp->row_lower_[iRow] = lb * rlen;
            agglp->row_upper_[iRow] = lb * rlen;
        }
        else{
            agglp->row_lower_[iRow] = lb * rlen;
            agglp->row_upper_[iRow] = ub * rlen;
        }
    }
}

void HighsOCAggregate::buildRhsFromSolutionExtended(){
    int iRow, pr, rf, rpf, rlen, prlen, rep;
    double pv, lb, ub, lbDiff, ubDiff, tol = 1e-6;
    for (iRow = 0; iRow < rowCnt; ++iRow){
        rep = rowrep[iRow];
        rf = ep->front[rep];
        rpf = epMinusOne->front[rep];
        rlen = ep->len[rf] + 1;
        pr = pFrontRow[rpf];
        prlen = epMinusOne->len[rpf] + 1;
        pv = (double)solution->row_value[pr]/prlen;
        lb = olp->row_lower_[rep - numCol];
        ub = olp->row_upper_[rep - numCol];
        lbDiff = std::fabs(pv - lb);
        ubDiff = std::fabs(pv - ub);
        if (ubDiff < tol){
            elp->row_lower_[iRow] = ub * rlen;
            elp->row_upper_[iRow] = ub * rlen;
        }
        else if (lbDiff < tol){
            elp->row_lower_[iRow] = lb * rlen;
            elp->row_upper_[iRow] = lb * rlen;
        }
        else{
            elp->row_lower_[iRow] = lb * rlen;
            elp->row_upper_[iRow] = ub * rlen;
        }
    }
    for (iRow = rowCnt; iRow < elp->num_row_; ++iRow){
        elp->row_lower_[iRow] = 0;
        elp->row_upper_[iRow] = 0;
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

void HighsOCAggregate::buildBndsFromScratch(){
    int iCol, crep;
    agglp->col_lower_.resize(colCnt);
    agglp->col_upper_.resize(colCnt);
    for (iCol = 0; iCol < colCnt; ++iCol){
        crep = colrep[iCol];
        agglp->col_lower_[iCol] = 
            olp->col_lower_[crep];
        agglp->col_upper_[iCol] = 
            olp->col_upper_[crep];
    }
}

void HighsOCAggregate::buildBndsFromScratchExtended(){
    int iCol, crep;
    for (iCol = 0; iCol < colCnt; ++iCol){
        crep = colrep[iCol];
        elp->col_lower_[iCol] = 
            olp->col_lower_[crep];
        elp->col_upper_[iCol] = 
            olp->col_upper_[crep];
    }
}

void HighsOCAggregate::buildBndsFromSolution(){
    int iCol, pc, cf, pcf, rep;
    double pv, ubDiff, lbDiff, lb, ub, tol = 1e-6;
    agglp->col_lower_.resize(colCnt);
    agglp->col_upper_.resize(colCnt);
    for (iCol = 0; iCol < colCnt; ++iCol){
        rep = colrep[iCol];
        cf = ep->front[rep];
        pcf = epMinusOne->front[rep];
        pc = pFrontCol[pcf];
        pv = solution->col_value[pc];
        ub = olp->col_upper_[rep];
        lb = olp->col_lower_[rep];
        ubDiff = std::fabs(pv - ub);
        lbDiff = std::fabs(pv - lb);
        if (ubDiff < tol){
            agglp->col_lower_[iCol] = ub;
            agglp->col_upper_[iCol] = ub;
        }
        else if (lbDiff < tol){
            agglp->col_lower_[iCol] = lb;
            agglp->col_upper_[iCol] = lb;
        }
        else{
            agglp->col_lower_[iCol] = lb;
            agglp->col_upper_[iCol] = ub;
        }
    }
}

void HighsOCAggregate::buildBndsFromSolutionExtended(){
    int iCol, pc, cf, pcf, rep;
    double pv, ubDiff, lbDiff, lb, ub, tol = 1e-6;
    for (iCol = 0; iCol < colCnt; ++iCol){
        rep = colrep[iCol];
        cf = ep->front[rep];
        pcf = epMinusOne->front[rep];
        pc = pFrontCol[pcf];
        pv = solution->col_value[pc];
        ub = olp->col_upper_[rep];
        lb = olp->col_lower_[rep];
        ubDiff = std::fabs(pv - ub);
        lbDiff = std::fabs(pv - lb);
        if (ubDiff < tol){
            elp->col_lower_[iCol] = ub;
            elp->col_upper_[iCol] = ub;
        }
        else if (lbDiff < tol){
            elp->col_lower_[iCol] = lb;
            elp->col_upper_[iCol] = lb;
        }
        else{
            elp->col_lower_[iCol] = lb;
            elp->col_upper_[iCol] = ub;
        }
    }
    for (iCol = colCnt; iCol < elp->num_col_; ++iCol){
        elp->col_lower_[iCol] = 0;
        elp->col_upper_[iCol] = 0;
    }
}   

void HighsOCAggregate::buildResiduals(){
    if (!ep->level) return;
    buildResidualLinks();
    buildResidualCols();
    buildResidualRows();
}

void HighsOCAggregate::buildResidualLinks(){
    int i, x1, x2, x0, nf, of, newNumRow = elp->num_row_;
    int rep, xOld, xNew, pcf, minParent = numCol, oldRep, newRep;
    HighsBasisStatus basic = HighsBasisStatus::kBasic;
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
    int elpNumRow = elp->num_row_;
    int elpNumCol = elp->num_col_;
    for (int iCol = pcolCnt; iCol < colCnt; ++iCol){
        xNew = iCol;
        rep = colrep[iCol];
        newRep = colrep[iCol];
        pcf = epMinusOne->front[rep];
        xOld = pFrontCol[pcf];
        oldRep = colrep[xOld];
        if (basis->col_status[xOld] != basic) continue;
        if (xOld < minParent) minParent = xOld;
        isParent[xOld] = 1;
        isChild[xNew] = 1;
        splitCells[xOld].push_back(xNew);
    }
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

void HighsOCAggregate::buildResidualCols(){
    int i, p, c, idx = elp->num_col_;
    for (i = 0; i < numResiduals; ++i){
        p = parent[i];
        c = child[i];
        elp->col_lower_[idx] = 0;
        elp->col_upper_[idx++] = 0;
    }
}

void HighsOCAggregate::buildResidualRows(){
    int i, idx = elp->num_row_;
    for (i = 0; i < numResiduals; ++i){
        elp->row_lower_[idx] = elp->row_upper_[idx] = 0;
        idx++;
    }
}

void HighsOCAggregate::buildBasis(bool finish, bool extended){
    buildColBasis();
    buildRowBasis();
}

void HighsOCAggregate::buildColBasis(){
    int iCol, pCol, pf, pc, crep, pcrep;
    HighsBasisStatus basic = HighsBasisStatus::kBasic, status;
    std::set<int> cols;
    int numOldBasic = 0;
    int numOldBasicToSplit = 0;
    int numNewBasic = 0;
    int numNonSingle = 0;
    std::fill(elpBasis->col_status.begin(), elpBasis->col_status.end(), basic);
    for (iCol = 0; iCol < colCnt; ++iCol){
        crep = colrep[iCol];
        pf = epMinusOne->front[crep];
        pCol = pFrontCol[pf];
        status = basis->col_status[pCol];
        elpBasis->col_status[iCol] = status;
    }
    for (iCol = colCnt; iCol < colCnt + numResiduals; ++iCol)
        elpBasis->col_status[iCol] = HighsBasisStatus::kLower;
    for (iCol = pcolCnt; iCol < colCnt; ++iCol){
        if (elpBasis->col_status[iCol] == HighsBasisStatus::kBasic)
            numNewBasic++;
    }
}

void HighsOCAggregate::buildRowBasis(){
    int iRow, r, pr, pf, rrep, rlen;
    int of, nf;
    HighsBasisStatus basic = HighsBasisStatus::kBasic, status;
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
        elpBasis->row_status[iRow] = HighsBasisStatus::kLower;
    }
}

void HighsOCAggregate::buildResidualColBasis(){
    int i;
    for (i = elp->num_col_ - elp->num_residual_cols_; i < elp->num_col_; ++i)
        elpBasis->col_status[i] = HighsBasisStatus::kLower;
    elpBasis->num_col_ = elp->num_col_;
}

void HighsOCAggregate::buildResidualRowBasis(){
    int i;
    for (i = elp->num_row_ - elp->num_residual_cols_; i < elp->num_row_; ++i)
        elpBasis->row_status[i] = HighsBasisStatus::kLower;
    elpBasis->num_row_ = elp->num_row_;
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
            fronts.insert(ep->front[i]);
            frontCol[ep->front[i]] = col[i];
            colFront[col[i]] = ep->front[i];
            continue;
        }
        newFront = fronts.insert(ep->front[i]).second;
        if (newFront){
            frontCol[ep->front[i]] = colCnt;
            colFront[colCnt] = ep->front[i];
            colrep[colCnt] = i;
            col[i] = colCnt++;
        }
    }
    elp->num_col_ = ep->ncsplits;
    elp->num_aggregate_cols_ = ep->ncsplits;
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
            fronts.insert(ep->front[i]);
            frontRow[ep->front[i]] = row[i - numCol];
            rowFront[row[i - numCol]] = ep->front[i];
            continue;
        }
        newFront = fronts.insert(ep->front[i]).second;
        if (newFront){
            frontRow[ep->front[i]] = rowCnt;
            rowFront[rowCnt] = ep->front[i];
            rowrep[rowCnt] = i;
            row[i - numCol] = rowCnt++;
        }
    }
    elp->num_row_ = ep->nrsplits;
    elp->num_aggregate_rows_ = ep->nrsplits;
}

void HighsOCAggregate::copyPartition(){
    epMinusOne->target = ep->target;
    epMinusOne->level = ep->level;
    epMinusOne->nsplits = ep->nsplits;
    epMinusOne->ncsplits = ep->ncsplits;
    epMinusOne->nrsplits = ep->nrsplits;
    // std::copy(ep->front.begin(), ep->front.end(), epMinusOne->front.begin());
    // std::copy(ep->label.begin(), ep->label.end(), epMinusOne->label.begin());
    // std::copy(ep->unlabel.begin(), ep->unlabel.end(), epMinusOne->unlabel.begin());
    // std::copy(ep->parent.begin(), ep->parent.end(), epMinusOne->parent.begin());
    // std::copy(ep->len.begin(), ep->len.end(), epMinusOne->len.begin());
    epMinusOne->front = ep->front;
    epMinusOne->label = ep->label;
    epMinusOne->unlabel = ep->unlabel;
    epMinusOne->parent = ep->parent;
    epMinusOne->len = ep->len;
}

HighsLp* HighsOCAggregate::getLp(){
    return elp;
}

HighsLp* HighsOCAggregate::getAggLp(){
    return agglp;
}

HighsBasis* HighsOCAggregate::getBasis(){
    return elpBasis;
}
