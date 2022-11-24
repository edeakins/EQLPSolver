#include "OCAggregate.h"
// #include "presolve/ICrashX.h"
#include "ipm/IpxWrapper.h"
#include "ipm/ipx/src/lp_solver.h"

ipx::Control ipx_control_;
ipx::Info ipx_info_;
ipx::Model ipx_model_;
std::unique_ptr<ipx::Basis> ipx_basis_;
std::vector<int> compare_index;

void HighsOCAggregate::passLpAndPartition(HighsLp& lp, OCPartition& partition, HighsOptions& options, HighsInfo& info){
    olp = lp;
    ep = partition;
    test_options = options;
    test_info = info;
    numCol = olp.num_col_;
    numRow = olp.num_row_;
    numTot = numCol + numRow;
    nnz = olp.a_matrix_.value_.size();
    numTotResiduals = numCol - ep.ncsplits;

    // Reserve memory for lp containers
    agglp.col_cost_.reserve(numCol);
    agglp.col_upper_.reserve(numCol);
    agglp.col_lower_.reserve(numCol);
    agglp.row_upper_.reserve(numRow);
    agglp.row_lower_.reserve(numRow);
    agglp.a_matrix_.start_.reserve(numCol);
    agglp.a_matrix_.index_.reserve(nnz);
    agglp.a_matrix_.value_.reserve(nnz);
    elp.col_cost_.reserve(numCol + numTotResiduals);
    elp.col_upper_.reserve(numCol + numTotResiduals);
    elp.col_lower_.reserve(numCol + numTotResiduals);
    elp.row_upper_.reserve(numRow + numTotResiduals);
    elp.row_lower_.reserve(numRow + numTotResiduals);
    elp.a_matrix_.start_.reserve(numCol + numTotResiduals +1);
    elp.a_matrix_.index_.reserve(nnz + 3 * numTotResiduals);
    elp.a_matrix_.value_.reserve(nnz + 3 * numTotResiduals);
    splitFrom.resize(numCol);
    splitFromNonbasicCount.resize(numCol);
    splitSize.resize(numCol);

    // Allocate col and row pointers
    col.assign(numCol, -1);
    colrep.assign(numCol, -1);
    row.assign(numRow, -1);
    rowrep.assign(numRow, -1);
    frontCol.assign(numCol, -1);
    colFront.assign(numCol, -1);
    frontRow.assign(numRow + numCol, -1);
    rowFront.assign(numRow + numCol, -1);
    frontMin.assign(numCol + numRow, -1);
    // Allocate temp column storage
    columnI.resize(nnz + numTotResiduals * 3);
    columnX.resize(numRow);
    columnF.resize(numRow);
    // Allocate link storage
    isParent.resize(numCol);
    parentRow.resize(numCol + 1);
    isChild.resize(numCol);
    parent.resize(numCol);
    childRow.resize(numCol);
    residualCol.resize(numCol);
    residualRow.resize(numCol);
    mark_degenerate.resize(numCol);
    zero_step_pivots.resize(numCol + numTotResiduals);
}

void HighsOCAggregate::resizeAlpContainers(){
    agglp.col_cost_.resize(colCnt);
    agglp.col_upper_.resize(colCnt);
    agglp.col_lower_.resize(colCnt);
    agglp.row_upper_.resize(rowCnt);
    agglp.row_lower_.resize(rowCnt);
    agglp.a_matrix_.start_.resize(0);
    agglp.a_matrix_.index_.resize(0);
    agglp.a_matrix_.value_.resize(0);
}

void HighsOCAggregate::resizeElpContainers(){
    // Resize Lp containers
    elp.col_cost_.resize(colCnt + numResiduals);
    elp.col_upper_.resize(colCnt + numResiduals);
    elp.col_lower_.resize(colCnt + numResiduals);
    elp.row_upper_.resize(rowCnt + numResiduals);
    elp.row_lower_.resize(rowCnt + numResiduals);
    elp.a_matrix_.start_.resize(0);
    elp.a_matrix_.index_.resize(0);
    elp.a_matrix_.value_.resize(0);
}

void HighsOCAggregate::resizeBasisTestContainers(){
    test_basic_start.assign(colCnt + rowCnt + 2 * numResiduals, 0);
    test_basic_finish.assign(colCnt + rowCnt + 2 * numResiduals, 0);
    test_basic_index.resize(0);
    col_needs_res_pivot.resize(colCnt + rowCnt);
    std::fill(col_needs_res_pivot.begin(), col_needs_res_pivot.end(), 0);
}

void HighsOCAggregate::resizeGramSchmidtMatrixContainers(){
    presolvelp.col_cost_.resize(colCnt);
    presolvelp.col_upper_.resize(colCnt);
    presolvelp.col_lower_.resize(colCnt);
    presolvelp.row_upper_.resize(rowCnt + numResiduals);
    presolvelp.row_lower_.resize(rowCnt + numResiduals);
    presolvelp.a_matrix_.start_.resize(0);
    presolvelp.a_matrix_.index_.resize(0);
    presolvelp.a_matrix_.value_.resize(0);
}

void HighsOCAggregate::buildLp(){
    findFrontMins();
    buildColPointers();
    buildRowPointers();
    resizeAlpContainers();
    buildObj();
    buildAmatrix();
    buildRhs();
    buildBnds();
    // buildResiduals();
    copyPartition();
    // pcol = col;
    // pcolrep = colrep;
    // prow = row;
    // prowrep = rowrep;
    // pFrontCol = frontCol;
    // pcolCnt = colCnt;
    // pFrontRow = frontRow;
    // prowCnt = rowCnt;
    agglp.level = level;
    ++level;
}

void HighsOCAggregate::buildLp(OCPartition& partition, HighsBasis& b,
                               HighsSolution& s, std::vector<HighsInt>& basic_index, HEkk& ekk_instance,
                               HighsOptions& options, HighsInfo& info, HighsTimer& timer){
    ep = partition;
    basis = b;
    solution = s;
    test_options = options;
    test_info = info;
    pcol = col;
    pcolrep = colrep;
    prow = row;
    prowrep = rowrep;
    pFrontCol = frontCol;
    pcolCnt = colCnt;
    pFrontRow = frontRow;
    prowCnt = rowCnt;
    findFrontMins();
    buildColPointers();
    buildRowPointers();
    // trackAndCountSplits();
    markDegenerate();
    buildResidualLinks();
    // countResiduals();
    // if (numResiduals < 1000)
    // resizeAlpContainers();
    // resizeBasisTestContainers();
    buildBasis(false, false);
    // buildLpForTestFactor();
    // printAMatrixToMatlabFormat(agglp.a_matrix_);
    // buildSolution();
    // timer.resetHighsTimer();
    // timer.startRunHighsClock();
    // double init_time = timer.readRunHighsClock();
    // // buildTestBasicIndex();
    // // fillTestBasicStart();
    // buildTestBasisLU();
    // double final_time = timer.readRunHighsClock();
    // build_test_factor_time = final_time - init_time;
    // std::cout << "\n factor test time: " << build_test_factor_time << "\n" << std::endl;
    // fillTestBasicFinish();
    // markDependentColumns();
    // buildResidualLinks();
    // buildResidualSubMatrix();
    resizeElpContainers();
    // changeBasis();
    buildObjExtended();
    buildAmatrixExtended();
    buildRhsExtended();
    buildBndsExtended();
    buildRowNames();
    buildColNames();
    resizeBasisTestContainers();
    // buildTestBasicIndex();
    ekk_instance.clear();
    elp.is_moved_ = false;
    elp.is_scaled_ = false;
    HighsLpSolverObject solver_object(elp, elpBasis,
                            lift_solution, info, ekk_instance, options, timer);
    buildTestBasisLU(solver_object);
    // updateTestBasisLU();
    copyPartition();
    agglp.level = level;
    elp.level = level; 
    presolvelp.level = level;
    // elp.pairs = pairs;
    ++level;
}


void HighsOCAggregate::buildLpForIPM(OCPartition& partition, HighsBasis& b,
                               HighsSolution& s, std::vector<HighsInt>& basic_index){
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
    findFrontMins();
    buildColPointers();
    buildRowPointers();
    // trackAndCountSplits();
    // markDegenerate();
    // buildResidualLinks();
    // countResiduals();
    // if (numResiduals < 1000)
    resizeAlpContainers();
    // buildBasis(false, false);
    buildObj();
    buildAmatrix();
    buildRhs();
    buildBnds();
    buildRowNames();
    buildColNames();
    // copyPartition();
    // agglp.level = level;
    // elp.level = level; 
    // presolvelp.level = level;
    // // elp.pairs = pairs;
    // ++level;
}

void HighsOCAggregate::buildLpForTestFactor(){
    buildObj();
    buildAmatrix();
    buildRhs();
    buildBnds();
    buildRowNames();
}

void HighsOCAggregate::buildElp(){
    // ep = partition;
    // basis = b;
    // solution = s;
    // pcol = col;
    // pcolrep = colrep;
    // prow = row;
    // prowrep = rowrep;
    // pFrontCol = frontCol;
    // pcolCnt = colCnt;
    // pFrontRow = frontRow;
    // prowCnt = rowCnt;
    // findFrontMins();
    // buildColPointers();
    // buildRowPointers();
    // trackAndCountSplits();
    // markDegenerate();
    // buildResidualLinks();
    // // countResiduals();
    // // if (numResiduals < 1000)
    resizeElpContainers();
    buildBasis(false, false);
    buildObjExtended();
    buildAmatrixExtended();
    buildRhsExtended();
    buildBndsExtended();
    buildRowNames();
    buildColNames();
    copyPartition();
    pcol = col;
    pcolrep = colrep;
    prow = row;
    prowrep = rowrep;
    pFrontCol = frontCol;
    pcolCnt = colCnt;
    pFrontRow = frontRow;
    prowCnt = rowCnt;
    agglp.level = level;
    elp.level = level; 
    presolvelp.level = level;
    // elp.pairs = pairs;
    ++level;
}

int HighsOCAggregate::countResiduals(OCPartition& partition, HighsBasis& b,
                               HighsSolution& s, std::vector<HighsInt>& basic_index){
    // temp_basis = basis;
    // temp_solution = solution;
    // temp_pcol = pcol;
    // temp_pcolrep = pcolrep;
    // temp_prow = prow;
    // temp_prowrep = prowrep;
    // temp_pFrontCol = pFrontCol;
    // temp_pFrontRow = pFrontRow;
    // temp_pcolCnt = pcolCnt;
    // temp_prowCnt = prowCnt;
    ep = partition;
    basis = b;
    solution = s;
    findFrontMins();
    buildColPointers();
    buildRowPointers();
    // trackAndCountSplits();
    markDegenerate();
    buildResidualLinks();
    return numResiduals;
}

void HighsOCAggregate::resetPreviousContainers(){
    basis = temp_basis;
    solution = temp_solution;
    pcol = temp_pcol;
    pcolrep = temp_pcolrep;
    prow = temp_prow;
    prowrep = temp_prowrep;
    pFrontCol = temp_pFrontCol;
    pFrontRow = temp_pFrontRow;
    pcolCnt = temp_pcolCnt;
    prowCnt = temp_prowCnt;
}

HighsSolution HighsOCAggregate::buildSolution(OCPartition& partition, HighsSolution& s){
    // ep = partition;
    // solution = s;
    // pcol = col;
    // pcolrep = colrep;
    // prow = row;
    // prowrep = rowrep;
    // pFrontCol = frontCol;
    // pcolCnt = colCnt;
    // pFrontRow = frontRow;
    // prowCnt = rowCnt;
    // findFrontMins();
    // buildColPointers();
    // buildRowPointers();
    HighsInt iCol, iRow, rep, cf, pcf, pc,
    rf, prf, pr, prlen; 
    double pv, pd;
    lift_solution.col_value.resize(numCol);
    lift_solution.row_value.resize(numRow);
    lift_solution.col_dual.resize(numCol);
    lift_solution.row_dual.resize(numRow);
    for (iCol = 0; iCol < colCnt; ++iCol){
        rep = colrep[iCol];
        cf = ep.front[rep];
        pcf = epMinusOne.front[rep];
        pc = pFrontCol[pcf];
        pv = solution.col_value[pc];
        pd = solution.col_dual[pc];
        lift_solution.col_value.at(iCol) = 
            pv;
        lift_solution.col_dual.at(iCol) = 
            pd;
    }
    for (iRow = 0; iRow < rowCnt; ++iRow){
        rep = rowrep[iRow];
        rf = ep.front[rep];
        prf = epMinusOne.front[rep];
        prlen = epMinusOne.len[prf] + 1;
        pr = pFrontRow[prf];
        pv = solution.row_value.at(pr);
        pd = solution.row_dual.at(pr);
        lift_solution.row_value.at(iRow) = 
            (double)pv/prlen;
        lift_solution.row_dual.at(iRow) = 
            pd;
    }
    copyPartition();
    agglp.level = level;
    elp.level = level; 
    presolvelp.level = level;
    // elp.pairs = pairs;
    ++level;
    return lift_solution;
}

void HighsOCAggregate::buildSolution(){
    // ep = partition;
    // solution = s;
    // pcol = col;
    // pcolrep = colrep;
    // prow = row;
    // prowrep = rowrep;
    // pFrontCol = frontCol;
    // pcolCnt = colCnt;
    // pFrontRow = frontRow;
    // prowCnt = rowCnt;
    // findFrontMins();
    // buildColPointers();
    // buildRowPointers();
    HighsInt iCol, iRow, rep, cf, pcf, pc,
    rf, prf, pr, prlen; 
    double pv, pd;
    lift_solution.col_value.resize(numCol);
    lift_solution.row_value.resize(numRow);
    lift_solution.col_dual.resize(numCol);
    lift_solution.row_dual.resize(numRow);
    for (iCol = 0; iCol < colCnt; ++iCol){
        rep = colrep[iCol];
        cf = ep.front[rep];
        pcf = epMinusOne.front[rep];
        pc = pFrontCol[pcf];
        pv = solution.col_value[pc];
        pd = solution.col_dual[pc];
        lift_solution.col_value.at(iCol) = 
            pv;
        lift_solution.col_dual.at(iCol) = 
            pd;
    }
    for (iRow = 0; iRow < rowCnt; ++iRow){
        rep = rowrep[iRow];
        rf = ep.front[rep];
        prf = epMinusOne.front[rep];
        prlen = epMinusOne.len[prf] + 1;
        pr = pFrontRow[prf];
        pv = solution.row_value.at(pr);
        pd = solution.row_dual.at(pr);
        lift_solution.row_value.at(iRow) = 
            (double)pv/prlen;
        lift_solution.row_dual.at(iRow) = 
            pd;
    }
    // copyPartition();
    agglp.level = level;
    elp.level = level; 
    presolvelp.level = level;
    // elp.pairs = pairs;
    ++level;
}

// Build A matrix without linkers if they exist
void HighsOCAggregate::buildAmatrix(){
    int iCol, cf, crep, clen, c;
    int mIdx, ro, rf, r, rlen;
    int nnzStan = 0, start = 0, startStan = 0;
    nnz = 0;
    double xlen, xv;
    std::set<int>::iterator it;
    agglp.a_matrix_.start_.push_back(0);
    for (iCol = 0; iCol < colCnt; ++iCol){
        crep = colrep[iCol];
        cf = colFront[iCol];
        clen = ep.len[cf] + 1;
        c = iCol;
        for (mIdx = olp.a_matrix_.start_[crep]; mIdx < olp.a_matrix_.start_[crep + 1]; ++mIdx){
            ro = olp.a_matrix_.index_[mIdx];
            rf = ep.front[ro + numCol];
            r = frontRow[rf];
            rlen = ep.len[rf] + 1;
            if (!columnF[r]++){ 
                if (nnz >= columnI.size())
                    std::cin.get();
                columnI[nnz++] = r;
                // nnzStan++;
            }
            columnX[r] += olp.a_matrix_.value_[mIdx] * clen;
        }
        for (mIdx = start; mIdx < nnz; ++mIdx){
            agglp.a_matrix_.index_.push_back(columnI[mIdx]);
            // AdegenRIndex[j] = columnI[j];
            agglp.a_matrix_.value_.push_back(columnX[columnI[mIdx]]);
            // AdegenRValue[j] = columnX[columnI[j]];
            columnF[columnI[mIdx]] = 0;
            columnX[columnI[mIdx]] = 0;
            columnI[mIdx] = 0;
        }
        agglp.a_matrix_.start_.push_back(nnz);
        // AdegenRStart[xi + 1] = nnzStan;
        start = nnz;
    }
    agglp.a_matrix_.format_ = MatrixFormat::kColwise;
    agglp.a_matrix_.num_col_ = colCnt;
    agglp.a_matrix_.num_row_ = rowCnt;
    agglp.num_col_ = agglp.a_matrix_.num_col_ = colCnt;
    agglp.num_row_ = agglp.a_matrix_.num_row_ = rowCnt;
    agglp.num_aggregate_cols_ = agglp.num_col_;
    agglp.num_aggregate_rows_ = agglp.num_row_;
    agglp.num_residual_cols_ = 0;
    agglp.num_residual_rows_ = 0;
}

// Build A matrix with linkers if they exist
void HighsOCAggregate::buildAmatrixExtended(){
    int iCol, cf, crep, clen, c;
    int mIdx, ro, rf, r, rlen;
    int nnzStan = 0, start = 0, startStan = 0;
    nnz = 0;
    double xlen, xv;
    std::set<int>::iterator it;
    elp.a_matrix_.start_.push_back(0);
    for (iCol = 0; iCol < colCnt; ++iCol){
        crep = colrep[iCol];
        cf = colFront[iCol];
        clen = ep.len[cf] + 1;
        c = iCol;
        for (mIdx = olp.a_matrix_.start_[crep]; mIdx < olp.a_matrix_.start_[crep + 1]; ++mIdx){
            ro = olp.a_matrix_.index_[mIdx];
            xv = olp.a_matrix_.value_[mIdx];
            // if (std::fabs(xv) < 1e-9) continue;
            rf = ep.front[ro + numCol];
            r = frontRow[rf];
            rlen = ep.len[rf] + 1;
            if (!columnF[r]++){ 
                columnI[nnz++] = r;
                // nnzStan++;
            }
            columnX[r] += xv * clen;
        }
        // int new_nnz = start;
        for (mIdx = start; mIdx < nnz; mIdx++){
            // if (std::fabs(columnX[columnI[mIdx]]) < 1e-9)
            elp.a_matrix_.index_.push_back(columnI[mIdx]);
            // AdegenRIndex[j] = columnI[j];
            elp.a_matrix_.value_.push_back(columnX[columnI[mIdx]]);
            // AdegenRValue[j] = columnX[columnI[j]];
            columnF[columnI[mIdx]] = 0;
            columnX[columnI[mIdx]] = 0;
            columnI[mIdx] = 0;
        }
        // Add linking row entry for parent columns
        if (isParent[c]){
            for (mIdx = parentRow[c]; mIdx < parentRow[c + 1]; ++mIdx){
                elp.a_matrix_.index_.push_back(mIdx);
                // AdegenRIndex[nnzStan] = j;
                elp.a_matrix_.value_.push_back(1);
                ++nnz;
                // AdegenRValue[nnzStan++] = 1;
            }
        } 
        // Add linking row entry for child columns
        if (isChild[c]){
            elp.a_matrix_.index_.push_back(childRow[c]);
            // AdegenRIndex[nnzStan] = childRow[xi];
            elp.a_matrix_.value_.push_back(-1);
            // AdegenRValue[nnzStan++] = -1;
            ++nnz;
        }
        elp.a_matrix_.start_.push_back(nnz);
        // AdegenRStart[xi + 1] = nnzStan;
        start = nnz;
    }
    // Add linking row entry for residual columns
    for (iCol = 0; iCol < numResiduals; ++iCol){
        elp.a_matrix_.index_.push_back(residualRow[iCol]);
        // AdegenRIndex[nnzStan] = residualRow[i];
        elp.a_matrix_.value_.push_back(-1);
        // AdegenRValue[nnzStan++] = -1;
        elp.a_matrix_.start_.push_back(++nnz);
        // AdegenRStart[elp.num_col_ + i + 1] = nnzStan;
    }
    elp.a_matrix_.format_ = MatrixFormat::kColwise;
    elp.a_matrix_.num_col_ = colCnt + numResiduals;
    // elp.a_matrix_.num_col_ = colCnt;
    elp.a_matrix_.num_row_ = rowCnt + numResiduals;
    elp.num_col_ = colCnt + numResiduals;
    // elp.num_col_ = colCnt;
    elp.num_row_ = rowCnt + numResiduals;
    elp.num_aggregate_cols_ = colCnt;
    elp.num_aggregate_rows_ = rowCnt;
    // elp.num_residual_cols_ = 0;
    elp.num_residual_cols_ = numResiduals;
    elp.num_residual_rows_ = numResiduals;
    elp.is_degenerate_residual = zero_step_pivots;
    // elp.residual_cols_ = residualCol;
}

// Build A matrix with linkers if they exist
void HighsOCAggregate::buildAmatrixExtendedNoResiduals(){
    int iCol, cf, crep, clen, c;
    int mIdx, ro, rf, r, rlen;
    int nnzStan = 0, start = 0, startStan = 0;
    nnz = 0;
    double xlen, xv;
    std::set<int>::iterator it;
    presolvelp.a_matrix_.start_.push_back(0);
    for (iCol = 0; iCol < colCnt; ++iCol){
        crep = colrep[iCol];
        cf = colFront[iCol];
        clen = ep.len[cf] + 1;
        c = iCol;
        for (mIdx = olp.a_matrix_.start_[crep]; mIdx < olp.a_matrix_.start_[crep + 1]; ++mIdx){
            ro = olp.a_matrix_.index_[mIdx];
            rf = ep.front[ro + numCol];
            r = frontRow[rf];
            rlen = ep.len[rf] + 1;
            if (!columnF[r]++){ 
                columnI[nnz++] = r;
                // nnzStan++;
            }
            columnX[r] += olp.a_matrix_.value_[mIdx] * clen;
        }
        for (mIdx = start; mIdx < nnz; ++mIdx){
            presolvelp.a_matrix_.index_.push_back(columnI[mIdx]);
            // AdegenRIndex[j] = columnI[j];
            presolvelp.a_matrix_.value_.push_back(columnX[columnI[mIdx]]);
            // AdegenRValue[j] = columnX[columnI[j]];
            columnF[columnI[mIdx]] = 0;
            columnX[columnI[mIdx]] = 0;
            columnI[mIdx] = 0;
        }
        // Add linking row entry for parent columns
        if (isParent[c]){
            for (mIdx = parentRow[c]; mIdx < parentRow[c + 1]; ++mIdx){
                presolvelp.a_matrix_.index_.push_back(mIdx);
                // AdegenRIndex[nnzStan] = j;
                presolvelp.a_matrix_.value_.push_back(1);
                ++nnz;
                // AdegenRValue[nnzStan++] = 1;
            }
        } 
        // Add linking row entry for child columns
        if (isChild[c]){
            presolvelp.a_matrix_.index_.push_back(childRow[c]);
            // AdegenRIndex[nnzStan] = childRow[xi];
            presolvelp.a_matrix_.value_.push_back(-1);
            // AdegenRValue[nnzStan++] = -1;
            ++nnz;
        }
        presolvelp.a_matrix_.start_.push_back(nnz);
        // AdegenRStart[xi + 1] = nnzStan;
        start = nnz;
    }
    presolvelp.a_matrix_.format_ = MatrixFormat::kColwise;
    presolvelp.a_matrix_.num_col_ = colCnt;
    presolvelp.a_matrix_.num_row_ = rowCnt + numResiduals;
    presolvelp.num_col_ = colCnt;
    presolvelp.num_row_ = rowCnt + numResiduals;
    presolvelp.num_aggregate_cols_ = colCnt;
    presolvelp.num_aggregate_rows_ = rowCnt;
    presolvelp.num_residual_cols_ = 0;
    presolvelp.num_residual_rows_ = numResiduals;
}

void HighsOCAggregate::checkAMatrix(){
    int iCol, crep, iRow, rrep, aIndex;
    std::vector<int> compareElp;
    std::vector<int> compareOlp(numRow);
    for (iCol = 0; iCol < elp.num_aggregate_cols_; ++iCol){
        int numColEntries = 0;
        compareElp.clear();
        crep = colrep[iCol];
        for (aIndex = elp.a_matrix_.start_[iCol]; aIndex < elp.a_matrix_.start_[iCol + 1]; ++aIndex){
            iRow = elp.a_matrix_.index_[aIndex];
            if (iRow >= elp.num_aggregate_rows_) continue;
            rrep = rowrep[iRow] - numCol;
            compareElp.push_back(rrep);
            numColEntries++;
        }
        for (aIndex = olp.a_matrix_.start_[crep]; aIndex < olp.a_matrix_.start_[crep + 1]; ++aIndex){
            iRow = olp.a_matrix_.index_[aIndex];
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

void HighsOCAggregate::printAMatrixToMatlabFormat(HighsSparseMatrix& matrix){
    // A Matrix
    std::vector<std::vector<double> > printMat;
    for (int iCol = 0; iCol < matrix.num_col_; ++iCol){
        std::vector<double> column(matrix.num_row_, 0);
        if (iCol < matrix.num_col_){
            for (int idx = matrix.start_[iCol]; idx < matrix.start_[iCol + 1]; ++idx){
                int iRow = matrix.index_[idx];
                int iRowValue = matrix.value_[idx];
                column[iRow] = iRowValue;
            }
        }
        else{
            int iRow = iCol - matrix.num_col_;
            column[iRow] = -1;
        }
        printMat.push_back(column);
    }
    std::cout << "Standard Form A Matrix" << std::endl;
    std::cout << "[ ";
    for (int iCol = 0; iCol < matrix.num_col_; ++iCol){
        for (int iRow = 0; iRow < matrix.num_row_; ++iRow){
            std::cout << printMat[iCol][iRow] << " ";
        }
        std::cout << ";" << std::endl;
    }
    std::cout << "]" << std::endl;
}

void HighsOCAggregate::buildObj(){
    int iCol, rep, cf;
    double c, clen;
    agglp.col_cost_.resize(colCnt);
    for (iCol = 0; iCol < colCnt; ++iCol){
        rep = colrep[iCol];
        cf = colFront[iCol];
        clen = ep.len[cf] + 1;
        c = olp.col_cost_[rep];
        agglp.col_cost_[iCol] = clen * c;
    }
    agglp.sense_ = olp.sense_;
}

void HighsOCAggregate::buildObjExtended(){
    int iCol, rep, cf;
    double c, clen;
    for (iCol = 0; iCol < colCnt; ++iCol){
        rep = colrep[iCol];
        cf = colFront[iCol];
        clen = ep.len[cf] + 1;
        c = olp.col_cost_[rep];
        elp.col_cost_[iCol] = clen * c;
    }
    elp.sense_ = olp.sense_;
}

void HighsOCAggregate::buildObjExtendedNoResiduals(){
    int iCol, rep, cf;
    double c, clen;
    presolvelp.col_cost_.resize(colCnt);
    for (iCol = 0; iCol < colCnt; ++iCol){
        rep = colrep[iCol];
        cf = colFront[iCol];
        clen = ep.len[cf] + 1;
        c = olp.col_cost_[rep];
        presolvelp.col_cost_[iCol] = clen * c;
    }
    presolvelp.sense_ = olp.sense_;
}

void HighsOCAggregate::buildRhs(){
    if (ep.level)
        buildRhsFromSolution();
    else
        buildRhsFromScratch();
}

void HighsOCAggregate::buildRhsExtended(){
    if (ep.level)
        buildRhsFromSolutionExtended();
    else
        buildRhsFromScratchExtended();
}

void HighsOCAggregate::buildRhsExtendedNoResiduals(){
    if (ep.level)
        buildRhsFromSolutionExtendedNoResiduals();
}

void HighsOCAggregate::buildRhsFromScratch(){
    int iRow, rf, rrep;
    double rlen;
    agglp.row_upper_.resize(rowCnt);
    agglp.row_lower_.resize(rowCnt);
    for (iRow = 0; iRow < rowCnt; ++iRow){
        rrep = rowrep[iRow];
        rf = rowFront[iRow];
        rlen = ep.len[rf] + 1;
        agglp.row_lower_[iRow] = 
            olp.row_lower_[rrep - numCol] * rlen;
        agglp.row_upper_[iRow] = 
            olp.row_upper_[rrep - numCol] * rlen;
    }
}

void HighsOCAggregate::buildRhsFromScratchExtended(){
    int iRow, rf, rrep;
    double rlen;
    for (iRow = 0; iRow < rowCnt; ++iRow){
        rrep = rowrep[iRow];
        rf = rowFront[iRow];
        rlen = ep.len[rf] + 1;
        elp.row_lower_[iRow] = 
            olp.row_lower_[rrep - numCol] * rlen;
        elp.row_upper_[iRow] = 
            olp.row_upper_[rrep - numCol] * rlen;
    }
}

void HighsOCAggregate::buildRhsFromSolution(){
    int iRow, pr, rf, rpf, rlen, prlen, rep;
    double pv, lb, ub, lbDiff, ubDiff, tol = 1e-6;
    agglp.row_upper_.resize(rowCnt);
    agglp.row_lower_.resize(rowCnt);
    for (iRow = 0; iRow < rowCnt; ++iRow){
        rep = rowrep[iRow];
        rf = ep.front[rep];
        rpf = epMinusOne.front[rep];
        rlen = ep.len[rf] + 1;
        pr = pFrontRow[rpf];
        prlen = epMinusOne.len[rpf] + 1;
        pv = (double)solution.row_value[pr]/prlen;
        lb = olp.row_lower_[rep - numCol];
        ub = olp.row_upper_[rep - numCol];
        lbDiff = std::fabs(pv - lb);
        ubDiff = std::fabs(pv - ub);
        if (ubDiff < tol){
            agglp.row_lower_[iRow] = ub * rlen;
            agglp.row_upper_[iRow] = ub * rlen;
        }
        else if (lbDiff < tol){
            agglp.row_lower_[iRow] = lb * rlen;
            agglp.row_upper_[iRow] = lb * rlen;
        }
        else{
            agglp.row_lower_[iRow] = lb * rlen;
            agglp.row_upper_[iRow] = ub * rlen;
        }
    }
}

void HighsOCAggregate::buildRhsFromSolutionExtended(){
    int iRow, pr, rf, rpf, rlen, prlen, rep;
    double pv, lb, ub, lbDiff, ubDiff, tol = 1e-6;
    for (iRow = 0; iRow < rowCnt; ++iRow){
        rep = rowrep[iRow];
        rf = ep.front[rep];
        rpf = epMinusOne.front[rep];
        rlen = ep.len[rf] + 1;
        pr = pFrontRow[rpf];
        prlen = epMinusOne.len[rpf] + 1;
        pv = (double)solution.row_value[pr]/prlen;
        lb = olp.row_lower_[rep - numCol];
        ub = olp.row_upper_[rep - numCol];
        lbDiff = std::fabs(pv - lb);
        ubDiff = std::fabs(pv - ub);
        if (ubDiff < tol){
            elp.row_lower_[iRow] = ub * rlen;
            elp.row_upper_[iRow] = ub * rlen;
        }
        else if (lbDiff < tol){
            elp.row_lower_[iRow] = lb * rlen;
            elp.row_upper_[iRow] = lb * rlen;
        }
        else{
            elp.row_lower_[iRow] = lb * rlen;
            elp.row_upper_[iRow] = ub * rlen;
        }
    }
    for (iRow = rowCnt; iRow < elp.num_row_; ++iRow){
        elp.row_lower_[iRow] = 0;
        elp.row_upper_[iRow] = 0;
    }
}

void HighsOCAggregate::buildRhsFromSolutionExtendedNoResiduals(){
    int iRow, pr, rf, rpf, rlen, prlen, rep;
    double pv, lb, ub, lbDiff, ubDiff, tol = 1e-6;
    for (iRow = 0; iRow < rowCnt; ++iRow){
        rep = rowrep[iRow];
        rf = ep.front[rep];
        rpf = epMinusOne.front[rep];
        rlen = ep.len[rf] + 1;
        pr = pFrontRow[rpf];
        prlen = epMinusOne.len[rpf] + 1;
        pv = (double)solution.row_value[pr]/prlen;
        lb = olp.row_lower_[rep - numCol];
        ub = olp.row_upper_[rep - numCol];
        lbDiff = std::fabs(pv - lb);
        ubDiff = std::fabs(pv - ub);
        if (ubDiff < tol){
            presolvelp.row_lower_[iRow] = ub * rlen;
            presolvelp.row_upper_[iRow] = ub * rlen;
        }
        else if (lbDiff < tol){
            presolvelp.row_lower_[iRow] = lb * rlen;
            presolvelp.row_upper_[iRow] = lb * rlen;
        }
        else{
            presolvelp.row_lower_[iRow] = lb * rlen;
            presolvelp.row_upper_[iRow] = ub * rlen;
        }
    }
    for (iRow = rowCnt; iRow < presolvelp.num_row_; ++iRow){
        presolvelp.row_lower_[iRow] = 0;
        presolvelp.row_upper_[iRow] = 0;
    }
}

void HighsOCAggregate::buildRowNames(){
    HighsInt min_rep, row_front, i_row;
    std::stringstream name;
    elp.row_names_.clear();
    for (i_row = 0; i_row < rowCnt; ++i_row){
        row_front = rowFront.at(i_row);
        min_rep = frontMin.at(row_front);
        name << "R" << min_rep;
        elp.row_names_.push_back(name.str());
        name.str(""); 
    }
}

void HighsOCAggregate::buildBnds(){
    if (ep.level)
        buildBndsFromSolution();
    else
        buildBndsFromScratch();
}

void HighsOCAggregate::buildBndsExtended(){
    if (ep.level)
        buildBndsFromSolutionExtended();
    else
        buildBndsFromScratchExtended();
}

void HighsOCAggregate::buildBndsExtendedNoResiduals(){
    if (ep.level)
        buildBndsFromSolutionExtendedNoResiduals();
}

void HighsOCAggregate::buildBndsFromScratch(){
    int iCol, crep;
    agglp.col_lower_.resize(colCnt);
    agglp.col_upper_.resize(colCnt);
    for (iCol = 0; iCol < colCnt; ++iCol){
        crep = colrep[iCol];
        agglp.col_lower_[iCol] = 
            olp.col_lower_[crep];
        agglp.col_upper_[iCol] = 
            olp.col_upper_[crep];
    }
}

void HighsOCAggregate::buildBndsFromScratchExtended(){
    int iCol, crep;
    for (iCol = 0; iCol < colCnt; ++iCol){
        crep = colrep[iCol];
        elp.col_lower_[iCol] = 
            olp.col_lower_[crep];
        elp.col_upper_[iCol] = 
            olp.col_upper_[crep];
    }
}

void HighsOCAggregate::buildBndsFromSolution(){
    int iCol, pc, cf, pcf, rep;
    double pv, ubDiff, lbDiff, lb, ub, tol = 1e-6;
    agglp.col_lower_.resize(colCnt);
    agglp.col_upper_.resize(colCnt);
    for (iCol = 0; iCol < colCnt; ++iCol){
        rep = colrep[iCol];
        cf = ep.front[rep];
        pcf = epMinusOne.front[rep];
        pc = pFrontCol[pcf];
        pv = solution.col_value[pc];
        ub = olp.col_upper_[rep];
        lb = olp.col_lower_[rep];
        ubDiff = std::fabs(pv - ub);
        lbDiff = std::fabs(pv - lb);
        if (ubDiff < tol){
            agglp.col_lower_[iCol] = ub;
            agglp.col_upper_[iCol] = ub;
        }
        else if (lbDiff < tol){
            agglp.col_lower_[iCol] = lb;
            agglp.col_upper_[iCol] = lb;
        }
        else{
            agglp.col_lower_[iCol] = lb;
            agglp.col_upper_[iCol] = ub;
        }
    }
}

void HighsOCAggregate::buildBndsFromSolutionExtended(){
    int iCol, pc, cf, pcf, rep;
    double pv, ubDiff, lbDiff, lb, ub, tol = 1e-6;
    for (iCol = 0; iCol < colCnt; ++iCol){
        rep = colrep[iCol];
        cf = ep.front[rep];
        pcf = epMinusOne.front[rep];
        pc = pFrontCol[pcf];
        pv = solution.col_value[pc];
        ub = olp.col_upper_[rep];
        lb = olp.col_lower_[rep];
        ubDiff = std::fabs(pv - ub);
        lbDiff = std::fabs(pv - lb);
        if (ubDiff < tol){
            elp.col_lower_[iCol] = ub;
            elp.col_upper_[iCol] = ub;
        }
        else if (lbDiff < tol){
            elp.col_lower_[iCol] = lb;
            elp.col_upper_[iCol] = lb;
        }
        else{
            elp.col_lower_[iCol] = lb;
            elp.col_upper_[iCol] = ub;
        }
    }
    for (iCol = colCnt; iCol < elp.num_col_; ++iCol){
        if (elpBasis.col_status.at(iCol) == HighsBasisStatus::kBasic){
            elp.col_lower_[iCol] = -kHighsInf;
            elp.col_upper_[iCol] = kHighsInf;
            continue;
        }
        elp.col_lower_[iCol] = 0;
        elp.col_upper_[iCol] = 0;
    }
}  

void HighsOCAggregate::buildBndsFromSolutionExtendedNoResiduals(){
    int iCol, pc, cf, pcf, rep;
    double pv, ubDiff, lbDiff, lb, ub, tol = 1e-6;
    for (iCol = 0; iCol < colCnt; ++iCol){
        rep = colrep[iCol];
        cf = ep.front[rep];
        pcf = epMinusOne.front[rep];
        pc = pFrontCol[pcf];
        pv = solution.col_value[pc];
        ub = olp.col_upper_[rep];
        lb = olp.col_lower_[rep];
        ubDiff = std::fabs(pv - ub);
        lbDiff = std::fabs(pv - lb);
        if (ubDiff < tol){
            presolvelp.col_lower_[iCol] = ub;
            presolvelp.col_upper_[iCol] = ub;
        }
        else if (lbDiff < tol){
            presolvelp.col_lower_[iCol] = lb;
            presolvelp.col_upper_[iCol] = lb;
        }
        else{
            presolvelp.col_lower_[iCol] = lb;
            presolvelp.col_upper_[iCol] = ub;
        }
    }
} 

void HighsOCAggregate::buildColNames(){
    HighsInt min_rep, col_front, i_col;
    std::stringstream name;
    elp.col_names_.clear();
    for (i_col = 0; i_col < colCnt; ++i_col){
        col_front = colFront.at(i_col);
        min_rep = frontMin.at(col_front);
        name << "C" << min_rep;
        elp.col_names_.push_back(name.str());
        name.str(""); 
    }
}
void HighsOCAggregate::buildResiduals(){
    if (!ep.level) return;
    buildResidualLinks();
    buildResidualCols();
    buildResidualRows();
}

void HighsOCAggregate::buildResidualLinks(){
    int i, x1, x2, x0, nf, of, newNumRow = elp.num_row_;
    int rep, xOld, xNew, pcf, minParent = numCol, oldRep, newRep;
    HighsBasisStatus basic = HighsBasisStatus::kBasic;
    std::pair<std::set<std::pair<int, int> >::iterator, bool> ret;
    numResiduals = 0;
    numParent = 0;
    numChild = 0;
    std::vector<HighsInt> temp_piv_need(col_needs_res_pivot.begin(), col_needs_res_pivot.end());
    std::fill(isParent.begin(), isParent.end(), 0);
    std::fill(parentRow.begin(), parentRow.end(), 0);
    std::fill(isChild.begin(), isChild.end(), 0);
    std::fill(parent.begin(), parent.end(), 0);
    std::fill(childRow.begin(), childRow.end(), 0);
    child.resize(0);
    parent_to_child.clear();
    child_to_res.assign(colCnt, -1);
    // res_to_child.assign(2 * colCnt - 1, -1);
    std::fill(residualRow.begin(), residualRow.end(), 0);
    std::fill(residualCol.begin(), residualCol.end(), 0);
    // std::fill(degenerate_basic_rows.begin(), degenerate_basic_rows.end(), 0);
    // std::fill(degenerate_basic_index.begin(), degenerate_basic_index.end(), 0);
    // std::fill(degenerate_basic_residuals.begin(), degenerate_basic_residuals.end(), 0);
    std::fill(zero_step_pivots.begin(), zero_step_pivots.end(), 0);
    std::vector<int> residual_cols_;
    std::vector<int> residual_to_old(numCol + numTotResiduals);
    int elpNumRow = rowCnt;
    int elpNumCol = colCnt;
    int numSkipped = 0;
    // pairs.clear();
    elp.residual_cols_.clear();
    elp.num_degenerate_cols_ = 0;
    for (int iCol = pcolCnt; iCol < colCnt; ++iCol){
        xNew = iCol;
        rep = colrep[iCol];
        newRep = colrep[iCol];
        pcf = epMinusOne.front[rep];
        xOld = pFrontCol[pcf];
        oldRep = colrep[xOld];
        if (basis.col_status[xOld] != basic){ numSkipped++; continue; }
        // if (!temp_piv_need.at(xNew) && !temp_piv_need.at(xOld)){ numSkipped++; continue;}
        // if (!temp_piv_need.at(xNew) && temp_piv_need.at(xOld)) temp_piv_need.at(xOld) = 0;
        // if (mark_degenerate.at(xOld)) { numSkipped++; continue; }
        if (xOld < minParent) minParent = xOld;
        splitCells[xOld].push_back(xNew);
    }
    elpNumCol = colCnt;
    for (const auto split : splitCells){
        isParent.at(split.first) = 1;
        numParent++;
        parentRow[split.first] = elpNumRow;
        parentRow[split.first + 1] = elpNumRow + split.second.size();
        for (int idx = 0; idx < split.second.size(); ++idx){
            int iCol = split.second.at(idx);
            // if (mark_degenerate.at(split.first)){
            //     zero_step_pivots.at(elpNumCol) = 1;
            //     elp.num_degenerate_cols_++;
            // }
            parent.at(iCol) = split.first;
            parent_to_child[split.first].push_back(iCol);
            isChild.at(iCol) = 1;
            numChild++;
            child.push_back(iCol);
            childRow[iCol] = elpNumRow;
            residual_cols_.push_back(elpNumCol);
            residual_to_old.at(elpNumCol) = split.first;
            child_to_res.at(iCol) = elpNumCol;
            // res_to_child.at(elpNumCol) = iCol;
            residualCol[numResiduals] = elpNumCol++;
            residualRow[numResiduals++] = elpNumRow++;
            // pairs.push_back(std::pair<int, int>(split.first, iCol));
        }
    }
    for (int iCol = 0; iCol < numResiduals; ++iCol){
        int r_col = residual_cols_.at(iCol);
        // int x_col = residual_to_old.at(r_col);
        // if (mark_degenerate.at(x_col))
        elp.residual_cols_.push_back(r_col);
    }
    // for (int iCol = 0; iCol < numResiduals; ++iCol){
    //     int r_col = residual_cols_.at(iCol);
    //     int x_col = residual_to_old.at(r_col);
    //     if (!mark_degenerate.at(x_col))
    //         elp.residual_cols_.push_back(r_col);
    // }
    splitCells.clear();
}

void HighsOCAggregate::buildResidualSubMatrix(){
    for (HighsInt i_row = 0; i_row < numResiduals; ++i_row){
        HighsInt i_child = child.at(i_row);
        test_basic_index.push_back(i_child);
        HighsInt i_parent = parent.at(i_child);
        HighsInt i_residual = residualCol.at(i_row);
        residual_ar.index_.push_back(i_parent);
        residual_ar.index_.push_back(i_child);
        residual_ar.index_.push_back(i_residual);
        residual_ar.value_.push_back(1);
        residual_ar.value_.push_back(-1);
        residual_ar.value_.push_back(-1);
        residual_ar.start_.push_back(residual_ar.index_.size());
    }
    residual_ar.num_col_ = numParent + numChild + numResiduals;
    residual_ar.num_row_ = numResiduals;
    residual_ar.format_ = MatrixFormat::kRowwise;
}

void HighsOCAggregate::buildResidualCols(){
    int i, p, c, idx = elp.num_col_;
    for (i = 0; i < numResiduals; ++i){
        p = parent[i];
        c = child[i];
        // elp.col_lower_[idx] = 0;
        // elp.col_upper_[idx++] = 0;
        elp.col_lower_[idx] = -kHighsInf;
        elp.col_upper_[idx++] = kHighsInf;
    }
}

void HighsOCAggregate::buildResidualRows(){
    int i, idx = elp.num_row_;
    for (i = 0; i < numResiduals; ++i){
        elp.row_lower_[idx] = elp.row_upper_[idx] = 0;
        idx++;
    }
}

void HighsOCAggregate::markDegenerate(){
    std::fill_n(mark_degenerate.begin(), mark_degenerate.size(), 0);
    HighsInt i_col;
    // elp.num_degenerate_cols_ = elp.num_aggregate_cols_;
    for (i_col = 0; i_col < pcolCnt; ++i_col){
        HighsInt crep = colrep.at(i_col);
        double lb = olp.col_lower_.at(crep);
        double ub = olp.col_upper_.at(crep);
        double value = solution.col_value.at(i_col);
        HighsInt ub_test = std::fabs(value - ub) < kHighsTiny ? 1 : 0;
        HighsInt lb_test = std::fabs(value - lb) < kHighsTiny ? 1 : 0;
        // std::cout << std::fabs(value) << std::endl;
        HighsInt basis_test = basis.col_status.at(i_col) == HighsBasisStatus::kBasic ? 1 : 0;
        if ((ub_test || lb_test) && basis_test){
            mark_degenerate.at(i_col) = 1;
            // elp.num_degenerate_cols_++;
        }
        // if ((ub_test) && basis_test)
        //     mark_degenerate.at(i_col) = 1;
    }
}

void HighsOCAggregate::buildTestBasisLU(HighsLpSolverObject& solver_object){
    //// Highs LU ///////
    // test_options.factor_pivot_tolerance = 1e-9;
    // HighsLpSolverObject solver_object(elp, elpBasis, lift_solution, test_info, 
    //                                 test_ekk, test_options, test_timer);
    HighsLp& so_lp = solver_object.lp_;
    HighsBasis& so_basis = solver_object.basis_;
    HighsOptions& so_options = solver_object.options_;
    HEkk& ekk_instance = solver_object.ekk_instance_;
    HighsLp& ekk_lp = ekk_instance.lp_;
    HighsSimplexStatus& ekk_status = ekk_instance.status_;
    ekk_instance.moveLp(solver_object);
    ekk_lp.residual_cols_.resize(0);
    ekk_instance.setBasis(elpBasis);
    ekk_instance.dual_edge_weight_.assign(elp.num_row_, 1.0);
    ekk_instance.scattered_dual_edge_weight_.resize(elp.num_col_ + elp.num_row_);
    ekk_instance.info_.backtracking_basis_edge_weight_.resize(elp.num_col_ + elp.num_row_);
    std::vector<HighsInt>& ekk_basic_index_start = ekk_instance.basis_.basicIndex_;
    std::vector<HighsInt> ekk_col_basic_row(colCnt + rowCnt + 2 * numResiduals);
    for (HighsInt i_basic_idx = 0; i_basic_idx < elp.num_row_; ++i_basic_idx){
        HighsInt i_var = ekk_basic_index_start.at(i_basic_idx);
        test_basic_start.at(i_var) = 1;
    }
    HighsInt rank_deficiency = ekk_instance.initialiseForAggregator();
    for (HighsInt i_basic_idx = 0; i_basic_idx < elp.num_row_; ++i_basic_idx){
        HighsInt i_var = ekk_basic_index_start.at(i_basic_idx);
        ekk_col_basic_row.at(i_var) = i_basic_idx;
    }
    std::vector<HighsInt> res_nonbasic(colCnt + numResiduals, -1);
    if (rank_deficiency){
        HVector col_aq;
        HVector row_ep;
        HighsInt rebuild = RebuildReason::kRebuildReasonNo;
        col_aq.setup(elp.num_row_);
        row_ep.setup(elp.num_row_);
        for (HighsInt i_basic_idx = 0; i_basic_idx < elp.num_row_; ++i_basic_idx){
            HighsInt i_var = ekk_basic_index_start.at(i_basic_idx);
            test_basic_finish.at(i_var) = 1;
        }
        for (HighsInt i_var = 0; i_var < elp.num_col_ + elp.num_row_; ++i_var){
            HighsInt was_basic = test_basic_start.at(i_var);
            HighsInt is_basic = test_basic_finish.at(i_var);
            if (was_basic && !is_basic){
                HighsInt i_residual;
                HighsInt i_parent;
                HighsInt i_child;
                HighsInt i_basic_row;
                HighsInt i_basic_slack;
                HighsInt i_pivot;
                HighsInt is_par = isParent.at(i_var);
                if (is_par){
                    continue;
                    // i_parent = i_var;
                    // std::vector<HighsInt>& children = parent_to_child.at(i_parent);
                    // for (auto chi : children){
                    //     i_residual = child_to_res.at(chi);
                    //     i_child = chi;
                    //     i_basic_row = ekk_col_basic_row.at(i_child);
                    //     i_basic_slack = i_basic_row + colCnt + numResiduals;
                    //     // Swap in slack var to basis so x child of x parent can become basic in residual row
                    //     ekk_instance.pivotColumnFtranForAggregator(i_basic_slack, col_aq); 
                    //     ekk_instance.unitBtranForAggregator(i_basic_row, row_ep);
                    //     ekk_instance.transformForUpdateForAggregator(&col_aq, &row_ep, i_basic_slack, &i_basic_row);
                    //     ekk_instance.updatePivotsForAggregator(i_basic_slack, i_basic_row, 1);
                    //     ekk_instance.updateFactorForAggregator(&col_aq, &row_ep, &i_basic_row, &rebuild);
                    //     ekk_instance.updateMatrixForAggregator(i_basic_slack, i_child);
                    //     test_basic_finish.at(i_child) = 1;
                    //     // Swap in x that we just removed to residual rows and take out corresponding r's
                    //     i_basic_row = residualRow.at(i_residual - colCnt);
                    //     ekk_instance.pivotColumnFtranForAggregator(i_child, col_aq); 
                    //     ekk_instance.unitBtranForAggregator(i_basic_row, row_ep);
                    //     ekk_instance.transformForUpdateForAggregator(&col_aq, &row_ep, i_child, &i_basic_row);
                    //     ekk_instance.updatePivotsForAggregator(i_child, i_basic_row, 1);
                    //     ekk_instance.updateFactorForAggregator(&col_aq, &row_ep, &i_basic_row, &rebuild);
                    //     ekk_instance.updateMatrixForAggregator(i_child, i_residual);
                    //     ekk_lp.residual_cols_.push_back(i_residual);
                    //     elpBasis.col_status.at(i_residual) = HighsBasisStatus::kLower;
                    //     ekk_lp.col_lower_.at(i_residual) = 0;
                    //     ekk_lp.col_upper_.at(i_residual) = 0;
                    //     res_nonbasic.at(i_residual) = 1;
                    // }
                }
                else{
                    i_child = i_var;
                    i_residual = child_to_res.at(i_child);
                    HighsInt i_row = residualRow.at(i_residual - colCnt);
                    ekk_instance.pivotColumnFtranForAggregator(i_child, col_aq); 
                    ekk_instance.unitBtranForAggregator(i_row, row_ep);
                    ekk_instance.transformForUpdateForAggregator(&col_aq, &row_ep, i_child, &i_row);
                    ekk_instance.updatePivotsForAggregator(i_child, i_row, 1);
                    ekk_instance.updateFactorForAggregator(&col_aq, &row_ep, &i_row, &rebuild);
                    ekk_instance.updateMatrixForAggregator(i_child, i_residual);
                    ekk_lp.residual_cols_.push_back(i_residual);
                    elpBasis.col_status.at(i_residual) = HighsBasisStatus::kLower;
                    ekk_lp.col_lower_.at(i_residual) = 0;
                    ekk_lp.col_upper_.at(i_residual) = 0;
                    res_nonbasic.at(i_residual) = 1;
                }
            }
        }
        ekk_instance.status_.has_invert = true;
        ekk_instance.status_.has_fresh_invert = true;
        bool is_now_full_rank = ekk_instance.getNonsingularInverseForAggregator();
        ekk_instance.computePrimalForAggregator();
        ekk_instance.computeSimplexPrimalInfeasibleForAggregator();
        std::cout << "is now full rank: " << is_now_full_rank << std::endl;
    }
    for (HighsInt i_res = colCnt; i_res < colCnt + numResiduals; ++i_res){
        if (elpBasis.col_status.at(i_res) != HighsBasisStatus::kBasic &&
            res_nonbasic.at(i_res) != 1)
            ekk_lp.residual_cols_.push_back(i_res);
    }
}

void HighsOCAggregate::updateTestBasisLU(){
    test_factor.addCols(numResiduals);
    test_factor.addRows(&residual_ar);
    test_factor.setup(elp.a_matrix_, test_basic_index);
    test_factor.build();
}

void HighsOCAggregate::buildTestBasicIndex(){
    HighsInt num_required = elp.num_row_;
    HighsInt num_have = 0;
    test_basic_index.resize(0);
    for (HighsInt i_col = 0; i_col < elp.num_col_; ++i_col){
        if (num_have == num_required) break;
        if (i_col < colCnt && elpBasis.col_status.at(i_col) == HighsBasisStatus::kBasic){
            test_basic_index.push_back(i_col);
            test_basic_start.at(i_col) = 1;
            num_have++;
        }
        else if (i_col >= colCnt){
            test_basic_index.push_back(i_col);
            test_basic_start.at(i_col) = 1;
            num_have++;
        }
    }
    for (HighsInt i_row = 0; i_row < prowCnt; ++i_row){
        if (num_have == num_required) break;
        if (elpBasis.row_status.at(i_row) == HighsBasisStatus::kBasic){
            test_basic_index.push_back(i_row + elp.num_col_);
            test_basic_start.at(i_row + elp.num_col_) = 1;
            num_have++;
        }
    }
}

std::vector<HighsInt> HighsOCAggregate::sortColWeights(std::vector<double>& weights){
    std::vector<HighsInt> perm(numCol + numRow);
    for (HighsInt i = 0; i < numCol + numRow; i++) perm[i] = i;

    pdqsort(perm.begin(), perm.end(), [&](HighsInt i, HighsInt j) {
        return std::make_pair(weights[i], i) > std::make_pair(weights[j], j);
    });

    return perm;
}

void HighsOCAggregate::fillTestBasicStart(){
    for (HighsInt idx = 0; idx < test_basic_index.size(); ++idx){
        HighsInt i_col = test_basic_index.at(idx);
        test_basic_start.at(i_col) = true;
    }
    // std::vector<HighsInt> swapped_in(agglp.a_matrix_.num_col_, 0);
    // for (HighsInt idx = 0; idx < test_basic_index.size(); ++idx){
    //     HighsInt i_col = test_basic_index.at(idx);
    //     if (i_col >= agglp.a_matrix_.num_col_) continue;
    //     HighsInt rep = colrep.at(i_col);
    //     HighsInt pcf = epMinusOne.front.at(rep);
    //     HighsInt i_old = pFrontCol.at(pcf);
    //     if (i_col == i_old) continue;
    //     if (!test_basic_start.at(i_old) && !swapped_in.at(i_old)++){
    //         test_basic_start.at(i_old) = 1;
    //         test_basic_start.at(i_col) = 0;
    //         test_basic_index.at(idx) = i_old;
    //     }
    // }
}

void HighsOCAggregate::fillTestBasicFinish(){
    for (HighsInt idx = 0; idx < test_basic_index.size(); ++idx){
        HighsInt i_col = test_basic_index.at(idx);
        if (i_col < colCnt + rowCnt)
            test_basic_finish.at(i_col) = true;
    }
}

void HighsOCAggregate::markDependentColumns(){
    for (HighsInt i_col = 0; i_col < colCnt; ++i_col){
        HighsInt start = test_basic_start.at(i_col);
        HighsInt finish = test_basic_finish.at(i_col);
        HighsInt rep = colrep.at(i_col);
        double lb = olp.col_lower_.at(rep);
        double ub = olp.col_upper_.at(rep);
        double primal_val = lift_solution.col_value.at(i_col);
        HighsInt lb_test = std::abs(primal_val - lb) < 1e-6;
        HighsInt ub_test = std::abs(primal_val - ub) < 1e-6;
        HighsInt bounded = lb_test || ub_test;
        if (start == 1 && finish != 1){
            col_needs_res_pivot.at(i_col) = 1;
            continue;
        }
        if (!start && !bounded){
            col_needs_res_pivot.at(i_col) = 1;
            continue;
        }
    }
}

void HighsOCAggregate::buildBasis(bool finish, bool extended){
    num_basic = 0;
    // buildRowBasis();
    // buildColBasis();
    buildTestColBasis();
    buildTestRowBasis();
    elpBasis.alien = false;
}

void HighsOCAggregate::changeBasis(){
    num_basic = 0;
    HighsInt new_num_col = agglp.a_matrix_.num_col_ + numResiduals;
    HighsInt new_num_row = agglp.a_matrix_.num_row_ + numResiduals;
    HighsInt base_size = agglp.a_matrix_.num_row_;
    elpBasis.col_status.resize(new_num_col);
    elpBasis.row_status.resize(new_num_row);
    std::fill(elpBasis.col_status.begin(), elpBasis.col_status.end(), HighsBasisStatus::kLower);
    std::fill(elpBasis.row_status.begin(), elpBasis.row_status.end(), HighsBasisStatus::kLower);
    HighsInt i_col, i_row, new_i_col, new_i_row, base_idx;
    // for (base_idx = 0; base_idx < base_size; ++base_idx){
    //     i_col = test_basic_index.at(base_idx);
    //     if (i_col < agglp.a_matrix_.num_col_){
    //         elpBasis.col_status.at(i_col) = HighsBasisStatus::kBasic;
    //         num_basic++;
    //     }
    //     else{
    //         i_row = i_col - agglp.a_matrix_.num_col_;
    //         // new_i_col = i_row + new_num_col;
    //         // new_i_row = new_i_col - new_num_col;
    //         elpBasis.row_status.at(i_row) = HighsBasisStatus::kBasic;
    //         num_basic++;
    //     }
    // }
    // for (i_col = 0; i_col < agglp.a_matrix_.num_col_; ++i_col){
    //     if (col_needs_res_pivot.at(i_col)){
    //         elpBasis.col_status.at(i_col) = HighsBasisStatus::kBasic;
    //         num_basic++;
    //     }
    // }
    // int num_bas = 0;
    // int num_set_bas = 0;
    // int num_sbas = 0;
    // int num_set_sbas = 0;
    // for (HighsInt i_col = 0; i_col < agglp.num_col_ + agglp.num_row_; ++i_col){
    //     bool basic = ipx_basis_->IsBasic(i_col);
    //     if (basic && i_col < agglp.num_col_) ++num_bas;
    //     else if (basic) ++num_sbas;
    //     HighsInt super_basic = col_needs_res_pivot.at(i_col);
    //     if ((basic || super_basic) && i_col < agglp.num_col_){
    //         elpBasis.col_status.at(i_col) = HighsBasisStatus::kBasic;
    //         num_set_bas++;
    //         num_basic++;
    //     }
    //     else if (basic){
    //         HighsInt i_row = i_col - agglp.num_col_;
    //         elpBasis.row_status.at(i_row) = HighsBasisStatus::kBasic;
    //         num_set_sbas++;
    //         num_basic++;
    //     }
    // }
    for (HighsInt i_base = 0; i_base < test_basic_index.size(); ++i_base){
        HighsInt i_col = test_basic_index.at(i_base);
        if (i_col < colCnt)
            elpBasis.col_status.at(i_col) = HighsBasisStatus::kBasic;
        else
            elpBasis.row_status.at(i_col - colCnt) = HighsBasisStatus::kBasic;
        num_basic++;
    }
    for (HighsInt i_col = 0; i_col < colCnt; ++i_col){
        HighsInt needs_pivot = col_needs_res_pivot.at(i_col);
        if (needs_pivot){
            elpBasis.col_status.at(i_col) = HighsBasisStatus::kBasic;
            num_basic++;
        }
    }
    test_basic_index.resize(0);
    for (HighsInt i_col = 0; i_col < colCnt; ++i_col){
        if (elpBasis.col_status.at(i_col) == HighsBasisStatus::kBasic)
            test_basic_index.push_back(i_col);
    }
    for (HighsInt i_row = 0; i_row < rowCnt; ++i_row){
        if (elpBasis.row_status.at(i_row) == HighsBasisStatus::kBasic){
            test_basic_index.push_back(i_row + colCnt + numResiduals);
        }
    }
}

void HighsOCAggregate::buildColBasis(){
    int iCol, pCol, pf, pc, crep, pcrep;
    HighsBasisStatus basic = HighsBasisStatus::kBasic, status;
    HighsBasisStatus nonbasic = HighsBasisStatus::kNonbasic;
    elpBasis.col_status.resize(colCnt + numResiduals);
    // elpBasis.col_status.resize(colCnt);
    std::fill_n(elpBasis.col_status.begin(), elpBasis.col_status.size(), basic);
    for (iCol = 0; iCol < colCnt; ++iCol){
        crep = colrep[iCol];
        pf = epMinusOne.front[crep];
        pCol = pFrontCol[pf];
        status = basis.col_status[pCol];
        // if (mark_degenerate.at(iCol)){
        //     elpBasis.col_status.at(iCol) = basic;
        //     splitFromNonbasicCount.at(pCol)++;
        //     continue;
        // }
        // if (mark_degenerate.at(pCol)){
        //     elpBasis.col_status.at(iCol) = basic;
        //     continue;
        // }
        // if (mark_degenerate.at(pCol) && splitFromNonbasicCount.at(pCol) < splitSize.at(pCol) - 1){
        //     elpBasis.col_status.at(iCol) = nonbasic;
        //     continue;
        // }
        // if (mark_degenerate.at(pCol)){
        //     elpBasis.col_status.at(iCol) = basic;
        //     continue;
        // }
        elpBasis.col_status[iCol] = status;
        // if (status == basic and num_basic < rowCnt){ 
        //     test_basic_index.push_back(iCol);
        //     test_basic_start.at(iCol) = 1;
        //     num_basic++;
        // }
    }
    HighsInt col_index = colCnt;
    for (iCol = colCnt; iCol < colCnt + numResiduals; ++iCol){
        // if (degenerate_basic_residuals.at(iCol)){
        //     elpBasis.col_status.at(iCol) = HighsBasisStatus::kBasic;
        //     continue;
        // }
        elpBasis.col_status[iCol] = HighsBasisStatus::kLower;
    }
    // for (iCol = 0; iCol < colCnt + numResiduals; ++iCol){
    //     if (elpBasis.col_status.at(iCol) == HighsBasisStatus::kBasic)
    //         num_basic++;
    // }
}

void HighsOCAggregate::buildTestColBasis(){
    elpBasis.col_status.resize(colCnt + numResiduals);
    std::fill_n(elpBasis.col_status.begin(), elpBasis.col_status.size(), HighsBasisStatus::kLower);
    for (HighsInt i_col = 0; i_col < colCnt; ++i_col){
        if (num_basic == rowCnt + numResiduals) return;
        HighsInt rep = colrep.at(i_col);
        HighsInt pf = epMinusOne.front.at(rep);
        HighsInt p_col = pFrontCol.at(pf);
        HighsBasisStatus status = basis.col_status.at(p_col);
        if (status == HighsBasisStatus::kBasic){
            elpBasis.col_status.at(i_col) = status;
            num_basic++;
        }
    }
    for (HighsInt i_col = colCnt; i_col < colCnt + numResiduals; ++i_col){
        if (num_basic == rowCnt + numResiduals) return;
        // HighsInt i_child = res_to_child.at(i_col);
        // if (elpBasis.col_status.at(i_child) == HighsBasisStatus::kBasic) continue;
        elpBasis.col_status.at(i_col) = HighsBasisStatus::kBasic;
        // elp.col_lower_.at(i_col) = -kHighsInf;
        // elp.col_upper_.at(i_col) = kHighsInf;
        num_basic++;
    }
}

void HighsOCAggregate::buildRowBasis(){
    int iRow, r, pr, pf, rrep, rlen;
    int of, nf;
    HighsBasisStatus status, basic = HighsBasisStatus::kBasic;
    // std::fill(elpBasis.row_status.begin(), elpBasis.row_status.end(), basic);
    elpBasis.row_status.resize(rowCnt + numResiduals);
    std::fill_n(elpBasis.row_status.begin(), elpBasis.row_status.size(), basic);
    int numNewBasic = 0;
    int numBasicToSplit = 0;
    int numBasicSplits = 0;
    int numNonBasicToSplit = 0;
    int numNonBasicSplits = 0;
    int numNonBasic = 0;
    for (iRow = 0; iRow < prowCnt; ++iRow){
        status = basis.row_status[iRow];
        elpBasis.row_status[iRow] = status;
        // if (status == basic && num_basic < rowCnt){ 
        //     test_basic_index.push_back(iRow + colCnt);
        //     test_basic_start.at(iRow + colCnt) = 1;
        // }
    }
    for (iRow = rowCnt; iRow < rowCnt + numResiduals; ++iRow){
        elpBasis.row_status[iRow] = HighsBasisStatus::kLower;
    }
    // for (iRow = 0; iRow < rowCnt + numResiduals; ++iRow){
    //     if (elpBasis.row_status.at(iRow) == HighsBasisStatus::kBasic)
    //         num_basic++;
    // }
}

void HighsOCAggregate::buildTestRowBasis(){
    elpBasis.row_status.resize(rowCnt + numResiduals);
    std::fill_n(elpBasis.row_status.begin(), elpBasis.row_status.size(), HighsBasisStatus::kLower);
    if (num_basic == rowCnt + numResiduals) return;
    for (HighsInt p_row = 0; p_row < prowCnt; ++p_row){
        if (num_basic == rowCnt + numResiduals) return;
        if (basis.row_status.at(p_row) == HighsBasisStatus::kBasic){
            elpBasis.row_status.at(p_row) = HighsBasisStatus::kBasic;
            num_basic++;
        }
    }
}

void HighsOCAggregate::buildResidualColBasis(){
    int i;
    for (i = elp.num_col_ - elp.num_residual_cols_; i < elp.num_col_; ++i)
        elpBasis.col_status[i] = HighsBasisStatus::kLower;
}

void HighsOCAggregate::buildResidualRowBasis(){
    int i;
    for (i = elp.num_row_ - elp.num_residual_cols_; i < elp.num_row_; ++i)
        elpBasis.row_status[i] = HighsBasisStatus::kLower;
}

void HighsOCAggregate::findFrontMins(){
    HighsInt v_min, i_front, i_part;
    for (i_front = 0; i_front < numTot; i_front += ep.len.at(i_front) + 1){
        v_min = kHighsIInf;
        for (i_part = i_front; i_part < i_front + ep.len.at(i_front) + 1; ++i_part){
            if (ep.label.at(i_part) < v_min) v_min = ep.label.at(i_part);
        }
        frontMin.at(i_front) = v_min;
    }
}

void HighsOCAggregate::buildColPointers(){
    int i, j, min_rep;
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
            fronts.insert(ep.front[i]);
            frontCol[ep.front[i]] = col[i];
            colFront[col[i]] = ep.front[i];
            continue;
        }
        newFront = fronts.insert(ep.front[i]).second;
        if (newFront){
            frontCol[ep.front[i]] = colCnt;
            colFront[colCnt] = ep.front[i];
            min_rep = frontMin.at(ep.front.at(i));
            colrep[colCnt] = min_rep;
            col[i] = colCnt++;
        }
    }
}

void HighsOCAggregate::buildRowPointers(){
    int i, j, min_rep;
    std::set<int> fronts;
    std::set<int> testRows;
    std::vector<int> rowFreq(numRow, 0);
    bool newFront;
    std::pair<std::set<int>::iterator,bool> newRep;
    newRowReps.clear();
    for (i = numCol; i < numTot; ++i){
        if (row[i - numCol] > -1){
            fronts.insert(ep.front[i]);
            frontRow[ep.front[i]] = row[i - numCol];
            rowFront[row[i - numCol]] = ep.front[i];
            continue;
        }
        newFront = fronts.insert(ep.front[i]).second;
        if (newFront){
            frontRow[ep.front[i]] = rowCnt;
            rowFront[rowCnt] = ep.front[i];
            min_rep = frontMin.at(ep.front.at(i));
            rowrep[rowCnt] = min_rep;
            row[i - numCol] = rowCnt++;
        }
    }
}

// void HighsOCAggregate::buildDegenerateAMatrix(){
//     int cnt = 0;
//     std::vector<HighsInt> deg_rows(numRow, 0);
//     for (int i_col = 0; i_col < colCnt; ++i_col){
//         if (mark_degenerate.at(i_col)){
//             for (int i_row = 0; i_row < basic_index_.size(); ++i_row){
//                 if (basic_index_.at(i_row) == i_col){
//                     cnt++;
//                     deg_rows.at(i_row) = 1;
//                 }
//             }
//         }
//     }
//     HighsInt i_vec, i_col, i_mat, i_row, crep, cf, clen, c,
//         ro, rf, r, rlen, nnz = 0, start = 0;
//     degenerate_col_map_.clear();
//     degenerate_matrix.clear();
//     // degenerate_matrix.start_.push_back(0);
//     std::set<HighsInt> degenerate_rows;
//     for (i_vec = 0; i_vec < basic_index_.size(); ++i_vec){
//         i_col = basic_index_.at(i_vec);
//         if (i_col >= numCol) continue;
//         if (!mark_degenerate.at(i_col)) continue;
//         degenerate_col_map_.push_back(i_col);
//         crep = colrep[i_col];
//         cf = colFront[i_col];
//         clen = ep.len[cf] + 1;
//         c = i_col;
//         int num = 0;
//         for (i_mat = olp.a_matrix_.start_[crep]; i_mat < olp.a_matrix_.start_[crep + 1]; ++i_mat){
//             ro = olp.a_matrix_.index_[i_mat];
//             rf = ep.front[ro + numCol];
//             r = frontRow[rf];
//             if (!deg_rows.at(r)) continue;
//             // if (basic_index_.at(r) != i_col) continue;
//             num++;
//             std::cout << "i_col: " << i_col << std::endl;
//             std::cout << "i_row_basic: " << r << std::endl;
//             std::cout << "num: " << num << std::endl;
//             rlen = ep.len[rf] + 1;
//             if (!columnF[r]++){ 
//                 if (nnz >= columnI.size())
//                     std::cin.get();
//                 columnI[nnz++] = r;
//                 // nnzStan++;
//             }
//             columnX[r] += olp.a_matrix_.value_[i_mat] * clen;
//         }
//         if (num == 0){
//             std::cin.get();
//         }
//         for (i_mat = start; i_mat < nnz; ++i_mat){
//             degenerate_matrix.index_.push_back(columnI[i_mat]);
//             degenerate_rows.insert(columnI.at(i_mat));
//             // AdegenRIndex[j] = columnI[j];
//             degenerate_matrix.value_.push_back(columnX[columnI[i_mat]]);
//             // AdegenRValue[j] = columnX[columnI[j]];
//             columnF[columnI[i_mat]] = 0;
//             columnX[columnI[i_mat]] = 0;
//             columnI[i_mat] = 0;
//         }
//         degenerate_matrix.start_.push_back(nnz);
//         // AdegenRStart[xi + 1] = nnzStan;
//         start = nnz;
//     }
//     degenerate_matrix.num_col_ = degenerate_col_map_.size();
//     degenerate_matrix.num_row_ = degenerate_rows.size();
// }

// void HighsOCAggregate::buildDegenerateLU(){
//     if (!degenerate_col_map_.size()) return;
//     std::vector<HighsInt> basic_degenerates(degenerate_col_map_.size());
//     std::iota(basic_degenerates.begin(), basic_degenerates.end(), 0);
//     degenerate_factor.setup(degenerate_matrix, basic_degenerates);
//     HighsInt rank_deficiency = degenerate_factor.build();
//     std::cout << "build LU" << std::endl;
// }

void HighsOCAggregate::trackAndCountSplits(){
    std::fill_n(splitFrom.begin(), numCol, -1);
    std::fill_n(splitSize.begin(), numCol, 0);
    std::fill_n(splitFromNonbasicCount.begin(), numCol, 0);
    for (int i_col = 0; i_col < colCnt; ++i_col){
        int crep = colrep[i_col];
        int pf = epMinusOne.front[crep];
        int p_col = pFrontCol[pf];
        if (p_col != i_col)
            splitFrom.at(i_col) = p_col;
        splitSize.at(p_col)++;
    }
}

void HighsOCAggregate::copyPartition(){
    epMinusOne.target = ep.target;
    epMinusOne.level = ep.level;
    epMinusOne.nsplits = ep.nsplits;
    epMinusOne.ncsplits = ep.ncsplits;
    epMinusOne.nrsplits = ep.nrsplits;
    // std::copy(ep.front.begin(), ep.front.end(), epMinusOne.front.begin());
    // std::copy(ep.label.begin(), ep.label.end(), epMinusOne.label.begin());
    // std::copy(ep.unlabel.begin(), ep.unlabel.end(), epMinusOne.unlabel.begin());
    // std::copy(ep.parent.begin(), ep.parent.end(), epMinusOne.parent.begin());
    // std::copy(ep.len.begin(), ep.len.end(), epMinusOne.len.begin());
    epMinusOne.front = ep.front;
    epMinusOne.label = ep.label;
    epMinusOne.unlabel = ep.unlabel;
    epMinusOne.parent = ep.parent;
    epMinusOne.len = ep.len;
}

HighsLp HighsOCAggregate::getLp(){
    return elp;
}

HighsLp HighsOCAggregate::getAggLp(){
    return agglp;
}

HighsLp HighsOCAggregate::getLpNoResiduals(){
    return presolvelp;
}

HighsBasis HighsOCAggregate::getBasis(){
    // elpBasis.alien = false;
    // elpBasis.was_alien = false;
    return elpBasis;
}

HEkk& HighsOCAggregate::getEkkInstance(){
    return test_ekk;
}

HighsSolution& HighsOCAggregate::getElpSolution(){
    return lift_solution;
}

HighsOptions& HighsOCAggregate::getElpOptions(){
    return test_options;
}

HighsInfo& HighsOCAggregate::getElpInfo(){
    return test_info;
}
