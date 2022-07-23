#ifndef PRESOLVE_OCAGGREGATOR_H
#define PRESOLVE_OCAGGREGATOR_H

// Includes
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <set>
#include <numeric>
#include "HighsLp.h"
#include "OCEquitable.h"
#include "HFactor.h"

/* Aggregator class for LPs based off a equitable partition */
class HighsOCAggregate{
public:
    void passLpAndPartition(HighsLp& lp, OCPartition& partition);
    void passNewPartitionBasisSolution(OCPartition& partition, HighsBasis& b,
                                        HighsSolution& s);
    void resizeAlpContainers();
    void resizeElpContainers();
    void resizeGramSchmidtMatrixContainers();
    void buildLp();
    void buildLpExtended();
    // Build solution for aggregate ipm solve to hand off to highs basic 
    // crossover
    HighsSolution buildSolution(OCPartition& partition, HighsSolution& s);
    // Build aggregate A matrices for ALP and EALP
    void buildAmatrix();
    void buildAmatrixExtended();
    void buildAmatrixExtendedNoResiduals();
    // Build aggregate b vectors for ALP and EALP 
    void buildRhs();
    void buildRhsFromScratch();
    void buildRhsFromSolution();
    void buildRhsExtended();
    void buildRhsFromScratchExtended();
    void buildRhsFromSolutionExtended();
    void buildRhsExtendedNoResiduals();
    void buildRhsFromSolutionExtendedNoResiduals();
    void buildRowNames();
    // Build aggregate LB/UB vectors for ALP and EALP
    void buildBnds();
    void buildBndsFromScratch();
    void buildBndsFromSolution();
    void buildBndsExtended();
    void buildBndsFromScratchExtended();
    void buildBndsFromSolutionExtended();
    void buildBndsExtendedNoResiduals();
    void buildBndsFromSolutionExtendedNoResiduals();
    void buildColNames();
    // Build aggregate c vectors for ALP and EALP
    void buildObj();
    void buildObjExtended();
    void buildObjExtendedNoResiduals();
    // Build and handle residual cols/rows for ALP and EALP
    int buildResiduals();
    void buildResidualLinks();
    void buildResidualCols();
    void buildResidualRows();
    // Build vectors containing the representative of an aggregate column/row
    void findFrontMins();
    void buildRowPointers();
    void buildRowFrontsAndReps();
    void buildColPointers();
    void buildColFrontsAndReps();
    HighsInt trackAndCountSplits();
    // Build Highs basis vectors for EALP
    void markDegenerate();
    void buildBasis(bool finish, bool extended);
    void buildColBasis();
    void buildRowBasis();
    void buildResidualColBasis();
    void buildResidualRowBasis();
    // Build LU initial factor check for good basis;
    void buildDegenerateLU();
    void buildDegenerateAMatrix();
    // Sanity checking for A matrices
    void checkAMatrix();
    void printAMatrixToMatlabFormat(HighsSparseMatrix& matrix);
    // Copy partition from level k when going to level k + 1
    void copyPartition();
    // Perform dependent residuals test using gram-schmidt process to
    // remove linking variables and rows that are not needed
    void gramSchmidt();
    void buildGramSchmidtExtendedMatrix();
    void updateGramSchmidtMatrix(HighsSparseMatrix& matrix, HighsInt i_col, HVector& column_v);
    void markLinkerDeleted(HighsInt i_link);
    void markRowIndependent(HighsInt i_row);
    void clearDeleteLinker();
    void clearIndependentRow();
    void updateHVectorIndex(HVector& column_v, HighsInt insert, HighsInt idx);
    void divideSparseVectorByScalar(HVector& column_q, double scalar, HVector& column_temp);
    void subtractSparseVector(HVector& column_v, HVector& column_q);
    double sparseDotProduct(HVector& column_v, HVector& column_q);
    // Get parts of agg
    HighsLp getLp();
    HighsLp getAggLp();
    HighsLp getLpNoResiduals();
    HighsBasis getBasis();
    HighsTimer timer_;

    // dev test functions
    void checkForBadNonBasics(HighsInt col);
    

    HighsLp elp; 
    HighsLp olp;
    HighsLp agglp;
    HighsLp presolvelp;
    OCPartition ep;
    OCPartition epMinusOne;
    HighsBasis basis;
    HighsBasis elpBasis;
    HighsSolution solution;
    HighsSolution lift_solution;
    HFactor degenerate_factor;
    HighsSparseMatrix degenerate_matrix;
    int level = 0;
    int numTot;
    int numCol; 
    int numRow; 
    int nnz;
    int numResiduals = 0;
    int prev_num_residuals = 0;
    int numTotResiduals;
    int degenCnt = 0;
    int colCnt = 0;
    int pcolCnt = 0;
    int rowCnt = 0;
    int prowCnt = 0;
    int num_basic;
    HighsInt num_deleted_links;
    std::vector<int> agg_col;
    std::vector<int> col;
    std::vector<int> colrep;
    std::vector<int> pcol;
    std::vector<int> pcolrep;
    std::vector<int> row;
    std::vector<int> rowrep;
    std::vector<int> prow;
    std::vector<int> prowrep;
    std::vector<int> splitFrom;
    std::vector<int> splitFromNonbasicCount;
    std::vector<int> splitSize;
    std::vector<double> columnX;
    std::vector<int> columnI;
    std::vector<int> columnF;
    std::vector<int> rowF;
    std::vector<int> rowNonbasic;
    std::vector<int> parent;
    std::vector<int> isParent;
    std::vector<int> parentStart;
    std::vector<int> parentFreq;
    std::vector<int> parentRow;
    std::vector<int> child;
    std::vector<int> isChild;
    std::vector<int> childRow;
    std::vector<int> residualCol;
    std::vector<int> residualRow;
    std::vector<int> degenSlack;
    std::vector<bool> degenRow;
    std::vector<HighsInt> mark_degenerate;
    std::vector<bool> finalRowRep;
    // For computing B^-1 * A_r submat
    std::vector<int> basicIndex;
    std::vector<int> nonbasicFlag;
    std::vector<int> indepRows;
    std::vector<int> rowPerm;
    std::vector<int> rowUnperm;
    std::vector<int> colPerm;
    std::vector<int> colUnperm;
    std::vector<int> frontCol;
    std::vector<int> colFront;
    std::vector<int> pFrontCol;
    std::vector<int> pColFront;
    std::vector<int> frontRow;
    std::vector<int> rowFront;
    std::vector<int> pFrontRow;
    std::vector<int> frontMin;
    std::vector<int> repFront;
    std::set<int> colReps;
    std::set<int> newColReps;
    std::set<int> rowReps;
    std::set<int> newRowReps;
    std::map<int, std::vector<int> > splitCells;
    // std::vector<std::pair<int, int> > splitCells;
    std::vector<int> zero_step_pivots;
};

#endif