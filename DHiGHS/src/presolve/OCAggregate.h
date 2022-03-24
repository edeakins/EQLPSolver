#ifndef PRESOLVE_OCAGGREGATOR_H
#define PRESOLVE_OCAGGREGATOR_H

// Includes
#include <string>
#include <set>
#include "HighsLp.h"
#include "OCEquitable.h"
#include "sparseMat.h"
#include "HFactor.h"
#include "HMatrix.h"
#include "HVector.h"

/* Aggregator class for LPs based off a equitable partition */
class HighsOCAggregate{
public:
    void allocate(HighsLp* lp, OCpartition* partition);
    void buildLp();
    void buildLp(OCpartition* partition, HighsBasis* b,
                HighsSolution* s, bool finish, bool extended);
    void buildAmatrix();
    void buildFinalAmatrix();
    void buildStandardMatrix();
    void buildRhs();
    void buildRhsFromScratch();
    void buildRhsFromSolution();
    void buildFinalRhs();
    void buildFinalRhsFromScratch();
    void buildFinalRhsFromSolution();
    void buildBnds();
    void buildBndsFromScratch();
    void buildBndsFromSolution();
    void buildFinalBnds();
    void buildFinalBndsFromScratch();
    void buildFinalBndsFromSolution();
    void buildObj();
    void buildFinalObj();
    void buildResiduals();
    void buildFinalResiduals();
    void buildResidualLinks();
    void buildFinalResidualLinks();
    void buildResidualSubMatrix();
    void buildResidualCols();
    void buildResidualRows();
    void buildRowPointers();
    void buildColPointers();
    void buildFinalPointers();
    void buildSolution(bool finish, bool extended);
    void buildColSolution();
    void buildRowSolution();
    void buildResidualColSolution();
    void buildResidualRowSolution();
    void buildBasis(bool finish, bool extended);
    void buildColBasis();
    void buildRowBasis();
    void buildResidualColBasis();
    void buildResidualRowBasis();
    void buildFinalColSolution();
    void buildFinalRowSolution();
    void buildFinalColBasis();
    void buildFinalRowBasis();
    // New sparse mat building
    void buildSparseMat();
    void buildLUFactor();
    // New method of trying to pivot out degenerate slacks
    void setUpPreBasicIndex();
    void setUpPreNonbasicFlag();
    void setUpColAq();
    void setUpPreLU();
    void setUpPreMatrix();
    void factor();
    void reduceColumn(int iCol);
    void buildSubLp();
    void buildSubLpBasis();
    void checkAMatrix();
    void printAMatrixToMatlabFormat();
    // Get parts of agg
    HighsLp* getLp();
    HighsSolution* getSolution();
    HighsBasis* getBasis();
    HighsLp* getSubLp();
    HighsBasis* getSubBasis();
    std::vector<int>& getUnPerm();

    HighsLp *elp; 
    HighsLp* olp;
    HighsLp* sublp;
    OCpartition* ep;
    OCpartition epMinusOne;
    HighsBasis* basis;
    HighsBasis* elpBasis;
    HighsBasis* sublpBasis;
    HighsSolution* solution;
    HighsSolution* elpSolution;
    bool buildFinalLp = false;
    int numTot;
    int numCol; 
    int numRow; 
    int nnz;
    int numResiduals = 0;
    int numTotResiduals;
    int degenCnt = 0;
    int colCnt = 0;
    int pcolCnt = 0;
    int rowCnt = 0;
    int prowCnt = 0;
    std::vector<int> AdegenRStart;
    std::vector<int> AdegenRIndex;
    std::vector<double> AdegenRValue;
    std::vector<int> col;
    std::vector<int> colrep;
    std::vector<int> pcol;
    std::vector<int> pcolrep;
    std::vector<int> row;
    std::vector<int> rowrep;
    std::vector<int> prow;
    std::vector<int> prowrep;
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
    std::vector<int> frontRow;
    std::vector<int> rowFront;
    std::vector<int> pFrontRow;
    std::set<int> colReps;
    std::set<int> newColReps;
    std::set<int> rowReps;
    std::set<int> newRowReps;
    std::map<int, std::vector<int> > splitCells;

    // New sparse mat class for initial LU to remove some r columns hopefull
    SparseMatrix A;
    SparseMatrix AT;
    HFactor LU;
    HMatrix matrix;
    HVector colAq;
};

#endif