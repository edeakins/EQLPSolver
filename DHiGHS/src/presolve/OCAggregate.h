#ifndef PRESOLVE_OCAGGREGATOR_H
#define PRESOLVE_OCAGGREGATOR_H

// Includes
#include <string>
#include "HighsLp.h"
#include "OCEquitable.h"

/* Aggregator class for LPs based off a equitable partition */
class HighsOCAggregate{
public:
    void allocate(HighsLp* lp, OCpartition* partition);
    void buildLp();
    void buildLp(OCpartition* partition, HighsBasis* b,
                HighsSolution* s, bool finish, bool extended);
    void buildAmatrix();
    void buildFinalAmatrix();
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
    void buildBasis(bool finish);
    void buildColBasis();
    void buildRowBasis();
    void buildResidualColBasis();
    void buildResidualRowBasis();
    void buildFinalColSolution();
    void buildFinalRowSolution();
    void buildFinalColBasis();
    void buildFinalRowBasis();
    // Get parts of agg
    HighsLp* getLp();
    HighsSolution* getSolution();
    HighsBasis* getBasis();

    HighsLp *elp; 
    HighsLp* olp;
    OCpartition* ep;
    HighsBasis* basis;
    HighsBasis* elpBasis;
    HighsSolution* solution;
    HighsSolution* elpSolution;
    bool buildFinalLp = false;
    int numTot;
    int numCol; 
    int numRow; 
    int nnz;
    int numResiduals;
    int numTotResiduals;
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
    std::vector<int> parentStart;
    std::vector<int> child;
};

#endif