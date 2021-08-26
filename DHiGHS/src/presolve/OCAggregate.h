#ifndef PRESOLVE_OCAGGREGATOR_H
#define PRESOLVE_OCAGGREGATOR_H

// Includes
#include <string>
#include "HighsLp.h"
#include "OCEquitable.h"

/* Aggregator class for LPs based off a equitable partition */
class HighsOCAggregate{
public:
    HighsLp* allocate(HighsLp* lp, OCpartition* partition, HighsBasis* b, HighsSolution* s);
    HighsLp* buildLp();
    HighsLp* buildLp(OCpartition* partition, HighsBasis* b, HighsSolution* s);
    void buildAmatrix();
    void buildRhs();
    void buildRhsFromScratch();
    void buildRhsFromSolution();
    void buildBnds();
    void buildBndsFromScratch();
    void buildBndsFromSolution();
    void buildObj();
    void buildResidualCols();
    void buildResidualRows();
    void buildRowPointers();
    void buildColPointers();
    HighsBasis* buildBasis();
    void buildColBasis();
    void buildRowBasis();
    void buildResidualsBasis();
    void buildResidualColBasis();
    void buildResidualRowBasis();

    HighsLp *elp, *olp;
    OCpartition* ep;
    HighsBasis* basis;
    HighsSolution* solution;
    int numTot;
    int numCol; 
    int numRow; 
    int nnz;
    int numResiduals;
    std::vector<int> col;
    std::vector<int> pcol;
    std::vector<int> row;
    std::vector<int> prow;
    std::vector<double> columnX;
    std::vector<int> columnI;
    std::vector<int> columnF;
};

#endif