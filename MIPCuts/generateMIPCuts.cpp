#include "generateMIPCuts.hpp"
using namespace std;

int main(int argc, const char *argv[]){
    // File name for problem being read in 
    string mpsFileName = argv[1];

    // Grab info for equitable partition scheme
    OsiClpSolverInterface si;
    si.readMps(mpsFileName.c_str(), "mps");
    const int nRows = si.getNumRows();
    const int nCols = si.getNumCols();
    const int nnz = si.getMatrixByCol()->getVectorStarts()[nCols];
    const vector<double> colCost(si.getObjCoefficients(), 
                                    si.getObjCoefficients() + nCols);
    const vector<double> colLower(si.getColLower(), si.getColLower() + nCols);
    const vector<double> colUpper(si.getColUpper(), si.getColUpper() + nCols);
    const vector<double> rowLower(si.getRowLower(), si.getRowLower() + nRows);
    const vector<double> rowUpper(si.getRowUpper(), si.getRowUpper() + nRows);
    const vector<double> Avalue(si.getMatrixByCol()->getElements(), 
                                    si.getMatrixByCol()->getElements() + nnz);
    const vector<int> Aindex(si.getMatrixByCol()->getIndices(), 
                                    si.getMatrixByCol()->getIndices() + nnz);
    const vector<int> Astart(si.getMatrixByCol()->getVectorStarts(), 
                                    si.getMatrixByCol()->getVectorStarts() + nCols + 1);
    const vector<double> ARvalue(si.getMatrixByRow()->getElements(), 
                                    si.getMatrixByRow()->getElements() + nnz);
    const vector<int> ARindex(si.getMatrixByRow()->getIndices(), 
                                    si.getMatrixByRow()->getIndices() + nnz);
    const vector<int> ARstart(si.getMatrixByRow()->getVectorStarts(), 
                                    si.getMatrixByRow()->getVectorStarts() + nCols + 1);

    // Pass LP info off to EP algorithm to compute an
    // EP of the LP
    EquitablePartition ep(nRows, nCols, colCost, colLower,
                            colUpper, rowLower, rowUpper,
                            Avalue, Aindex, Astart,
                            ARvalue, ARindex, ARstart);
    // Pass EP to aggregator
    AggregateLp aLp(ep);
    
    si.initialSolve();
    CglLiftAndProject cg1;
}