#include "generateMIPCuts.hpp"
using namespace std;

int main(int argc, const char *argv[]){
    // File name for problem being read in 
    string mpsFileName = argv[1];
    // Intialize solvers and cut generators
    OsiClpSolverInterface si;
    OsiClpSolverInterface aggSi;
    CglSimpleRounding liftProjectCuts;
    OsiCuts cuts;
    OsiSolverInterface::ApplyCutsReturnCode acRc;

    // Grab info for equitable partition scheme;
    si.readMps(mpsFileName.c_str(), "mps");
    si.initialSolve();
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
                                    si.getMatrixByRow()->getVectorStarts() + nRows + 1);

    // Pass LP info off to EP algorithm to compute an
    // EP of the LP
    EquitablePartition ep(nRows, nCols, colCost, colLower,
                            colUpper, rowLower, rowUpper,
                            Avalue, Aindex, Astart,
                            ARvalue, ARindex, ARstart);
    // Pass EP to aggregator for first aggregation
    AggregateLp aLp(ep);
    aLp.aggregate();
    // Grab aggregate Lp info
    int nCols_ = aLp.nCols_;
    int nRows_ = aLp.nRows_;
    double* colCost_ = &aLp.colCost_[0];
    double* colUpper_ = &aLp.colUpper_[0];
    double* colLower_ = &aLp.colLower_[0];
    double* rowUpper_ = &aLp.rowUpper_[0];
    double* rowLower_ = &aLp.rowLower_[0];
    double* Avalue_ = &aLp.Avalue_[0];
    int* Aindex_ = &aLp.Aindex_[0];
    int* Astart_ = &aLp.Astart_[0];
    // Load problem into a native Coin clpmodel type
    aggSi.loadProblem(nCols_, nRows_, Astart_, Aindex_,
                        Avalue_, colLower_, colUpper_,
                        colCost_, rowLower_, rowUpper_);
    // Solve initial aggregate
    aggSi.initialSolve();
    // Generate l & p cuts for this model
    liftProjectCuts.generateCuts(si, cuts);
    acRc = si.applyCuts(cuts,0.0);
    cout <<endl <<endl;
      cout <<cuts.sizeCuts() <<" cuts were generated" <<endl;
      cout <<"  " <<acRc.getNumInconsistent() <<" were inconsistent" <<endl;
      cout <<"  " <<acRc.getNumInconsistentWrtIntegerModel() 
           <<" were inconsistent for this problem" <<endl;
      cout <<"  " <<acRc.getNumInfeasible() <<" were infeasible" <<endl;
      cout <<"  " <<acRc.getNumIneffective() <<" were ineffective" <<endl;
      cout <<"  " <<acRc.getNumApplied() <<" were applied" <<endl;
      cout <<endl <<endl;
    aggSi.resolve();
    // Load problem into a native Coin clpmodel type
    /* Build first Coin aggregate model and pass off to another 
    solver interface (may not need to solver interfaces, but for
    safe keeping early on shall have two) */
    // while (!ep.isDiscrete()){
    //     ep.refine();
    //     aLp.updateEP(ep);
    //     aLp.aggregate();
    // }
    
    // si.initialSolve();
    // CglLiftAndProject cg1;
}