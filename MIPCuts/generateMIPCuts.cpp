#include "generateMIPCuts.hpp"
using namespace std;

void update(int& nCols, int& nRows, int nnz, vector<int>& Astart,
            vector<int>& Aindex, vector<double>& Avalue, vector<int>& ARstart,
            vector<int>& ARindex, vector<double>& ARvalue, vector<double>& colLower,
            vector<double>& colUpper, vector<double>& rowLower, vector<double>& rowUpper,
            EquitablePartition& ep, AggregateLp& aLp, OsiClpSolverInterface& aggSi, OsiClpSolverInterface& si){
  // Grab info we will need for update
  int nAggRowsWCuts = aggSi.getNumRows();
  int nAggColsWCuts = aggSi.getNumCols();
  int nnzAggWCuts = aggSi.getMatrixByCol()->getVectorStarts()[nAggColsWCuts];
  int nAggRows = aLp.nRows_;
  int nAggCols = aLp.nCols_;
  vector<double> rowLowerAgg(aggSi.getRowLower(), aggSi.getRowLower() + nAggRowsWCuts);
  vector<double> rowUpperAgg(aggSi.getRowUpper(), aggSi.getRowUpper() + nAggRowsWCuts);
  vector<double> ARvalueAgg(aggSi.getMatrixByRow()->getElements(), 
                                aggSi.getMatrixByRow()->getElements() + nnzAggWCuts);
  vector<int> ARindexAgg(aggSi.getMatrixByRow()->getIndices(), 
                                aggSi.getMatrixByRow()->getIndices() + nnzAggWCuts);
  vector<int> ARstartAgg(aggSi.getMatrixByRow()->getVectorStarts(), 
                                aggSi.getMatrixByRow()->getVectorStarts() + nAggRowsWCuts + 1);
  vector<vector<int> > C = ep.C; 
  // Update info for original LP
  int nRowsAdded = 0;
  for (int i = nAggRows; i < nAggRowsWCuts; ++i){
    double lb = rowLowerAgg[i];
    double ub = rowUpperAgg[i];
    vector<int> indices;
    vector<double> values;
    // Do row wise update (This is easy because cuts are in row format)
    for (int j = ARstartAgg[i]; j < ARstartAgg[i + 1]; ++j){
      int color = ARindexAgg[j];
      double alpha = ARvalueAgg[j];
      indices.insert(indices.end(), C[color].begin(), C[color].end());
      vector<double> temp(indices.size(), alpha);
      values.insert(values.end(), temp.begin(), temp.end());
    }
    int numEl = indices.size();
    si.addRow(numEl, &indices[0], &values[0], lb, ub);
  } 
  // Columns will be automagically updated using clp model class in background
  // Pull new LP info with added cuts in the highest space
  int nRowsWCuts = si.getNumRows();
  int nColsWCuts = si.getNumCols();
  int nnzWCuts = si.getMatrixByCol()->getVectorStarts()[nColsWCuts];
  rowLower.insert(rowLower.end(), si.getRowLower() + nRows, si.getRowLower() + nRowsWCuts);
  rowUpper.insert(rowUpper.end(), si.getRowUpper() + nRows, si.getRowUpper() + nRowsWCuts);
  //Avalue.clear();
  Avalue.assign(si.getMatrixByCol()->getElements(), 
                si.getMatrixByCol()->getElements() + nnzWCuts);
  //Aindex.clear();
  Aindex.assign(si.getMatrixByCol()->getIndices(), 
                si.getMatrixByCol()->getIndices() + nnzWCuts);
  //Astart.clear();
  Astart.assign(si.getMatrixByCol()->getVectorStarts(),
                si.getMatrixByCol()->getVectorStarts() + nColsWCuts + 1);
  ARvalue.insert(ARvalue.end(), si.getMatrixByRow()->getElements() + nnz, 
                 si.getMatrixByRow()->getElements() + nnzWCuts);
  ARindex.insert(ARindex.end(), si.getMatrixByRow()->getIndices() + nnz, 
                 si.getMatrixByRow()->getIndices() + nnzWCuts);
  ARstart.insert(ARstart.end(), si.getMatrixByRow()->getVectorStarts() + nRows + 1, 
                 si.getMatrixByRow()->getVectorStarts() + nRowsWCuts + 1);
  nRows = nRowsWCuts;
  nCols = nColsWCuts;  
  nnz = nnzWCuts;   
}

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
    int nRows = si.getNumRows();
    int _nRows = si.getNumRows();
    int nCols = si.getNumCols();
    int _nCols = si.getNumCols();
    int nnz = si.getMatrixByCol()->getVectorStarts()[nCols];
    int _nnz = si.getMatrixByCol()->getVectorStarts()[nCols];
    vector<double> colCost(si.getObjCoefficients(), 
                                  si.getObjCoefficients() + nCols);
    vector<double> colLower(si.getColLower(), si.getColLower() + nCols);
    vector<double> colUpper(si.getColUpper(), si.getColUpper() + nCols);
    vector<double> rowLower(si.getRowLower(), si.getRowLower() + nRows);
    vector<double> rowUpper(si.getRowUpper(), si.getRowUpper() + nRows);
    vector<double> Avalue(si.getMatrixByCol()->getElements(), 
                                  si.getMatrixByCol()->getElements() + nnz);
    vector<int> Aindex(si.getMatrixByCol()->getIndices(), 
                                  si.getMatrixByCol()->getIndices() + nnz);
    vector<int> Astart(si.getMatrixByCol()->getVectorStarts(), 
                                  si.getMatrixByCol()->getVectorStarts() + nCols + 1);
    vector<double> ARvalue(si.getMatrixByRow()->getElements(), 
                                  si.getMatrixByRow()->getElements() + nnz);
    vector<int> ARindex(si.getMatrixByRow()->getIndices(), 
                                  si.getMatrixByRow()->getIndices() + nnz);
    vector<int> ARstart(si.getMatrixByRow()->getVectorStarts(), 
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
    for (int i = 0; i < nCols; ++i){
      aggSi.setInteger(i);
    }
    // Solve initial aggregate
    aggSi.initialSolve();
    // Generate l & p cuts for this model
    liftProjectCuts.generateCuts(aggSi, cuts);
    acRc = aggSi.applyCuts(cuts,0.0);
    update(nCols, nRows, nnz, Astart, Aindex,
            Avalue, ARstart, ARindex, ARvalue,
            colLower, colUpper, rowLower, rowUpper,
            ep, aLp, aggSi, si);
    ep.refine();
    aLp.updateMasterLpAndEp(ep, _nCols, _nRows, 
                        _nnz, Astart, Aindex,
                        Avalue, rowLower, rowUpper);
    aLp.aggregate();
    cout << "finish" << endl;
    // Establish number of cuts added and new matrix data
    // cout <<endl <<endl;
    //   cout <<cuts.sizeCuts() <<" cuts were generated" <<endl;
    //   cout <<"  " <<acRc.getNumInconsistent() <<" were inconsistent" <<endl;
    //   cout <<"  " <<acRc.getNumInconsistentWrtIntegerModel() 
    //        <<" were inconsistent for this problem" <<endl;
    //   cout <<"  " <<acRc.getNumInfeasible() <<" were infeasible" <<endl;
    //   cout <<"  " <<acRc.getNumIneffective() <<" were ineffective" <<endl;
    //   cout <<"  " <<acRc.getNumApplied() <<" were applied" <<endl;
    // cout <<endl <<endl;
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
