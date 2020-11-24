#ifndef AGGREGATE_H_
#define AGGREGATE_H_

#include "EquitablePartition.hpp"

class AggregateLp{
public:
    AggregateLp(EquitablePartition& ep);
    void updateEP(EquitablePartition& ep);
    void clear();
    void findDimensions();
    void aggregateColBounds();
    void aggregateRowBounds();
    void aggregateAMatrix();
    void aggregateCostVector();
    void aggregate();
    vector<double>& getColUpper();
    vector<double>& getColLower();
    vector<double>& getColCost();
    vector<double>& getRowUpper();
    vector<double>& getRowLower();
    vector<double>& getAvalue();
    vector<int>& getAindex();
    vector<int>& getAstart();

    // Original LP info
    int nRows;
    int nCols; 
    int numTot;
    vector<double> colCost;
    vector<double> colLower;
    vector<double> colUpper; 
    vector<double> rowLower;
    vector<double> rowUpper;
    vector<double> Avalue; 
    vector<int> Aindex; 
    vector<int> Astart;
    vector<int> AindexP;

    // Reduced LP 
    int nRows_ = 0;
    int tempNRows_ = 0;
    int nCols_ = 0; 
    int tempNCols_ = 0;
    int numTot_ = 0;
    vector<double> colCost_;
    vector<double> colLower_;
    vector<double> colUpper_; 
    vector<double> rowLower_;
    vector<double> rowUpper_;
    vector<double> Avalue_; 
    vector<int> Aindex_; 
    vector<int> Astart_;

    // EP info
    vector<vector<int> > C;
    vector<int> Csize;
    vector<int> color;

};

#endif 