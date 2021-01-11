#ifndef HIGHS_EQUITABLE_H
#define HIGHS_EQUITABLE_H

#include "HighsLp.h"
#include "HighsSimpleDec.h"
#include "HighsTimer.h"

#include <string>
#include <vector>
#include <algorithm>
#include <stack>
#include <set>
#include <list>
#include <map>
#include <tuple>
#include <numeric>
#include <functional>
#include <forward_list>
using namespace std;

class HighsEquitable {
public:
	// Setup for equitable ptn
	HighsEquitable(const HighsLp& lp);
	void transpose();
    void initRefinement();
    void handleNegatives();
	void refine();
    void splitColor(int s);
    void isolate(int i);
    void findTarget();
    bool isDiscrete();
    void packVectors();

	/// Original LP info
    int nRows;
    int nCols; 
    int nTot;
    vector<double> colCost;
    vector<double> colLower;
    vector<double> colUpper; 
    vector<double> rowLower;
    vector<double> rowUpper;
    vector<double> Avalue; 
    vector<int> Aindex; 
    vector<int> Astart;
    vector<double> AvaluePos;
    vector<double> ARvalue;
    vector<int> ARindex;
    vector<int> ARstart;
    vector<double> ARvaluePos;
    vector<int> AindexP;

    // EP info
    int vCol, cCol, refinements = 0;
    vector<vector<int> > C;
    vector<vector<int> > A;
    vector<int> color;
    vector<bool> SCheck;
    vector<double> cdeg;
	vector<double> mincdeg;
	vector<double> maxcdeg;
    vector<bool> isAdj;
    vector<int> colorsAdj;
    vector<int> colorsToSplit;
    vector<bool> isolates;
    stack<int> S;
    vector<int> Csize;
    vector<int> Asize;
    vector<int> initialParts;
    vector<int> numEdges;
    vector<int> parentPartition;

};

#endif