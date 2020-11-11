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
	void clear();
	void setup(const HighsLp& lp);
	void handleNegatives();
	void createRowCopy();
	void initialRefinement();
	void splitColor(int color);
	void split(int s);
	void refine();
	void findTarget();
	void isolate(int i);
	void collectLinkingPairs();
	bool isPartitionDiscrete();
	void packVectors();
	void shiftVectors(int& start, int& finish, int& shift);

	// Some storage for info about lp and ptn
	string model_name = "";
	string lp_name = "";
	// (Scalars)
	int numRow;
	int numCol;
	int numTot;
	int sense;
	int nnz;
	double offset;
	int vCol;
	int cCol;
	int numParts;
	int r;
	int rep;
	int isolated = -1;
	int refinements = 0;
	
	// timers
	double handleNegativesTime = 0;
	double transposeTime = 0;
	double initialRefinementTime = 0;
	double refineTime = 0;
	double splitColorTime = 0;
	double findTargetTime = 0;
	double isolateTime = 0;
	double collectLinkingPairsTime = 0;
	double isPartitionDiscreteTime = 0;
	double packVectorsTime = 0;
	double allocateStorageTime = 0;
	double loopForNewColorsTime = 0;
	double setNewStackTime = 0;
	double removeAndAddColorsTime = 0;
	double init_time;
	double init_time2;
	double init_time3;

	// (Sparse vectors)
	forward_list<int> colorsToSplit;
	vector<int> initialParts;
	vector<int> color;
	vector<int> colReps;
	vector<int> rowReps;
	vector<bool> SCheck;
	vector<double> cdeg;
	vector<double> mincdeg;
	vector<double> maxcdeg;
	vector<int> numEdges;
	vector<bool> isAdj;
	list<int> colorsAdj;
	vector<bool> isolates;
	vector<bool> targets;
	vector<int> Astart;
    vector<int> Aindex;
    vector<double> Avalue;
	vector<double> AvalueCopy;
	vector<double> ARvalueCopy;
	vector<int> ARindex;
	vector<int> ARstart;
	vector<int> AR_Nend;
	vector<int> Xstart_;
	vector<int> Xindex_;
	vector<double> Xvalue_;
    vector<double> colCost;
    vector<double> colLower;
    vector<double> colUpper;
    vector<double> rowLower;
    vector<double> rowUpper;
    stack<int> S;
	vector<vector<double> > coeff;

	// (Dense vectors)
	vector<vector<int>> adjListLab;
	vector<vector<double>> adjListWeight;
	vector<vector<double>> adjListWeightReal;
	vector<list<int>* > C;
	vector<int> Csize;
	vector<list<int> > prevC;
	vector<forward_list<int>* > A;
	vector<int> Asize;

	// Partition storage
	vector<int> label;
	vector<int> labelTemp;
	vector<int> ptn;
	vector<int> ptnTemp;
	vector<int> colorStart;
	vector<int> colorStartTemp;

	// Contains linked variables from the splitting up of the ptns
	vector<pair<int, int> > linkingPairs;
	map<int, vector<int> > commonLinkers;

	// Contains previous ptns coloring so we know how to link variables
	vector<int> previousColumnColoring;
	vector<int> previousRowColoring;

	// Contains the representatives of each color class within a ptn at each major refinement
	vector<int> columnColorReps;

	// Contains ptn size information for tableau scaling
	vector<int> partSize;
	vector<int> previousPartSize;

	// List iterator pointers
	list<int>::iterator vPointer;
	forward_list<int>::iterator wPointer;
	forward_list<int>::iterator sPointer;
	list<int>::iterator cPointer;

	// Map to store the nodes with the same degree sum
	map<double, vector<int> > degSumNode;

	// List iterator 
	int v;
	int w;
	int c;
	int s;

	// Timer
	HighsTimer timer;

};

#endif