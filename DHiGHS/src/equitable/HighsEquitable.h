#ifndef HIGHS_EQUITABLE_H
#define HIGHS_EQUITABLE_H

#include "HighsLp.h"
#include "HighsSimpleDec.h"

#include <string>
#include <vector>
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
	// Setup for equitable partition
	void clear();
	void setup(const HighsLp& lp);
	void lp2Graph();
	void initialRefinement();
	void splitColor(int color);
	void refine();
	void findTarget();
	void isolate(int i);
	void collectLinkingPairs();
	bool isPartitionDiscrete();

	// Some storage for info about lp and partition
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
	int isolated = -1;
	int refinements = 0;

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
	vector<int> Xstart_;
	vector<int> Xindex_;
	vector<double> Xvalue_;
    vector<double> colCost;
    vector<double> colLower;
    vector<double> colUpper;
    vector<double> rowLower;
    vector<double> rowUpper;
    stack<int> S;

	// (Dense vectors)
	vector<vector<int>> adjListLab;
	vector<vector<double>> adjListWeight;
	vector<vector<double>> adjListWeightReal;
	vector<list<int>* > C;
	vector<vector<int>> prevC;
	vector<forward_list<int>* > A;

	// Contains linked variables from the splitting up of the partitions
	vector<pair<int, int> > linkingPairs;
	map<int, vector<int> > commonLinkers;
	// Contains previous partitions coloring so we know how to link variables
	vector<int> previousColumnColoring;
	vector<int> previousRowColoring;

	// Contains the representatives of each color class within a partition at each major refinement
	vector<int> columnColorReps;

	// Contains partition size information for tableau scaling
	vector<int> partSize;
	vector<int> previousPartSize;

	// List iterator pointers
	list<int>::iterator vPointer;
	forward_list<int>::iterator wPointer;
	forward_list<int>::iterator sPointer;
	list<int>::iterator cPointer;

	// List iterator 
	int v;
	int w;
	int c;
	int s;

};

#endif