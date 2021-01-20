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

// struct HighsColoring{
//     int *lab;
//     int *unlab;
//     int *clen;
//     int *cfront;
// };

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
    /* mimicing saucy */
    // void setLabel(int index, int value);
    // void addInduce(int who);
    // void fixFronts(int cf, int ff);
    // void colorAlloc();
   
    // int refineSaucy();
    // bool atTerminal();
    // void clearRefine();
    // void swapLabels(int a, int b);
    // void moveToBack(int k);
    // void dataMark(int k);
    // void lp2Graph();
    // void fixAdj();
    // int refineCell();
    // int refineCells();
    // int refSingleton(int* adj, int* edg, int cf);
    // int refNonsingle(int* adj, int* edg, int cf);
    // int refSingletonUndirected(int cf);
    // int refNonsingleUndirected(int cf);
    // int refSingleCell(int cf);
    // int refNonsingleCell(int cf);
    // // int refNonsingle();
    // // int refine2(int cf);
    // void introsort(int* a, int n);
    // void introsortLoop(int* a, int n, int lim);
    // void insertionSort(int* a, int n);
    // int logBase2(int n);
    // int median(int a, int b, int c);
    // void heapSort(int* a, int n);
    // int maybeSplit(int cf, int ff);
    // int partition(int* a, int n, int m);
    // void swap(int* a, int x, int y);
    // void siftUp(int* a, int k);
    // void siftDown(int* a, int n);

    // // Saucy allocation functions
    // int* ints(int n) {return (int*)malloc(n * sizeof(int));}
    // double* doubles(int n) {return (double*)malloc(n * sizeof(double));}
    // int* zeros(int n) {return(int*)malloc(n * sizeof(int));}
    // char* bits(int n) {return(char*)malloc(n * sizeof(char));}

	/// Original LP info
    int nRows;
    int nCols; 
    int nTot;
    int nnz;
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
    vector<vector<int> > prevC;
    vector<vector<int> > A;
    vector<int> color;
    vector<int> previousRowColoring;
    vector<int> previousColumnColoring;
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

    // // Saucy like storage (trying to optimize ep algorithm)
    // struct HighsColoring coloring;

    // // Saucy miscelaneous used in their code
    // int* junk;

    // /* graph storage */
    // int* adj;
    // int* edg;
    // double* wt;

    // /* Refinement: inducers */
	// char* indmark;   /* Induce marks */
	// int* ninduce;    /* Nonsingletons that might induce refinement */
	// int* sinduce;    /* Singletons that might induce refinement */
	// int nninduce = 0;    /* Size of ninduce stack */
	// int nsinduce = 0;    /* Size of sinduce stack */

	// /* Refinement: marked cells */
	// int* clist;      /* List of cells marked for refining */
	// int csize = 0;  
    // int lev = 0;     /* Number of cells in clist */

	// /* Refinement: workspace */
	// char* stuff;     /* Bit vector, but one char per bit */
	// int* ccount;     /* Number of connections to refining cell */
	// int* bucket;     /* Workspace */
	// int* count;      /* Num vertices with same adj count to ref cell */
	// int* junk;       /* More workspace */
	// int* gamma;      /* Working permutation */
	// int* conncnts;   /* Connection counts for cell fronts */    /* Backward next-nonsingleton pointers */
    // int* nextnon;
    // int* prevnon;

    // /* Search: split record */
	// int* splitwho;   /* List of where splits occurred */
	// int* splitfrom;  /* List of cells which were split */
	// int* splitlev;   /* Where splitwho/from begins for each level */
	// int nsplits = 0;     /* Number of splits at this point */

};

#endif