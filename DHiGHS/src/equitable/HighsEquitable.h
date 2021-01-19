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

struct HighsColoring{
    std::vector<int> lab;
    std::vector<int> unlab;
    std::vector<int> clen;
    std::vector<int> cfront;
};

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
    void setLabel(int index, int value);
    void addInduce(int who);
    void fixFronts(int cf, int ff);
    void colorAlloc();
    int *ints(int n) { return (int*)malloc(n * sizeof(int)); }

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

    // Saucy like storage (trying to optimize ep algorithm)
    struct HighsColoring coloring;

    /* Refinement: inducers */
	vector<char> indmark;   /* Induce marks */
	vector<int> ninduce;    /* Nonsingletons that might induce refinement */
	vector<int> sinduce;    /* Singletons that might induce refinement */
	int nninduce = 0;    /* Size of ninduce stack */
	int nsinduce = 0;    /* Size of sinduce stack */

	/* Refinement: marked cells */
	vector<int> clist;      /* List of cells marked for refining */
	int csize;       /* Number of cells in clist */

	/* Refinement: workspace */
	vector<char> stuff;     /* Bit vector, but one char per bit */
	vector<int> ccount;     /* Number of connections to refining cell */
	vector<int> bucket;     /* Workspace */
	vector<int> count;      /* Num vertices with same adj count to ref cell */
	vector<int> junk;       /* More workspace */
	vector<int> gamma;      /* Working permutation */
	vector<int> conncnts;   /* Connection counts for cell fronts */
    int *nextnon;    /* Forward next-nonsingleton pointers */
	int *prevnon;    /* Backward next-nonsingleton pointers */

    /* Search: split record */
	vector<int> splitwho;   /* List of where splits occurred */
	vector<int> splitfrom;  /* List of cells which were split */
	vector<int> splitlev;   /* Where splitwho/from begins for each level */
	int nsplits;     /* Number of splits at this point */

};

#endif