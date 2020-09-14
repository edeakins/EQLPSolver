#ifndef AGGREGATE_H_
#define AGGREGATE_H_

#include "HighsLp.h"
#include "HighsEquitable.h"
#include "HighsQRmodule.h"
#include "HighsTimer.h"
//#include "Highs.h"

class node{
public:
	node(int lab, int idx){
		label = lab;
		index = idx;
	}
	int label;
	int index;
};

// graph class to be used to determine connected components for GS
class linkGraph{
public:
	void initializeMat(int numNodes){
		adjMat.assign(numNodes, vector<int>(numNodes));
	}
	void addLinkNode(int label, int index){
		linkNodes.push_back(node(label, index));
		nodes.push_back(node(label, index));
	}
	void addRowNode(int label, int index){
		oldRowNodes.push_back(node(label, index));
		nodes.push_back(node(label, index));
	}
	void addNode(int label){
		nodes.push_back(node(label, numNodes));
		labelIndex[label] = numNodes;
		indexLabel[numNodes] = label;
		numNodes++;
	}
	void addEdge(int u, int v){
		int uIdx = labelIndex[u];
		int vIdx = labelIndex[v];
		adjMat[uIdx][vIdx] = 1;
		adjMat[vIdx][uIdx] = 1;
	}
	void connectedComponents(){
		vector<bool> visited(nodes.size());
		for (int i = 0; i < nodes.size(); ++i)
			if (!visited[i])
				components.push_back(DFS(i, visited));
	}
	vector<int> DFS(int i, vector<bool>& visited){
		vector<int> comp;
		stack<int> nodeStack;
		nodeStack.push(i);
		while (!nodeStack.empty()){
			i = nodeStack.top();
			nodeStack.pop();
			if (!visited[i]){
				comp.push_back(nodes[i].label);
				visited[i] = true;
			}
			for (int j = 0; j < adjMat[i].size(); ++j)
				if (adjMat[i][j] && !visited[j])
					nodeStack.push(j);
		}
		return comp;
	}
	map<int, int> labelIndex;
	map<int, int> indexLabel;
	vector<vector<int> > adjMat;
	vector<vector<int> > components;
	vector<node> linkNodes;
	vector<node> oldRowNodes;
	vector<node> nodes;
	int numNodes = 0;
};

class HighsAggregate{
public:
	HighsAggregate(HighsLp& lp, const HighsEquitable& ep, HighsSolution& solution, HighsBasis& basis, HighsTableau& tableau, bool flag);
	//virtual ~HighsAggregate(){};
	void clear();
	void aggregateAMatrix();
	void aggregateCVector();
	void collectColumns();
	void findLpBasis();
	void setInitialRhsAndBounds();
	void setRhsAndBounds();
	void appendLinkersToAMatrix(vector<double>& row, int rowIdx, int rIdx);
	void appendLinkersToRowRhs(int rowIdx);
	void appendLinkersToColBounds(int rIdx);
	void createRowWiseAMatrix();
	void setAggregateRealRowsRhs();
	void setAggregateRealColsBounds();
	void setAggregateLinkerColsBounds();
	void setAggregateLinkerRowsRhs();
	void findLinks();
	void initGSMatricesAndGraphs();
	void findMissingBasicColumns();
	void doGramSchmidt(int oldPart, int idx);
	void setAlpBasis();
	HighsLp& getAlp();
	HighsBasis& getAlpBasis();
	bool dependanceCheck(vector<double> &v); // int impliedRowsIdx); // int linkIdx = -1);
	void findPreviousBasisForRows();
	void findPreviousBasisForColumns();
	vector<double> rowCoeff(int column);
	void eraseLinkersIfNotNeeded();
	bool varIsBounded(pair<int, int> link);
	void editRowWiseMatrix(int domLink, int slavLink);
	void createImpliedRows(HighsLp& lp);
	vector<double> aggregateImpliedRow(int impliedRow);
	void getAggImpliedRows();
	vector<double> createImpliedLinkRows(vector<double> &linkRow, int linkIdx);
	void liftTableau();
	void transposeMatrix();
	void collectRowsForGS();
	void collectPartsForGS();
	void examinePartition();
	void initialAggregateAMatrix();
	void trackRowColors(HighsLp& ep);
	void findLinkerComponents();

	// Lp to store the aggregate LP into
	HighsLp alp;
	HighsBasis alpBasis;
	bool flag_;
	int iterations = 0;
	int pivots = 0;
	int masterIter = 0;
	int masterIter_ = 0;

	// Copy original lp data from equitable partition
	// scalars
	int numRow;
	int impliedNumRow;
	int numCol;
	int numTot;
	int realNumCol;
	int realNumRow;
	string model_name_;
	string lp_name_;

	// (sparse storage) Original lp
	vector<int> Astart;
	vector<int> impliedARstart;
    vector<int> Aindex;
    vector<int> impliedARindex;
    vector<double> Avalue;
    vector<double> impliedARvalue;
    vector<int> color;
    vector<int> linkerColor;
    vector<double> colCost;
    vector<double> colLower;
    vector<double> colUpper;
    vector<double> rowLower;
    vector<double> rowUpper;

    // Linking pairs for new aggregate lp
    vector<pair<int, int> > linkingPairs;

   	// Storage to obtain "warm start" for new small lp
    vector<double> col_value;
    vector<double> row_value;
    vector<HighsBasisStatus> col_status;
    vector<HighsBasisStatus> row_status;
    HighsQR QR;
    vector<vector<vector<double> > > QRStorage;
    vector<int> QRIndexUpdate;

    // (dense storage)
	vector<vector<int> > C;
	vector<vector<int> > prevC;
	vector<vector<int> > adjListLab;
	vector<vector<double> > adjListWeight;

	// For the new "smaller" lp
	int numRowAfterImp_ = 0;
	int numLiftedRow_ = 0;
	int numRow_ = 0;
	int numCol_ = 0;
	int numTot_ = 0;
	int realNumRow_ = 0;
	int realNumCol_ = 0;
	int realNumTot_ = 0;
	int numActiveRows_ = 0;
	int numActiveBounds_ = 0;
	int numLinkers_ = 0;
	int numLinkers;
	int originalNumLinkers;
	vector<vector<int> > adjMatLinks_;
	vector<int> startingBasicRows_;
	vector<int> startingBasicColumns_;
	vector<int> potentialBasicRows_;
	vector<int> potentialBasicColumns_;
	vector<int> Astart_;
	vector<int> ARstart_;
	vector<int> ARtableauStart;
    vector<int> Aindex_;
    vector<int> ARindex_;
    vector<int> ARtableauIndex;
    vector<int> A_Nend_;
    vector<int> GSRstart_;
    vector<int> GSRindex_;
    vector<int> linkers;
	vector<int> artificialVariables;
    vector<bool> activeConstraints_;
    vector<bool> activeBounds_;
    vector<double> Avalue_;
    vector<double> ARvalue_;
    vector<double> ARtableauValue;
    vector<double> GSRvalue_;
    vector<double> colCost_;
    vector<double> colLower_;
    vector<double> colUpper_;
    vector<double> rowLower_;
    vector<double> rowUpper_;
    vector<double> scale;
    vector<double> ARreducedRHS;

    // Previous row coloring
    map<int, HighsBasisStatus> previousRowInfo;
	map<int, double> previousRowValue;
	vector<int> previousRowColoring;

    // Previous column coloring
    map<int, HighsBasisStatus> previousColumnInfo;
	map<int, double> previousColumnValue;
	vector<int> previousColumnColoring;

	map<int, vector<int> > commonLinkers;

	// Collect implied linkers for constraint manipulation
	vector<pair<int, int> > equalColors;

	// Holds all aggregated implied rows
	vector<vector<double> > aggImpliedRows;
	map<int, vector<vector<double> > > impliedRowsByColor;

	// Contains partition size information for tableau scaling
	vector<int> partSize;
	vector<int> previousPartSize;

	// Stores the actual color of rows
	vector<int> rowColor;
	vector<int> rowColor_;

	// Keeps track of the history of constraints and their colors so we know when 
	// a constraint become active (constraint became active while in this color class)
	vector<bool> activeColorHistory;
	vector<bool> activeColorHistory_;

	/* Gram schmidt storage and indices.  Mostly storage for Gram schmidt matrices
	and their index information which we use to determine whether a "linker" variable 
	requires pivoting or not. */
	int numRowGS_ = 0;
	vector<int> AstartGS_;
	vector<int> ARstartGS_;
	vector<int> AindexGS_;
	vector<int> ARindexGS_;
	vector<int> AR_NendGS_;
	vector<double> AvalueGS_;
	vector<double> ARvalueGS_;
	vector<int> partsForGS;
	vector<bool> colorsForGS;
	vector<int> numSplits;
	vector<bool> linkIsNeeded;
	vector<bool> linkIsErased;
	map<int, int> impliedIdx;
	map<int, vector<int> > linksInMatrix;
	map<int, vector<int> > varToLink;
	map<int, vector<int> > connectedColorsAndLinks;
	map<int, linkGraph> cLGraphs;
	map<int, vector<vector<double> > > QRMatrices;
	};

#endif

