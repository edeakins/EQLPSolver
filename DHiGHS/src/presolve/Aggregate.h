#ifndef AGGREGATE_H_
#define AGGREGATE_H_

#include "HighsLp.h"
#include "HighsEquitable.h"
#include "HighsQRmodule.h"
#include "HighsTimer.h"
//#include "Highs.h"

class HighsAggregate{
public:
	HighsAggregate(HighsLp& lp, const HighsEquitable& ep, HighsSolution& solution, HighsBasis& basis, HighsTableau& tableau, bool flag);
	//virtual ~HighsAggregate(){};
	void clear();
	void aggregate();
	void aggregateAMatrix();
	void aggregateColBounds();
	void aggregateRowBounds();
	void aggregateCVector();
	void appendColsToLpVectors(const int num_new_col,
                                  const double* XcolCost,
                                  const double* XcolLower,
                                  const double* XcolUpper);
	void appendColsToLpMatrix(const int num_new_col,
                                 const int num_new_nz, const int* XAstart,
                                 const int* XAindex, const double* XAvalue);
	void appendRowsToLpMatrix(const int num_new_row,
                                 const int num_new_nz, const int* XARstart,
                                 const int* XARindex, const double* XARvalue);
	void appendRowsToLpVectors(const int num_new_row,
                                  const double* XrowLower,
                                  const double* XrowUpper);							   
	void transpose(vector<int>& xAstart, vector<int>& xAindex, vector<double>& xAvalue,
					vector<int>& xARstart, vector<int>& xARindex, vector<double> &xARvalue);
	void addLinkingRows();
	void collectColumns();
	void findLpBasis();
	void setInitialRhsAndBounds();
	void setRhsAndBounds();
	void appendLinkersToAMatrix(vector<double>& row, vector<int> idx, int rowIdx, int rIdx);
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
	void rowCoeff(int column);
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
    vector<int> Xindex;
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

    // (dense storage)
	vector<list<int>* > C;
	vector<list<int>* > prevC;
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
	vector<double> coeff;
	vector<int> coeffIdx;

    // Previous row coloring
	vector<int> previousRowColoring;

    // Previous column coloring
	vector<int> previousColumnColoring;

	// Timer
	
	};

#endif

