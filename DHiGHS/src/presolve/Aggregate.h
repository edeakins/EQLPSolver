#ifndef AGGREGATE_H_
#define AGGREGATE_H_

#include "HighsLp.h"
#include "HighsEquitable.h"
#include "HighsQRmodule.h"
#include "HighsTimer.h"
//#include "Highs.h"

class HighsAggregate{
public:
	HighsAggregate(HighsLp& lp, const struct eq_part* ep, HighsSolution& solution, HighsBasis& basis,
	int numRefinements);
	int update(const HighsSolution& solution, const HighsBasis& basis);
	void translateFrontsToColors();
	void findNonbasicRows();
	void findNonbasicCols();
	void packVectors();
	void foldObj();
	void foldMatrix();
	void fixMatrix();
	void findRowRepsToFix();
	void findColRepsToFix();
	void foldRhsInit();
	void foldRhs();
	void foldBndsInit();
	void foldBnds();
	void setRowBasis();
	void setColBasis();
	void addRows();
	void appendRowsToMatrix();
	void addCols();
	void identifyLinks();
	void createLinkRows();
	void reset();
	void countNumRefinements();
	void savePartition();
	void saveRowsAndColsFromLastSolve();
	void clearLp();
	void clearLinks();
	HighsLp* getAlp();
	HighsBasis* getBasis();
	void appendLinkersToLp();
	void appendColsToLpVectors(const int num_new_col,
                                  vector<double>& XcolCost,
                                  vector<double>& XcolLower,
                                  vector<double>& XcolUpper);
	void appendColsToLpMatrix(const int num_new_col,
                                 const int num_new_nz, vector<int>& XAstart,
                                 vector<int>& XAindex, vector<double>& XAvalue);
	void appendRowsToLpMatrix(const int num_new_row,
                                 const int num_new_nz, vector<int>& XAstart,
                                 vector<int>& XAindex, vector<double>& XAvalue);
	void appendRowsToLpVectors(const int num_new_row,
                                  vector<double>& XrowLower,
                                  vector<double>& XrowUpper);							   
	void transpose(vector<int>& xAstart, vector<int>& xAindex, vector<double>& xAvalue,
					vector<int>& xARstart, vector<int>& xARindex, vector<double> &xARvalue);
	// HighsLp& getAlp();
	// HighsBasis& getAlpBasis();

	/* Data structs to house the aggregated lp 
	and the aggregated lp basis information */
	HighsLp* alp;
	HighsBasis* alpBasis;
	HighsBasis* prevBasis;
	HighsSolution* prevSol;
	bool solve;
	bool solved;

	/* Scalar values that contain dimensional and
	name data about the original lp */
	int numRef;
	int iter;
	int numRow;
	int numCol;
	int numTot;
	int nnz;
	int numPairs;
	string model_name_;
	string lp_name_;

	/* Sparse vectors that store information about the original lp
	such as the A matrix, cost vector, and right hand sides.  Data is 
	stored in csc data format */
	vector<int> Astart;
    vector<int> Aindex;
    vector<double> Avalue;
    vector<double> colCost;
    vector<double> colLower;
    vector<double> colUpper;
    vector<double> rowLower;
    vector<double> rowUpper;

    // Linking pairs for new aggregate lp
	vector<int> parentPartition;
	vector<pair<int, int> > linkingPairs;

	/* Sparse storage for the previous iterations basis
	information.  We use this to set the basis of upcoming 
	aggregated lp so that no phase 1 simplex is needed. */
	vector<double> col_value;
    vector<double> row_value;
    vector<HighsBasisStatus> col_status;
	vector<HighsBasisStatus> row_status;
	vector<HighsBasisStatus> col_status_;
	vector<HighsBasisStatus> row_status_;
	vector<bool> nonBasicCol;
	vector<bool> nonBasicRow;


	/* New scalars for the current aggregated lp */
	int numRow_ = 0;
	int numCol_ = 0;
	int numTot_ = 0;
	int oldNumCol_ = 0;
	int newNumCol_ = 0;
	int numLinkers_ = 0;
	int previousNumCol_ = 0;
	int previousNumRow_ = 0;
	int previousNumSolveCol_ = 0;
	int previousNumSolveRow_ = 0;

	/* This is the sparse storage for the current 
	 aggregate lp.  This data will be uploaded to highs
	 with the current basis to obtain the lifted solution
	 and then used to pivot on linkers to a solution only in the 
	 lifted space. */
	vector<int> Astart_;
    vector<int> Aindex_;
    vector<int> linkers;
	vector<double> Avalue_;
    vector<double> colCost_;
    vector<double> colLower_;
    vector<double> colUpper_;
    vector<double> rowLower_;
    vector<double> rowUpper_;
	vector<double> AvaluePacked_;
	vector<int> AindexPacked_;
	vector<bool> rowRepsToFix;
	vector<double> rowRepsValue;
	vector<double> rowRepsScale;
	vector<bool> colRepsToFix;
	vector<double> colRepsValue;
	vector<bool> inMat;
	vector<int> parent;
	vector<int> child;
	vector<double> coeff;
	// For linker additions
	int maxLinkCols;
	int maxLinkSpace;
	vector<int> linkARstart;
	vector<int> linkARindex;
	vector<double> linkARvalue;
	vector<int> linkAlength;
	vector<bool> linked;
	vector<double> linkLB;
	vector<double> linkUB;
	vector<bool> skipLink;

	// For equitable partitions
	struct eq_part* partition;
	struct eq_part* previousPartition;
	vector<int> labels;
	vector<int> cell;
	vector<int> cellFront;
	vector<int> cellSize;
	vector<int> previousCell;
	vector<int> previousCellSize;
	vector<int> previousCellFront;
	vector<int> previousLabels;
	// Mapping for rows and cols from partition cells
	vector<int> colsToReps;
	vector<int> prevColsToReps;
	vector<int> repsToCols;
	vector<int> prevRepsToCols;
	vector<int> repsToRows;
	vector<int> prevRepsToRows;
	vector<int> rowsToReps;
	vector<int> prevRowsToReps;
	vector<int> prevCol;
	vector<int> col;
	vector<int> prevRow;
	vector<int> row;
	vector<int> lastSolveCol;
	vector<int> lastSolveRow;
	vector<bool> cellMarked;
	vector<bool> isRep;
	vector<int> cellToRow;
	vector<int> cellToCol;
	vector<int> newCellForRep;
	vector<int> cellReps;
	vector<bool> cellContainsOldRep;
};

#endif

