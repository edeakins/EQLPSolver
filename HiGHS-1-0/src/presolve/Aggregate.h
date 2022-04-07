#ifndef AGGREGATE_H_
#define AGGREGATE_H_

#include <cstring>

#include "HighsEquitable.h"
#include "HighsTimer.h"
//#include "Highs.h"

class HighsAggregate{
public:
	// Create and allocate
	void setUp(const struct lpPartition* ep, HighsLp* lp);
	void allocateAlp();
	void resizeElp();
	void resizeLpSym();
	void copyPartition();
	// Fold to alp
	void fold();
	void translateFrontsToColors();
	void packVectors();
	void foldObj();
	void foldMatrix();
	void foldRebuild();
	void foldRetract();
	void foldRhs();
	void foldBnd();
	// lift to elp
	void lift(HighsSolution& solution, HighsBasis& basis);
	void liftAMatrix();
	void liftObjective();
	void liftBnd();
	void liftRhs();
	void liftColBasis();
	void liftColBasis(HighsBasis& aBasis);
	void liftRowBasis();
	void liftRowBasis(HighsBasis& aBasis);
	void liftSolution(HighsSolution& aSolution);
	// utility functions
	HighsLp* getAlp();
	HighsBasis* getAlpBasis();
	HighsLp* getElp();
	HighsBasis* getElpBasis();
	HighsBasis* getLpBasis();
	HighsSolution* getLpSolution();
	void countNonbasicSplits();
	void makeLinks();
	void createLinkRows();
	void addRows();
	void appendRowsToMatrix();
	void addCols();
	void fixAstart();
	void appendLinkersToLp();
	/* Data structs to house the aggregated lp 
	and the aggregated lp basis information */
	HighsLp alpRetract; // Testing 
	HighsLp *alp;
	HighsLp *elp;
	HighsBasis* alpBasis;
	HighsBasis* elpBasis;
	HighsBasis* lpBasis;
	HighsSolution* lpSolution;
	HighsBasis prevBasis;
	HighsSolution prevSol;
	/* Scalar values that contain dimensional and
	name data about the original lp */
	int numRef;
	int iter;
	int numRow;
	int numCol;
	int numTot;
	int nnz;
	int numPairs;
	int elpNumRow_;
	int elpNumCol_;
	int elpNumTot_;
	int elpNnz_;
	int elpNumResNnz_;
	int elpNumResCols_;
	int elpNumResRows_;
	string model_name_;
	string lp_name_;
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
	int alpNumCol_ = 0;
	int alpNumRow_ = 0;
	int alpNumTot_ = 0;
	int alpNnz_ = 0;

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
	vector<bool> fixed;
	vector<bool> skip;
	vector<bool> rowRepsToFix;
	vector<double> rowRepsValue;
	vector<double> rowRepsScale;
	vector<bool> colRepsToFix;
	vector<double> colRepsValue;
	vector<bool> inMat;
	vector<int> parent;
	vector<int> parents;
	vector<int> child;
	vector<double> coeff;
	vector<bool> nonzero;
	vector<int> index;
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
	const struct lpPartition* partition;
	const struct eq_part* previousPartition;
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

