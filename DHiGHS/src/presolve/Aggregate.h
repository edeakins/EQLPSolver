#ifndef AGGREGATE_H_
#define AGGREGATE_H_

#include "HighsLp.h"
#include "HighsEquitable.h"
#include "HighsQRmodule.h"
#include "HighsTimer.h"
//#include "Highs.h"

class HighsAggregate{
public:
	HighsAggregate(HighsLp& lp, const HighsEquitable& ep, HighsSolution& solution, HighsBasis& basis);
	//virtual ~HighsAggregate(){};
	void clear();
	void aggregate();
	void aggregateAMatrix();
	void aggregateColBounds();
	void aggregateRowBounds();
	void aggregateCVector();
	void setColBasis();
	void setRowBasis();
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
	HighsLp& getAlp();
	HighsBasis& getAlpBasis();

	/* Data structs to house the aggregated lp 
	and the aggregated lp basis information */
	HighsLp alp;
	HighsBasis alpBasis;

	/* Scalar values that contain dimensional and
	name data about the original lp */
	int numRow;
	int numCol;
	int numTot;
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
    vector<pair<int, int> > linkingPairs;

	/* Sparse storage for the previous iterations basis
	information.  We use this to set the basis of upcoming 
	aggregated lp so that no phase 1 simplex is needed. */
	vector<double> col_value;
    vector<double> row_value;
    vector<HighsBasisStatus> col_status;
	vector<HighsBasisStatus> row_status;

   	/* Sparse storage for the new aggregate lp's "warm start" 
	basis.  Note that we do not pre load the solutions to this
	lp as they are clear from the pre loaded basis. */
	vector<HighsBasisStatus> col_status_;
	vector<HighsBasisStatus> row_status_;

    /* Dense storage for equitable partition information.
	 These are doubly linked lists of elements by their color class.
	 We have the previous color class information and current color class
	 information. */
	vector<list<int>* > C;
	vector<list<int> > prevC;
	vector<int> previousRowColoring;
	vector<int> previousColumnColoring;
	vector<int> AindexSub;
	vector<int> color;

	/* New scalars for the current aggregated lp */
	int numRow_ = 0;
	int numCol_ = 0;
	int numTot_ = 0;
	int numLinkers_ = 0;

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
};

#endif

