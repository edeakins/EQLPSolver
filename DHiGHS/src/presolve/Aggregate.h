#ifndef AGGREGATE_H_
#define AGGREGATE_H_

#include "HighsLp.h"
#include "HighsEquitable.h"
#include "HighsQRmodule.h"
//#include "Highs.h"

class HighsAggregate{
public:
	HighsAggregate(){
	}
	HighsAggregate(const HighsLp& lp, const HighsEquitable& ep, HighsSolution* solution, HighsBasis* basis, bool flag);
	~HighsAggregate(){
		printf("\ndtor HighsAggregate called\n");
		cin.get();
	}
	void clear();
	void build(const HighsEquitable& ep, HighsSolution& solution, HighsBasis& basis, bool flag);
	void setup(HighsLp& lp);
	void aggregateAMatrix();
	void setLinkerMatrix();
	void createRowWiseAMatrix();
	void setAggregateRealRowsRhs();
	void setAggregateRealColsBounds();
	void setAggregateLinkerColsBounds();
	void setAggregateLinkerRowsRhs();
	void findLinks();
	void findMissingBasicColumns();
	void setAlpBasis();
	HighsLp& getAlp();
	HighsBasis& getAlpBasis();
	bool dependanceCheck(vector<double> &v);
	void findPreviousBasisForRows();
	void findPreviousBasisForColumns();
	vector<double> rowCoeff(int column);

	// Lp to store the aggregate LP into
	HighsLp alp;
	HighsBasis alpBasis;
	bool flag_;
	int iterations = 0;
	int pivots = 0;

	// Copy original lp data from equitable partition
	// scalars
	int numRow;
	int numCol;
	int numTot;
	string model_name_;
	string lp_name_;

	// (sparse storage) Original lp
	vector<int> Astart;
    vector<int> Aindex;
    vector<int> color;
    vector<int> linkerColor;
    vector<double> Avalue;
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

    // (dense storage)
	vector<vector<int>> C;
	vector<vector<int>> prevC;
	vector<vector<int>> adjListLab;
	vector<vector<double>> adjListWeight;

	// For the new "smaller" lp
	int numRow_;
	int numCol_;
	int numTot_;
	int numActiveRows_;
	int numActiveBounds_;
	int numLinkers_;
	vector<int> startingBasicRows_;
	vector<int> startingBasicColumns_;
	vector<int> potentialBasicRows_;
	vector<int> potentialBasicColumns_;
	vector<int> Astart_;
	vector<int> ARstart_;
    vector<int> Aindex_;
    vector<int> ARindex_;
    vector<int> AR_Nend_;
    vector<bool> activeConstraints_;
    vector<bool> activeBounds_;
    vector<double> Avalue_;
    vector<double> ARvalue_;
    vector<double> colCost_;
    vector<double> colLower_;
    vector<double> colUpper_;
    vector<double> rowLower_;
    vector<double> rowUpper_;

    // Basis for warm starting the lp
    // HighsBasis warmStart_;

    // Previous row coloring
    map<int, HighsBasisStatus> previousRowInfo;
	map<int, double> previousRowValue;
	vector<int> previousRowColoring;

    // Previous column coloring
    map<int, HighsBasisStatus> previousColumnInfo;
	map<int, double> previousColumnValue;
	vector<int> previousColumnColoring;
};

#endif

