#ifndef HMATRIX_H_
#define HMATRIX_H_

#include "HVector.h"

#include <vector>
#include <cmath> 
#include <iostream> 
using namespace std;

class HMatrix {
public:
    const int *getAstart() const {
        return &Astart[0];
    }
    const int *getAindex() const {
        return &Aindex[0];
    }
    const double *getAvalue() const {
        return &Avalue[0];
    }
    void setup(int numCol, int numRow, const int *Astart, const int *Aindex,
            const double *Avalue);
    void setupOC(int numCol, int numRow, const int *Astart, const int *Aindex,
            const double *Avlaue, const vector<bool> &startingBasis);
    void getActiveConstraints(vector<int> &basicIndex, vector<double> &baseValue, 
                              vector<bool> &activeConstraints, vector<double> &rhs, vector<double> &workValue);
    void price_by_col(HVector& row_ap, HVector& row_ep) const;
    void price_by_row(HVector& row_ap, HVector& row_ep) const;
    void update(int columnIn, int columnOut);
    double compute_dot(HVector& vector, int iCol) const;
    void collect_aj(HVector& vector, int iCol, double multi) const;

    void compute_vecT_matB(const double *vec, const int *base, HVector *res);
    void compute_matB_vec(const double *vec, const int *base, HVector *res);
public:
    int numCol;
    int numRow;
    vector<int> Astart;
    vector<int> Aindex;
    vector<double> Avalue;

    vector<int> ARstart;
    vector<int> AR_Nend;
    vector<int> ARindex;
    vector<double> ARvalue;

};

#endif /* HMATRIX_H_ */
