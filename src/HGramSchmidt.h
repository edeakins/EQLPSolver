
#ifndef HGRAM_H_
#define HGRAM_H_
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
using namespace std;

class HGramSchmidt{
public:
    /* Deakins - Classes */
    //HGramSchmidt();
    void gramSchmidt (vector<vector<double> > &activeSet, 
                                       vector<int> &startingBasis, 
                                       vector<int> &activeIdx,
                                       int &aggNumCol,
                                       int &aggNumRow,
                                       int &aggNumTot,
                                       int &rCnt,
                                       vector<int> &residuals);
    double norm (vector<double> &x, int length);
    void vec_copy (double * x, double * y, int length);
    void partialvec_copy (double * x, double * y, int length, int index);
    void scalar_div (double * x, double r, int length, double * y);
    void scalar_sub (double * x, double r, int length, double * y);
    void partialscalar_sub (double * x, double r, int length, int index, double * y); 
    double dotProduct (vector<double> &u, vector<double> &v, int length);
    double partialdot_product (double * x, double * y, int length, int index); 
    double subdot_product (double * x, double * y, int length, int index);   
    vector<double> project(vector<double> &v, vector<double> &q);
    void vecMinus(vector<double> &v, vector<double> &projection);
    bool indepCheck(vector<double> &u);
    vector<double> vecDivide(vector<double> &v, double scale);
};
#endif /* HGRAM_H_ */