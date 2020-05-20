#ifndef HIGHS_QR
#define HIGHS_QR

#include <algorithm>
#include <cmath>
#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <vector>

using namespace std;
class HighsQR{
public:
	double norm(vector<double>& Aj);
	void normalize(vector<double>& Aj, vector<double>& temp);
	void scalarMultiply(vector<double>& Aj, double scale, vector<double>& temp);
	void scalarDiv(vector<double>& Aj, double r, vector<double>& temp);
	void vectorSub(vector<double>& Aj, vector<double>& projection);
	double dotProduct(vector<double>& Aj, vector<double>& Ak);
	void project(vector<double>& u, vector<double>& v, vector<double>& temp);
	void gramSchmidt(vector<vector<double> >& Amatrix, int idx);
	int linkStartIdx;
};

#endif
