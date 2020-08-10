#include "HighsQRmodule.h"

double HighsQR::norm(vector<double>& Aj){
	double sum = 0;
	for (int i = 0; i < linkStartIdx; ++i)
		sum += Aj[i] * Aj[i];
	return sqrt(sum);
}

void HighsQR::normalize(vector<double>& Aj, vector<double>& temp){
	double scale = norm(Aj);
	if (isnan(scale)) scale = 0;
	scalarDiv(Aj, scale, temp);
}

void HighsQR::scalarMultiply(vector<double>& Aj, double scale, vector<double>& temp){
	if (isnan(scale))
		scale = 0;
	for (int i = 0; i < Aj.size(); ++i)
		temp[i] = Aj[i] * scale;
}

void HighsQR::scalarDiv(vector<double>& Aj, double r, vector<double>& temp){
	for (int i = 0; i < Aj.size(); ++i){
		double val = Aj[i] / r;
		temp[i] = (isnan(val) || isinf(val)) ? 0 : val;
	}
}

void HighsQR::vectorSub(vector<double>& Aj, vector<double>& projection){
	for (int i = 0; i <	Aj.size(); ++i){
		Aj[i] -= projection[i];
		if (fabs(Aj[i]) < 1e-5) Aj[i] = 0;
	}
}

double HighsQR::dotProduct(vector<double>& Aj, vector<double>& Ak){
	double sum = 0;
	for (int i = 0; i < linkStartIdx; ++i){
		sum += Aj[i] * Ak[i];
	}
	return sum;
}

void HighsQR::project(vector<double>& u, vector<double>& v, vector<double>& temp){
	double numer = dotProduct(v, u);
	double denom = dotProduct(u, u);
	double scale = numer/denom;
	scalarMultiply(u, scale, temp);
}

void HighsQR::gramSchmidt(vector<vector<double> >& Amatrix, int colIdx, int rowIdx){
	int i, j; 
	double scale = 0;
	int n = Amatrix.size();
	int m = Amatrix[0].size();
	linkStartIdx = colIdx;
	vector<double> temp_i(m, 0);
	vector<double> temp_j(m, 0);
	for (i = rowIdx; i < n; ++i){
		normalize(Amatrix[i], temp_i);
		for (j = i + 1; j < n; ++j){
			double scale = dotProduct(temp_i, Amatrix[j]);
			scalarMultiply(temp_i, scale, temp_j);
			vectorSub(Amatrix[j], temp_j);
		}
	}
}