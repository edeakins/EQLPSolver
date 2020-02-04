#include "HighsQRmodule.h"

double HighsQR::norm(vector<double>& Aj){
	double sum = 0; 
	for (int i = 0; i < Aj.size(); ++i)
		sum += Aj[i] * Aj[i];
	return sqrt(sum);
}

void HighsQR::scalarMultiply(vector<double>& Aj, double scale, vector<double>& temp){
	if (isnan(scale))
		scale = 0;
	for (int i = 0; i < Aj.size(); ++i)
		temp[i] = Aj[i] * scale;
}

void HighsQR::scalarDiv(vector<double>& Aj, double r){
	for (int i = 0; i < Aj.size(); ++i)
		Aj[i] /= r;
}

void HighsQR::vectorSub(vector<double>& Aj, vector<double>& projection){
	for (int i = 0; i < Aj.size(); ++i){
		Aj[i] -= projection[i];
		if (fabs(Aj[i]) < 1e-5) Aj[i] = 0;
	}
}

double HighsQR::dotProduct(vector<double>& Aj, vector<double>& Ak){
	double sum = 0; 
	for (int i = 0; i < Aj.size(); ++i){
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

void HighsQR::gramSchmidt(vector<vector<double> >& Amatrix){
	int i, j;
	int n = Amatrix.size();
	int m = Amatrix[0].size();
	vector<double> temp(m, 0);
	for (i = 0; i < n; ++i){
		for (j = i + 1; j < n; ++j){
			project(Amatrix[i], Amatrix[j], temp);
			vectorSub(Amatrix[j], temp);
		}
	}
}