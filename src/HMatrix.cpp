#include "HMatrix.h"
#include "HConst.h"
#include <cmath>
#include <iostream>
#include <cstdio>

using namespace std;

void HMatrix::setup(int numCol_, int numRow_, const int *Astart_,
        const int *Aindex_, const double *Avalue_) {
    // Copy A
    numCol = numCol_;
    numRow = numRow_;
    Astart.assign(Astart_, Astart_ + numCol_ + 1);

    int AcountX = Astart_[numCol_];
    Aindex.assign(Aindex_, Aindex_ + AcountX);
    Avalue.assign(Avalue_, Avalue_ + AcountX);

    // Build row copy - pointers
    ARstart.resize(numRow + 1);
    AR_Nend.assign(numRow, 0);
    for (int k = 0; k < AcountX; k++)
        AR_Nend[Aindex[k]]++;
    ARstart[0] = 0;
    for (int i = 1; i <= numRow; i++)
        ARstart[i] = ARstart[i - 1] + AR_Nend[i - 1];
    for (int i = 0; i < numRow; i++)
        AR_Nend[i] = ARstart[i];
    // Build row copy - elements
    ARindex.resize(AcountX);
    ARvalue.resize(AcountX);
    for (int iCol = 0; iCol < numCol; iCol++) {
        for (int k = Astart[iCol]; k < Astart[iCol + 1]; k++) {
            int iRow = Aindex[k];
            int iPut = AR_Nend[iRow]++;
            ARindex[iPut] = iCol;
            ARvalue[iPut] = Avalue[k];
        }
        // cout << "AR_Nend: { ";
        // for (int k = 0; k < AR_Nend.size(); ++k){
        //     cout << AR_Nend[k] << " ";
        // }
        // cout << "}" << endl;
        // cin.get();
    }

    // for (int i = 0; i < numRow; ++i){
    //     cout << "row: " << i << " [" ;
    //     for (int j = ARstart[i]; j < ARstart[i + 1]; ++j){
    //         cout << ARvalue[j] << " ";
    //     }
    //     cout << "]" << endl;
    // }
    // cin.get();
}

void HMatrix::setupOC(int numCol_, int numRow_, const int *Astart_,
        const int *Aindex_, const double *Avalue_, const vector<bool> &Basis){
    // Copy A
    numCol = numCol_;
    numRow = numRow_;
    Astart.assign(Astart_, Astart_ + numCol_ + 1);

    int AcountX = Astart_[numCol_];
    Aindex.assign(Aindex_, Aindex_ + AcountX);
    Avalue.assign(Avalue_, Avalue_ + AcountX);

    // Build row copy - pointers
    ARstart.resize(numRow + 1);
    AR_Nend.assign(numRow, 0);

    // Build row copy - pointers
    ARstart.resize(numRow + 1);
    AR_Nend.assign(numRow, 0);
    for (int k = 0; k < AcountX; k++){
            AR_Nend[Aindex[k]]++; 
    }
    ARstart[0] = 0;
    for (int i = 1; i <= numRow; i++)
        ARstart[i] = ARstart[i - 1] + AR_Nend[i - 1];
    for (int i = 0; i < numRow; i++)
        AR_Nend[i] = ARstart[i];
    // Build row copy - elements
    ARindex.resize(AcountX);
    ARvalue.resize(AcountX);
    for (int iCol = 0; iCol < numCol; iCol++) {
        if (Basis[iCol])
            continue;
        for (int k = Astart[iCol]; k < Astart[iCol + 1]; k++){
            int iRow = Aindex[k];
            int iPut = AR_Nend[iRow]++;
            ARindex[iPut] = iCol;
            ARvalue[iPut] = Avalue[k];
        }
    }
    for (int iCol = 0; iCol < numCol; iCol++) {
        if (!Basis[iCol])
            continue;
        for (int k = Astart[iCol]; k < Astart[iCol + 1]; k++){
            int iRow = Aindex[k];
            int iPut = AR_Nend[iRow]++;
            ARindex[iPut] = iCol;
            ARvalue[iPut] = Avalue[k];
        }
    }
    // cout << "AR_Nend: { ";
    //     for (int k = 0; k < AR_Nend.size(); ++k){
    //         cout << AR_Nend[k] << " ";
    //     }
    //     cout << "}" << endl;
    //     cin.get();
    
    for (int i = 0; i < numRow; ++i){
        for (int j = ARstart[i]; j < ARstart[i + 1]; ++j){
            if (Basis[ARindex[j]]){
                AR_Nend[i]--;
            }
        }
    }
    // cout << "AR_Nend: { ";
    //     for (int k = 0; k < AR_Nend.size(); ++k){
    //         cout << AR_Nend[k] << " ";
    //     }
    //     cout << "}" << endl;
    //     cin.get();
    vector<vector<int> > mat;
    vector<int> ins(numCol, 0);
    for (int i = 0; i < numRow; ++i){
        for (int j = ARstart[i]; j < ARstart[i + 1]; ++j){
            ins[ARindex[j]] = ARvalue[j];
        }
        mat.push_back(ins);
        ins.assign(numCol,0);
    }
    // cout << "[ " << endl;
    // for (int i = 0; i < mat.size(); ++i){
    //     for (int j = 0; j < mat[i].size(); ++j){
    //         if (j == mat[i].size() - 1){
    //             cout << mat[i][j] << "";
    //             continue;
    //         }
    //         cout << mat[i][j] << ", ";
    //     }
    //     cout << ";" << endl;
    // }
    // cout << " ]" << endl;
    // cin.get();
    
}

void HMatrix::getActiveConstraints(vector<int> &basicIndex, vector<double> &baseValue, 
                                   vector<bool> &activeConstraints, vector<double> &rhs, vector<double> &workValue){
    int i, j, k;
    vector<double> trueRhs(rhs.size(), 0);
    vector<bool> nonBasic(numCol, true);
    // for (i = 0; i < rhs.size(); ++i){
    //     cout << i << " is " << rhs[i] << endl;
    // }
    // cout << rhs.size() << endl;
    // cin.get();
    for (i = 0; i < basicIndex.size(); ++i){
        if (basicIndex[i] < numCol)
            nonBasic[basicIndex[i]] = false;
    }
    for (i = 0; i < nonBasic.size(); ++i){
        //cout << "var: " << basicIndex[i] << " = " << baseValue[i] << endl;
        //cin.get();
        for (j = 0; j < rhs.size(); ++j){
            for (k = ARstart[j]; k < ARstart[j + 1]; ++k){
                if (nonBasic[i] && i == ARindex[k])
                    trueRhs[j] += workValue[i] * ARvalue[k];
            }
        }
    }
    for (i = 0; i < basicIndex.size(); ++i){
        if (basicIndex[i] >= numCol)
            continue;
        for (j = 0; j < rhs.size(); ++j){
            for (k = ARstart[j]; k < ARstart[j + 1]; ++k){
                if (basicIndex[i] == ARindex[k])
                    trueRhs[j] += baseValue[i] * ARvalue[k];
            }
        }
    }
    for (i = 0; i < trueRhs.size(); ++i){
        // cout << "i: " << i << endl;
        if (fabs(trueRhs[i]) <= 1e-6)
            trueRhs[i] = 0;
        // cout << "row: " << i << " true: " << trueRhs[i] << "  given: " << rhs[i] << endl; 
        // cin.get();
        if (fabs(trueRhs[i] - rhs[i]) <= 1e-6){
            activeConstraints[i] = true;
            //cout << "row: " << i << endl;
        }
        // else{
        //     cout << "weird row: " << i << endl;
        // }
    }
    // for (i = 0; i < activeConstraints.size(); ++i){
    //     cout << i << " is " << activeConstraints[i] << endl;
    // }
    // cin.get();
}

void HMatrix::update(int columnIn, int columnOut) {
    if (columnIn < numCol) {
        for (int k = Astart[columnIn]; k < Astart[columnIn + 1]; k++) {
            int iRow = Aindex[k];
            int iFind = ARstart[iRow];
            int iSwap = --AR_Nend[iRow];
            while (ARindex[iFind] != columnIn)
                iFind++;
            swap(ARindex[iFind], ARindex[iSwap]);
            swap(ARvalue[iFind], ARvalue[iSwap]);
        }
    }

    if (columnOut < numCol) {
        for (int k = Astart[columnOut]; k < Astart[columnOut + 1]; k++) {
            int iRow = Aindex[k];
            int iFind = AR_Nend[iRow];
            // cout << "iRow: " << iRow << endl;
            // cout << "ARstart: " << ARstart[iRow] << endl;
            // cout << "AR_Nend: " << AR_Nend[iRow] << endl;
            int iSwap = AR_Nend[iRow]++;
            // cout << "AR_Nend: " << AR_Nend[iRow] << endl;
            // cin.get();
            while (ARindex[iFind] != columnOut)
                iFind++;
            // cout << "ARindex of col out: " << ARindex[iFind] << endl;
            // cout << "ARindex of swap: " << ARindex[iSwap] << endl;
            // cout << "ARvalue of col out: " << ARvalue[iFind] << endl;
            // cout << "ARvalue of swap: " << ARvalue[iSwap] << endl;
            // cout << "ARindex: { ";
            // for (int j = ARstart[iRow]; j < ARstart[iRow + 1]; ++j){
            //     cout << ARindex[j] << " ";
            // }
            // cout << "}" << endl;
            // cin.get();
            swap(ARindex[iFind], ARindex[iSwap]);
            swap(ARvalue[iFind], ARvalue[iSwap]);
            // cout << "ARindex: { ";
            // for (int j = ARstart[iRow]; j < ARstart[iRow + 1]; ++j){
            //     cout << ARindex[j] << " ";
            // }
            // cout << "}" << endl;
            // cin.get();
        } 
    }
    // cout << "AR_Nend: { ";
    //     for (int k = 0; k < AR_Nend.size(); ++k){
    //         cout << AR_Nend[k] << " ";
    //     }
    //     cout << "}" << endl;
    //     cin.get();
    
    // for (int i = 0; i < numRow; ++i){
    //     cout << "row: " << i << " { ";
    //     for (int j = ARstart[i]; j < ARstart[i + 1]; ++j){
    //         cout << ARvalue[j] << " ";
    //     }
    //     cout << "}" << endl;
    // }
    // cin.get();
    // cout << "ARindex: { ";
    //     for (int k = 0; k < ARindex.size(); ++k){
    //         cout << ARindex[k] << " ";
    //     }
    //     cout << "}" << endl;
    //     cin.get();
}

double HMatrix::compute_dot(HVector& vector, int iCol) const {
    double result = 0;
    if (iCol < numCol) {
        for (int k = Astart[iCol]; k < Astart[iCol + 1]; k++)
            result += vector.array[Aindex[k]] * Avalue[k];
    } else {
        result = vector.array[iCol - numCol];
    }
    return result;
}

void HMatrix::collect_aj(HVector& vector, int iCol, double multi) const {
	
    //cout << "iCol: " << iCol << endl;
    if (iCol < numCol) {
        for (int k = Astart[iCol]; k < Astart[iCol + 1]; k++) {
            int index = Aindex[k];
            //cout << "index: " << index << endl;
            double value0 = vector.array[index];
            double value1 = value0 + multi * Avalue[k];
            if (value0 == 0)
                vector.index[vector.count++] = index;
            vector.array[index] =
                    (fabs(value1) < HSOL_CONST_TINY) ? HSOL_CONST_ZERO : value1;
        }
    } else {
        int index = iCol - numCol;
        //cout << "index: " << index << endl;
        double value0 = vector.array[index];
        double value1 = value0 + multi;
        if (value0 == 0)
            vector.index[vector.count++] = index;
        vector.array[index] =
                (fabs(value1) < HSOL_CONST_TINY) ? HSOL_CONST_ZERO : value1;
    }
}

void HMatrix::price_by_col(HVector& row_ap, HVector& row_ep) const {
    // Alias
    int ap_count = 0;
    int *ap_index = &row_ap.index[0];
    double *ap_array = &row_ap.array[0];
    const double *ep_array = &row_ep.array[0];

    // Computation
    for (int iCol = 0; iCol < numCol; iCol++) {
        double value = 0;
        for (int k = Astart[iCol]; k < Astart[iCol + 1]; k++) {
            value += ep_array[Aindex[k]] * Avalue[k];
        }
        if (fabs(value) > HSOL_CONST_TINY) {
            ap_array[iCol] = value;
            ap_index[ap_count++] = iCol;
        }
    }
    row_ap.count = ap_count;
    // for (int i = 0; i < row_ap.array.size(); ++i){
    // 	cout << row_ap.array[i] << endl;
    // 	cout << "\n" << endl;
    // }
    // for (int i = 0; i < numCol; ++i)
    // 	cout << ap_index[i] << endl;

}

void HMatrix::price_by_row(HVector& row_ap, HVector& row_ep) const {
    // Alias
    int ap_count = 0;
    int *ap_index = &row_ap.index[0];
    double *ap_array = &row_ap.array[0];
    const int ep_count = row_ep.count;
    const int *ep_index = &row_ep.index[0];
    const double *ep_array = &row_ep.array[0];

    // Computation
    for (int i = 0; i < ep_count; i++) {
        int iRow = ep_index[i];
        double multi = ep_array[iRow];
        for (int k = ARstart[iRow]; k < AR_Nend[iRow]; k++) {
            int index = ARindex[k];
            double value0 = ap_array[index];
            double value1 = value0 + multi * ARvalue[k];
            if (value0 == 0)
                ap_index[ap_count++] = index;
            ap_array[index] =
                    (fabs(value1) < HSOL_CONST_TINY) ? HSOL_CONST_ZERO : value1;
        }
    }

    // Try to remove cancellation
    const int apcount1 = ap_count;
    ap_count = 0;
    for (int i = 0; i < apcount1; i++) {
        const int index = ap_index[i];
        const double value = ap_array[index];
        if (fabs(value) > HSOL_CONST_TINY) {
            ap_index[ap_count++] = index;
        } else {
            ap_array[index] = 0;
        }
    }
    row_ap.count = ap_count;
}

void HMatrix::compute_vecT_matB(const double *vec, const int *base,
        HVector *result) {
    result->clear();
    int resultCount = 0;
    int *resultIndex = &result->index[0];
    double *resultArray = &result->array[0];
    for (int i = 0; i < numRow; i++) {
        int iCol = base[i];
        double value = 0;
        if (iCol < numCol) {
            for (int k = Astart[iCol]; k < Astart[iCol + 1]; k++)
                value += vec[Aindex[k]] * Avalue[k];
        } else {
            value = vec[iCol - numCol];
        }
        if (fabs(value) > HSOL_CONST_TINY) {
            resultArray[i] = value;
            resultIndex[resultCount++] = i;
        }
    }
    result->count = resultCount;
}

void HMatrix::compute_matB_vec(const double *vec, const int *base,
        HVector *result) {
    result->clear();
    int resultCount = 0;
    int *resultIndex = &result->index[0];
    double *resultArray = &result->array[0];

    for (int i = 0; i < numRow; i++) {
        int iCol = base[i];
        double value = vec[i];
        if (fabs(value) > HSOL_CONST_TINY) {
            if (iCol < numCol) {
                for (int k = Astart[iCol]; k < Astart[iCol + 1]; k++)
                    resultArray[Aindex[k]] += value * Avalue[k];
            } else {
                resultArray[iCol - numCol] += value;
            }
        }
    }

    for (int i = 0; i < numRow; i++) {
        if (fabs(resultArray[i]) > HSOL_CONST_TINY) {
            resultIndex[resultCount++] = i;
        } else {
            resultArray[i] = 0;
        }
    }
    result->count = resultCount;
}
