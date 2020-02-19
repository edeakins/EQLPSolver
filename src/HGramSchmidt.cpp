
#include "HGramSchmidt.h"
using namespace std;

void HGramSchmidt::gramSchmidt (vector<vector<double> > &activeSet, 
                                       vector<int> &startingBasis, 
                                       vector<int> &activeIdx,
                                       int &aggNumCol,
                                       int &aggNumRow,
                                       int &aggNumTot, 
                                       int &rCnt,
                                       vector<int> &residuals){
    // cout << "aggNumTot: " << aggNumTot << endl;
    // cout << "rCnt: " << rCnt << endl;
    vector<vector<double> > basis;
    vector<int> basicVars;
    vector<double> u;
    double l2Norm = 0;
    vector<vector<double> > q;
    vector<vector<double> > v; 
    vector<double> projection;
    int rIdx = -1;
    int cnt = 1;
    // basis.push_back(activeSet[0]);
    int rank = 1;
    int i, j, k;
    cout << "[ ";
    for (i = 0; i < activeSet.size(); ++i){
        v.push_back(activeSet[i]);
        for (j = 0; j < v[i].size(); ++j){
            cout << v[i][j] << " ";
        }
        cout << ";" << endl;
    }
    cout << " ]" << endl;
    cin.get();
    for (i = 0; i < activeSet.size(); ++i){
        // cout << "original vector" << endl;
        // for (j = 0; j < v[i].size(); ++j){
        //     cout << v[i][j] << endl;
        // }
        // cin.get();
        // cout << "row: " << activeIdx[i] << endl;
        // cin.get();
        l2Norm = norm(v[i], v[i].size());
        q.push_back(vecDivide(v[i], l2Norm));
        for (j = i + 1; j < activeSet.size(); ++j){
            projection = project(v[j], q[i]);
            // cout << "projection: { ";
            // for (k = 0; k < projection.size(); ++k)
            //     cout << projection[k] << " ";
            // cout << "}" << endl;
            // cin.get();
            vecMinus(v[j], projection);
        }
        cout << "column of gramSchmidt: " << cnt << endl;
        cin.get();
        for (int j = 0; j < q[i].size(); ++j)
            cout << q[i][j] << endl;
        cin.get();
        cnt++;
        if (indepCheck(q[i])){
            //basis.push_back(u);
            rank++;
            // cout << "rank: " << rank << endl;
            // cin.get();
            if (activeIdx[i] >= aggNumTot){
                rCnt++;
                rIdx++;
                // cout << "rIdx: " << activeIdx[i] - aggNumTot << endl;
                // cin.get();
                residuals.push_back(rIdx);
            }
        }
        else if (!indepCheck(q[i]) && activeIdx[i] < aggNumTot){
            cout << "linDependent row: " << activeIdx[i] << endl;
            startingBasis.push_back(activeIdx[i]);
        }
        // if(activeIdx[i] == aggNumTot - 1){
        // cout << "rank: " << rank << endl;
        // cout << "rCnt: " << rCnt << endl;
        // cin.get();
        // }
    }
    // cout << "q" << endl;
    // for (int i = 0; i < q.size(); ++i){
    //     cout << "row: " << i << " { ";
    //     for (int j = 0; j < q[i].size(); ++j){
    //         cout << q[i][j] << " ";
    //     }
    //     cout << "}" << endl;
    // }
    // cin.get();
    // if(activeIdx[i] < aggNumTot){
    // cout << "rank: " << rank << endl;
    // cout << "rCnt: " << rCnt << endl;
    // cin.get();
    // }
    // for (int i = 0; i < startingBasis.size(); ++i)
    //     cout << startingBasis[i] << endl;
    cout << "startingBasis" << endl;
    for (int i = 0; i < startingBasis.size(); ++i){
        if (startingBasis[i] >= aggNumCol){
            startingBasis[i] += rCnt;
        }
        cout << startingBasis[i] << endl;
    }
    cout << "rCnt: " << rCnt << endl;
    cin.get();
    aggNumRow += rCnt;
    aggNumCol += rCnt;
    aggNumTot = aggNumCol + aggNumRow;
}

bool HGramSchmidt::indepCheck(vector<double> &u){
    bool condition = false;
    for (int i = 0; i < u.size(); ++i){
        if (fabs(u[i]) > 1e-5){
            condition = true;
        }
        // else{
        //     u[i] = 0;
        // }
    }
    return condition;
}

vector<double> HGramSchmidt::vecDivide(vector<double> &v, double scale){
    vector<double> temp(v.size(), 0.0);
    for (int i = 0; i < v.size(); ++i)
        if (fabs(scale) < 1e-6){
            temp[i] = 0;
        }
        else{
            temp[i] = v[i]/scale;
        }
    return temp;
}

void HGramSchmidt::vecMinus(vector<double> &v, vector<double> &projection){
    for (int i = 0; i < v.size(); ++i){
        v[i] -= projection[i];
    }

}

vector<double> HGramSchmidt::project(vector<double> &v, vector<double> &q){
    int length = v.size();
    vector<double> projection(v.size());
    double scale = 0;
    scale = dotProduct(v, q, v.size());
    double factor = 1.0/scale;
    projection = vecDivide(q, factor);
    return projection;
    
}

double HGramSchmidt::norm (vector<double> &x, int length) {
    int i, length5;
    double a, sum = 0;

    length5 = length % 5;

    for(i = 0; i < length5; i++) {
        sum += x[i] * x[i];
    }
    for(; i < length; i += 5) {
        sum += x[i] * x[i] + x[i + 1] * x[i + 1] + x[i + 2] * x[i + 2]
                           + x[i + 3] * x[i + 3] + x[i + 4] * x[i + 4];
    }

    return sqrt(sum);
}

void HGramSchmidt::vec_copy (double * x, double * y, int length) {
    int i, length5;

    length5 = length % 5;

    for(i = 0; i < length5; i++) {
        y[i] = x[i];
    }
    for(; i < length; i += 5) {
        y[i] = x[i];
        y[i + 1] = x[i + 1];
        y[i + 2] = x[i + 2];
        y[i + 3] = x[i + 3];
        y[i + 4] = x[i + 4];
    }
}

void HGramSchmidt::partialvec_copy (double * x, double * y, int length, int index) {
    int i, length5;

    length5 = length % 5;

    for(i = 0; i < length5; i++) {
        y[i] = x[i + index];
    }
    for(; i < length; i += 5) {
        y[i] = x[i + index];
        y[i + 1] = x[i + index + 1];
        y[i + 2] = x[i + index + 2];
        y[i + 3] = x[i + index + 3];
        y[i + 4] = x[i + index + 4];
    }
}

void HGramSchmidt::scalar_div (double * x, double r, int length, double * y) {
    int i, length5;

    length5 = length % 5;

    for(i = 0; i < length5; i++) {
        y[i] = x[i]/r;
    }
    for(; i < length; i += 5) {
        y[i] = x[i]/r;
        y[i + 1] = x[i + 1]/r;
        y[i + 2] = x[i + 2]/r;
        y[i + 3] = x[i + 3]/r;
        y[i + 4] = x[i + 4]/r;
    }
}

void HGramSchmidt::scalar_sub (double * x, double r, int length, double * y) {
    int i, length5;

    length5 = length % 5;

    for(i = 0; i < length5; i++) {
        y[i] -= r * x[i];
    }
    for(; i < length; i += 5) {
        y[i] -= r * x[i];
        y[i + 1] -= r * x[i + 1];
        y[i + 2] -= r * x[i + 2];
        y[i + 3] -= r * x[i + 3];
        y[i + 4] -= r * x[i + 4];
    }
}

void HGramSchmidt::partialscalar_sub (double * x, double r, int length, 
                                              int index, double * y) 
{
    int i, length5;

    length5 = length % 5;

    for(i = 0; i < length5; i++) {
        y[i + index] -= r * x[i];
    }
    for(; i < length; i += 5) {
        y[i + index] -= r * x[i];
        y[i + index + 1] -= r * x[i + 1];
        y[i + index + 2] -= r * x[i + 2];
        y[i + index + 3] -= r * x[i + 3];
        y[i + index + 4] -= r * x[i + 4];
    }
}

double HGramSchmidt::dotProduct (vector<double> &u, vector<double> &v, int length) {
    int i, length5;
    double sum = 0;

    length5 = length % 5;

    for(i = 0; i < length5; i++) {
        sum += u[i] * v[i];
    }
    for(; i < length; i += 5) {
        sum += u[i] * v[i] + u[i + 1] * v[i + 1] + u[i + 2] * v[i + 2]
                           + u[i + 3] * v[i + 3] + u[i + 4] * v[i + 4];
    }

    return sum;
}

double HGramSchmidt::partialdot_product (double * x, double * y, int length, int index) {
    int i, length5;
    double sum = 0;

    length5 = length % 5;

    for(i = index; i < length5; i++) {
        sum += x[i] * y[i];
    }
    for(; i < length; i += 5) {
        sum += x[i] * y[i] + x[i + 1] * y[i + 1] + x[i + 2] * y[i + 2]
                           + x[i + 3] * y[i + 3] + x[i + 4] * y[i + 4];
    }

    return sum;
}

double HGramSchmidt::subdot_product (double * x, double * y, int length, int index) {
    int i, length5;
    double sum = 0;

    length5 = length % 5;

    for(i = 0; i < length5; i++) {
        sum += x[i + index] * y[i];
    }
    for(; i < length; i += 5) {
        sum += x[i + index] * y[i] + x[i + index + 1] * y[i + 1] 
                                   + x[i + index + 2] * y[i + 2]
                                   + x[i + index + 3] * y[i + 3]
                                   + x[i + index + 4] * y[i + 4];
    }

    return sum;
}

// int main () {
//     int i, j, n, m, q_n, r_m;
//     bool full;
//     double x;

//     /* let user set the dimension of matrix A */
//     std::cout << "Enter the dimension m (where A is a m by n matrix): ";
//     std::cin >> m;
//     std::cout << "Enter the dimension n (where A is a m by n matrix): ";
//     std::cin >> n;

//     if(m != n) {
//         /* check if m < n */
//         if(m < n) {
//             printf("For a successful factorization, this implementation "
//                    "requires n <= m.\nTerminating program.\n");
//             return 0;
//         }
//         /* let user choose either full or thin QR factorization */
//         std::cout << "Enter either 0 to compute a thin QR factorization"
//                   << std::endl;
//         std::cout << "          or 1 to compute a full QR factorization: ";
//         std::cin >> full;
//     }
//     else { // else m == n so full and thin QR factorization are identical */
//         full = 1;
//     }

//     /* set dimensions of matrices Q and R based on full or thin QR */
//     if(full) { // Q is m by m and R is m by n
//         q_n = m;
//         r_m = m;
//     }
//     else { // Q is m by n and R is n by n
//         q_n = n;
//         r_m = n;
//     }

//     /* allocate memory for the matrices A and R */
//     double ** a = new double*[q_n];
//     double ** r = new double*[n];
//     for(i = 0; i < n; i++) {
//         a[i] = new double[m];
//         r[i] = new double[r_m];
//     }
//     for(; i < q_n; i++) {
//         a[i] = new double[m];
//     }

//     /* initialize the values in matrix A (only n columns regardless of
//        thin QR or full QR) */
//     for(i = 0; i < n; i++) {
//         for(j = i; j < m; j++) {
//             a[i][j] = j - i + 1; // this choice of values was arbitrary
//         }
//     }

//     /* print the matrix A before calling gramSchmidt */
//     std::cout << "A = " << std::endl;
//     for(i = 0; i < m; i++) {
//         for(j = 0; j < n; j++) {
//             printf("%9.6lg ", a[j][i]);
//         }
//         std::cout << std::endl;
//     }
//     std::cout << std::endl;

//     /* execute gramSchmidt to compute QR factorization */
//     gramSchmidt(a, r, m, n, full);

//     /* print the matrix Q resulting from gramSchmidt */
//     std::cout << "Q = " << std::endl;
//     for(i = 0; i < m; i++) {
//         for(j = 0; j < q_n; j++) {
//             if(a[j][i] >= 0) {
//                 std::cout << " ";
//             }
//             printf("%9.6lg ", a[j][i]);
//         }
//         std::cout << std::endl;
//     }
//     std::cout << std::endl;

//     /* print the matrix R resulting from gramSchmidt */
//     std::cout << "R = " << std::endl;
//     for(i = 0; i < r_m; i++) {
//         for(j = 0; j < n; j++) {
//             printf("%9.6lg ", r[j][i]);
//         }
//         std::cout << std::endl;
//     }
//     std::cout << std::endl;

//     /* print numerical evidence that columns of Q are orthonormal */
//     printf("Numerical verification that {q_1, ..., q_%i} is an "
//            "orthonormal set:\n", q_n);
//     for(i = 0; i < q_n; i++) {
//         for(j = i; j < q_n; j++) {
//             x = dot_product(a[i], a[j], m);
//             printf("q_%i * q_%i = %lg\n", i + 1, j + 1, x);
//         }
//     }
    
//     /* free memory */
//     for(i = 0; i < n; i++) {
//         delete[] a[i];
//         delete[] r[i];
//     }
//     for(; i < q_n; i++) {
//         delete[] a[i];
//     }
//     delete[] a;  
//     delete[] r;

//     return 0;       // exit main
// }
