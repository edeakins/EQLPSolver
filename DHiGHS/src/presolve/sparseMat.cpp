#include <sparseMat.h>

SparseMatrix::SparseMatrix() {
    resize(0,0);
}

SparseMatrix::SparseMatrix(int nrow, int ncol) {
    resize(nrow, ncol);
}

SparseMatrix::SparseMatrix(int nrow, int ncol, int min_capacity) {
    resize(nrow, ncol, min_capacity);
}

void SparseMatrix::reserve(int min_capacity) {
    if (min_capacity > capacity()) {
        rowidx_.resize(min_capacity);
        values_.resize(min_capacity);
    }
}

void SparseMatrix::resize(int nrow, int ncol, int min_capacity) {
    assert(nrow >= 0);
    assert(ncol >= 0);
    assert(min_capacity >= 0);
    nrow_ = nrow;
    colptr_.resize(ncol+1);
    colptr_.shrink_to_fit();
    std::fill(colptr_.begin(), colptr_.end(), 0);
    rowidx_.resize(min_capacity);
    rowidx_.shrink_to_fit();
    values_.resize(min_capacity);
    values_.shrink_to_fit();
}

void SparseMatrix::clear() {
    resize(0,0);
}

void SparseMatrix::LoadFromArrays(int nrow, int ncol, const int* Abegin,
                                  const int* Aend, const int* Ai,
                                  const double* Ax) {
    int nz = 0;
    for (int j = 0; j < ncol; j++)
        nz += Aend[j]-Abegin[j];
    // ncol > nrow ? resize(ncol, ncol, nz) : resize(nrow, nrow, nz);
    resize(nrow, ncol, nz);
    int put = 0;
    for (int j = 0; j < ncol; j++) {
        colptr_[j] = put;
        for (int p = Abegin[j]; p < Aend[j]; p++) {
            if (Ax[p] != 0.0) {
                rowidx_[put] = Ai[p];
                values_[put] = Ax[p];
                put++;
            }
        }
    }
    colptr_[ncol] = put;
    int end = nrow > ncol ? nrow + 1 : 0;
    for (int i = ncol + 1; i < end; ++i)
        colptr_[i] = put;
    SortIndices();
}

void SparseMatrix::SortIndices() {
    if (IsSorted())
        return;
    std::vector<std::pair<int,double>> work(rows());
    for (int j = 0; j < cols(); j++) {
        int nz = 0;             // # entries in column j
        for (int p = begin(j); p < end(j); p++) {
            work[nz].first = index(p);
            work[nz].second = value(p);
            nz++;
        }
        std::sort(work.begin(), work.begin() + nz);
        for (int k = 0, p = begin(j); p < end(j); k++, p++) {
            index(p) = work[k].first;
            value(p) = work[k].second;
        }
    }
}

void SparseMatrix::add_column() {
    int nz = entries();
    int nznew = nz + queue_size();
    reserve(nznew);
    std::copy(rowidx_queue_.begin(), rowidx_queue_.end(), rowidx_.begin() + nz);
    std::copy(values_queue_.begin(), values_queue_.end(), values_.begin() + nz);
    colptr_.push_back(nznew);
    clear_queue();
}

void SparseMatrix::clear_queue() {
    rowidx_queue_.clear();
    values_queue_.clear();
}

bool SparseMatrix::IsSorted() const {
    for (int j = 0; j < cols(); j++) {
        for (int p = begin(j); p < end(j)-1; p++)
            if (index(p) > index(p+1))
                return false;
    }
    return true;
}

SparseMatrix SparseMatrix::Transpose(const SparseMatrix& A) {
    SparseMatrix AT;
    Transpose(A, AT);
    return AT;
}

void SparseMatrix::Transpose(const SparseMatrix& A, SparseMatrix& AT) {
    const int m = A.rows();
    const int n = A.cols();
    const int nz = A.entries();
    AT.resize(n, m, nz);

    // Compute row counts of A in workspace.
    std::vector<int> work(m);
    for (int p = 0; p < nz; p++)
        work[A.index(p)]++;

    // Set column pointers for AT.
    int* ATp = AT.colptr();
    int sum = 0;
    for (int i = 0; i < m; i++) {
        ATp[i] = sum;
        sum += work[i];
        work[i] = ATp[i];
    }
    assert(sum == nz);
    ATp[m] = sum;

    // Fill AT with one column of A at a time.
    // work[i] is the next free slot in column i of AT.
    for (int j = 0; j < n; j++) {
        for (int p = A.begin(j); p < A.end(j); p++) {
            int put = work[A.index(p)]++;
            AT.index(put) = j;
            AT.value(put) = A.value(p);
        }
    }
}

SparseMatrix SparseMatrix::CopyColumns(const SparseMatrix& A, const std::vector<int>& cols) {
    SparseMatrix A2(A.rows(), 0);
    for (int j : cols) {
        for (int p = A.begin(j); p < A.end(j); p++)
            A2.push_back(A.index(p), A.value(p));
        A2.add_column();
    }
    return A2;
}