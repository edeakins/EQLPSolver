#ifndef SPARSE_MATRIX_H_
#define SPARSE_MATRIX_H_

#include <vector>
#include <algorithm>
#include <cassert>

class SparseMatrix {
public:
    SparseMatrix();
    SparseMatrix(int nrow, int ncol);
    SparseMatrix(int nrow, int ncol, int min_capacity);

    int rows() const { return nrow_; }
    int cols() const { return colptr_.size()-1; }
    int entries() const { return colptr_.back(); }

    // # entries in column j
    int entries(int j) const { return end(j)-begin(j); }

    // Maximum # entries that can be stored in the matrix.
    int capacity() const { return rowidx_.size(); }

    // Increases capacity if necessary such that capacity() >= min_capacity.
    // Matrix remains unchanged, pointers are invalidated.
    void reserve(int min_capacity);

    // Changes matrix dimensions. Matrix becomes empty, pointers are
    // invalidated.
    void resize(int nrow, int ncol, int min_capacity = 0);

    // Identical to resize(0,0).
    void clear();

    // Builds matrix from data in compressed column format. The matrix data must
    // be valid (no duplicates, no out of range indices, no inf/nan values);
    // this is not checked. The row indices in the input matrix need not be
    // sorted, but those in the output matrix will be.
    void LoadFromArrays(int nrow, int ncol, const int* Abegin, const int* Aend,
                        const int* Ai, const double* Ax);

    int begin(int j) const { return colptr_[j]; }
    int end(int j) const { return colptr_[j+1]; }

    // Accesses entry at position @pos by value.
    int index(int pos) const { return rowidx_[pos]; }
    double value(int pos) const { return values_[pos]; }

    // Accesses entry at position @pos by reference.
    int& index(int pos) { return rowidx_[pos]; }
    double& value(int pos) { return values_[pos]; }

    // Accesses underlying arrays.
    const int *colptr() const { return colptr_.data(); }
    const int *rowidx() const { return rowidx_.data(); }
    const double *values() const { return values_.data(); }
    int *colptr() { return colptr_.data(); }
    int *rowidx() { return rowidx_.data(); }
    double *values() { return values_.data(); }

    // Stores the entries in each column in increasing order of index.
    void SortIndices();

    // The following methods provide a queue for adding new columns to the
    // matrix. Entries in the queue are not part of the matrix (so do not
    // contribute to entries()).

    // Appends an entry to the end of the queue.
    void push_back(int i, double x) {
        rowidx_queue_.push_back(i);
        values_queue_.push_back(x);
    }

    // Returns # entries in the queue.
    int queue_size() const { return rowidx_queue_.size(); }

    // Accesses entry at position @pos in the queue by value.
    int qindex(int pos) const { return rowidx_queue_[pos]; }
    double qvalue(int pos) const { return values_queue_[pos]; }

    // Accesses entry at position @pos in the queue by reference.
    int& qindex(int pos) { return rowidx_queue_[pos]; }
    double& qvalue(int pos) { return values_queue_[pos]; }

    // Makes new column from queue. The queue becomes empty and cols()
    // increases by 1.
    void add_column();

    // Discards queue.
    void clear_queue();

    // Returns true if row indices are sorted.
    bool IsSorted() const;

    // Builds transpose of matrix.
    SparseMatrix Transpose(const SparseMatrix& A);

    // Resizes @AT as necessary and fills with the tranpose of A.
    void Transpose(const SparseMatrix& A, SparseMatrix& AT);

    // Returns a copy of A[:,cols].
    SparseMatrix CopyColumns(const SparseMatrix& A, const std::vector<int>& cols);

    int nrow_;
    std::vector<int> colptr_;
    std::vector<int> rowidx_;
    std::vector<double> values_;
    std::vector<int> rowidx_queue_;
    std::vector<double> values_queue_;
};
#endif