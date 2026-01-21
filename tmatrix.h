#pragma once
#include <vector>
#include <iostream>
#include <stdexcept>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include "common.h"

class TVector; // forward

class TMatrix {
protected:
    int n_ = 0, m_ = 0;
    std::vector<ld> data_;

    size_t idx(int i, int j) const {
        return static_cast<size_t>(i) * static_cast<size_t>(m_) + static_cast<size_t>(j);
    }

public:
    TMatrix() = default;
    TMatrix(int n, int m, ld value = 0.0L)
        : n_(n), m_(m), data_(static_cast<size_t>(n) * m, value) {
        if (n < 0 || m < 0) throw std::invalid_argument("Negative size");
    }

    int rowCount() const { return n_; }
    int colCount() const { return m_; }
    int rowHigh()  const { return n_ - 1; }
    int colHigh()  const { return m_ - 1; }

    void resize(int n, int m, ld value = 0.0L) {
        if (n < 0 || m < 0) throw std::invalid_argument("Negative size");
        n_ = n; m_ = m;
        data_.assign(static_cast<size_t>(n) * m, value);
    }

    ld& operator()(int i, int j) {
        if (i < 0 || i >= n_ || j < 0 || j >= m_) throw std::out_of_range("Matrix index out of range");
        return data_[idx(i, j)];
    }
    const ld& operator()(int i, int j) const {
        if (i < 0 || i >= n_ || j < 0 || j >= m_) throw std::out_of_range("Matrix index out of range");
        return data_[idx(i, j)];
    }

    TMatrix operator-() const {
        TMatrix r(n_, m_);
        for (size_t k = 0; k < data_.size(); ++k) r.data_[k] = -data_[k];
        return r;
    }

    TMatrix& operator+=(const TMatrix& rhs) {
        if (n_ != rhs.n_ || m_ != rhs.m_) throw std::invalid_argument("Different size of matrices");
        for (size_t k = 0; k < data_.size(); ++k) data_[k] += rhs.data_[k];
        return *this;
    }
    TMatrix& operator-=(const TMatrix& rhs) {
        if (n_ != rhs.n_ || m_ != rhs.m_) throw std::invalid_argument("Different size of matrices");
        for (size_t k = 0; k < data_.size(); ++k) data_[k] -= rhs.data_[k];
        return *this;
    }

    friend TMatrix operator+(TMatrix lhs, const TMatrix& rhs) { return lhs += rhs; }
    friend TMatrix operator-(TMatrix lhs, const TMatrix& rhs) { return lhs -= rhs; }

    TMatrix operator*(ld k) const {
        TMatrix r(n_, m_);
        for (size_t t = 0; t < data_.size(); ++t) r.data_[t] = data_[t] * k;
        return r;
    }
    friend TMatrix operator*(ld k, const TMatrix& A) { return A * k; }

    TMatrix operator*(const TMatrix& B) const {
        if (m_ != B.n_) throw std::invalid_argument("Matrix multiply: incompatible sizes");
        TMatrix C(n_, B.m_, 0.0L);
        for (int i = 0; i < n_; ++i) {
            for (int k = 0; k < m_; ++k) {
                ld aik = (*this)(i, k);
                for (int j = 0; j < B.m_; ++j) C(i, j) += aik * B(k, j);
            }
        }
        return C;
    }

    // Объявление есть, определение будет в Operators.h
    TVector operator*(const TVector& v) const;

    TMatrix t() const {
        TMatrix R(m_, n_);
        for (int i = 0; i < n_; ++i)
            for (int j = 0; j < m_; ++j)
                R(j, i) = (*this)(i, j);
        return R;
    }

    static TMatrix E(int n) {
        TMatrix I(n, n, 0.0L);
        for (int i = 0; i < n; ++i) I(i, i) = 1.0L;
        return I;
    }

    TMatrix& swapRows(int i, int j) {
        if (i < 0 || i >= n_ || j < 0 || j >= n_) throw std::out_of_range("swapRows out of range");
        if (i == j) return *this;
        for (int col = 0; col < m_; ++col) std::swap((*this)(i, col), (*this)(j, col));
        return *this;
    }

    ld det() const {
        if (n_ != m_) {
            throw std::invalid_argument("det: matrix must be square");
        }

        const int size = n_;
        TMatrix work_matrix = *this;

        ld determinant_value = 1.0L;
        int determinant_sign = 1;

        for (int pivot_col = 0; pivot_col < size; ++pivot_col) {

            int pivot_row = pivot_col;
            ld max_abs_value = fabsl(work_matrix(pivot_col, pivot_col));

            for (int row = pivot_col + 1; row < size; ++row) {
                ld candidate = fabsl(work_matrix(row, pivot_col));
                if (candidate > max_abs_value) {
                    max_abs_value = candidate;
                    pivot_row = row;
                }
            }

            if (max_abs_value < EPS) {
                return 0.0L;
            }

            if (pivot_row != pivot_col) {
                work_matrix.swapRows(pivot_row, pivot_col);
                determinant_sign = -determinant_sign;
            }

            const ld leading_element = work_matrix(pivot_col, pivot_col);
            determinant_value *= leading_element;

            for (int row = pivot_col + 1; row < size; ++row) {
                const ld elimination_factor =
                        work_matrix(row, pivot_col) / leading_element;

                work_matrix(row, pivot_col) = 0.0L;

                for (int col = pivot_col + 1; col < size; ++col) {
                    work_matrix(row, col) -=
                            elimination_factor * work_matrix(pivot_col, col);
                }
            }
        }

        return determinant_value * static_cast<ld>(determinant_sign);
    }


    TMatrix operator!() const {
        if (n_ != m_) {
            throw std::invalid_argument("inverse: matrix must be square");
        }

        const int size = n_;
        TMatrix left_matrix = *this;
        TMatrix right_matrix = TMatrix::E(size);

        for (int pivot_col = 0; pivot_col < size; ++pivot_col) {

            int pivot_row = pivot_col;
            ld max_abs_value = fabsl(left_matrix(pivot_col, pivot_col)); // |a_ii|

            for (int row = pivot_col + 1; row < size; ++row) {
                ld candidate = fabsl(left_matrix(row, pivot_col));
                if (candidate > max_abs_value) {
                    max_abs_value = candidate;
                    pivot_row = row;
                }
            }
            if (max_abs_value < EPS) {
                throw std::runtime_error(
                            "Gauss-Jordan inverse: matrix is degenerate"
                            );
            }
            if (pivot_row != pivot_col) {
                left_matrix.swapRows(pivot_row, pivot_col);
                right_matrix.swapRows(pivot_row, pivot_col);
            }

            ld leading_element = left_matrix(pivot_col, pivot_col);

            for (int col = 0; col < size; ++col) {
                left_matrix(pivot_col, col)  /= leading_element;
                right_matrix(pivot_col, col) /= leading_element;
            }

            for (int row = 0; row < size; ++row) {

                if (row == pivot_col) continue;

                ld elimination_factor = left_matrix(row, pivot_col);

                if (fabsl(elimination_factor) < EPS) continue;

                for (int col = 0; col < size; ++col) {
                    left_matrix(row, col)  -= elimination_factor * left_matrix(pivot_col, col);
                    right_matrix(row, col) -= elimination_factor * right_matrix(pivot_col, col);
                }
            }
        }
        return right_matrix;
    }

    friend std::ostream& operator<<(std::ostream& os, const TMatrix& A) {
        os << "[\n";
        for (int i = 0; i < A.n_; ++i) {
            os << "  ";
            for (int j = 0; j < A.m_; ++j)
                os << std::setw(12) << static_cast<double>(A(i, j)) << " ";
            os << "\n";
        }
        os << "]";
        return os;
    }
};
