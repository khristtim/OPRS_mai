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
        if (n_ != m_) throw std::invalid_argument("det: matrix must be square");
        TMatrix A = *this;
        ld detv = 1.0L;
        int sign = 1;

        for (int col = 0; col < n_; ++col) {
            int piv = col;
            ld best = fabsl(A(col, col));
            for (int r = col + 1; r < n_; ++r) {
                ld v = fabsl(A(r, col));
                if (v > best) { best = v; piv = r; }
            }
            if (best < EPS) return 0.0L;

            if (piv != col) { A.swapRows(piv, col); sign = -sign; }

            ld pivot = A(col, col);
            detv *= pivot;

            for (int r = col + 1; r < n_; ++r) {
                ld f = A(r, col) / pivot;
                A(r, col) = 0.0L;
                for (int c = col + 1; c < n_; ++c) A(r, c) -= f * A(col, c);
            }
        }
        return detv * static_cast<ld>(sign);
    }

    TMatrix operator!() const {
        if (n_ != m_) throw std::invalid_argument("inverse: matrix must be square");
        int n = n_;
        TMatrix A = *this;
        TMatrix Inv = TMatrix::E(n);

        for (int col = 0; col < n; ++col) {
            int piv = col;
            ld best = fabsl(A(col, col));
            for (int r = col + 1; r < n; ++r) {
                ld v = fabsl(A(r, col));
                if (v > best) { best = v; piv = r; }
            }
            if (best < EPS) throw std::runtime_error("Gauss inverse: matrix is singular/degenerate");

            if (piv != col) { A.swapRows(piv, col); Inv.swapRows(piv, col); }

            ld pivot = A(col, col);
            for (int j = 0; j < n; ++j) { A(col, j) /= pivot; Inv(col, j) /= pivot; }

            for (int r = 0; r < n; ++r) {
                if (r == col) continue;
                ld f = A(r, col);
                if (fabsl(f) < EPS) continue;
                for (int j = 0; j < n; ++j) {
                    A(r, j) -= f * A(col, j);
                    Inv(r, j) -= f * Inv(col, j);
                }
            }
        }
        return Inv;
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
