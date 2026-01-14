#pragma once
#include "tmatrix.h"
#include <vector>
#include <cmath>
#include <stdexcept>

class TSymmetricMatrix : public TMatrix {
public:
    TSymmetricMatrix() = default;
    explicit TSymmetricMatrix(int n, ld value = 0.0L) : TMatrix(n, n, value) {}

    bool isSymmetric(ld eps = 1e-10L) const {
        if (rowCount() != colCount()) return false;
        int n = rowCount();
        for (int i = 0; i < n; ++i)
            for (int j = i + 1; j < n; ++j)
                if (fabsl((*this)(i, j) - (*this)(j, i)) > eps) return false;
        return true;
    }

    TMatrix choleskyL() const {
        if (rowCount() != colCount()) throw std::invalid_argument("Cholesky: matrix must be square");
        if (!isSymmetric()) throw std::runtime_error("Cholesky: matrix is not symmetric");

        int n = rowCount();
        TMatrix L(n, n, 0.0L);

        for (int i = 0; i < n; ++i) {
            for (int j = 0; j <= i; ++j) {
                ld sum = (*this)(i, j);
                for (int k = 0; k < j; ++k) sum -= L(i, k) * L(j, k);

                if (i == j) {
                    if (sum <= EPS) throw std::runtime_error("Cholesky: matrix is not positive definite");
                    L(i, j) = std::sqrt(sum);
                } else {
                    if (fabsl(L(j, j)) < EPS) throw std::runtime_error("Cholesky: zero pivot");
                    L(i, j) = sum / L(j, j);
                }
            }
        }
        return L;
    }

    TSymmetricMatrix inverseCholesky() const {
        int n = rowCount();
        TMatrix L = choleskyL();

        TMatrix X(n, n, 0.0L);
        for (int col = 0; col < n; ++col) {
            // forward: L y = e_col
            std::vector<ld> y(static_cast<size_t>(n), 0.0L);
            for (int i = 0; i < n; ++i) {
                ld sum = (i == col ? 1.0L : 0.0L);
                for (int k = 0; k < i; ++k) sum -= L(i, k) * y[static_cast<size_t>(k)];
                ld diag = L(i, i);
                if (fabsl(diag) < EPS) throw std::runtime_error("Cholesky solve: zero diagonal");
                y[static_cast<size_t>(i)] = sum / diag;
            }

            // backward: L^T x = y
            std::vector<ld> x(static_cast<size_t>(n), 0.0L);
            for (int i = n - 1; i >= 0; --i) {
                ld sum = y[static_cast<size_t>(i)];
                for (int k = i + 1; k < n; ++k) sum -= L(k, i) * x[static_cast<size_t>(k)];
                ld diag = L(i, i);
                if (fabsl(diag) < EPS) throw std::runtime_error("Cholesky solve: zero diagonal");
                x[static_cast<size_t>(i)] = sum / diag;
            }

            for (int i = 0; i < n; ++i) X(i, col) = x[static_cast<size_t>(i)];
        }

        TSymmetricMatrix Inv(n, 0.0L);
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                Inv(i, j) = (X(i, j) + X(j, i)) / 2.0L;

        return Inv;
    }
};
