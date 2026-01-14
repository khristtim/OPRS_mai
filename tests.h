#pragma once
#include "operators.h"
#include "tsymmetricMatrix.h"
#include <iostream>
#include <string>
#include <vector>
#include <stdexcept>
#include <algorithm>

struct TestFailure : std::runtime_error {
    using std::runtime_error::runtime_error;
};

inline bool almostEqual(ld a, ld b, ld eps = 1e-9L) {
    ld diff = fabsl(a - b);
    ld scale = std::max<ld>(1.0L, std::max(fabsl(a), fabsl(b)));
    return diff <= eps * scale;
}

inline void requireTrue(bool cond, const std::string& msg) {
    if (!cond) throw TestFailure(msg);
}

inline void requireNear(ld a, ld b, const std::string& msg, ld eps = 1e-9L) {
    if (!almostEqual(a, b, eps)) {
        throw TestFailure(msg + " (got " + std::to_string((double)a) +
                          ", expected " + std::to_string((double)b) + ")");
    }
}

inline void requireVectorNear(const TVector& a, const TVector& b, const std::string& msg, ld eps = 1e-9L) {
    requireTrue(a.size() == b.size(), msg + " (vector size mismatch)");
    for (size_t i = 0; i < a.size(); ++i)
        requireNear(a[i], b[i], msg + " at i=" + std::to_string(i), eps);
}

inline void requireMatrixNear(const TMatrix& A, const TMatrix& B, const std::string& msg, ld eps = 1e-9L) {
    requireTrue(A.rowCount() == B.rowCount() && A.colCount() == B.colCount(), msg + " (matrix size mismatch)");
    for (int i = 0; i < A.rowCount(); ++i)
        for (int j = 0; j < A.colCount(); ++j)
            requireNear(A(i, j), B(i, j), msg + " at (" + std::to_string(i) + "," + std::to_string(j) + ")", eps);
}

// ==== Tests ====
inline void testVectorBasics() {
    TVector a{1, 2, 3};
    TVector b{4, 5, 6};

    requireVectorNear(a + b, TVector{5, 7, 9}, "Vector add");
    requireVectorNear(b - a, TVector{3, 3, 3}, "Vector sub");
    requireNear(a * b, 32.0L, "Dot product");
    requireVectorNear(a ^ b, TVector{-3, 6, -3}, "Cross product");
    requireNear(a.length(), std::sqrt(14.0L), "Vector length");
    requireVectorNear(a * 2.0L, TVector{2, 4, 6}, "Vector*scalar");
    requireVectorNear(2.0L * a, TVector{2, 4, 6}, "scalar*Vector");
}

inline void testMatrixBasics() {
    TMatrix A(2, 3);
    A(0,0)=1; A(0,1)=2; A(0,2)=3;
    A(1,0)=4; A(1,1)=5; A(1,2)=6;

    TMatrix B = A;
    requireMatrixNear(A + B, A * 2.0L, "Matrix add");
    requireMatrixNear(A - B, TMatrix(2,3,0.0L), "Matrix sub");
    TMatrix At = A.t();
    requireTrue(At.rowCount() == 3 && At.colCount() == 2, "Transpose size");
    requireNear(At(2,1), 6.0L, "Transpose content");
}

inline void testMatrixVectorMultiplication() {
    TMatrix A(2, 3);
    A(0,0)=1; A(0,1)=2; A(0,2)=3;
    A(1,0)=4; A(1,1)=5; A(1,2)=6;

    TVector v{10, 20, 30};
    requireVectorNear(A * v, TVector{140, 320}, "Matrix*Vector");

    TVector u{7, 8};
    requireVectorNear(u * A, TVector{39, 54, 69}, "Vector*Matrix");
}

inline void testMatrixMatrixMultiplication() {
    TMatrix A(2, 2);
    A(0,0)=1; A(0,1)=2;
    A(1,0)=3; A(1,1)=4;

    TMatrix B(2, 2);
    B(0,0)=5; B(0,1)=6;
    B(1,0)=7; B(1,1)=8;

    TMatrix C = A * B;
    TMatrix Exp(2,2);
    Exp(0,0)=19; Exp(0,1)=22;
    Exp(1,0)=43; Exp(1,1)=50;

    requireMatrixNear(C, Exp, "Matrix*Matrix");
}

inline void testDetAndInverseGauss() {
    TMatrix A(2,2);
    A(0,0)=4; A(0,1)=7;
    A(1,0)=2; A(1,1)=6;

    requireNear(A.det(), 10.0L, "Determinant");

    TMatrix Inv = !A;
    TMatrix Exp(2,2);
    Exp(0,0)= 0.6L; Exp(0,1)= -0.7L;
    Exp(1,0)= -0.2L; Exp(1,1)= 0.4L;

    requireMatrixNear(Inv, Exp, "Inverse (Gauss)", 1e-9L);
    requireMatrixNear(A * Inv, TMatrix::E(2), "A*inv(A)=I", 1e-9L);

    TMatrix S(2,2);
    S(0,0)=1; S(0,1)=2;
    S(1,0)=2; S(1,1)=4;
    bool thrown = false;
    try { (void)!S; } catch (...) { thrown = true; }
    requireTrue(thrown, "Inverse must throw on singular matrix");
}

inline void testCholeskyInverse() {
    TSymmetricMatrix A(2);
    A(0,0)=4; A(0,1)=2;
    A(1,0)=2; A(1,1)=3;

    requireTrue(A.isSymmetric(), "Symmetry check");
    TSymmetricMatrix Inv = A.inverseCholesky();

    TMatrix Exp(2,2);
    Exp(0,0)=0.375L; Exp(0,1)=-0.25L;
    Exp(1,0)=-0.25L; Exp(1,1)=0.5L;

    requireMatrixNear(Inv, Exp, "Inverse (Cholesky)", 1e-9L);
    requireMatrixNear(static_cast<TMatrix>(A) * static_cast<TMatrix>(Inv), TMatrix::E(2), "A*inv=I", 1e-9L);

    TSymmetricMatrix B(2);
    B(0,0)=1; B(0,1)=2;
    B(1,0)=2; B(1,1)=1;
    bool thrown = false;
    try { (void)B.inverseCholesky(); } catch (...) { thrown = true; }
    requireTrue(thrown, "Cholesky must throw on non-PD");
}

inline void testSymmetryCheck() {
    TSymmetricMatrix S(3);
    S(0,0)=1; S(0,1)=2; S(0,2)=3;
    S(1,0)=2; S(1,1)=4; S(1,2)=5;
    S(2,0)=3; S(2,1)=5; S(2,2)=6;
    requireTrue(S.isSymmetric(), "Symmetric matrix must be symmetric");

    TSymmetricMatrix N(2);
    N(0,0)=1; N(0,1)=2;
    N(1,0)=3; N(1,1)=4;
    requireTrue(!N.isSymmetric(), "Non-symmetric matrix must be NOT symmetric");
}

inline void testQuaternionAndRotations() {
    // Поворот (1,0,0) вокруг оси Z на 90° => (0,1,0)
    TVector v{1,0,0};
    TVector axis{0,0,1};

    const ld deg = 90.0L;
    const ld phi = deg * acosl(-1.0L) / 180.0L;

    // Rodrigues
    TVector r1 = v.rotateByRodrigFormula(axis, (double)phi);
    requireVectorNear(r1, TVector{0,1,0}, "Rodrigues rotation 90deg about Z", 1e-9L);

    // Quaternion from axis-angle
    TQuaternion q(axis, phi);
    TVector r2 = v.rotateByQuaternion(q);
    requireVectorNear(r2, TVector{0,1,0}, "Quaternion rotation 90deg about Z", 1e-9L);

    // Quaternion computed inside method
    TVector r3 = v.rotateByQuaternion(axis, (double)phi);
    requireVectorNear(r3, TVector{0,1,0}, "Quaternion(axis,phi) rotation", 1e-9L);

    // Проверим базовые операции кватернионов
    TQuaternion a(1,2,3,4);
    TQuaternion b(5,6,7,8);

    TQuaternion sum = a + b;
    requireNear(sum[0], 6.0L, "Quat add q0");
    requireNear(sum[1], 8.0L, "Quat add q1");
    requireNear(sum[2],10.0L, "Quat add q2");
    requireNear(sum[3],12.0L, "Quat add q3");

    // conjugate
    TQuaternion ac = a.conjugate();
    requireNear(ac[0], 1.0L, "Quat conjugate q0");
    requireNear(ac[1],-2.0L, "Quat conjugate q1");
    requireNear(ac[2],-3.0L, "Quat conjugate q2");
    requireNear(ac[3],-4.0L, "Quat conjugate q3");

    // norm / normalized
    ld n = a.norm();
    requireTrue(n > 0, "Quat norm positive");
    TQuaternion an = a.normalized();
    requireNear(an.norm(), 1.0L, "Quat normalized norm==1", 1e-9L);
}

inline void runAllTests() {
    struct TestCase { const char* name; void(*fn)(); };
    std::vector<TestCase> tests = {
        {"Vector basics", testVectorBasics},
        {"Matrix basics", testMatrixBasics},
        {"Matrix-Vector multiplication", testMatrixVectorMultiplication},
        {"Matrix-Matrix multiplication", testMatrixMatrixMultiplication},
        {"Det and inverse (Gauss)", testDetAndInverseGauss},
        {"Inverse (Cholesky) for symmetric", testCholeskyInverse},
        {"Symmetry check", testSymmetryCheck},
        {"Quaternion + rotations", testQuaternionAndRotations},
    };

    int passed = 0;
    std::cout << "=== Running tests (" << tests.size() << ") ===\n";
    for (const auto& t : tests) {
        try {
            t.fn();
            std::cout << "[PASS] " << t.name << "\n";
            ++passed;
        } catch (const TestFailure& e) {
            std::cout << "[FAIL] " << t.name << " : " << e.what() << "\n";
        } catch (const std::exception& e) {
            std::cout << "[FAIL] " << t.name << " : unexpected exception: " << e.what() << "\n";
        } catch (...) {
            std::cout << "[FAIL] " << t.name << " : unknown error\n";
        }
    }

    std::cout << "=== Result: " << passed << "/" << tests.size() << " passed ===\n";
    if (passed != static_cast<int>(tests.size())) throw std::runtime_error("Some tests failed");
}


