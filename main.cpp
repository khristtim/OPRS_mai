#include "tests.h"
#include <iostream>
#include <iomanip>
#include <limits>
#include <string>

// ---------- helpers: safe input ----------
static void clearBadInput() {
    std::cin.clear();
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
}

static int readInt(const std::string& prompt, int minVal = std::numeric_limits<int>::min()) {
    while (true) {
        std::cout << prompt;
        int x;
        if (std::cin >> x) {
            if (x >= minVal) return x;
            std::cout << "Value must be >= " << minVal << "\n";
        } else {
            std::cout << "Bad input. Try again.\n";
            clearBadInput();
        }
    }
}

static ld readLD(const std::string& prompt) {
    while (true) {
        std::cout << prompt;
        ld x;
        if (std::cin >> x) return x;
        std::cout << "Bad input. Try again.\n";
        clearBadInput();
    }
}

static TVector readVector(int n, const std::string& name = "v") {
    TVector v(n);
    std::cout << "Enter vector " << name << " elements (" << n << " numbers):\n";
    for (int i = 0; i < n; ++i) {
        v[static_cast<size_t>(i)] = readLD("  " + name + "[" + std::to_string(i) + "] = ");
    }
    return v;
}

static TMatrix readMatrix(int n, int m, const std::string& name = "A") {
    TMatrix A(n, m);
    std::cout << "Enter matrix " << name << " elements (" << n << "x" << m << "), row by row:\n";
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            A(i, j) = readLD("  " + name + "(" + std::to_string(i) + "," + std::to_string(j) + ") = ");
        }
    }
    return A;
}

static TSymmetricMatrix readSymmetricMatrix(int n, const std::string& name = "S") {
    TSymmetricMatrix S(n);
    std::cout << "Enter symmetric matrix " << name << " (" << n << "x" << n << ")\n";
    std::cout << "Tip: you can enter all elements; program will check symmetry.\n";
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            S(i, j) = readLD("  " + name + "(" + std::to_string(i) + "," + std::to_string(j) + ") = ");
        }
    }
    return S;
}

static ld degToRad(ld deg) {
    return deg * acosl(-1.0L) / 180.0L;
}

static TQuaternion readQuaternion(const std::string& name = "Q") {
    std::cout << "Enter quaternion " << name << " components (q0 q1 q2 q3):\n";
    ld q0 = readLD("  q0 = ");
    ld q1 = readLD("  q1 = ");
    ld q2 = readLD("  q2 = ");
    ld q3 = readLD("  q3 = ");
    return TQuaternion(q0, q1, q2, q3);
}


// ---------- menu ----------
static void printMenu() {
    std::cout << "\n=== Linear Algebra Console ===\n"
              << "Vectors:\n"
              << "  1) v + w\n"
              << "  2) v - w\n"
              << "  3) k * v\n"
              << "  4) v * k\n"
              << "  5) dot(v, w)\n"
              << "  6) cross(v, w)  (only 3D)\n"
              << "  7) |v|\n"
              << "  8) normalize(v)\n"
              << "\nMatrices:\n"
              << "  9) A + B\n"
              << " 10) A - B\n"
              << " 11) k * A\n"
              << " 12) A * k\n"
              << " 13) A * B\n"
              << " 14) A * v\n"
              << " 15) v * A\n"
              << " 16) transpose(A)\n"
              << " 17) det(A) (square)\n"
              << " 18) inv(A) by Gauss  (square)\n"
              << "\nSymmetric matrices:\n"
              << " 19) check symmetry (S)\n"
              << " 20) inv(S) by Cholesky (SPD)\n"
              << "\nQuaternions & rotations:\n"
              << " 21) build Q from axis e + angle (deg)\n"
              << " 22) Q1 + Q2\n"
              << " 23) Q1 * Q2\n"
              << " 24) conjugate(Q)\n"
              << " 25) norm(Q) and normalized(Q)\n"
              << " 26) rotate v by Rodrigues formula (axis+angle)\n"
              << " 27) rotate v by quaternion Q (q0..q3)\n"
              << " 28) rotate v by quaternion built from axis+angle\n"

              << "\nOther:\n"
              << "  0) Exit\n";
}

static void runConsole() {
    while (true) {
        printMenu();
        int op = readInt("Select operation: ", 0);

        try {
            switch (op) {
                case 0:
                    std::cout << "Bye!\n";
                    return;

                // ----- vectors -----
                case 1: { // v+w
                    int n = readInt("Vector size n: ", 1);
                    auto v = readVector(n, "v");
                    auto w = readVector(n, "w");
                    std::cout << "v + w = " << (v + w) << "\n";
                    break;
                }
                case 2: { // v-w
                    int n = readInt("Vector size n: ", 1);
                    auto v = readVector(n, "v");
                    auto w = readVector(n, "w");
                    std::cout << "v - w = " << (v - w) << "\n";
                    break;
                }
                case 3: { // k*v
                    int n = readInt("Vector size n: ", 1);
                    ld k = readLD("k: ");
                    auto v = readVector(n, "v");
                    std::cout << "k * v = " << (k * v) << "\n";
                    break;
                }
                case 4: { // v*k
                    int n = readInt("Vector size n: ", 1);
                    auto v = readVector(n, "v");
                    ld k = readLD("k: ");
                    std::cout << "v * k = " << (v * k) << "\n";
                    break;
                }
                case 5: { // dot
                    int n = readInt("Vector size n: ", 1);
                    auto v = readVector(n, "v");
                    auto w = readVector(n, "w");
                    std::cout << "dot(v, w) = " << static_cast<double>(v * w) << "\n";
                    break;
                }
                case 6: { // cross
                    std::cout << "Cross product is only for 3D.\n";
                    auto v = readVector(3, "v");
                    auto w = readVector(3, "w");
                    std::cout << "v x w = " << (v ^ w) << "\n";
                    break;
                }
                case 7: { // length
                    int n = readInt("Vector size n: ", 1);
                    auto v = readVector(n, "v");
                    std::cout << "|v| = " << static_cast<double>(v.length()) << "\n";
                    break;
                }
                case 8: { // normalize
                    int n = readInt("Vector size n: ", 1);
                    auto v = readVector(n, "v");
                    v.norm();
                    std::cout << "normalized v = " << v << "\n";
                    break;
                }

                // ----- matrices -----
                case 9: { // A+B
                    int n = readInt("Rows n: ", 1);
                    int m = readInt("Cols m: ", 1);
                    auto A = readMatrix(n, m, "A");
                    auto B = readMatrix(n, m, "B");
                    std::cout << "A + B =\n" << (A + B) << "\n";
                    break;
                }
                case 10: { // A-B
                    int n = readInt("Rows n: ", 1);
                    int m = readInt("Cols m: ", 1);
                    auto A = readMatrix(n, m, "A");
                    auto B = readMatrix(n, m, "B");
                    std::cout << "A - B =\n" << (A - B) << "\n";
                    break;
                }
                case 11: { // k*A
                    int n = readInt("Rows n: ", 1);
                    int m = readInt("Cols m: ", 1);
                    ld k = readLD("k: ");
                    auto A = readMatrix(n, m, "A");
                    std::cout << "k * A =\n" << (k * A) << "\n";
                    break;
                }
                case 12: { // A*k
                    int n = readInt("Rows n: ", 1);
                    int m = readInt("Cols m: ", 1);
                    auto A = readMatrix(n, m, "A");
                    ld k = readLD("k: ");
                    std::cout << "A * k =\n" << (A * k) << "\n";
                    break;
                }
                case 13: { // A*B
                    int n = readInt("A rows n: ", 1);
                    int m = readInt("A cols m: ", 1);
                    int p = readInt("B cols p: ", 1);
                    std::cout << "B rows must be m = " << m << "\n";
                    auto A = readMatrix(n, m, "A");
                    auto B = readMatrix(m, p, "B");
                    std::cout << "A * B =\n" << (A * B) << "\n";
                    break;
                }
                case 14: { // A*v
                    int n = readInt("A rows n: ", 1);
                    int m = readInt("A cols m: ", 1);
                    auto A = readMatrix(n, m, "A");
                    std::cout << "Vector size must be m = " << m << "\n";
                    auto v = readVector(m, "v");
                    std::cout << "A * v = " << (A * v) << "\n";
                    break;
                }
                case 15: { // v*A
                    int n = readInt("A rows n: ", 1);
                    int m = readInt("A cols m: ", 1);
                    auto A = readMatrix(n, m, "A");
                    std::cout << "Vector size must be n = " << n << " (row-vector)\n";
                    auto v = readVector(n, "v");
                    std::cout << "v * A = " << (v * A) << "\n";
                    break;
                }
                case 16: { // transpose
                    int n = readInt("Rows n: ", 1);
                    int m = readInt("Cols m: ", 1);
                    auto A = readMatrix(n, m, "A");
                    std::cout << "transpose(A) =\n" << A.t() << "\n";
                    break;
                }
                case 17: { // det
                    int n = readInt("Square size n: ", 1);
                    auto A = readMatrix(n, n, "A");
                    std::cout << "det(A) = " << static_cast<double>(A.det()) << "\n";
                    break;
                }
                case 18: { // inverse by Gauss
                    int n = readInt("Square size n: ", 1);
                    auto A = readMatrix(n, n, "A");
                    auto Inv = !A;
                    std::cout << "inv(A) =\n" << Inv << "\n";
                    break;
                }

                // ----- symmetric -----
                case 19: { // check symmetry
                    int n = readInt("Square size n: ", 1);
                    auto S = readSymmetricMatrix(n, "S");
                    std::cout << "S.isSymmetric() = " << (S.isSymmetric() ? "true" : "false") << "\n";
                    break;
                }
                case 20: { // inverse cholesky
                    int n = readInt("Square size n: ", 1);
                    auto S = readSymmetricMatrix(n, "S");
                    std::cout << "Checking symmetry...\n";
                    if (!S.isSymmetric()) {
                        std::cout << "Matrix is NOT symmetric. Cholesky inverse cannot be used.\n";
                        break;
                    }
                    auto Inv = S.inverseCholesky();
                    std::cout << "inv(S) by Cholesky =\n" << static_cast<TMatrix>(Inv) << "\n";
                    break;
                }
                // ----- quaternions & rotations -----
                case 21: { // build Q from axis+angle (deg)
                    std::cout << "Axis e (3D):\n";
                    TVector e = readVector(3, "e");
                    ld deg = readLD("Angle (deg): ");
                    ld phi = degToRad(deg);
                    TQuaternion Q(e, phi);
                    std::cout << "Q = " << Q << "\n";
                    break;
                }
                case 22: { // Q1 + Q2
                    auto Q1 = readQuaternion("Q1");
                    auto Q2 = readQuaternion("Q2");
                    std::cout << "Q1 + Q2 = " << (Q1 + Q2) << "\n";
                    break;
                }
                case 23: { // Q1 * Q2
                    auto Q1 = readQuaternion("Q1");
                    auto Q2 = readQuaternion("Q2");
                    std::cout << "Q1 * Q2 = " << (Q1 * Q2) << "\n";
                    break;
                }
                case 24: { // conjugate(Q)
                    auto Q = readQuaternion("Q");
                    std::cout << "conjugate(Q) = " << Q.conjugate() << "\n";
                    break;
                }
                case 25: { // norm(Q) and normalized(Q)
                    auto Q = readQuaternion("Q");
                    std::cout << "norm(Q) = " << static_cast<double>(Q.norm()) << "\n";
                    std::cout << "normalized(Q) = " << Q.normalized() << "\n";
                    break;
                }
                case 26: { // rotate Rodrigues
                    std::cout << "Vector v (3D):\n";
                    TVector v = readVector(3, "v");
                    std::cout << "Axis e (3D):\n";
                    TVector e = readVector(3, "e");
                    ld deg = readLD("Angle (deg): ");
                    ld phi = degToRad(deg);

                    TVector vr = v.rotateByRodrigFormula(e, (double)phi);
                    std::cout << "v' = " << vr << "\n";
                    break;
                }
                case 27: { // rotate by given quaternion
                    std::cout << "Vector v (3D):\n";
                    TVector v = readVector(3, "v");
                    auto Q = readQuaternion("Q");
                    TVector vr = v.rotateByQuaternion(Q);
                    std::cout << "v' = " << vr << "\n";
                    break;
                }
                case 28: { // rotate by quaternion built from axis+angle
                    std::cout << "Vector v (3D):\n";
                    TVector v = readVector(3, "v");
                    std::cout << "Axis e (3D):\n";
                    TVector e = readVector(3, "e");
                    ld deg = readLD("Angle (deg): ");
                    ld phi = degToRad(deg);

                    TVector vr = v.rotateByQuaternion(e, (double)phi);
                    std::cout << "v' = " << vr << "\n";
                    break;
                }


                default:
                    std::cout << "Unknown operation.\n";
                    break;
            }
        } catch (const std::exception& e) {
            std::cout << "ERROR: " << e.what() << "\n";
        }
    }
}

int main() {
    std::cout << std::fixed << std::setprecision(6);

    // 1) run tests
    try {
        runAllTests();
        std::cout << "All tests OK.\n";
    } catch (const std::exception& e) {
        std::cout << "TESTS FAILED: " << e.what() << "\n";
        return 1;
    }

    // 2) run interactive console
    runConsole();
    return 0;
}
