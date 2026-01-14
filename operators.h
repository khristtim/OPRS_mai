#pragma once
#include "tvector.h"
#include "tmatrix.h"
#include "tquaternion.h"
#include <cmath>

// ===== Matrix * Vector =====
inline TVector TMatrix::operator*(const TVector& v) const {
    if (m_ != static_cast<int>(v.size()))
        throw std::invalid_argument("Matrix*Vector: incompatible sizes");
    TVector r(n_);
    for (int i = 0; i < n_; ++i) {
        ld s = 0.0L;
        for (int j = 0; j < m_; ++j) s += (*this)(i, j) * v[static_cast<size_t>(j)];
        r[static_cast<size_t>(i)] = s;
    }
    return r;
}

// ===== Vector * Matrix =====
inline TVector TVector::operator*(const TMatrix& A) const {
    if (static_cast<int>(size()) != A.rowCount())
        throw std::invalid_argument("Vector*Matrix: incompatible sizes");
    TVector r(A.colCount());
    for (int j = 0; j < A.colCount(); ++j) {
        ld s = 0.0L;
        for (int i = 0; i < A.rowCount(); ++i) s += vec_[static_cast<size_t>(i)] * A(i, j);
        r[static_cast<size_t>(j)] = s;
    }
    return r;
}

// ===== Quaternion ctor from axis+angle =====
inline TQuaternion::TQuaternion(const TVector& e, ld phi) {
    if (e.size() != 3) throw std::invalid_argument("Quaternion(axis,angle): axis must be 3D");
    ld ex = e[0], ey = e[1], ez = e[2];
    ld len = std::sqrt(ex*ex + ey*ey + ez*ez);
    if (len < EPS) throw std::invalid_argument("Quaternion(axis,angle): axis must be non-zero");

    ld kx = ex/len, ky = ey/len, kz = ez/len;
    ld half = phi / 2.0L;
    ld c = std::cos(half);
    ld s = std::sin(half);

    (*this)[0] = c;
    (*this)[1] = kx*s;
    (*this)[2] = ky*s;
    (*this)[3] = kz*s;
}

// ===== TVector rotations =====
inline TVector TVector::rotateByRodrigFormula(const TVector& e, double phi) const {
    if (size() != 3 || e.size() != 3) throw std::invalid_argument("Rodrigues rotation: vectors must be 3D");

    // k = e/|e|
    ld ex=e[0], ey=e[1], ez=e[2];
    ld elen = std::sqrt(ex*ex + ey*ey + ez*ez);
    if (elen < EPS) throw std::invalid_argument("Rodrigues rotation: axis must be non-zero");
    ld kx=ex/elen, ky=ey/elen, kz=ez/elen;

    ld vx=(*this)[0], vy=(*this)[1], vz=(*this)[2];

    ld c = std::cos((ld)phi);
    ld s = std::sin((ld)phi);

    // k x v
    ld cx = ky*vz - kz*vy;
    ld cy = kz*vx - kx*vz;
    ld cz = kx*vy - ky*vx;

    // k · v
    ld dot = kx*vx + ky*vy + kz*vz;

    TVector r(3);
    r[0] = vx*c + cx*s + kx*dot*(1.0L - c);
    r[1] = vy*c + cy*s + ky*dot*(1.0L - c);
    r[2] = vz*c + cz*s + kz*dot*(1.0L - c);
    return r;
}

inline TVector TVector::rotateByQuaternion(const TQuaternion& Q) const {
    if (size() != 3) throw std::invalid_argument("Quaternion rotation: vector must be 3D");

    // Нормируем кватернион поворота
    TQuaternion q = Q.normalized();

    // p = (0, v)
    TQuaternion p(0.0L, (*this)[0], (*this)[1], (*this)[2]);

    // p' = q * p * q_conj
    TQuaternion qp = q * p;
    TQuaternion res = qp * q.conjugate();

    // Векторная часть результата
    TVector r(3);
    r[0] = res[1];
    r[1] = res[2];
    r[2] = res[3];
    return r;
}

inline TVector TVector::rotateByQuaternion(const TVector& e, double phi) const {
    TQuaternion q(e, (ld)phi);
    return rotateByQuaternion(q);
}
