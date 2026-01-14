#pragma once
#include <cmath>
#include <stdexcept>
#include <initializer_list>
#include <iostream>
#include "common.h"

class TVector; // forward

class TQuaternion {
    // Родрига-Гамильтона: q = (q0, q1, q2, q3)
    // где q0 — скалярная часть, (q1,q2,q3) — векторная часть
    ld q_[4]{0,0,0,0};

public:
    // 1) Конструктор по 4 параметрам Родрига-Гамильтона
    TQuaternion(ld q0, ld q1, ld q2, ld q3) {
        q_[0]=q0; q_[1]=q1; q_[2]=q2; q_[3]=q3;
    }

    // Удобно: default
    TQuaternion() = default;

    // initializer_list (чтобы можно было TQuaternion{1,0,0,0})
    TQuaternion(std::initializer_list<ld> init) {
        if (init.size() != 4) throw std::invalid_argument("Quaternion init list must have 4 elements");
        size_t i=0;
        for (auto v: init) q_[i++] = v;
    }

    // 2) Конструктор по оси вращения e (3D) и углу phi (радианы)
    // q = (cos(phi/2), k*sin(phi/2)), где k = e/|e|
    TQuaternion(const TVector& e, ld phi);

    // 3) Конструктор копии (по умолчанию норм)
    TQuaternion(const TQuaternion&) = default;
    TQuaternion& operator=(const TQuaternion&) = default;

    // Свойство доступа по индексу
    ld& operator[](size_t i) {
        if (i >= 4) throw std::out_of_range("Quaternion index out of range");
        return q_[i];
    }
    const ld& operator[](size_t i) const {
        if (i >= 4) throw std::out_of_range("Quaternion index out of range");
        return q_[i];
    }

    // 4) Сложение
    TQuaternion operator+(const TQuaternion& rhs) const {
        return TQuaternion(q_[0]+rhs.q_[0], q_[1]+rhs.q_[1], q_[2]+rhs.q_[2], q_[3]+rhs.q_[3]);
    }

    // 5) Умножение (гамильтоново произведение)
    TQuaternion operator*(const TQuaternion& r) const {
        const ld a0=q_[0], a1=q_[1], a2=q_[2], a3=q_[3];
        const ld b0=r.q_[0], b1=r.q_[1], b2=r.q_[2], b3=r.q_[3];
        return TQuaternion(
            a0*b0 - a1*b1 - a2*b2 - a3*b3,
            a0*b1 + a1*b0 + a2*b3 - a3*b2,
            a0*b2 - a1*b3 + a2*b0 + a3*b1,
            a0*b3 + a1*b2 - a2*b1 + a3*b0
        );
    }

    // 6) Норма (длина) кватерниона
    ld norm() const {
        return std::sqrt(q_[0]*q_[0] + q_[1]*q_[1] + q_[2]*q_[2] + q_[3]*q_[3]);
    }

    // 7) Нормированный кватернион (возвращает новый)
    TQuaternion normalized() const {
        ld n = norm();
        if (n < EPS) throw std::runtime_error("Cannot normalize zero quaternion");
        return TQuaternion(q_[0]/n, q_[1]/n, q_[2]/n, q_[3]/n);
    }

    // 8) Сопряжённый: q* = (q0, -qv)
    TQuaternion conjugate() const {
        return TQuaternion(q_[0], -q_[1], -q_[2], -q_[3]);
    }

    friend std::ostream& operator<<(std::ostream& os, const TQuaternion& q) {
        os << "{" << (double)q.q_[0] << ", " << (double)q.q_[1] << ", "
           << (double)q.q_[2] << ", " << (double)q.q_[3] << "}";
        return os;
    }
};
