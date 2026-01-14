#pragma once
#include <vector>
#include <initializer_list>
#include <iostream>
#include <stdexcept>
#include <cmath>
#include "common.h"

class TMatrix; // forward

class TQuaternion;

class TVector {
protected:
    std::vector<ld> vec_;

public:
    TVector() = default;
    explicit TVector(int n, ld value = 0.0L) : vec_(static_cast<size_t>(n), value) {}
    explicit TVector(const std::vector<ld>& v) : vec_(v) {}
    TVector(std::initializer_list<ld> init) : vec_(init) {}

    size_t size() const { return vec_.size(); }
    size_t high() const {
        if (vec_.empty()) throw std::runtime_error("Vector is empty");
        return vec_.size() - 1;
    }
    void resize(int n, ld value = 0.0L) { vec_.assign(static_cast<size_t>(n), value); }

    ld& operator[](size_t i) { return vec_.at(i); }
    const ld& operator[](size_t i) const { return vec_.at(i); }

    TVector operator-() const {
        TVector r(static_cast<int>(size()));
        for (size_t i = 0; i < size(); ++i) r.vec_[i] = -vec_[i];
        return r;
    }

    TVector& operator+=(const TVector& rhs) {
        if (size() != rhs.size()) throw std::invalid_argument("Different size of vectors");
        for (size_t i = 0; i < size(); ++i) vec_[i] += rhs.vec_[i];
        return *this;
    }
    TVector& operator-=(const TVector& rhs) {
        if (size() != rhs.size()) throw std::invalid_argument("Different size of vectors");
        for (size_t i = 0; i < size(); ++i) vec_[i] -= rhs.vec_[i];
        return *this;
    }

    friend TVector operator+(TVector lhs, const TVector& rhs) { return lhs += rhs; }
    friend TVector operator-(TVector lhs, const TVector& rhs) { return lhs -= rhs; }

    TVector operator*(ld k) const {
        TVector r(static_cast<int>(size()));
        for (size_t i = 0; i < size(); ++i) r.vec_[i] = vec_[i] * k;
        return r;
    }
    friend TVector operator*(ld k, const TVector& v) { return v * k; }

    // dot
    ld operator*(const TVector& rhs) const {
        if (size() != rhs.size()) throw std::invalid_argument("Different size of vectors");
        ld s = 0.0L;
        for (size_t i = 0; i < size(); ++i) s += vec_[i] * rhs.vec_[i];
        return s;
    }

    // cross (3D)
    TVector operator^(const TVector& rhs) const {
        if (size() != rhs.size()) throw std::invalid_argument("Different size of vectors");
        if (size() != 3) throw std::invalid_argument("Cross product is defined only for 3D vectors");
        TVector r(3);
        r[0] = vec_[1] * rhs[2] - vec_[2] * rhs[1];
        r[1] = vec_[2] * rhs[0] - vec_[0] * rhs[2];
        r[2] = vec_[0] * rhs[1] - vec_[1] * rhs[0];
        return r;
    }

    ld norm2() const {
        ld s = 0.0L;
        for (ld x : vec_) s += x * x;
        return s;
    }
    ld length() const { return std::sqrt(norm2()); }

    TVector& norm() {
        ld len = length();
        if (len < EPS) throw std::runtime_error("Cannot normalize zero vector");
        for (auto& x : vec_) x /= len;
        return *this;
    }

    // Объявление есть, определение будет в Operators.h
    TVector operator*(const TMatrix& A) const;

    friend std::ostream& operator<<(std::ostream& os, const TVector& v) {
        os << "{";
        for (size_t i = 0; i < v.size(); ++i) {
            os << static_cast<double>(v.vec_[i]);
            if (i + 1 != v.size()) os << ", ";
        }
        os << "}";
        return os;
    }

    TVector rotateByRodrigFormula(const TVector& e, double phi) const;
    TVector rotateByQuaternion(const TQuaternion& Q) const;
    TVector rotateByQuaternion(const TVector& e, double phi) const;
};
