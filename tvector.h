#pragma once
#include <vector>
#include <iostream>
#include <cmath>

class TVector {
protected:
    // Размерность вектора
    // Элементы вектора
    std::vector<long double> vec_;
public:
    // Конструктор по умолчанию
    TVector() = default;
    TVector(const std::vector<long double>& vector) : vec_(vector)
    {}
    // Конструктор с заданным кол-вом элементов
    TVector(int n)
    {
        vec_.resize(n);
    }
    // Конструктор копий
    TVector(const TVector& rvalue)
    {
        vec_ = rvalue.vec_;
    }
    // Оператор присваивания
    TVector& operator=(const TVector& rvalue)
    {
        vec_ = rvalue.vec_;
        return *this;
    }
    virtual ~TVector() = default;
    // Функция получение кол-ва элементов вектора
    inline size_t size() const
    {
        return vec_.size();
    }
    // Функция получения индекса последнего элемента
    inline size_t high() const
    {
        if(!vec_.empty())
            return vec_.size() - 1;
        std::cout << "Vector is empty" << std::endl;
        return -1;
    }
    // Функция задания кол-ва элементов вектора
    void resize(int n)
    {
        vec_.resize(n);
    }
    // Оператор доступа к элементам вектора
    // inline double& operator[](int i) { return data[i]; } //REALIZE LATER
    // Оператор константного доступа к элементам вектора
    // inline const double& operator[](int i) const { return data[i]; } //REALIZE LATER
    // Оператор - унарный минус
    TVector operator-() const
    {
        TVector result;
        result.vec_.resize(vec_.size());
        for (size_t i = 0; i < vec_.size(); ++i)
        {
            result.vec_[i] = -vec_[i];
        }
        return result;
    }
    // Оператор вычитания векторов
    TVector& operator-=(const TVector& rhs)
    {
        if (vec_.size() != rhs.size())
        {
            throw std::invalid_argument("Different size of vectors");
        }
        for (size_t i = 0; i < vec_.size(); ++i)
        {
            vec_[i] -= rhs.vec_[i];
        }
        return *this;
    }
    friend TVector operator-(TVector lhs, const TVector& rhs)
    {
        lhs -= rhs;
        return lhs;
    }
    TVector& operator+=(const TVector& rhs)
    {
        if (vec_.size() != rhs.size())
        {
            throw std::invalid_argument("Different size of vectors");
        }
        for (size_t i = 0; i < vec_.size(); ++i)
        {
            vec_[i] += rhs.vec_[i];
        }
        return *this;
    }
    friend TVector operator+(TVector lhs, const TVector& rhs)
    {
        lhs += rhs;
        return lhs;
    }
    // Оператор умножения вектора на число
    TVector operator*(const long double& arg) const
    {
        if(vec_.empty())
        {
            throw std::invalid_argument("Vector is empty");
        }
        TVector result;
        result.vec_.resize(vec_.size());
        for(size_t i = 0; i < vec_.size(); ++i)
        {
            result.vec_[i] = vec_[i] * arg;
        }
        return result;
    }
    // Оператор скалярного умножения векторов
    long double operator*(const TVector& rhs) const
    {
        if (vec_.size() != rhs.size())
        {
            throw std::invalid_argument("Different size of vectors");
        }
        long double result = 0.0;
        for(size_t i = 0; i < vec_.size(); ++i)
        {
            result += vec_[i] * rhs.vec_[i];
        }
        return result;
    }


    // Оператор умножения вектора на матрицу
    //    TVector operator * (const TMatrix& arg) const; //REALIZE LATER
    // Оператор умножения вектора на кватернион
    //    TQuaternion operator * (const TQuaternion& arg) const; //REALIZE LATER
    // Оператор векторного умножения векторов
    TVector operator^(const TVector& arg) const
    {
        if (vec_.size() != rhs.size())
        {
            throw std::invalid_argument("Different size of vectors");
        }

    }
    // Дружественная функция - оператор умножения числа на вектор
    friend TVector operator*(double lvalue, const TVector& rvalue);
    // Функция получения модуля вектора
    inline long double norm() const
    {
        long double norm_ = 0.0;
        for(auto i : vec_)
        {
            norm_ += i * i;
        }
        return norm_;
    }

    long double length() const
    {
        return std::sqrt(norm());
    }

    friend std::ostream& operator<<(std::ostream& os, const TVector& vec)
    {
        os << '{';
        bool l_comma = false;
        for(size_t i = 0; i < vec.size(); ++i)
        {
            l_comma = i == vec.high() ? true : false;
            if(l_comma)
            {
                os << vec.vec_[i];
                break;
            }
            os << vec.vec_[i] << ", ";

        }
        os << '}';
        return os;
    }
};
