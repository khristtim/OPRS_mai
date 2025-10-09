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
        TVector vec;
        for(auto i : vec_)
        {
            vec.vec_.push_back(-i);
        }
        return vec;
    }
    // Оператор вычитания векторов
    TVector operator-(const TVector& arg) const;
    // Оператор сложения векторов
    TVector operator+(const TVector& arg) const;
    // Оператор умножения вектора на число
    TVector operator*(long double arg) const;
    // Оператор скалярного умножения векторов
    double operator*(const TVector& arg) const;
    // Оператор умножения вектора на матрицу
//    TVector operator * (const TMatrix& arg) const; //REALIZE LATER
    // Оператор умножения вектора на кватернион
//    TQuaternion operator * (const TQuaternion& arg) const; //REALIZE LATER
    // Оператор векторного умножения векторов
    TVector operator ^ (const TVector& arg) const;
    // Дружественная функция - оператор умножения числа на вектор
    friend TVector operator * (double lvalue, const TVector& rvalue);
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
