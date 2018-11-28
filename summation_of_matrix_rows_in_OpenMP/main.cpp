//
//  main.cpp
//  summation_of_matrix_rows_in_OpenMP
//
//  Created by Артём Семёнов on 31/10/2018.
//  Copyright © 2018 Артём Семёнов. All rights reserved.
//

#include <iostream>
#include <vector>
#include <cmath>
#include <omp.h>
#include "mkl.h"

const int rows = 100000, cols = 960; // размеры будущей матрицы

void matrixGeneration(std::vector<float>& vec); // функция генерирует случайную линеризованную матрицу
void sumInCP(float *mat, const int r, const int c, std::vector<float>& res); // функция суммирует строки на процессоре
void sumInCPOMP(float *mat, const int r, const int c, std::vector<float>& res);
void sumInCPOMKL(float *mat, const int r, const int c, float *res);
void sumInCPOMKL2(float *mat, const int r, const int c, float *res);

int main(int argc, const char * argv[]) {
    std::vector<float> matrix, res1, res2, res3, res4; // вектор матрицы и векторы результатов
    matrix.resize(cols * rows); // выделяем память для матрицы
    res1.resize(rows); // выделяем память для первого массива результатов.
    res2.resize(rows); // выделяем память для второго массива результатов
    res3.resize(rows);
    res4.resize(rows);
    matrixGeneration(matrix); // генирируем матрицу случайных чисел
    sumInCP(matrix.data(), rows, cols, res1); // выполняем суммирование на процессоре.
    sumInCPOMP(matrix.data(), rows, cols, res2); // выполняем суммирование на openmp
    sumInCPOMKL(matrix.data(), rows, cols, res3.data());
    sumInCPOMKL2(matrix.data(), rows, cols, res4.data());
    // сверяемся
    long good = 0;
    for (int i = 0; i < res2.size(); i++) {
        good = res2[i] == res1[i] ? good+1: good;
    }
    if (good == rows) {
        std::cout << "Всё хорошо\n";
    } else std::cout << "Что-то пошло не так\n";
    good = 0;
    for (int i = 0; i < res3.size(); i++) {
        good = res3[i] == res1[i] ? good+1: good;
    }
    if (good == rows) {
        std::cout << "Всё хорошо\n";
    } else std::cout << "Что-то пошло не так\n";
    good = 0;
    for (int i = 0; i < res4.size(); i++) {
        good = res4[i] == res1[i] ? good+1: good;
    }
    if (good == rows) {
        std::cout << "Всё хорошо\n";
    } else std::cout << "Что-то пошло не так\n";
    return 0;
}

void matrixGeneration(std::vector<float>& vec) {
    srand(static_cast<int>(100));
    for (int i = 0; i < vec.size(); i++) { // идём по массиву
        float r = rand() % 3 - 1;
        vec[i] = r; // заполняем его случайными числами
    }
}

void sumInCP(float *mat, const int r, const int c, std::vector<float>& res) {
    double t1 = omp_get_wtime();
    for (int i = 0; i < r; i++) {
        float sum = 0;
        int cr = i * c;
        for (int j = 0; j < c; j++) {
            sum += mat[cr + j];
        }
        res[i] = sum;
    }
    double t2 = omp_get_wtime();
    std::cout << "Время последовательно: " << t2-t1 << std::endl; // выводим время в милисекундах
}

void sumInCPOMP(float *mat, const int r, const int c, std::vector<float>& res) {
    srand(static_cast<int>(100));
    double t1 = omp_get_wtime();
#pragma omp parallel for num_threads(16)
    for (int i = 0; i < 1000; i++) { // гнездуем циклы для повышения производительности
        for (int k = 0; k < 100; k++) {
            float sum = 0;
            int cr = (i * 100 + k) * 960;
            for (int j = 0; j < 960; j++) {
                    sum += mat[cr + j];
            }
            res[i * 100 + k] = sum;
        }
    }
    double t2 = omp_get_wtime();
    std::cout << res[rand() % res.size()] << " Время параллельно: " << t2-t1 << std::endl; // выводим время в милисекундах
}

void sumInCPOMKL(float *mat, const int r, const int c, float *res) {
    srand(static_cast<int>(100));
    std::vector<float> ones;
    ones.resize(c);
    for (int i = 0; i < c; i++) {
        ones[i] = 1;
    }
    double t1 = omp_get_wtime();
    int one = 1, zero = 0;
    float onef = 1, zerof = 0;
    sgemv("t", &c, &r, &onef, mat, &c, ones.data(), &one, &zerof, res, &one);
    double t2 = omp_get_wtime();
    std::cout << res[rand() % r] << " Время на МКЛ: " << t2-t1 << std::endl; // выводим время в милисекундах
}

void sumInCPOMKL2(float *mat, const int r, const int c, float *res) {
    srand(static_cast<int>(100));
    std::vector<float> ones;
    ones.resize(c);
    for (int i = 0; i < c; i++) {
        ones[i] = 1;
    }
    int one = 1, zero = 0, part = 100000/250;
    float onef = 1, zerof = 0;
    double t1 = omp_get_wtime();
#pragma omp parallel for num_threads(240)
    for (int i = 0; i < r/part; i++) {
        sgemv("t", &c, &part, &onef, mat + (i * c) * part, &c, ones.data(), &one, &zerof, res + i * part, &one);
    }
    double t2 = omp_get_wtime();
    std::cout << res[rand() % r] << " Время на МКЛ по кускам: " << t2-t1 << std::endl; // выводим время в милисекундах
}

/*
void dgemv_(char* TRANS, const int* M, const int* N,
            double* alpha, double* A, const int* LDA, double* X,
            const int* INCX, double* beta, double* C, const int* INCY);
}
*/
