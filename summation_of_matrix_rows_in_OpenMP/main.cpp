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

const int rows = 100000, cols = 960; // размеры будущей матрицы

void matrixGeneration(std::vector<float>& vec); // функция генерирует случайную линеризованную матрицу
void sumInCP(std::vector<float>& mat, const int r, const int c, std::vector<float>& res); // функция суммирует строки на процессоре
void sumInCPOMP(std::vector<float>& mat, const int r, const int c, std::vector<float>& res);
int main(int argc, const char * argv[]) {
    std::vector<float> matrix, res1, res2; // вектор матрицы и векторы результатов
    matrix.resize(cols * rows); // выделяем память для матрицы
    res1.resize(rows); // выделяем память для первого массива результатов.
    res2.resize(rows); // выделяем память для второго массива результатов
    matrixGeneration(matrix); // генирируем матрицу случайных чисел
    sumInCP(matrix, rows, cols, res1); // выполняем суммирование на процессоре.
    sumInCPOMP(matrix, rows, cols, res2); // выполняем суммирование на openmp
    // сверяемся
    long good = 0;
    for (int i = 0; i < res2.size(); i++) {
        good = res2[i] == res1[i] ? good+1: good;
    }
    if (good == rows) {
        std::cout << "Всё хорошо\n";
    } else std::cout << "Что-то пошло не так\n";
    return 0;
}

void matrixGeneration(std::vector<float>& vec) {
    srand(static_cast<int>(100));
    for (int i = 0; i < vec.size(); i++) { // идём по массиву
        double r = double((rand() % 200)/ 100 - 1);
        vec[i] = r; // заполняем его случайными числами
    }
}

void sumInCP(std::vector<float>& mat, const int r, const int c, std::vector<float>& res) {
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

void sumInCPOMP(std::vector<float>& mat, const int r, const int c, std::vector<float>& res) {
    srand(static_cast<int>(100));
    double t1 = omp_get_wtime();
#pragma omp parallel for
    for (int i = 0; i < 1000; i++) {
        for (int k = 0; k < 100; k++) {
            float sum = 0;
            int cr = (i * 100 + k) * c;
            for (int j = 0; j < c; j++) {
                sum += mat[cr + j];
            }
            res[i * 100 + k] = sum;
        }
    }
    double t2 = omp_get_wtime();
    std::cout << res[rand() % res.size()] << " Время параллельно: " << t2-t1 << std::endl; // выводим время в милисекундах
}
