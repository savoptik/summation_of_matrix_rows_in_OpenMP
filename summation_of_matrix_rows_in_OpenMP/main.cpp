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
void sumInCP(float *mat, const int r, const int c, std::vector<float>& res); // функция суммирует строки на процессоре
void sumInCPOMP(float *mat, const int r, const int c, std::vector<float>& res);

int main(int argc, const char * argv[]) {
    std::vector<float> matrix, res1, res2; // вектор матрицы и векторы результатов
    matrix.resize(cols * rows); // выделяем память для матрицы
    res1.resize(rows); // выделяем память для первого массива результатов.
    res2.resize(rows); // выделяем память для второго массива результатов
    matrixGeneration(matrix); // генирируем матрицу случайных чисел
    sumInCP(matrix.data(), rows, cols, res1); // выполняем суммирование на процессоре.
    sumInCPOMP(matrix.data(), rows, cols, res2); // выполняем суммирование на openmp
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
            float *sum = new float [16];
            int cr = (i * 100 + k) * 960;
            for (int j = 0; j < 60; j++) {
                for (int l = 0; l < 16; l++) {
                    sum[l] += mat[cr + j * 16 + l];
                }
            }
            res[i * 100 + k] = sum[0] + sum[1] + sum[2] + sum[3] + sum[4] + sum[5] + sum[6] + sum[7] + sum[8] + sum[9] + sum[10] + sum[11] + sum[12] + sum[13] + sum[14] + sum[15];
        }
    }
    double t2 = omp_get_wtime();
    std::cout << res[rand() % res.size()] << " Время параллельно: " << t2-t1 << std::endl; // выводим время в милисекундах
}
