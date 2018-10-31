//
//  main.cpp
//  summation_of_matrix_rows_in_OpenMP
//
//  Created by Артём Семёнов on 31/10/2018.
//  Copyright © 2018 Артём Семёнов. All rights reserved.
//

#include <iostream>
#include <vector>
#include <random>
#include <chrono>

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
    for (uint i = 0; i < res2.size(); i++) {
        good = res2[i] == res1[i] ? good+1: good;
    }
    if (good == rows) {
        std::cout << "Всё хорошо\n";
    } else std::cout << "Что-то пошло не так\n";
    return 0;
}

void matrixGeneration(std::vector<float>& vec) {
    std::mt19937 gen(static_cast<int>(time(0))); // создаём генератор.
    std::uniform_real_distribution<float> urd(-1.0, 1.0); // задаём диапазон
    for (long i = 0; i < vec.size(); i++) { // идём по массиву
        vec[i] = urd(gen); // заполняем его случайными числами
    }
}

void sumInCP(std::vector<float>& mat, const int r, const int c, std::vector<float>& res) {
    std::chrono::high_resolution_clock::time_point t0 = std::chrono::high_resolution_clock::now(); // фиксируем время начала вычияления
    for (int i = 0; i < r; i++) {
        float sum = 0;
        int cr = i * c;
        for (int j = 0; j < c; j++) {
            sum += mat[cr + j];
        }
        res[i] = sum;
    }
    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now(); // фиксируем время конца вычисления
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count(); // получаем время выполнения в микросекундах
    std::cout << "Время на ЦП: " << (double)duration / 1000 << std::endl; // выводим время в милисекундах
}

void sumInCPOMP(std::vector<float>& mat, const int r, const int c, std::vector<float>& res) {
    std::chrono::high_resolution_clock::time_point t0 = std::chrono::high_resolution_clock::now(); // фиксируем время начала вычияления
    for (int i = 0; i < r; i++) {
        float sum = 0;
        int cr = i * c;
        for (int j = 0; j < c; j++) {
            sum += mat[cr + j];
        }
        res[i] = sum;
    }
    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now(); // фиксируем время конца вычисления
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count(); // получаем время выполнения в микросекундах
    std::cout << "Время на ЦП openmp: " << (double)duration / 1000 << std::endl; // выводим время в милисекундах
}
