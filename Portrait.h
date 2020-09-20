#pragma once
#define _CRT_SECURE_NO_WARNINGS
#include <fstream>
#include <iostream>
#include <string>
using namespace std;

class Portrait {
private:
    int countOfNodes;
    int countOfLinks;

    int countOfRows;
    int countOfColumns;

    // Позиция начала списка столбцов данной строки
    int* IA;
    // Для всех строк матрицы хранятся номера столбцов с ненулевыми значениями
    int* JA;
    // Коеффициенты матрицы
    int* A;

    // Кол-во неделящихся клеток
    int k1;
    // Кол-во делящихся клеток
    int k2;

public:
#pragma region Свойства
    int getCountOfNodes()
    {
        return this->countOfNodes;
    }

    void setCountOfNodes(int value)
    {
        this->countOfNodes = value;
    }

    int getCountOfLinks()
    {
        return this->countOfLinks;
    }

    void setCountOfLinks(int value)
    {
        this->countOfLinks = value;
    }

    int getCountOfRows()
    {
        return this->countOfRows;
    }

    void setCountOfRows(int value)
    {
        this->countOfRows = value;
    }

    int getCountOfColumns()
    {
        return this->countOfColumns;
    }

    void setCountOfColumns(int value)
    {
        this->countOfColumns = value;
    }

    int getK1()
    {
        return this->k1;
    }

    void setK1(int value)
    {
        this->k1 = value;
    }

    int getK2()
    {
        return this->k2;
    }

    void setK2(int value)
    {
        this->k2 = value;
    }

    int* getIA()
    {
        return this->IA;
    }

    void setIA(int* value)
    {
        this->IA = value;
    }

    int* getJA()
    {
        return this->JA;
    }

    void setJA(int* value)
    {
        this->JA = value;
    }

    int* getA()
    {
        return this->A;
    }

    void setCountOfNodes(int* value)
    {
        this->A = value;
    }
#pragma endregion

    Portrait()
    {
        this->countOfNodes = 0;
        this->countOfLinks = 0;
        this->IA = new int[this->countOfNodes + 1];
        this->JA = new int[this->countOfLinks * 2];
        this->A = new int[this->countOfLinks * 2];
        this->countOfColumns = 0;
        this->countOfRows = 0;
    }
};

Portrait createPortraitFromFile();

Portrait readPortrainFromStream(ifstream& input);

void makeDefaultPortraitFile();