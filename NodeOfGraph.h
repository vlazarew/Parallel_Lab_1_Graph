#pragma once
#define _CRT_SECURE_NO_WARNINGS
#include "stdio.h"
#include <vector>

class NodeOfGraph {
private:
    int id;
    //int value;
    int indexOfColumn;
    int indexOfRow;
    std::vector<NodeOfGraph> neighbors;

public:
    void setIndexOfColumn(int value) { this->indexOfColumn = value; };
    int getIndexOfColumn() { return this->indexOfColumn; };

    void setIndexOfRow(int value) { this->indexOfRow = value; };
    int getIndexOfRow() { return this->indexOfRow; };

    int getID() { return this->id; };

    /*void setValue(int value) { this->value = value; };
    int getValue() { return this->value; };*/

    void setNeighbors(std::vector<NodeOfGraph> value) { this->neighbors = value; };
    std::vector<NodeOfGraph> getNeightbors() { return this->neighbors; };

    NodeOfGraph()
    {
        this->indexOfColumn = 0;
        this->indexOfRow = 0;
        this->neighbors;
        this->id = 0;
        //this->value = 0;
    }

    NodeOfGraph(int indexOfColumn, int indexOfRow, int countOfRows, int countOfColumns)
    {
        this->indexOfColumn = indexOfColumn;
        this->indexOfRow = indexOfRow;
        this->neighbors;
        this->id = countOfColumns * indexOfRow + indexOfColumn;
        //this->value = id;
    }

};

void printResult(std::vector<NodeOfGraph> *graph, double seconds);