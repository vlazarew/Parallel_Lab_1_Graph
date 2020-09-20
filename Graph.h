#pragma once
#include "NodeOfGraph.h"
using namespace std;

class Graph {
public:
    Graph();

private:
    int countOfNodes;
    int countOfLinks;
    NodeOfGraph* nodes;

    // ������� ������ ������ �������� ������ ������
    int* IA;
    // ��� ���� ����� ������� �������� ������ �������� � ���������� ����������
    int* JA;
    // ������������ �������
    int* A;

    int countOfColums;
    int countOfRows;
};

Graph::Graph()
{
    this->countOfNodes = 0;
    this->countOfLinks = 0;
    this->nodes = new NodeOfGraph[this->countOfNodes];
    this->IA = new int[this->countOfNodes + 1];
    this->JA = new int[this->countOfLinks * 2];
    this->A = new int[this->countOfLinks * 2];
}
