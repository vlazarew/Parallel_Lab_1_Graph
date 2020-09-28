using namespace std;
#include "NodeOfGraph.h"
#include <iostream>
#include <time.h>
#include <vector>

void handleNode(NodeOfGraph* node, int countOfColumns, int countOfRows, int nodesToThroughLink, vector<NodeOfGraph>* resultGraph, int i, int connectNodes, int k1, int k2)
{
    int nodeIndexOfColumn = node->getIndexOfColumn();
    int nodeIndexOfRow = node->getIndexOfRow();

    bool throughLink = false;

    bool rightNodeExist = (nodeIndexOfColumn < countOfColumns - 1);
    bool leftNodeExist = (nodeIndexOfColumn > 0);
    bool downNodeExist = (nodeIndexOfRow < countOfRows - 1);
    bool upNodeExist = (nodeIndexOfRow > 0);

    if (nodesToThroughLink > 0 && rightNodeExist && downNodeExist) {
        nodesToThroughLink--;
    } else {
        throughLink = rightNodeExist && downNodeExist;
    }

    vector<NodeOfGraph>* nodesNeighbors = &node->getNeightbors();

    // Добавляем верхнюю вершину в соседи
    if (upNodeExist) {
        NodeOfGraph* neededNode = &resultGraph->at(i - countOfColumns);
        nodesNeighbors->push_back(*neededNode);
    }

    // Добавляем левую вершину в соседи
    if (leftNodeExist) {
        NodeOfGraph* neededNode = &resultGraph->at(i - 1);
        nodesNeighbors->push_back(*neededNode);
    }

    // Добавляем текущую вершину в соседи
    nodesNeighbors->push_back(*node);

    // Добавляем правую вершину в соседи
    if (rightNodeExist) {
        NodeOfGraph* neededNode = &resultGraph->at(i + 1);
        nodesNeighbors->push_back(*neededNode);
    }

    // Добавляем нижнюю вершину в соседи
    if (downNodeExist) {
        NodeOfGraph* neededNode = &resultGraph->at(i + countOfColumns);
        nodesNeighbors->push_back(*neededNode);
    }

    // Добавление боковых
    if (throughLink) {
        NodeOfGraph* throwingNode = &resultGraph->at(i + countOfColumns + 1);
        nodesNeighbors->push_back(*throwingNode);

        vector<NodeOfGraph> throwingNodeNeighbors = throwingNode->getNeightbors();
        throwingNodeNeighbors.push_back(*node);
        throwingNode->setNeighbors(throwingNodeNeighbors);
        resultGraph->at(throwingNode->getID()) = *throwingNode;

        // Сколько узлов еще надо соединить
        connectNodes--;
        if (connectNodes == 0) {
            nodesToThroughLink = k1;
            connectNodes = k2;
        }
    }

    node->setNeighbors(*nodesNeighbors);
}

void createNodesOfGraph(int Nx, int Ny, int k1, int k2, vector<NodeOfGraph>* resultGraph)
{
    // Поскольку узлы
    int countOfRows = Nx + 1;
    int countOfColumns = Ny + 1;
    //
    int loopSize = k1 + k2;

    int nodesToThroughLink = k1;
    int connectNodes = k2;

    for (int i = 0; i < countOfRows; i++) {
        for (int j = 0; j < countOfColumns; j++) {
            NodeOfGraph* node = &NodeOfGraph(j, i, countOfRows, countOfColumns);
            resultGraph->push_back(*node);
        }
    }

    int sizeOfResultGraph = resultGraph->size();
    for (int i = 0; i < sizeOfResultGraph; i++) {
        NodeOfGraph node = resultGraph->at(i);
        handleNode(&node, countOfColumns, countOfRows, nodesToThroughLink, resultGraph, i, connectNodes, k1, k2);
        resultGraph->at(i) = node;
    }

}

int main(int argc, char* argv[])
{
    clock_t start = clock();

    setlocale(LC_ALL, "Russian");

    int Nx, Ny, k1, k2;
    bool isPrint = false;
    vector<int> arguments;

    // Первый параметр - ссылка на сборку
    for (int i = 1; i < argc; i++) {
        arguments.push_back(atoi(argv[i]));
    }

    if (arguments.empty() || arguments.size() < 4) {
        cout << "Введите 4 параметра портрета графа (Nx, Ny, k1, k2)!" << endl;
        return 0;
    }

    Nx = arguments[0];
    Ny = arguments[1];
    k1 = arguments[2];
    k2 = arguments[3];

    if (arguments.size() > 4) {
        isPrint = (arguments[4] == 1) ? true : false;
    }

    vector<NodeOfGraph> resultGraph;
    createNodesOfGraph(Nx, Ny, k1, k2, &resultGraph);

    clock_t end = clock();
    double seconds = (double)(end - start) / CLOCKS_PER_SEC;

    if (isPrint) {
        printResult(&resultGraph, seconds);
    }

    return 0;
}
