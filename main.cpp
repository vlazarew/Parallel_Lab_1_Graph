using namespace std;
#include "NodeOfGraph.h"
#include <iostream>
#include <time.h>
#include <vector>

void createNodesOfGraph(int Nx, int Ny, int k1, int k2, vector<NodeOfGraph>* resultGraph, vector<int>* IA, vector<int>* JA)
{
    // Поскольку узлы
    int countOfRows = Nx;
    int countOfColumns = Ny;
    //
    int loopSize = k1 + k2;

    int nodesToThroughLink = k1;
    int connectNodes = k2;

    for (int i = 0; i < countOfRows; i++) {
        for (int j = 0; j < countOfColumns; j++) {
            int nodeId = countOfColumns * i + j;
            NodeOfGraph* node = &NodeOfGraph(j, i, nodeId);
            resultGraph->at(nodeId) = *node;
        }
    }

    int sizeOfResultGraph = resultGraph->size();
    vector<int> linkedNodes;
    for (int i = 0; i < sizeOfResultGraph; i++) {
        NodeOfGraph node = resultGraph->at(i);
        int nodeIndexOfColumn = node.getIndexOfColumn();
        int nodeIndexOfRow = node.getIndexOfRow();

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

        vector<int> nodesNeighbors = node.getNeightbors();
        for (int index = 0; index < linkedNodes.size(); index++) {
            if (linkedNodes[index] == i) {
                NodeOfGraph* neededNode = &resultGraph->at(i - countOfColumns -1);
                int nedeedNodeId = neededNode->getID();
                JA->push_back(nedeedNodeId);
            }
        }

        // Добавляем верхнюю вершину в соседи
        if (upNodeExist) {
            NodeOfGraph* neededNode = &resultGraph->at(i - countOfColumns);
            int nedeedNodeId = neededNode->getID();
            nodesNeighbors.push_back(nedeedNodeId);
            JA->push_back(nedeedNodeId);
        }

        // Добавляем левую вершину в соседи
        if (leftNodeExist) {
            NodeOfGraph* neededNode = &resultGraph->at(i - 1);
            int nedeedNodeId = neededNode->getID();
            nodesNeighbors.push_back(nedeedNodeId);
            JA->push_back(nedeedNodeId);
        }

        // Добавляем текущую вершину в соседи
        int nodeId = node.getID();
        nodesNeighbors.push_back(node.getID());
        JA->push_back(nodeId);

        // Добавляем правую вершину в соседи
        if (rightNodeExist) {
            NodeOfGraph* neededNode = &resultGraph->at(i + 1);
            int nedeedNodeId = neededNode->getID();
            nodesNeighbors.push_back(nedeedNodeId);
            JA->push_back(nedeedNodeId);
        }

        // Добавляем нижнюю вершину в соседи
        if (downNodeExist) {
            NodeOfGraph* neededNode = &resultGraph->at(i + countOfColumns);
            int nedeedNodeId = neededNode->getID();
            nodesNeighbors.push_back(nedeedNodeId);
            JA->push_back(nedeedNodeId);
        }

        // Добавление боковых
        if (throughLink) {
            NodeOfGraph* throwingNode = &resultGraph->at(i + countOfColumns + 1);
            int throwingNodeId = throwingNode->getID();
            nodesNeighbors.push_back(throwingNodeId);
            JA->push_back(throwingNodeId);
            linkedNodes.push_back(throwingNodeId);

            vector<int> throwingNodeNeighbors = throwingNode->getNeightbors();
            throwingNodeNeighbors.push_back(nodeId);
            throwingNode->setNeighbors(throwingNodeNeighbors);
            resultGraph->at(throwingNode->getID()) = *throwingNode;
            IA->at(throwingNodeId + 1) = IA->at(throwingNodeId + 1) + 1;

            // Сколько узлов еще надо соединить
            connectNodes--;
            if (connectNodes == 0) {
                nodesToThroughLink = k1;
                connectNodes = k2;
            }
        }

        node.setNeighbors(nodesNeighbors);
        resultGraph->at(i) = node;
        int currentIA = IA->at(i) + nodesNeighbors.size();
        IA->at(i + 1) = currentIA;
    }
}

int main(int argc, char* argv[])
{

    setlocale(LC_ALL, "Russian");

    int Nx, Ny, k1, k2;
    bool isPrint = false;
    vector<int> arguments, IA, JA;
    vector<NodeOfGraph> resultGraph;

    // Первый параметр - ссылка на сборку
    for (int i = 1; i < argc; i++) {
        arguments.push_back(atoi(argv[i]));
    }

    if (arguments.empty() || arguments.size() < 4) {
        cout << "Введите 4 параметра портрета графа (Nx, Ny, k1, k2)!" << endl;
        return 0;
    }

    Ny = arguments[0] + 1;
    Nx = arguments[1] + 1;
    k1 = arguments[2];
    k2 = arguments[3];

    if (arguments.size() > 4) {
        isPrint = (arguments[4] == 1) ? true : false;
    }

    int countOfNodes = Nx * Ny;
    resultGraph.resize(countOfNodes);
    IA.resize(countOfNodes);
    //JA.resize(countOfNodes);
    IA.push_back(0);

    clock_t start = clock();
    createNodesOfGraph(Nx, Ny, k1, k2, &resultGraph, &IA, &JA);
    clock_t end = clock();
    double seconds = (double)(end - start) / CLOCKS_PER_SEC;

    if (isPrint) {
        printResult(&resultGraph, &IA, &JA, seconds);
    }

    clock_t endPrint = clock();
    double secondsPrint = (double)(endPrint - start) / CLOCKS_PER_SEC;
    cout << "Время первого этапа с печатью: " << secondsPrint << " c." << endl;

    return 0;
}
