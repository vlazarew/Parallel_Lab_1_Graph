using namespace std;
#include "NodeOfGraph.h"
#include <iostream>
#include <time.h>
#include <vector>

void createNodesOfGraph(int Nx, int Ny, int k1, int k2, vector<NodeOfGraph>* resultGraph)
{
    // ��������� ����
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
        NodeOfGraph* node = &resultGraph->at(i);

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

        vector<int> nodesNeighbors = node->getNeightbors();

        // ��������� ������� ������� � ������
        if (upNodeExist) {
            NodeOfGraph* neededNode = &resultGraph->at(i - countOfColumns);
            nodesNeighbors.push_back(neededNode->getID());
        }

        // ��������� ����� ������� � ������
        if (leftNodeExist) {
            NodeOfGraph* neededNode = &resultGraph->at(i - 1);
            nodesNeighbors.push_back(neededNode->getID());
        }

        // ��������� ������� ������� � ������
        nodesNeighbors.push_back(node->getID());

        // ��������� ������ ������� � ������
        if (rightNodeExist) {
            NodeOfGraph* neededNode = &resultGraph->at(i + 1);
            nodesNeighbors.push_back(neededNode->getID());
        }

        // ��������� ������ ������� � ������
        if (downNodeExist) {
            NodeOfGraph* neededNode = &resultGraph->at(i + countOfColumns);
            nodesNeighbors.push_back(neededNode->getID());
        }

        // ���������� �������
        if (throughLink) {
            NodeOfGraph* throwingNode = &resultGraph->at(i + countOfColumns + 1);
            nodesNeighbors.push_back(throwingNode->getID());

            vector<int> throwingNodeNeighbors = throwingNode->getNeightbors();
            throwingNodeNeighbors.push_back(node->getID());
            throwingNode->setNeighbors(throwingNodeNeighbors);
            resultGraph->at(throwingNode->getID()) = *throwingNode;

            // ������� ����� ��� ���� ���������
            connectNodes--;
            if (connectNodes == 0) {
                nodesToThroughLink = k1;
                connectNodes = k2;
            }
        }

        node->setNeighbors(nodesNeighbors);
        resultGraph->at(i) = *node;
    }

    //return resultGraph;
}

int main(int argc, char* argv[])
{
    clock_t start = clock();

    setlocale(LC_ALL, "Russian");

    int Nx, Ny, k1, k2;
    bool isPrint = false;
    vector<int> arguments;

    // ������ �������� - ������ �� ������
    for (int i = 1; i < argc; i++) {
        arguments.push_back(atoi(argv[i]));
    }

    if (arguments.empty() || arguments.size() < 4) {
        cout << "������� 4 ��������� �������� ����� (Nx, Ny, k1, k2)!" << endl;
        return 0;
    }

    Nx = arguments[0];
    Ny = arguments[1];
    k1 = arguments[2];
    k2 = arguments[3];

    if (arguments.size() > 4) {
        isPrint = (arguments[4] == 1) ? true : false;
    }

    //Portrait portrait = createPortrait(Nx, Ny, k1, k2);
    vector<NodeOfGraph> resultGraph;
    createNodesOfGraph(Nx, Ny, k1, k2, &resultGraph);

    clock_t end = clock();
    double seconds = (double)(end - start) / CLOCKS_PER_SEC;

    if (isPrint) {
        printResult(&resultGraph, seconds);
    }

    return 0;
}
