#include "NodeOfGraph.h"
#include <iostream>
#include <sstream>
#include <vector>
using namespace std;

string vectorToString(vector<int>* intVector)
{
    stringstream result;
    result << "[";

    for (int i = 0; i < intVector->size(); i++) {

        if (i == intVector->size() - 1) {
            result << to_string(intVector->at(i));
        } else {
            result << to_string(intVector->at(i)) << ", ";
        }
    }

    result << "]";

    return result.str();
}

vector<string> createAdjacencyList(vector<NodeOfGraph>* graph)
{
    vector<string> result;

    for (NodeOfGraph& node : *graph) {
        stringstream resultString;
        resultString << "[" << node.getID() << "] �> {";

        vector<int> neighbors = node.getNeightbors();
        for (int i = 0; i < neighbors.size(); i++) {
            int neighborID = neighbors[i];

            if (i == neighbors.size() - 1) {
                resultString << to_string(neighborID);
            } else {
                resultString << to_string(neighborID) << ", ";
            }
        }

        resultString << "}";
        result.push_back(resultString.str());
    }

    return result;
}

void printResult(vector<NodeOfGraph>* graph, vector<int>* IA, vector<int>* JA, double seconds)
{
    cout << "N (������ �������) � " << graph->size() << " ������" << endl;

    string IAstring = vectorToString(IA);
    cout << "IA: " << IAstring << endl;

    string JAstring = vectorToString(JA);
    cout << "JA: " << JAstring << endl;

    vector<string> adjacencyList = createAdjacencyList(graph);
    cout << "������ ���������: " << endl;
    for (string adjancecyString : adjacencyList) {
        cout << adjancecyString << endl;
    }

    cout << "����� ����������: " << seconds << " �." << endl;
    cout << "����� �� 1 ����: " << seconds / graph->size() << " �." << endl;
}
