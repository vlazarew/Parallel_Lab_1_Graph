#include "NodeOfGraph.h"
#include <sstream>
#include <vector>
#include <iostream>
using namespace std;

string createIA(vector<NodeOfGraph> *graph)
{
    string result = "[";
    int countOfLinks = 0;

    result += to_string(countOfLinks);

    for (NodeOfGraph &node : *graph) {
        result += ",";
        countOfLinks += node.getNeightbors().size();
        result += to_string(countOfLinks);
    }

    return result + "]";
}

string createJA(vector<NodeOfGraph> *graph)
{
    string result = "[";

    for (int j = 0; j < graph->size(); j++) {
        NodeOfGraph node = graph->at(j);
        bool lastNode = (j == graph->size() - 1);

        vector<int> neighbors = node.getNeightbors();
        for (int i = 0; i < neighbors.size(); i++) {
            int neighborID = neighbors[i];

            if ((i == neighbors.size() - 1) && lastNode) {
                result += to_string(neighborID);
            } else {
                result += to_string(neighborID) + ", ";
            }
        }
    }

    return result + "]";
}

vector<string> createAdjacencyList(vector<NodeOfGraph> *graph)
{
    vector<string> result;

    for (NodeOfGraph &node : *graph) {
        stringstream resultString;
        resultString << "[" << node.getID() << "] Ц> {";

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

void printResult(vector<NodeOfGraph> *graph, double seconds)
{
    cout << "N (–азмер матрицы) Ц " << graph->size() << " вершин" << endl;

    string IA = createIA(graph);
    cout << "IA: " << IA << endl;

    string JA = createJA(graph);
    cout << "JA: " << JA << endl;

    vector<string> adjacencyList = createAdjacencyList(graph);
    cout << "—писок смежности: " << endl;
    for (string adjancecyString : adjacencyList) {
        cout << adjancecyString << endl;
    }

    cout << "¬рем€ выполнени€: " << seconds << " с." << endl;
    cout << "¬рем€ на 1 узел: " << seconds / graph->size() << " с." << endl;
}
