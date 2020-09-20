#pragma once
#include <list>
using namespace std;

class NodeOfGraph {

public:
    NodeOfGraph();

private:
    int indexOfColumn;
    int indexOfRow;
    list<NodeOfGraph> neighbors;
};

NodeOfGraph::NodeOfGraph()
{
    this->indexOfColumn = 0;
    this->indexOfRow = 0;
    this->neighbors;
}
