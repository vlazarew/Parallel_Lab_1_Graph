using namespace std;
#include "NodeOfGraph.h"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <omp.h>
#include <sstream>
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

	#pragma omp parallel
	{
		#pragma omp for
		for (int i = 0; i < countOfRows; i++)
		{
			for (int j = 0; j < countOfColumns; j++)
			{
				int nodeId = countOfColumns * i + j;
				NodeOfGraph* node = &NodeOfGraph(j, i, nodeId);
				resultGraph->at(nodeId) = *node;
			}
		}
	}

	int sizeOfResultGraph = resultGraph->size();
	map<int, int> linkedNodes;
	/*#pragma omp parallel shared(sizeOfResultGraph, countOfColumns, countOfRows, resultGraph, IA, JA)
	{
		#pragma omp for*/
	for (int i = 0; i < sizeOfResultGraph; i++)
	{
		NodeOfGraph node = resultGraph->at(i);
		int nodeIndexOfColumn = node.getIndexOfColumn();
		int nodeIndexOfRow = node.getIndexOfRow();

		bool throughLink = false;

		bool rightNodeExist = (nodeIndexOfColumn < countOfColumns - 1);
		bool leftNodeExist = (nodeIndexOfColumn > 0);
		bool downNodeExist = (nodeIndexOfRow < countOfRows - 1);
		bool upNodeExist = (nodeIndexOfRow > 0);

		if (nodesToThroughLink > 0 && rightNodeExist && downNodeExist)
		{
			nodesToThroughLink--;
		}
		else
		{
			throughLink = rightNodeExist && downNodeExist;
		}

		vector<int> nodesNeighbors = node.getNeightbors();
		auto linkedId = linkedNodes.find(i);
		if (linkedNodes.find(i) != linkedNodes.end())
		{
			//#pragma omp critical
			//{
			JA->push_back(linkedId->second);
			linkedNodes.erase(i);
		//}
		}

		// Добавляем верхнюю вершину в соседи
		if (upNodeExist)
		{
			NodeOfGraph* neededNode = &resultGraph->at(i - countOfColumns);
			int nedeedNodeId = neededNode->getID();
			nodesNeighbors.push_back(nedeedNodeId);
			JA->push_back(nedeedNodeId);
		}

		// Добавляем левую вершину в соседи
		if (leftNodeExist)
		{
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
		if (rightNodeExist)
		{
			NodeOfGraph* neededNode = &resultGraph->at(i + 1);
			int nedeedNodeId = neededNode->getID();
			nodesNeighbors.push_back(nedeedNodeId);
			JA->push_back(nedeedNodeId);
		}

		// Добавляем нижнюю вершину в соседи
		if (downNodeExist)
		{
			NodeOfGraph* neededNode = &resultGraph->at(i + countOfColumns);
			int nedeedNodeId = neededNode->getID();
			nodesNeighbors.push_back(nedeedNodeId);
			JA->push_back(nedeedNodeId);
		}

		// Добавление боковых
		if (throughLink)
		{
			NodeOfGraph* throwingNode = &resultGraph->at(i + countOfColumns + 1);
			int throwingNodeId = throwingNode->getID();
			nodesNeighbors.push_back(throwingNodeId);
			JA->push_back(throwingNodeId);

			//#pragma omp critical
			//{
			linkedNodes.insert(pair<int, int>(throwingNodeId, nodeId));
		//}

			vector<int> throwingNodeNeighbors = throwingNode->getNeightbors();
			throwingNodeNeighbors.push_back(nodeId);
			throwingNode->setNeighbors(throwingNodeNeighbors);
			resultGraph->at(throwingNode->getID()) = *throwingNode;
			IA->at(throwingNodeId + 1) = IA->at(throwingNodeId + 1) + 1;

			// Сколько узлов еще надо соединить
			connectNodes--;
			if (connectNodes == 0)
			{
				nodesToThroughLink = k1;
				connectNodes = k2;
			}
		}

		node.setNeighbors(nodesNeighbors);
		resultGraph->at(i) = node;
		int currentIA = IA->at(i) + nodesNeighbors.size();
		IA->at(i + 1) = currentIA;
	}
//}
}

void makeSLAE(vector<int>* IA, vector<int>* JA, vector<double>* A, vector<double>* b, int countOfNodes)
{
	// i - номер строки, j - номер столбца
	#pragma omp parallel
	{
		#pragma omp for
		for (int i = 0; i < countOfNodes; i++)
		{
			double rowSum = 0;
			int diagonalIndex = 0;
			int startIndex = IA->at(i);
			int endIndex = IA->at(i + 1);

			for (int j = startIndex; j < endIndex; ++j)
			{
				int indexOfNode = JA->at(j);
				if (i == indexOfNode)
				{
					diagonalIndex = j;
					continue;
				}
				else
				{
					double value = cos(i * indexOfNode + i + indexOfNode);
					A->at(j) = value;
					rowSum += abs(value);
				}
			}
			A->at(diagonalIndex) = 1.234 * rowSum;
			b->at(i) = sin(i);
		}
	}
}

void printSLAE(vector<double>* A, vector<double>* b, vector<int>* IA)
{
	cout << "Коэффициенты матрицы и правой части: " << endl;

	int countOfLinks = 0;
	int currentIndex = 0;
	for (int i = 1; i < IA->size(); i++)
	{
		countOfLinks = IA->at(i) - IA->at(i - 1);

		stringstream resultString;
		resultString << to_string(i - 1) << " -> [";
		for (int j = 0; j < countOfLinks; j++)
		{
			resultString << fixed << showpoint << setprecision(5) << A->at(currentIndex);
			currentIndex++;
			if (j != countOfLinks - 1)
			{
				resultString << ", ";
			}
		}

		resultString << "]";
		cout << resultString.str() << " = " << fixed << showpoint << setprecision(5) << b->at(i - 1) << endl;
	}
}

int main(int argc, char* argv[])
{

	#pragma region Подготовка данных
	setlocale(LC_ALL, "Russian");

	int Nx, Ny, k1, k2, T;
	bool isPrint = false;
	vector<int> arguments, IA, JA;
	vector<NodeOfGraph> resultGraph;

	// Первый параметр - ссылка на сборку
	for (int i = 1; i < argc; i++)
	{
		try
		{
			arguments.push_back(atoi(argv[i]));
		}
		catch (const std::exception&)
		{
			cout << "Некорректный ввод параметров запуска" << endl;
		}
	}

	if (arguments.empty() || arguments.size() < 4)
	{
		cout << "Введите 4 параметра портрета графа (Nx, Ny, k1, k2, T)!" << endl;
		return 0;
	}

	Ny = arguments[0] + 1;
	Nx = arguments[1] + 1;
	k1 = arguments[2];
	k2 = arguments[3];
	T = arguments[4];

	if (arguments.size() > 5)
	{
		isPrint = (arguments[5] == 1) ? true : false;
	}

	int countOfNodes = Nx * Ny;
	resultGraph.resize(countOfNodes);
	IA.resize(countOfNodes);
	IA.push_back(0);

	omp_set_num_threads(T);
	#pragma endregion

	#pragma region Первый этап
	double start = omp_get_wtime();
	createNodesOfGraph(Nx, Ny, k1, k2, &resultGraph, &IA, &JA);
	double end = omp_get_wtime();
	double seconds = end - start;
	cout << "Всего элементов: " << countOfNodes << " шт." << endl;
	cout << "Время первого этапа: " << seconds << " c." << endl;
	cout << "Время первого этапа на 1 элемент: " << seconds / countOfNodes << " c." << endl
		<< endl;
	#pragma endregion

	#pragma region Второй этап
		// Вектор ненулевых коэффициентов матрицы
	vector<double> A;
	A.resize(JA.size());
	// Вектор правой части
	vector<double> b;
	b.resize(countOfNodes);

	double startSecondStage = omp_get_wtime();
	makeSLAE(&IA, &JA, &A, &b, countOfNodes);
	double endSecondStage = omp_get_wtime();
	double endSecondStagePrint = endSecondStage - startSecondStage;
	cout << "Время второго этапа: " << endSecondStagePrint << " c." << endl;
	cout << "Время второго этапа на 1 элемент: " << endSecondStagePrint / countOfNodes << " c." << endl << endl;

	if (isPrint)
	{
		printResult(&resultGraph, &IA, &JA, seconds);
		printSLAE(&A, &b, &IA);
	}

	#pragma endregion

	return 0;
}
