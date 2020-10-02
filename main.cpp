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

	vector<int> tempIA;
	vector<vector<int>> tempJA;
	int sizeOfResultGraph = resultGraph->size();
	tempIA.resize(sizeOfResultGraph);

	vector<int> throughLinkNodesUp;
	vector<int> throughLinkNodesDown;
	for (int i = 0; i < sizeOfResultGraph; i++)
	{
		int nodeIndexOfColumn = i % (countOfRows + 1);
		int nodeIndexOfRow = i / (countOfRows + 1);

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

		if (throughLink)
		{
			throughLinkNodesUp.push_back(i);
			throughLinkNodesDown.push_back(i);

			// Сколько узлов еще надо соединить
			connectNodes--;
			if (connectNodes == 0)
			{
				nodesToThroughLink = k1;
				connectNodes = k2;
			}
		}
	}

	#pragma omp parallel
	{
		#pragma omp for
		for (int i = 0; i < countOfRows; i++)
		{
			for (int j = 0; j < countOfColumns; j++)
			{
				// Создание узла
				int nodeId = countOfColumns * i + j;
				NodeOfGraph* node = &NodeOfGraph(j, i, nodeId);

				// Обработка узлов
				int nodeIndexOfColumn = node->getIndexOfColumn();
				int nodeIndexOfRow = node->getIndexOfRow();

				auto throughLinkNodesIdDown = find(throughLinkNodesDown.begin(), throughLinkNodesDown.end(), nodeId);
				bool throughLinkDown = (throughLinkNodesIdDown != throughLinkNodesDown.end());
				auto throughLinkNodesIdUp = find(throughLinkNodesUp.begin(), throughLinkNodesUp.end(), nodeId - countOfColumns - 1);
				bool throughLinkUp = (throughLinkNodesIdUp != throughLinkNodesUp.end());

				bool rightNodeExist = (nodeIndexOfColumn < countOfColumns - 1);
				bool leftNodeExist = (nodeIndexOfColumn > 0);
				bool downNodeExist = (nodeIndexOfRow < countOfRows - 1);
				bool upNodeExist = (nodeIndexOfRow > 0);

				vector<int> nodesNeighbors = node->getNeightbors();

				// Добавление боковых сверху
				if (throughLinkUp)
				{
					int throwingNodeId = nodeId - countOfColumns - 1;
					nodesNeighbors.push_back(throwingNodeId);
					/*#pragma omp critical
					{
						throughLinkNodesUp.erase(remove(throughLinkNodesUp.begin(), throughLinkNodesUp.end(), nodeId), throughLinkNodesUp.end());
					}*/
				}

				// Добавляем верхнюю вершину в соседи
				if (upNodeExist)
				{
					nodesNeighbors.push_back(nodeId - countOfColumns);
				}

				// Добавляем левую вершину в соседи
				if (leftNodeExist)
				{
					nodesNeighbors.push_back(nodeId - 1);
				}

				// Добавляем текущую вершину в соседи
				nodesNeighbors.push_back(nodeId);

				// Добавляем правую вершину в соседи
				if (rightNodeExist)
				{
					nodesNeighbors.push_back(nodeId + 1);
				}

				// Добавляем нижнюю вершину в соседи
				if (downNodeExist)
				{
					nodesNeighbors.push_back(nodeId + countOfColumns);
				}

				// Добавление боковых снизу
				if (throughLinkDown)
				{
					int throwingNodeId = nodeId + countOfColumns + 1;
					nodesNeighbors.push_back(throwingNodeId);
					/*#pragma omp critical
					{
						throughLinkNodesDown.erase(remove(throughLinkNodesDown.begin(), throughLinkNodesDown.end(), nodeId), throughLinkNodesDown.end());
					}*/
				}
				//

				node->setNeighbors(nodesNeighbors);
				resultGraph->at(nodeId) = *node;

				tempIA.at(nodeId) = nodesNeighbors.size();
			}
		}
	}

	// Заполнение портретов
	for (int i = 0; i < sizeOfResultGraph; i++)
	{
		NodeOfGraph node = resultGraph->at(i);
		vector<int> nodesNeighbors = node.getNeightbors();
		IA->at(i + 1) = IA->at(i) + tempIA.at(i);

		for (int element : nodesNeighbors)
		{
			JA->push_back(element);
		}
	}

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

void doTask(int T, int Nx, int Ny, int k1, int k2, vector<NodeOfGraph> resultGraph, vector<int> IA, vector<int> JA, int countOfNodes, int isPrint)
{
	cout << "Количество потоков: " << to_string(T) << endl << endl;

	#pragma region Первый этап
	omp_set_num_threads(T);
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

	cout << "//////" << endl << endl;

	#pragma endregion
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
	#pragma endregion

	doTask(1, Nx, Ny, k1, k2, resultGraph, IA, JA, countOfNodes, isPrint);
	doTask(T, Nx, Ny, k1, k2, resultGraph, IA, JA, countOfNodes, isPrint);

	return 0;
}
