#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <omp.h>
#include <sstream>
#include <time.h>
#include <vector>

#pragma region Этап 1. Генерация
void createNodesOfGraph(int Nx, int Ny, int k1, int k2, std::map<int, std::vector<int>> resultGraph, std::vector<int>* IA, std::vector<int>* JA, int sizeOfResultGraph)
{
	// Поскольку узлы
	int countOfRows = Nx;
	int countOfColumns = Ny;
	//
	int loopSize = k1 + k2;

	int nodesToThroughLink = k1;
	int connectNodes = k2;

	std::vector<int> throughLinkNodesUp;
	std::vector<int> throughLinkNodesDown;

	for (int i = 0; i < sizeOfResultGraph; i++)
	{
		int nodeIndexOfColumn = i % (countOfRows + 1);
		int nodeIndexOfRow = i / (countOfRows + 1);

		std::vector<int> nodesNeighbors;
		resultGraph.insert(std::pair<int, std::vector<int>>(i, nodesNeighbors));

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
				std::vector<int> nodesNeighbors = resultGraph.at(nodeId);

				// Обработка узлов
				int nodeIndexOfColumn = j;
				int nodeIndexOfRow = i;

				auto throughLinkNodesIdDown = find(throughLinkNodesDown.begin(), throughLinkNodesDown.end(), nodeId);
				bool throughLinkDown = (throughLinkNodesIdDown != throughLinkNodesDown.end());
				auto throughLinkNodesIdUp = find(throughLinkNodesUp.begin(), throughLinkNodesUp.end(), nodeId - countOfColumns - 1);
				bool throughLinkUp = (throughLinkNodesIdUp != throughLinkNodesUp.end());

				bool rightNodeExist = (nodeIndexOfColumn < countOfColumns - 1);
				bool leftNodeExist = (nodeIndexOfColumn > 0);
				bool downNodeExist = (nodeIndexOfRow < countOfRows - 1);
				bool upNodeExist = (nodeIndexOfRow > 0);

				// Добавление боковых сверху
				if (throughLinkUp)
				{
					int throwingNodeId = nodeId - countOfColumns - 1;
					nodesNeighbors.push_back(throwingNodeId);
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
				}
				//

				resultGraph.at(nodeId) = nodesNeighbors;
			}
		}
	}

	// Заполнение портретов
	for (int i = 0; i < sizeOfResultGraph; i++)
	{
		std::vector<int> nodesNeighbors = resultGraph.at(i);
		IA->at(i + 1) = IA->at(i) + nodesNeighbors.size();

		for (int element : nodesNeighbors)
		{
			JA->push_back(element);
		}
	}

}

std::string vectorToString(std::vector<int>* intVector)
{
	std::stringstream result;
	result << "[";

	for (int i = 0; i < intVector->size(); i++)
	{

		if (i == intVector->size() - 1)
		{
			result << std::to_string(intVector->at(i));
		}
		else
		{
			result << std::to_string(intVector->at(i)) << ", ";
		}
	}

	result << "]";

	return result.str();
}

std::vector<std::string> createAdjacencyList(std::map<int, std::vector<int>> graph)
{
	std::vector<std::string> result;

	for (std::pair<int, std::vector<int>> nodePair : graph)
	{
		std::stringstream resultString;
		resultString << "[" << nodePair.first << "] –> {";

		std::vector<int> neighbors = nodePair.second;
		for (int i = 0; i < neighbors.size(); i++)
		{
			int neighborID = neighbors[i];

			if (i == neighbors.size() - 1)
			{
				resultString << std::to_string(neighborID);
			}
			else
			{
				resultString << std::to_string(neighborID) << ", ";
			}
		}

		resultString << "}";
		result.push_back(resultString.str());
	}

	return result;
}

void printResult(std::map<int, std::vector<int>> graph, std::vector<int>* IA, std::vector<int>* JA, double seconds, int countOfNode)
{
	std::cout << "N (Размер матрицы) – " << countOfNode << " вершин" << std::endl;

	std::string IAstring = vectorToString(IA);
	std::cout << "IA: " << IAstring << std::endl;

	std::string JAstring = vectorToString(JA);
	std::cout << "JA: " << JAstring << std::endl;

	std::vector<std::string> adjacencyList = createAdjacencyList(graph);
	std::cout << "Список смежности: " << std::endl;
	for (std::string adjancecyString : adjacencyList)
	{
		std::cout << adjancecyString << std::endl;
	}

	std::cout << "Время выполнения: " << seconds << " с." << std::endl;
	std::cout << "Время на 1 узел: " << seconds / countOfNode << " с." << std::endl;
}

#pragma endregion

#pragma region Этап 2. Построение СЛАУ
void makeSLAE(std::vector<int>* IA, std::vector<int>* JA, std::vector<double>* A, std::vector<double>* b, int countOfNodes)
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

void printSLAE(std::vector<double>* A, std::vector<double>* b, std::vector<int>* IA)
{
	std::cout << "Коэффициенты матрицы и правой части: " << std::endl;

	int countOfLinks = 0;
	int currentIndex = 0;
	for (int i = 1; i < IA->size(); i++)
	{
		countOfLinks = IA->at(i) - IA->at(i - 1);

		std::stringstream resultString;
		resultString << std::to_string(i - 1) << " -> [";
		for (int j = 0; j < countOfLinks; j++)
		{
			resultString << std::fixed << std::showpoint << std::setprecision(5) << A->at(currentIndex);
			currentIndex++;
			if (j != countOfLinks - 1)
			{
				resultString << ", ";
			}
		}

		resultString << "]";
		std::cout << resultString.str() << " = " << std::fixed << std::showpoint << std::setprecision(5) << b->at(i - 1) << std::endl;
	}
}

#pragma endregion

#pragma region Этап 3. Решение СЛАУ
void solveSLAE(std::vector<int>* IA, std::vector<int>* JA, std::vector<double>* A, std::vector<double>* b, int Nx, double tol, std::vector<double>* x, int n, double res)
{

}
#pragma endregion

#pragma region Главные методы
void doTask(int T, int Nx, int Ny, int k1, int k2, std::map<int, std::vector<int>> resultGraph, std::vector<int> IA, std::vector<int> JA, int countOfNodes, int isPrint, double tol)
{
	std::cout << "Количество потоков: " << std::to_string(T) << std::endl << std::endl;

	#pragma region Первый этап
	omp_set_num_threads(T);
	double start = omp_get_wtime();
	createNodesOfGraph(Nx, Ny, k1, k2, resultGraph, &IA, &JA, countOfNodes);
	double end = omp_get_wtime();
	double seconds = end - start;
	std::cout << "Всего элементов: " << countOfNodes << " шт." << std::endl;
	std::cout << "Время первого этапа: " << seconds << " c." << std::endl;
	std::cout << "Время первого этапа на 1 элемент: " << seconds / countOfNodes << " c." << std::endl
		<< std::endl;
	#pragma endregion

	#pragma region Второй этап
	// Вектор ненулевых коэффициентов матрицы
	std::vector<double> A;
	A.resize(JA.size());
	// Вектор правой части
	std::vector<double> b;
	b.resize(countOfNodes);

	double startSecondStage = omp_get_wtime();
	makeSLAE(&IA, &JA, &A, &b, countOfNodes);
	double endSecondStage = omp_get_wtime();
	double endSecondStagePrint = endSecondStage - startSecondStage;
	std::cout << "Время второго этапа: " << endSecondStagePrint << " c." << std::endl;
	std::cout << "Время второго этапа на 1 элемент: " << endSecondStagePrint / countOfNodes << " c." << std::endl << std::endl;

	#pragma region Третий этап
	// Вектор решений
	std::vector<double> x;
	x.resize(Nx);
	// Количество итераций
	int n;
	// L2 норма невзяки
	double res;

	solveSLAE(&IA, &JA, &A, &b, Nx, tol, &x, n, res);
	#pragma endregion

	if (isPrint)
	{
		printResult(resultGraph, &IA, &JA, seconds, countOfNodes);
		printSLAE(&A, &b, &IA);
	}

	std::cout << "//////" << std::endl << std::endl;

	#pragma endregion
}

int main(int argc, char* argv[])
{

	#pragma region Подготовка данных
	setlocale(LC_ALL, "Russian");

	int Nx, Ny, k1, k2, T;
	bool isPrint = false;
	std::vector<int> arguments, IA, JA;
	std::map<int, std::vector<int>> resultGraph;

	// Первый параметр - ссылка на сборку
	for (int i = 1; i < argc; i++)
	{
		try
		{
			arguments.push_back(atoi(argv[i]));
		}
		catch (const std::exception&)
		{
			std::cout << "Некорректный ввод параметров запуска" << std::endl;
		}
	}

	if (arguments.empty() || arguments.size() < 4)
	{
		std::cout << "Введите 4 параметра портрета графа (Nx, Ny, k1, k2, T)!" << std::endl;
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
	IA.resize(countOfNodes);
	IA.push_back(0);
	#pragma endregion

	doTask(1, Nx, Ny, k1, k2, resultGraph, IA, JA, countOfNodes, isPrint);
	doTask(T, Nx, Ny, k1, k2, resultGraph, IA, JA, countOfNodes, isPrint);

	return 0;
}

#pragma endregion
