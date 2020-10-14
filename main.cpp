#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <omp.h>
#include <sstream>
#include <time.h>
#include <vector>
#include <string>

#pragma region ���� 1. ���������
void createNodesOfGraph(int Nx, int Ny, int k1, int k2, std::map<int, std::vector<int>> resultGraph, std::vector<int>* IA, std::vector<int>* JA, int sizeOfResultGraph)
{
	// ��������� ����
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
		int nodeIndexOfColumn = i % (countOfRows);
		int nodeIndexOfRow = i / (countOfRows);

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

			// ������� ����� ��� ���� ���������
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
				// �������� ����
				int nodeId = countOfColumns * i + j;
				std::vector<int> nodesNeighbors = resultGraph.at(nodeId);

				// ��������� �����
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

				// ���������� ������� ������
				if (throughLinkUp)
				{
					int throwingNodeId = nodeId - countOfColumns - 1;
					nodesNeighbors.push_back(throwingNodeId);
				}

				// ��������� ������� ������� � ������
				if (upNodeExist)
				{
					nodesNeighbors.push_back(nodeId - countOfColumns);
				}

				// ��������� ����� ������� � ������
				if (leftNodeExist)
				{
					nodesNeighbors.push_back(nodeId - 1);
				}

				// ��������� ������� ������� � ������
				nodesNeighbors.push_back(nodeId);

				// ��������� ������ ������� � ������
				if (rightNodeExist)
				{
					nodesNeighbors.push_back(nodeId + 1);
				}

				// ��������� ������ ������� � ������
				if (downNodeExist)
				{
					nodesNeighbors.push_back(nodeId + countOfColumns);
				}

				// ���������� ������� �����
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

	// ���������� ���������
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
		resultString << "[" << nodePair.first << "] �> {";

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
	std::cout << "N (������ �������) � " << countOfNode << " ������" << std::endl;

	std::string IAstring = vectorToString(IA);
	std::cout << "IA: " << IAstring << std::endl;

	std::string JAstring = vectorToString(JA);
	std::cout << "JA: " << JAstring << std::endl;

	std::vector<std::string> adjacencyList = createAdjacencyList(graph);
	std::cout << "������ ���������: " << std::endl;
	for (std::string adjancecyString : adjacencyList)
	{
		std::cout << adjancecyString << std::endl;
	}

	std::cout << "����� ����������: " << seconds << " �." << std::endl;
	std::cout << "����� �� 1 ����: " << seconds / countOfNode << " �." << std::endl;
}

#pragma endregion

#pragma region ���� 2. ���������� ����
void makeSLAE(std::vector<int>* IA, std::vector<int>* JA, std::vector<double>* A, std::vector<double>* b, int countOfNodes)
{
// i - ����� ������, j - ����� �������
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
	std::cout << "������������ ������� � ������ �����: " << std::endl;

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

#pragma region ���� 3. ������� ����

double scalar(std::vector<double>* x1, std::vector<double>* x2)
{
	double result = 0;
	int vectorSize = x1->size();

	#pragma omp parallel
	{
		#pragma omp for
		for (int i = 0; i < vectorSize; i++)
		{
			result += x1->at(i) * x2->at(i);
		}
	}

	return result;
}

double normalizeVector(std::vector<double>* x)
{
	double result = 0;

	#pragma omp parallel
	{
		#pragma omp for
		for (int i = 0; i < x->size(); i++)
		{
			result += x->at(i) * x->at(i);
		}
	}

	return std::sqrt(result);
}

std::vector<double> spMV(int countOfNodes, std::vector<int>* IA, std::vector<int>* JA, std::vector<double>* A, std::vector<double>* x)
{
	std::vector<double> result;
	result.resize(countOfNodes);
	if (x->size() == 0)
	{
		x->resize(countOfNodes);
	}

	#pragma omp parallel
	{
		#pragma omp for
		for (int i = 0; i < countOfNodes; i++)
		{
			for (int k = IA->at(i); k < IA->at(i + 1); k++)
			{
				result.at(i) += A->at(k) * x->at(JA->at(k));
			}
		}
	}

	return result;
}

std::vector<double> linearCombination(std::vector<double>* x1, std::vector<double>* x2, double a1, double a2)
{
	std::vector<double> result;
	int vectorSize = x1->size();
	result.resize(vectorSize);


	#pragma omp parallel
	{
		#pragma omp for
		for (int i = 0; i < vectorSize; i++)
		{
			result.at(i) = x1->at(i) * a1 + x2->at(i) * a2;
		}
	}

	return result;
}

void createMMatrixFromA(int countOfNodes, std::vector<int>* IA, std::vector<int>* JA, std::vector<double>* A, std::vector<int>* IAM, std::vector<int>* JAM, std::vector<double>* AM)
{
	IAM->resize(countOfNodes + 1);
	AM->resize(countOfNodes);
	JAM->resize(countOfNodes);
	IAM->at(0) = 0;

	#pragma omp parallel
	{
		#pragma omp for
		for (int i = 0; i < countOfNodes; i++)
		{
			int startIndex = IA->at(i);
			int endIndex = IA->at(i + 1);

			for (int j = startIndex; j < endIndex; ++j)
			{
				int indexOfNode = JA->at(j);
				if (i == indexOfNode)
				{
					AM->at(i) = A->at(j);
					IAM->at(i + 1) = i + 1;
					JAM->at(i) = i;
				}
			}
		}
	}
}

void reverseMMatrix(std::vector<double>* AM)
{
	#pragma omp parallel
	{
		#pragma omp for
		for (int i = 0; i < AM->size(); i++)
		{
			AM->at(i) = 1 / AM->at(i);
		}
	}
}

void solveSLAE(std::vector<int>* IA, std::vector<int>* JA, std::vector<double>* A, std::vector<double>* b, int countOfNodes, double tol, std::vector<std::vector<double>>* x, std::vector<std::vector<double>>* r, int n, double res, std::vector<double>* xRes)
{
	double normB = normalizeVector(b);

	std::vector<std::vector<double>> z;
	std::vector<std::vector<double>> p;
	std::vector<std::vector<double>> q;
	std::vector<double> po;
	std::vector<double> betta;
	std::vector<double> alpha;
	// ���� ��� �������� ������� �������
	z.push_back(std::vector<double> {0});
	p.push_back(std::vector<double> {0});
	q.push_back(std::vector<double> {0});
	po.push_back(0);
	betta.push_back(0);
	alpha.push_back(0);

	std::vector<double> AX0 = spMV(countOfNodes, IA, JA, A, &x->at(0));
	r->at(0) = linearCombination(b, &AX0, 1, -1);

	bool convergence = false;
	int k = 1;
	std::vector<int> IAM;
	std::vector<int> JAM;
	std::vector<double> AM;
	createMMatrixFromA(countOfNodes, IA, JA, A, &IAM, &JAM, &AM);
	reverseMMatrix(&AM);

	do
	{
		// ���� ��������
		if (k + 1 > z.size())
		{
			z.resize(k + 1);
		}
		if (k + 1 > p.size())
		{
			p.resize(k + 1);
		}
		if (k + 1 > q.size())
		{
			q.resize(k + 1);
		}
		if (k + 1 > po.size())
		{
			po.resize(k + 1);
		}
		if (k + 1 > betta.size())
		{
			betta.resize(k + 1);
		}
		if (k + 1 > alpha.size())
		{
			alpha.resize(k + 1);
		}
		//

		/*createMMatrixFromA(countOfNodes, IA, JA, A, &IAM, &JAM, &AM);
		reverseMMatrix(&AM);*/
		z.at(k) = spMV(countOfNodes, &IAM, &JAM, &AM, &r->at(k - 1));
		po.at(k) = scalar(&r->at(k - 1), &z.at(k));

		if (k == 1)
		{
			p.at(k) = z.at(k);
		}
		else
		{
			betta.at(k) = po.at(k) / po.at(k - 1);
			p.at(k) = linearCombination(&z.at(k), &p.at(k - 1), 1, b->at(k));
		}

		q.at(k) = spMV(countOfNodes, IA, JA, A, &p.at(k));
		alpha.at(k) = po.at(k) / scalar(&p.at(k), &q.at(k));
		x->at(k) = linearCombination(&x->at(k - 1), &p.at(k), 1, alpha.at(k));
		r->at(k) = linearCombination(&r->at(k - 1), &q.at(k), 1, -alpha.at(k));

		std::cout << "�������� " << k << " ||b - Ax|| = " << normalizeVector(&r->at(k)) << " po = " << po.at(k) << std::endl;
		if (po.at(k) < tol || k >= 15000)
		{
			convergence = true;
		}
		else
		{
			k++;
		}
	}
	while (!convergence);

	res = normalizeVector(&r->at(k));
	n = k;
	xRes = &x->at(k);
}

void printSolveVector(double res, int n, std::vector<double>* x)
{
	std::cout << "x = [";
	for (int i = 0; i < x->size(); i++)
	{
		if (i == x->size() - 1)
		{
			std::cout << std::fixed << std::showpoint << std::setprecision(5) << x->at(i) << " ]" << std::endl;
		}
		else
		{
			std::cout << std::fixed << std::showpoint << std::setprecision(5) << x->at(i) << ", ";
		}
	}

	std::cout << "����� ������� = " << res << std::endl;
	std::cout << "���������� �������� = " << n << std::endl;
}

#pragma endregion

#pragma region ������� ������
void doTask(int T, int Nx, int Ny, int k1, int k2, std::map<int, std::vector<int>> resultGraph, std::vector<int> IA, std::vector<int> JA, int countOfNodes, int isPrint, double tol)
{
	std::cout << "���������� �������: " << std::to_string(T) << std::endl
		<< std::endl;

	#pragma region ������ ����
	omp_set_num_threads(T);
	double start = omp_get_wtime();
	createNodesOfGraph(Nx, Ny, k1, k2, resultGraph, &IA, &JA, countOfNodes);
	double end = omp_get_wtime();
	double seconds = end - start;
	std::cout << "����� ���������: " << countOfNodes << " ��." << std::endl;
	std::cout << "����� ������� �����: " << seconds << " c." << std::endl;
	std::cout << "����� ������� ����� �� 1 �������: " << seconds / countOfNodes << " c." << std::endl
		<< std::endl;
	#pragma endregion

	#pragma region ������ ����
		// ������ ��������� ������������� �������
	std::vector<double> A;
	A.resize(JA.size());
	// ������ ������ �����
	std::vector<double> b;
	b.resize(countOfNodes);

	double startSecondStage = omp_get_wtime();
	makeSLAE(&IA, &JA, &A, &b, countOfNodes);
	double endSecondStage = omp_get_wtime();
	double endSecondStagePrint = endSecondStage - startSecondStage;
	std::cout << "����� ������� �����: " << endSecondStagePrint << " c." << std::endl;
	std::cout << "����� ������� ����� �� 1 �������: " << endSecondStagePrint / countOfNodes << " c." << std::endl
		<< std::endl;

	#pragma region ������ ����
	// ������ �������� �������
	std::vector<std::vector<double>> x;
	std::vector<double> xRes;
	// ������ � ��������
	std::vector<std::vector<double>> r;
	x.resize(countOfNodes);
	r.resize(countOfNodes);
	// ���������� ��������
	int n = 0;
	// L2 ����� �������
	double res = 0;

	double startThirdStage = omp_get_wtime();
	solveSLAE(&IA, &JA, &A, &b, countOfNodes, tol, &x, &r, n, res, &xRes);
	double endThirdStage = omp_get_wtime();
	double endThirdStagePrint = endThirdStage - startThirdStage;
	std::cout << "����� �������� �����: " << endThirdStagePrint << " c." << std::endl;
	std::cout << "����� ������� ����� �� 1 �������: " << endThirdStagePrint / countOfNodes << " c." << std::endl
		<< std::endl;
	#pragma endregion

	if (isPrint)
	{
		printResult(resultGraph, &IA, &JA, seconds, countOfNodes);
		printSLAE(&A, &b, &IA);
		printSolveVector(res, n, &xRes);
	}

	std::cout << "//////" << std::endl
		<< std::endl;

	#pragma endregion
}

int main(int argc, char* argv[])
{

	#pragma region ���������� ������
	setlocale(LC_ALL, "Russian");

	int Nx, Ny, k1, k2, T;
	// �������� �������
	double tol;
	bool isPrint = false;
	std::vector<int> arguments, IA, JA;
	std::map<int, std::vector<int>> resultGraph;

	// ������ �������� - ������ �� ������
	for (int i = 1; i < argc; i++)
	{
		try
		{
			if (i == 6)
			{
				tol = std::stod(argv[i]);
			}
			else
			{
				arguments.push_back(atoi(argv[i]));
			}
		}
		catch (const std::exception&)
		{
			std::cout << "������������ ���� ���������� �������" << std::endl;
		}
	}

	if (arguments.empty() || arguments.size() < 5)
	{
		std::cout << "������� 6 ��������� �������� ����� (Nx, Ny, k1, k2, T, tol, isPrint)!" << std::endl;
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

		//doTask(1, Nx, Ny, k1, k2, resultGraph, IA, JA, countOfNodes, isPrint, tol);
	doTask(T, Nx, Ny, k1, k2, resultGraph, IA, JA, countOfNodes, isPrint, tol);

	return 0;
}

#pragma endregion
