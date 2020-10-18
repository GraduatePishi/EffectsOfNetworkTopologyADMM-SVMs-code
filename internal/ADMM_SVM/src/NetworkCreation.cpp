#include <vector>
#include <iostream>
#include <random>
#include <Eigen/Core>
#include <algorithm>
#include <cstdlib>
#include "utils.h"
#include "NetworkCreation.h"
#include <Eigen/Eigenvalues> 
#include "Printer.h"
#include "FileReader.h"
#include "ReadDataset.h"
#include <map>
using namespace std;
using ExceptionError = std::runtime_error;

 Eigen::VectorXd SecondSmallEigenVal(Eigen::MatrixXd& iMat) {
	Eigen::EigenSolver<Eigen::MatrixXd> es(iMat);
	Eigen::VectorXd tmpVec = es.eigenvalues().real();
	std::sort(begin(tmpVec), end(tmpVec), std::less<double>());
	return tmpVec;
}


//std::vector<int> Randomizing(std::vector<int>& vec) {
//	std::random_device rd;
//	std::mt19937 eng(rd());
//	std::shuffle(vec.begin(), vec.end(), eng);
//	return vec;
//}
 std::pair<Neighbors, Eigen::MatrixXd> createLineOrCircleNeigh(int NrOfNodes, Graph iGraph) {
	Neighbors retNeigh(NrOfNodes);
	Eigen::MatrixXd AdjacencyMat = Eigen::MatrixXd::Zero(NrOfNodes, NrOfNodes);
	for (int i = 0; i < NrOfNodes - 1; i++) {
		auto j = i + 1;
		retNeigh[i].emplace_back(j);
		retNeigh[j].emplace_back(i);
		AdjacencyMat(i, j) = 1;
		AdjacencyMat(j, i) = 1;
	}

	if ( iGraph=="circle") {
		retNeigh[0].emplace_back(NrOfNodes - 1);
		retNeigh[NrOfNodes - 1].emplace_back(0);
		AdjacencyMat(0, NrOfNodes - 1) = 1; //Connecting the first node to the last node to make a cycle
		AdjacencyMat(NrOfNodes - 1, 0) = 1;
	}
	return { retNeigh,AdjacencyMat };
}

std::pair<double, double> CalcEigenLaplacian(Eigen::MatrixXd const& iAdjacencyMat) {
	/*calculating the eignenvalues of LAplacian matrix*/
	int nrNodes = (int) iAdjacencyMat.rows();
	Eigen::MatrixXd LaplacianMat = iAdjacencyMat *(-1.0);
	Eigen::MatrixXd Lprime(nrNodes, nrNodes);
	for (int i = 0; i < nrNodes; i++) {
		LaplacianMat(i, i) = iAdjacencyMat.row(i).sum();
		Lprime.row(i) = LaplacianMat.row(i) / LaplacianMat(i, i);
	}

	auto Eigen = SecondSmallEigenVal(LaplacianMat); //returns the vector of all eigenvalues of laplacian
	auto EigenPrime = SecondSmallEigenVal(Lprime); //returns the vector of all eigenvalues of scaled laplacian

	if (Eigen[1] < 1e-10 && nrNodes != 1) {
		l::log("SecondSmallEigenVal for laplacian: {}, SecondSmallEigenVal for scaled laplacian: {}", Eigen[1], EigenPrime[1]);
		std::cout << "THE GRAPH IS NOT CONNECTED!!!!!" << "\n";
		throw ExceptionError("the graph is not connected!!!");
	}
	//if (iCrossVal == false) 
		std::cout << "SecondSmallEigenVal for laplacian : " << Eigen[1] << " ,  SecondSmallEigenVal for scaled laplacian:  " << EigenPrime[1] << "\n";

	return { Eigen[1], EigenPrime[1] };
}

Eigen::MatrixXd GetAdjacencyMat(int iNrOfNodes, std::string const& iGraphName, std::string iGraphFolder) {
	std::string graphFile;
	graphFile = ReadDataSetFile();
	graphFile.append(fmt::format("/graphs/{}/{}/", iNrOfNodes, iGraphFolder));
	graphFile.append(iGraphName);
	Eigen::MatrixXd AdjaMat = GraphReader(fs::u8path(graphFile));

	if (AdjaMat.rows()!= iNrOfNodes) {
		cout << "Graph nodes are " << AdjaMat.rows() << " while core numbers are !!!!!!!!!!!!!! " << iNrOfNodes << endl;
		throw ExceptionError("GRaph nodes NOT matches core numbers!!!!!!!!!!!!!!!!!!");
	}
	if (AdjaMat.rows() != AdjaMat.cols() ) {
		cout << "Graph is not symmetric " << AdjaMat.rows() << " Not " << AdjaMat.cols() << endl;
		throw ExceptionError("GRaph is not symmetric !!!!!!!!!!!!!!!!!!");

	}
	return AdjaMat;
}
Neighbors createRegularsNeigh(std::string const& iName, int NrOfNodes, std::string iGraphFolder) {
	Eigen::MatrixXd AdjacencyMat = GetAdjacencyMat(NrOfNodes, iName, iGraphFolder);
	CalcEigenLaplacian(AdjacencyMat);
	Neighbors retNeigh(NrOfNodes);
	for (int i = 0; i < NrOfNodes; i++) {
		for (int j = 0; j < NrOfNodes; j++) {
			if (AdjacencyMat(i, j) == 1.0) {
				retNeigh[i].push_back(j);
			}
		}
	}
	return retNeigh;
}

Neighbors Network(Graph iGraph, int NrOfNodes, std::string GraphFile, bool iCrossVal) {
	Neighbors NeighInfo(NrOfNodes);
	std::pair<Neighbors, Eigen::MatrixXd> NeighTemp;
	iGraph.append(".graph");
	if (iGraph=="Circle" || iGraph=="Line") {
		NeighTemp = createLineOrCircleNeigh(NrOfNodes, iGraph);
		NeighInfo = NeighTemp.first;
	}
	else {
		NeighInfo = createRegularsNeigh(iGraph, NrOfNodes, GraphFile);
	}


	//switch (iGraph)
	//{
	//case Graph::line:
	//case Graph::cycle:
	//	NeighTemp = createLineOrCircleNeigh(NrOfNodes,iGraph);
	//	NeighInfo = NeighTemp.first;
	//	break;
	//case Graph::TenReg:
	//	NeighInfo = createRegularsNeigh("10Regular.graph", NrOfNodes, GraphFile);
	//	break;
	//case Graph::TwentyReg:
	//	NeighInfo = createRegularsNeigh("20Regular.graph", NrOfNodes, GraphFile);
	//	break;
	//case Graph::TwentyFiveReg:
	//	NeighInfo = createRegularsNeigh("25Regular.graph", NrOfNodes, GraphFile);
	//	break;
	//case Graph::ThirtyReg:
	//	NeighInfo = createRegularsNeigh("30Regular.graph", NrOfNodes, GraphFile);
	//	break;
	//case Graph::FourtyReg:
	//	NeighInfo = createRegularsNeigh("40Regular.graph", NrOfNodes, GraphFile);
	//	break;
	//case Graph::FourtyFiveReg:
	//	NeighInfo = createRegularsNeigh("45Regular.graph", NrOfNodes, GraphFile);
	//	break;
	//case Graph::FourtySevenReg:
	//	NeighInfo = createRegularsNeigh("47Regular.graph", NrOfNodes, GraphFile);
	//	break;
	//case Graph::FiftyReg:
	//	NeighInfo = createRegularsNeigh("50Regular.graph", NrOfNodes, GraphFile);
	//	break;
	//case Graph::SixtyReg:
	//	NeighInfo = createRegularsNeigh("60Regular.graph", NrOfNodes, GraphFile);
	//	break;
	//case Graph::SeventyReg:
	//	NeighInfo = createRegularsNeigh("70Regular.graph", NrOfNodes, GraphFile);
	//	break;
	//case Graph::SeventyNineReg:
	//	NeighInfo = createRegularsNeigh("79Regular.graph", NrOfNodes, GraphFile);
	//	break;
	//case Graph::EightyReg:
	//	NeighInfo = createRegularsNeigh("80Regular.graph", NrOfNodes, GraphFile);
	//	break;
	//case Graph::NintyReg:
	//	NeighInfo = createRegularsNeigh("90Regular.graph", NrOfNodes, GraphFile);
	//	break;
	//case Graph::HundredReg:
	//	NeighInfo = createRegularsNeigh("100Regular.graph", NrOfNodes, GraphFile);
	//	break;
	//case Graph::HundredOneReg:
	//	NeighInfo = createRegularsNeigh("101Regular.graph", NrOfNodes, GraphFile);
	//	break;
	//case Graph::HundredTenReg:
	//	NeighInfo = createRegularsNeigh("110Regular.graph", NrOfNodes, GraphFile);
	//	break;
	//case Graph::HundredTwentyReg:
	//	NeighInfo = createRegularsNeigh("120Regular.graph", NrOfNodes, GraphFile);
	//	break;
	//case Graph::TwoHundredReg:
	//	NeighInfo = createRegularsNeigh("200Regular.graph", NrOfNodes, GraphFile);
	//	break;
	//case Graph::HundredNineteenReg:
	//	NeighInfo = createRegularsNeigh("119Regular.graph", NrOfNodes, GraphFile);
	//	break;
	//case Graph::TwoReg:
	//	NeighInfo = createRegularsNeigh("2Regular.graph", NrOfNodes, GraphFile);
	//	break;
	//case Graph::ThreeReg:
	//	NeighInfo = createRegularsNeigh("3Regular.graph", NrOfNodes, GraphFile);
	//	break;
	//case Graph::FiveReg:
	//	NeighInfo = createRegularsNeigh("5Regular.graph", NrOfNodes, GraphFile);
	//	break;
	//case Graph::SevenReg:
	//	NeighInfo = createRegularsNeigh("7Regular.graph", NrOfNodes, GraphFile);
	//	break;
	//case Graph::NineReg:
	//	NeighInfo = createRegularsNeigh("9Regular.graph", NrOfNodes, GraphFile);
	//	break;
	//case Graph::ElevenReg:
	//	NeighInfo = createRegularsNeigh("11Regular.graph", NrOfNodes, GraphFile);
	//	break;
	//case Graph::ThirteenReg:
	//	NeighInfo = createRegularsNeigh("13Regular.graph", NrOfNodes, GraphFile);
	//	break;
	//case Graph::FifteenReg:
	//	NeighInfo = createRegularsNeigh("15Regular.graph", NrOfNodes, GraphFile);
	//	break;
	//case Graph::NineteenReg:
	//	NeighInfo = createRegularsNeigh("19Regular.graph", NrOfNodes, GraphFile);
	//	break;
	//case Graph::TwentyOneReg:
	//	NeighInfo = createRegularsNeigh("21Regular.graph", NrOfNodes, GraphFile);
	//	break;
	//case Graph::TwentyThreeReg:
	//	NeighInfo = createRegularsNeigh("23Regular.graph", NrOfNodes, GraphFile);
	//	break;
	//case Graph::TwentySevenReg:
	//	NeighInfo = createRegularsNeigh("27Regular.graph", NrOfNodes, GraphFile);
	//	break;
	//case Graph::ThirtyOneReg:
	//	NeighInfo = createRegularsNeigh("31Regular.graph", NrOfNodes, GraphFile);
	//	break;
	//case Graph::ThirtyThreeReg:
	//	NeighInfo = createRegularsNeigh("33Regular.graph", NrOfNodes, GraphFile);
	//	break;
	//case Graph::ThirtyFiveReg:
	//	NeighInfo = createRegularsNeigh("35Regular.graph", NrOfNodes, GraphFile);
	//	break;
	//case Graph::ThirtyNineReg:
	//	NeighInfo = createRegularsNeigh("39Regular.graph", NrOfNodes, GraphFile);
	//	break;
	//case Graph::ThirtyFiveRandom:
	//	NeighInfo = createRegularsNeigh("35Random.graph", NrOfNodes, GraphFile);
	//	break;
	//case Graph::ThirtyOneRandom:
	//	NeighInfo = createRegularsNeigh("31Random.graph", NrOfNodes, GraphFile);
	//	break;
	//case Graph::TwentySevenRandom:
	//	NeighInfo = createRegularsNeigh("27Random.graph", NrOfNodes, GraphFile);
	//	break;
	//case Graph::TwentyThreeRandom:
	//	NeighInfo = createRegularsNeigh("25Random.graph", NrOfNodes, GraphFile);
	//	break;
	//case Graph::NineteenRandom:
	//	NeighInfo = createRegularsNeigh("19Random.graph", NrOfNodes, GraphFile);
	//	break;
	//case Graph::FifteenRandom:
	//	NeighInfo = createRegularsNeigh("15Random.graph", NrOfNodes, GraphFile);
	//	break;
	//case Graph::ElevenRandom:
	//	NeighInfo = createRegularsNeigh("11Random.graph", NrOfNodes, GraphFile);
	//	break;
	//case Graph::SevenRandom:
	//	NeighInfo = createRegularsNeigh("7Random.graph", NrOfNodes, GraphFile);
	//	break;
	//case Graph::hardCoded:
	//	std::cout << "running special case " << (int)iGraph;
	//	break;
	//default:	
	//	throw ExceptionError("The graph type does not exist!!");
	//	break;
	//}
	return NeighInfo;
}
