#pragma once
#include "Types.h"
//enum class Graph {
//	line,
//	cycle,
//	hardCoded,
//	unconnected,
//	unConnIslands,
//	self, //Each node is its own neighbor
//	Shirin_ThreeReg,
//	Shirin_FiveReg,
//	Comb,
//	Island,
//	Star,
//	Alone,
//	/************************/
//	Random,
//	/********Network for big datasets************/
//	TenReg,
//	TwentyReg,
//	TwentyFiveReg,
//	ThirtyReg,
//	FourtyReg,
//	FourtyFiveReg,
//	FourtySevenReg,
//	FiftyReg,
//	SixtyReg,
//	SeventyReg,
//	SeventyNineReg,
//	EightyReg,
//	NintyReg,
//	HundredReg,
//	HundredOneReg,
//	HundredTenReg,
//	HundredTwentyReg,
//	TwoHundredReg,
//	HundredNineteenReg,
//	/********Network for small datasets************/
//	TwoReg,
//	ThreeReg,
//	FiveReg,
//	SevenReg,
//	NineReg,
//	ElevenReg,
//	ThirteenReg,
//	FifteenReg,
//	NineteenReg,
//	TwentyOneReg,
//	TwentyThreeReg,
//	TwentySevenReg,
//	ThirtyOneReg,
//	ThirtyThreeReg,
//	ThirtyFiveReg,
//	ThirtyNineReg,
//	/**************** Random with max degree*/
//	ThirtyFiveRandom,
//	ThirtyOneRandom,
//	TwentySevenRandom,
//	TwentyThreeRandom,
//	NineteenRandom,
//	FifteenRandom,
//	ElevenRandom,
//	SevenRandom
//};
Neighbors Network(Graph iGraph, int NrOfNodes, std::string GraphFile, bool iCrossVal);
Neighbors createRegularsNeigh(std::string const& iName, int NrOfNodes, std::string iGraphFolder);
std::pair<Neighbors, Eigen::MatrixXd> createLineOrCircleNeigh(int NrOfNodes, Graph iGraph);
std::pair<double, double> CalcEigenLaplacian(Eigen::MatrixXd const& iAdjacencyMat);
Eigen::MatrixXd GetAdjacencyMat(int iNrOfNodes, std::string const& iGraphName, std::string iGraphFolder);
