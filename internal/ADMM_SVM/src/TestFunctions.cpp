#include "SequencesOfTests.h"
#include "Types.h"
#include "ChoosingT.h"
#include "Kernel.h"
#include "NodeSolver_linear.h"
#include "NodeSolver.h"
#include <iostream>
#include "Partitioning.h"
#include "FileReader.h"
#include "Solver.h"
#include "CrossValidation.h"
#include "InverseOfMat.h"
#include "NetworkCreation.h"
#include "ReadDataset.h"
#include "utils.h"
#include "TestPhase.h"
#include "Scaling.h"
#include <math.h>
#include "Randomizing.h"
#include "FoldData.h"
#include <algorithm>
#include "ProblemSettings.h"
#include "DistributedData.h"
#ifdef _WIN32
#include <limits.h>
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "C:\Program Files (x86)\Microsoft SDKs\MPI\Include\mpi.h"
#else
#include "mpi.h"
#endif
using namespace std;
#define RUNTESTS
#ifndef RUNTESTS
struct TestProblem
{
	TestProblem(int iRow, int iCol, double iSeed, int iWorldSize, bool iShuffle, int iWorldRank, string iScaling, bool iPreT, TypeForDesign iTtype, int iReducedSize, double iJcTest, double iEtaTest, double iGamma, bool iADMMParam, std::string igraphFile, std::string iDataname)
	{
		std::string iPath = ReadDataSetFile();
		iPath.append(iDataname);
		fs::path filePath = fs::u8path(iPath);
		mTrainfilePath = filePath;
		L_test = iReducedSize;
		//const ProblemStatment prob(0, 0);// = FileReader(filePath);
		//auto DataPart = DistributedData(iWorldSize, iWorldRank, prob, iShuffle, iScaling, mTrainfilePath, L_test, iPreT, iTtype, iSeed);
		//Data = std::get<0>(DataPart);
		//Label = std::get<1>(DataPart);
		//cout << "DataMat[ node: " << iWorldRank << "]\n" << Data << endl;
		//cout << "DataLbl[ node: " << iWorldRank << "]\n" << Label << endl;
		int rows = iRow;
		int cols = iCol;
		NoDataPoints = rows;
		p_test = cols + 1;
		colSize = cols;
		JC_test = iJcTest;
		EtaTest = iEtaTest;
		wj_result = WeighVector::Ones(iReducedSize)*(-1000);
		average_wi = WeighVector::Zero(iReducedSize);
		gamma = iGamma;
		ADMMParam = iADMMParam;
		graphFile = igraphFile;
		nrNodes = iWorldSize;
	}

	DataMatrix Data;
	Graph graph = Graph::line;
	std::string graphFile;
	fs::path mTrainfilePath;
	std::string mDataname;
	RowMajorMatirx Data_lin;
	LabelVector Label;
	Row NoDataPoints;
	Column colSize;
	int p_test{ -1 };
	int L_test{ -1 };
	int nrNodes{ -1 };
	double JC_test{ -1 };
	double EtaTest{ -1 };
	double gamma{ -1 };
	double mSeed{ 32 };
	int mWorld_rank{ -1 }, mWorld_size{ -1 };
	bool ADMMParam;
	bool mShuffle{ false };
	std::string mScaling{ "no" };
	WeighVector wj_result;
	Bias bj_result = -1000;
	Multipliers initLambda;
	WeighVector initWj;
	Bias init_bj = 5;
	WeighVector average_wi;
	Bias average_bi = 0.0;
	Neighbors neighborIndexes;
	IndexOfNeighbors localNeigVec;;
	Multipliers lambda;
	TypeForDesign Ttype;// = TypeForDesign::identity;
	DimReductionMat consensusT;
	TypeForKernel kerneltype;
	string name;
	string testTerminationCriteria = "AllThree";
};
enum class ProblemToTest
{
	SimpleLinear2_first6Values,
	Covtype,
	SimpleLinear2,
	Simple1,
	fourclass1,
	ForBroad,
	splice20GradTest,
	fourclass100GradTest,
	heart,
	linear,
	splice,
	spliceForHitRate,
	simpleLinear2_16,
	neighborsCheck,
	Higgs
};
auto createTestProblem(ProblemToTest iVal)//int type)
{
	if (ProblemToTest::SimpleLinear2_first6Values == iVal)
	{
		std::string name = "simpleLinear2_first6values.dat";
		double JC_test = 200.0;//Best value for simpleLinear2 is eta=1 and JC=10, 200
		double EtaTest = 1;
		double gamma = 0.0;

		bool ADMMPar = false, shuffle{ false };
		std::string GraphFolder = "Regular1";
		int Tsize{ 1000 }, world_size{ 2 }, world_rank{ 0 };
		double iSeed = 32;
		int row{ 6 }, col{ 2 };
		TestProblem retProblem(row,col,iSeed, world_size, shuffle, world_rank, "no", false, TypeForDesign::random, Tsize, JC_test, EtaTest, gamma, ADMMPar, GraphFolder, name);
		retProblem.name = "simpleLinear2_first6values";
		retProblem.graphFile = GraphFolder;
		retProblem.p_test = 3;
		retProblem.L_test = 2;
		retProblem.lambda = Multipliers::Ones(retProblem.NoDataPoints);
		retProblem.nrNodes = 1;
		retProblem.wj_result << -1, -1;
		retProblem.bj_result = 3;
		retProblem.mScaling = "no";
		retProblem.initLambda = Multipliers::Ones(retProblem.NoDataPoints)*5.0;
		retProblem.initWj = WeighVector::Ones(retProblem.L_test)*5.0;
		retProblem.init_bj = 5.0;
		retProblem.kerneltype = TypeForKernel::linear;
		retProblem.neighborIndexes = Network(Graph::hardCoded, retProblem.nrNodes, retProblem.graphFile, false);
		retProblem.localNeigVec= retProblem.neighborIndexes[0];
		MinVec MinVector;
		MaxVec MaxVector;
		retProblem.consensusT = ConsensusT(retProblem.L_test, retProblem.Data, TypeForDesign::identity, MaxVector,MinVector, 32);
		return retProblem;
	}
	if (ProblemToTest::Covtype == iVal)
	{
		std::string name = "covtype.dat";
		double JC_test = 1.0;//Best value for simpleLinear2 is eta=1 and JC=10, 200
		double EtaTest = 100;
		double gamma = 32.0;
		bool ADMMPar = false, shuffle{ false };
		std::string GraphFolder = "Regular1";
		int Tsize{ 1000 }, world_size{ 2 }, world_rank{ 0 };
		double iSeed = 32;
		int row{ 6 }, col{ 2 };
		TestProblem retProblem(row, col, iSeed, world_size, shuffle, world_rank, "no", false, TypeForDesign::random, Tsize, JC_test, EtaTest, gamma, ADMMPar, GraphFolder, name);
		retProblem.name = "covtype";
		retProblem.graphFile = GraphFolder;
		retProblem.L_test = 1000;
		retProblem.lambda = Multipliers::Ones(retProblem.NoDataPoints);
		retProblem.nrNodes = 1;
		retProblem.wj_result << -1, -1;
		retProblem.bj_result = 3;
		retProblem.mScaling = "no";
		retProblem.initLambda = Multipliers::Ones(retProblem.NoDataPoints)*5.0;
		retProblem.initWj = WeighVector::Ones(retProblem.L_test)*5.0;
		retProblem.init_bj = 5.0;
		retProblem.kerneltype = TypeForKernel::rbf;
		//retProblem.neighborIndexes = Network(Graph::hardCoded, retProblem.nrNodes, retProblem.graphFile, true);
		//retProblem.neighborIndexes = getNeighbors(retProblem.nrNodes);
		MinVec MinVector;
		MaxVec MaxVector;
		retProblem.consensusT = ConsensusT(retProblem.L_test, retProblem.Data, TypeForDesign::random, MaxVector, MinVector, 32);
		return retProblem;
	}
	if (ProblemToTest::Higgs == iVal)
	{
		std::string name = "higgs10k.dat";
		double JC_test = 1.0;//Best value for simpleLinear2 is eta=1 and JC=10, 200
		double EtaTest = 100;
		double gamma = 32.0;
		bool ADMMPar = false, shuffle{ false };
		std::string GraphFolder = "Regular1";
		int Tsize{ 1000 }, world_size{ 2 }, world_rank{ 0 };
		double iSeed = 32;
		int row{ 6 }, col{ 2 };
		TestProblem retProblem(row, col, iSeed, world_size, shuffle, world_rank, "no", false, TypeForDesign::random, Tsize, JC_test, EtaTest, gamma, ADMMPar, GraphFolder, name);

		retProblem.name = "higgs10k";
		retProblem.graphFile = GraphFolder;
		retProblem.L_test = 100;
		retProblem.lambda = Multipliers::Ones(retProblem.NoDataPoints);
		retProblem.nrNodes = 1;
		retProblem.wj_result << -1, -1;
		retProblem.bj_result = 3;
		retProblem.mScaling = "no";
		retProblem.initLambda = Multipliers::Ones(retProblem.NoDataPoints)*5.0;
		retProblem.initWj = WeighVector::Ones(retProblem.L_test)*5.0;
		retProblem.init_bj = 5.0;
		retProblem.kerneltype = TypeForKernel::rbf;
		//retProblem.neighborIndexes = Network(Graph::hardCoded, retProblem.nrNodes, retProblem.graphFile, true);
		//retProblem.neighborIndexes = getNeighbors(retProblem.nrNodes);
		MinVec MinVector;
		MaxVec MaxVector;
		retProblem.consensusT = ConsensusT(retProblem.L_test, retProblem.Data, TypeForDesign::random, MaxVector, MinVector, 32);
		return retProblem;
	}
	if (ProblemToTest::linear == iVal)
	{
		int iWorld_size{ -1 }, iWorld_rank{ -1 };
		MPI_Comm_size(MPI_COMM_WORLD, &iWorld_size);
		MPI_Comm_rank(MPI_COMM_WORLD, &iWorld_rank);
		std::string name = "linear.dat";
		double JC_test = 1.0;//Best value for simpleLinear2 is eta=1 and JC=10, 200
		double EtaTest = 100;
		double gamma = 32.0;
		bool ADMMPar = false, shuffle{ false };
		std::string GraphFolder = "Regular1";
		int Tsize{ 1000 };
		double iSeed = 32;
		int row{ 6 }, col{ 2 };
		TestProblem retProblem(row, col, iSeed, iWorld_size, shuffle, iWorld_rank, "no", false, TypeForDesign::random, Tsize, JC_test, EtaTest, gamma, ADMMPar, GraphFolder, name);
		retProblem.name = "linear";
		retProblem.graphFile = GraphFolder;
		retProblem.L_test = 100;
		retProblem.lambda = Multipliers::Ones(retProblem.NoDataPoints);
		retProblem.mScaling = "no";
		retProblem.kerneltype = TypeForKernel::rbf;
		retProblem.neighborIndexes = Network(Graph::line, retProblem.nrNodes, retProblem.graphFile, true);
		//retProblem.neighborIndexes = getNeighbors(retProblem.nrNodes);
		MinVec MinVector;
		MaxVec MaxVector;
		retProblem.consensusT = ConsensusT(retProblem.L_test, retProblem.Data, TypeForDesign::random, MaxVector, MinVector, 32);
		return retProblem;
	}
}
bool DistributedDataTest(std::string iDataName,  bool PositiveFlag) { // this is for 4 nodes

	int iWorld_size{ -1 }, iWorld_rank{ -1 };
	MPI_Comm_size(MPI_COMM_WORLD, &iWorld_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &iWorld_rank);


	ProblemStatment iProb(0, 0); // this is for reading from file directly
	ProblemStatment TestProb(0, 0);
	std::string scaling = "no";
	std::string Testset = ReadDataSetFile();
	Testset.append(iDataName);
	fs::path TestfilePath = fs::u8path(Testset);
	bool shuffle = false, finishTest = false;
	int LocalTest = 0;
	iProb = FileReader(Testset);
	int l = 10;
	TypeForDesign t_type = TypeForDesign::random;
	double seed = 32;
	auto DataPart = DistributedData(iWorld_size, iWorld_rank, TestProb, shuffle, scaling, TestfilePath, l, false, t_type, seed);
	int partSize = (int)(iProb.Data.rows() / iWorld_size);
	std::vector<RowMajorMatirx> A(iWorld_size);
	A[iWorld_rank] = iProb.Data.block(iWorld_rank*partSize, 0, partSize, (int)iProb.Data.cols());
	if (PositiveFlag == false) {
		if (iWorld_rank == 1) {
			A[1] = RowMajorMatirx::Identity(partSize, (int)iProb.Data.cols())*2.0;
		}
	}
	if (!A[iWorld_rank].isApprox(std::get<0>(DataPart))) {
		//cout << "**********************A[" << iWorld_rank << "] and Data[" << iWorld_rank << "] is not the same for " << iDataName <<"!!!!" << endl;
		LocalTest = 0;
	}
	else {
		LocalTest = 1;
	}
	int GlobalTest = 0;

	MPI_Allreduce(&LocalTest, &GlobalTest, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	if (GlobalTest != iWorld_size) {
		if (iWorld_rank == 0) {
			//cout << "Distributing data for " << iDataName << " is NOT passed!!!!" << endl;
		}
		finishTest = false;
	}
	else {
		finishTest = true;
		if (iWorld_rank == 0) {
			//cout << "Distributing data for " << iDataName << " is PASSED!!!!" << endl;
		}
	}

	return finishTest;
}

std::tuple<MinVec, MaxVec, RowMajorMatirx, RowMajorMatirx> TestScaling(std::string scaling, const fs::path orgData, const fs::path scaleData) {
	std::tuple<MinVec, MaxVec, RowMajorMatirx, RowMajorMatirx> Data;
	int dataLength = 3650;
	ProblemStatment Prob(0, 0);
	RowMajorMatirx Example1;
	std::get<2>(Data) = ReadT(orgData, dataLength, 1);//Original Data
	std::get<3>(Data) = ReadT(scaleData, dataLength, 1);//scaled Data
	MinVec minTemp(1);
	MaxVec maxTemp(1);
	ScalingTrain(minTemp, maxTemp, std::get<2>(Data), scaling);
	std::get<0>(Data) = minTemp;
	std::get<1>(Data) = maxTemp;
	return Data;
}
double CondNumMatrix(TestProblem& problem) {

	double cond = 0.0;
	std::vector<std::pair<std::unique_ptr<nodesolver>, NodeData>> NodeInfo(1);
	NodeInfo[0] = std::make_pair(std::make_unique<nodesolver>(problem.Data, problem.Label,problem.consensusT, problem.kerneltype, problem.localNeigVec,problem.JC_test,problem.EtaTest, problem.gamma,problem.mWorld_rank),
		NodeData(problem.L_test,problem.Data.rows(), problem.Data.cols())
	);
	Eigen::JacobiSVD<DataMatrix> svd_Fixed(NodeInfo[0].first->mFixedMatrix);
	double cond_Fixed = svd_Fixed.singularValues()(0) / svd_Fixed.singularValues()(svd_Fixed.singularValues().size() - 1);
	return cond;
}

double ConditionTesting(string name) { // this is for 4 nodes
	int iWorld_size{ -1 }, iWorld_rank{ -1 };
	MPI_Comm_size(MPI_COMM_WORLD, &iWorld_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &iWorld_rank);

	ProblemStatment iProb(0, 0); // this is for reading from file directly
	std::string Testset = ReadDataSetFile();
	Testset.append(name);
	fs::path TestfilePath = fs::u8path(Testset);
	int L = 500;
	auto DataPart = DistributedData(iWorld_size, iWorld_rank, iProb, false, "no", TestfilePath, L, false, TypeForDesign::random, 32);
	auto dataMat = std::get<0>(DataPart);
	auto dataLbl = std::get<1>(DataPart);
	int row = std::get<2>(DataPart);
	int col = std::get<3>(DataPart);
	auto TMat = std::get<8>(DataPart);
	//cout << "DataMat[node " << iWorld_rank <<"] :\n " << dataMat << endl;
	//cout << "DataLbl[node " << iWorld_rank << "] :\n " << dataLbl << endl;

	//cout << "problem.initLambda :\n " << problem.consensusT << endl;
	IndexOfNeighbors localNeigVec;
	Neighbors neighborsVectors;
	neighborsVectors = Network(Graph::line, iWorld_size, "", false);

	std::vector<std::pair<std::unique_ptr<nodesolver>, NodeData>> NodeInfo(1);
	NodeInfo[0] = std::make_pair(std::make_unique<nodesolver>(dataMat, dataLbl, TMat, TypeForKernel::rbf, neighborsVectors[iWorld_rank], 10, 10, 0.5, iWorld_rank),
		NodeData(L, row, col)
	);

	Eigen::JacobiSVD<DataMatrix> svd_Fixed(NodeInfo[0].first->mFixedMatrix);
	double cond_Fixed = svd_Fixed.singularValues()(0) / svd_Fixed.singularValues()(svd_Fixed.singularValues().size() - 1);
	double globalCond{ 0.0 };

	MPI_Allreduce(&cond_Fixed, &globalCond, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	//cout << "cond_Fixed[node " << iWorld_rank << "]: " << cond_Fixed << endl;

	return globalCond;
}
bool testNetwork(int iNrOfNodes,int iExpectedRegularity, std::string const& iGraphName, std::string iGraphFolder) {
	Eigen::MatrixXd Adja1 = GetAdjacencyMat(iNrOfNodes, iGraphName, iGraphFolder);
	if (Adja1.rows() < 1 || Adja1.cols() < 1) {
		return false;
	}
	for (int i = 0; i< Adja1.rows(); i++) {
		if (Adja1.row(i).sum() != iExpectedRegularity) {
			return false;
		}
		if (Adja1.col(i).sum() != iExpectedRegularity) {
			return false;
		}
	}
	if (!Adja1.isApprox(Adja1.transpose())) {
		return false;
	}
	return true;
}
bool testEigenvalLaplac(int iNrOfNodes, std::string const& iGraphName, std::string iGraphFolder) {
	Eigen::MatrixXd Adja = GetAdjacencyMat(iNrOfNodes, iGraphName, iGraphFolder);
	std::pair<double, double> EigLAplac;
	EigLAplac=CalcEigenLaplacian(Adja);
	if (EigLAplac.first < 0 || EigLAplac.second < 0) 
		return false;
	if (EigLAplac.first > 1e+15 || EigLAplac.second > 1e+15)
		return false;
	return true;
}
bool TestingManualKernlandFunc(const RowMajorMatirx& iTestData, const DataMatrix& PartialDataMat,TypeForKernel iKerneltype, const  DimReductionMat& iT, KrnlMatrx iNormT_squared, KrnlMatrx iNorm_Xsquared, double iGamma) {
	auto KrnlTest = getTestKernel(iTestData, PartialDataMat, iT, iNormT_squared, iNorm_Xsquared, iGamma);
	Eigen::MatrixXd Krnl_TestDatXj_1 = KrnlTest.first;
	Eigen::MatrixXd K_iTestAndT_1 = KrnlTest.second;

	Eigen::MatrixXd K_iTestAndT_2 = K(iTestData, iT, iKerneltype, iGamma);
	Eigen::MatrixXd Krnl_TestDatXj_2 = K(iTestData, PartialDataMat, iKerneltype, iGamma);
	if (!Krnl_TestDatXj_1.isApprox(Krnl_TestDatXj_2) || !K_iTestAndT_1.isApprox(K_iTestAndT_2))
		return false;
	
	return true;
}
bool TestKernel(std::string iName, int iL) {
	int world_rank{ -1 }, world_size{ -1 };
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	ProblemStatment prob(0, 0), probTest(0, 0);
	std::string Trainset;
	Trainset = ReadDataSetFile();
	std::string BoxPath = Trainset; //the address of all datasets
	Trainset.append(iName);
	fs::path TrainfilePath = fs::u8path(Trainset);
	auto const& DataPart = DistributedData(world_size, world_rank, prob, false, "no", TrainfilePath, iL, true, TypeForDesign::random, 0);
	auto localDataMat = std::get<0>(DataPart);
	auto localDataVec = std::get<1>(DataPart);
	auto ProblemRowSize = std::get<2>(DataPart);
	auto ProblemColSize = std::get<3>(DataPart);
	probTest = std::get<9>(DataPart);
	DimReductionMat T(iL, ProblemColSize);
	T = std::get<8>(DataPart);
	double jc = 1, eta = 5, gamma = 0.05;

	IndexOfNeighbors localNeigVec ;
	std::vector<std::pair<std::unique_ptr<nodesolver>, NodeData>> NodeInfo(1);
	localNeigVec = DistrNeighToAllNnodes(world_rank, Graph::cycle, world_size, "Regular1", false, jc, gamma, eta);

	NodeInfo[0] = std::make_pair(std::make_unique<nodesolver>(localDataMat, localDataVec, T, TypeForKernel::rbf, localNeigVec, jc, eta, gamma, world_rank),
		NodeData(iL, (int)localDataMat.rows(), (int)localDataMat.cols()));

	cout << "mFixedMatrix[node: " << world_rank << "]: \n" << NodeInfo[0].first->mFixedMatrix << endl;
	//if (iL<10) {
	//	cout << "T[node: " << world_rank << "]: \n" << T << endl;
	//	cout << "localDataMat[node: " << world_rank << "]: \n" << localDataMat << endl;
	//	cout << "probTest.Data[node: " << world_rank << "]: \n" << probTest.Data << endl;
	//	cout << "m_normT_squared[node: " << world_rank << "]: \n" << NodeInfo[0].first->m_normT_squared << endl;
	//	cout << "m_norm_Xsquared[node: " << world_rank << "]: \n" << NodeInfo[0].first->m_norm_Xsquared << endl;
	//}
	//else {
	//	cout << "T[node: " << world_rank << "]: \n" << T.block(0, 0,20, 10) << endl;
	//	cout << "localDataMat[node: " << world_rank << "]: \n" << localDataMat.block(0, 0, 20, localDataMat.cols()/2) << endl;
	//	cout << "probTest.Data[node: " << world_rank << "]: \n" << probTest.Data.block(0, 0, 20, 10) << endl;
	//	cout << "m_normT_squared[node: " << world_rank << "]: \n" << NodeInfo[0].first->m_normT_squared.block(0, 0, 20, 10) << endl;
	//	cout << "m_norm_Xsquared[node: " << world_rank << "]: \n" << NodeInfo[0].first->m_norm_Xsquared.block(0, 0, 20, 10) << endl;

	//}
	bool testResult=TestingManualKernlandFunc(probTest.Data, localDataMat, TypeForKernel::rbf,T, NodeInfo[0].first->m_normT_squared, NodeInfo[0].first->m_norm_Xsquared,gamma);
	return testResult;
}
bool TestKernelForNonLin() { // This test should be with only 2 nodes, i.e., "mpiexec -n 2"
	int world_rank{ -1 }, world_size{ -1 };
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	ProblemStatment prob(0, 0), probTest(0, 0);
	std::string Trainset;
	Trainset = ReadDataSetFile();
	std::string BoxPath = Trainset; //the address of all datasets
	Trainset.append("nonLin.dat");
	int iL = 2;
	fs::path TrainfilePath = fs::u8path(Trainset);
	auto const& DataPart = DistributedData(world_size, world_rank, prob, false, "no", TrainfilePath, iL, true, TypeForDesign::random, 0);
	auto localDataMat = std::get<0>(DataPart);
	auto localDataVec = std::get<1>(DataPart);
	auto ProblemRowSize = std::get<2>(DataPart);
	auto ProblemColSize = std::get<3>(DataPart);
	probTest = std::get<9>(DataPart);
	DimReductionMat T(iL, ProblemColSize);
	T = std::get<8>(DataPart);
	double jc = 1, eta = 5, gamma = 0.05;

	IndexOfNeighbors localNeigVec;
	std::vector<std::pair<std::unique_ptr<nodesolver>, NodeData>> NodeInfo(1);
	localNeigVec = DistrNeighToAllNnodes(world_rank, Graph::cycle, world_size, "Regular1", false, jc, gamma, eta);

	NodeInfo[0] = std::make_pair(std::make_unique<nodesolver>(localDataMat, localDataVec, T, TypeForKernel::rbf, localNeigVec, jc, eta, gamma, world_rank),
		NodeData(iL, (int)localDataMat.rows(), (int)localDataMat.cols()));

	double normMat=NodeInfo[0].first->mFixedMatrix.norm();
	double normMAtall = -1;
	MPI_Allreduce(&normMat, &normMAtall, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	return normMAtall / 2 == normMat;
	//if (iL<10) {
	//	cout << "T[node: " << world_rank << "]: \n" << T << endl;
	//	cout << "localDataMat[node: " << world_rank << "]: \n" << localDataMat << endl;
	//	cout << "probTest.Data[node: " << world_rank << "]: \n" << probTest.Data << endl;
	//	cout << "m_normT_squared[node: " << world_rank << "]: \n" << NodeInfo[0].first->m_normT_squared << endl;
	//	cout << "m_norm_Xsquared[node: " << world_rank << "]: \n" << NodeInfo[0].first->m_norm_Xsquared << endl;
	//}
	//else {
	//	cout << "T[node: " << world_rank << "]: \n" << T.block(0, 0,20, 10) << endl;
	//	cout << "localDataMat[node: " << world_rank << "]: \n" << localDataMat.block(0, 0, 20, localDataMat.cols()/2) << endl;
	//	cout << "probTest.Data[node: " << world_rank << "]: \n" << probTest.Data.block(0, 0, 20, 10) << endl;
	//	cout << "m_normT_squared[node: " << world_rank << "]: \n" << NodeInfo[0].first->m_normT_squared.block(0, 0, 20, 10) << endl;
	//	cout << "m_norm_Xsquared[node: " << world_rank << "]: \n" << NodeInfo[0].first->m_norm_Xsquared.block(0, 0, 20, 10) << endl;

	//}

}
namespace {
	///*************** The follwoing scaling tests run with a simple "mpiexec -n 2 MachineLearningExec.exe" script******/
	//TEST(DistributingData, positive_simpleLinear2) { //#1
	//	EXPECT_EQ(true, DistributedDataTest("simpleLinear2.dat",  true));
	//	MPI_Barrier(MPI_COMM_WORLD);
	//}
	//TEST(DistributingData, positive_svmguid) {//#2
	//	EXPECT_EQ(true, DistributedDataTest("svmguid.dat", true));
	//	MPI_Barrier(MPI_COMM_WORLD);
	//}
	//TEST(DistributingData, positive_higgs10k) {//#3
	//	EXPECT_EQ(true, DistributedDataTest("higgs10k.dat", true));
	//	MPI_Barrier(MPI_COMM_WORLD);
	//}
	//TEST(DistributingData, Negative_simpleLinear2) {//#4
	//	EXPECT_EQ(false, DistributedDataTest("simpleLinear2.dat", false));
	//	MPI_Barrier(MPI_COMM_WORLD);
	//}
	//TEST(DistributingData, Negative_svmguid) {//#5
	//	EXPECT_EQ(false, DistributedDataTest("svmguid.dat", false));
	//	MPI_Barrier(MPI_COMM_WORLD);
	//}
	//TEST(DistributingData, Negative_higgs10k){//#6
	//	EXPECT_EQ(false, DistributedDataTest("higgs10k.dat", false));
	//	MPI_Barrier(MPI_COMM_WORLD);
	//}

	/////*************** The follwoing scaling tests run with a simple "MachineLearningExec.exe" script******/
	//TEST(Scaling, NormalScaling) {//#7
	//	auto DataSet= TestScaling("normal", "Examp1.dat", "scaledExamp1.dat");
	//	EXPECT_TRUE(std::get<2>(DataSet).isApprox(std::get<3>(DataSet))) << "Wrong normal scaling!!!!!!!!!!!!!!";; //the expected results from Python and my actuall results from my caling function are the same when the expected results are written with high precision into the file 
	//}
	//TEST(Scaling, MinVecScaling) {//#8
	//	auto DataSet = TestScaling("normal", "Examp1.dat", "scaledExamp1.dat");
	//	EXPECT_TRUE(abs(0 - std::get<0>(DataSet)[0])<1e-14) << "Min value is wrong!!";
	//}
	//TEST(Scaling, MaxVecScaling) {//#9
	//	auto DataSet = TestScaling("normal", "Examp1.dat", "scaledExamp1.dat");
	//	EXPECT_TRUE(abs(26.3 - std::get<1>(DataSet)[0])< 1e-14) << "Max value is wrong!!";
	//}

	/////*************** The follwoing scaling tests run with a simple "MachineLearningExec.exe" script******/
	//TEST(Scaling, StandardScaling) {//#10 if I take away "-1" from the denominator in stndard scaling, this will not work
	//	auto DataSet = TestScaling("standard", "ExampStandard.dat", "scaledExampStandard.dat"); 
	//	//cout.precision(50);
	//	//cout << "std::get<2>(DataSet): \n" << std::get<2>(DataSet)<< endl;
	//	//cout << "std::get<3>(DataSet)" << std::get<3>(DataSet).transpose() << endl;
	//	EXPECT_TRUE(std::get<2>(DataSet).isApprox(std::get<3>(DataSet))) << "Wrong Standard scaling!!!!!!!!!!!!!!";; //the expected results from Python and my actuall results from my caling function are the same when the expected results are written with high precision into the file 
	//}
	//TEST(Scaling, MeanVecScaling) {//#11
	//	auto DataSet = TestScaling("standard", "ExampStandard.dat", "scaledExampStandard.dat");
	//	EXPECT_TRUE(abs(11.17775342465753- std::get<0>(DataSet)[0]) <1e-14) << "Min value is wrong!!";
	//}
	//TEST(Scaling, STDScaling) {//#12
	//	auto DataSet = TestScaling("standard", "ExampStandard.dat", "scaledExampStandard.dat");
	//	EXPECT_TRUE(abs(4.07127907531081- std::get<1>(DataSet)[0])<1e-14) << "Max value is wrong!!";
	//}

	///***************** Testing the condition number of fixed matrix ***********************/
	//TEST(ConditionNumber, HiggsDataSet) {
	//	double obtainedCond= ConditionTesting("higgs10k.dat");
	//	//cout << "condition number: " << obtainedCond << endl;
	//	EXPECT_TRUE(obtainedCond < 1e+10) << "The condition number is very large !!!";
	//}
	/****************** Testing Adjacency Matrices******************/
	///* positive testes*/
	//TEST(TestRegularity, 5Reg8Nodes_positive) {
	//	EXPECT_TRUE(testNetwork(8, 5,"5Regular.graph", "Regular1"));
	//}
	//TEST(TestRegularity, 7Reg8Nodes_positive) {
	//	EXPECT_TRUE(testNetwork(8, 7, "7Regular.graph", "Regular1"));
	//}
	///* Negative testes*/
	//TEST(TestRegularity, 5Reg8Nodes_negative) {
	//	EXPECT_FALSE(testNetwork(8, 3, "5Regular.graph", "Regular1"));
	//}
	//TEST(TestRegularity, 7Reg8Nodes_negative) {
	//	EXPECT_FALSE(testNetwork(8, 5, "7Regular.graph", "Regular1"));
	//}
	///********** Testing Eigenvalues of laplcaian ********/
	//TEST(testLaplacianEigenvalues, 7Reg8Nodes_positive) {
	//	EXPECT_TRUE(testEigenvalLaplac(8, "7Regular.graph", "Regular1"));
	//}
	//TEST(testLaplacianEigenvalues, 5Reg8Nodes_positive) {
	//	EXPECT_TRUE(testEigenvalLaplac(8, "5Regular.graph", "Regular1"));
	//}

	/****************** Testing Manual kernel with function kernels******************/
	//TEST(testKernel, nonlin) {
	//	EXPECT_TRUE(TestKernel("nonlin.dat",2));
	//}
	//TEST(testKernel, higgs) {
	//	EXPECT_TRUE(TestKernel("higgs10k.dat",100));
	//}
	/****************** Testing mFixedMatrices when two nodes have the same data******************/
	TEST(testFixedMatrix, nonlin) {
		EXPECT_TRUE(TestKernelForNonLin());
	}
	
}//END NAMESPACE
#endif
/**
2. test mfixedmatrix
3. teste kernel functions

**/