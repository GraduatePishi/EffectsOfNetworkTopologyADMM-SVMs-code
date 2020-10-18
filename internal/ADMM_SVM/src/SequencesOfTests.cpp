//#include "SequencesOfTests.h"
//#include "Types.h"
//#include "ChoosingT.h"
//#include "Kernel.h"
//#include "NodeSolver_linear.h"
//#include "NodeSolver.h"
//#include <iostream>
//#include "Partitioning.h"
//#include "FileReader.h"
//#include "Solver.h"
//#include "CrossValidation.h"
//#include "InverseOfMat.h"
//#include "NetworkCreation.h"
//#include "ReadDataset.h"
//#include "utils.h"
//#include "TestPhase.h"
//#include "Scaling.h"
//#include <math.h>
//#include "Randomizing.h"
//#include "FoldData.h"
//#include <algorithm>
//#include "ProblemSettings.h"
//#include "DistributedData.h"
//#ifndef _MSC_VER
//#ifndef MPIRUN
//#include "pstl/execution"
//#include "pstl/algorithm"
//namespace my_paral = pstl;
//#endif
//#else
//#include <execution>
//namespace my_paral = std;
//#endif
//using namespace std;
//
//struct TestProblem
//{
//	TestProblem(int iReducedSize , double iJcTest, double iEtaTest, double iGamma, bool iADMMParam, std::string igraphFile, std::string iDataname)
//	{
//		std::string iPath = ReadDataSetFile();
//		iPath.append(iDataname);
//		fs::path filePath = fs::u8path(iPath);
//		mTrainfilePath = filePath;
//		const auto prob = FileReader(filePath);
//		Data = prob.Data;
//		int rows = prob.Data.rows();
//		int cols = prob.Data.cols();
//		Data_lin = RowMajorMatirx(rows, cols + 1);
//		Label = prob.lables;
//		NoDataPoints = rows;
//		p_test = cols + 1;
//		L_test = iReducedSize;
//		JC_test = iJcTest;
//		EtaTest = iEtaTest;
//		wj_result= WeighVector::Ones(iReducedSize)*(-1000);
//		average_wi = WeighVector::Zero(iReducedSize);
//		initLambda = Multipliers::Ones(rows)*5.0;
//		initWj = WeighVector::Ones(iReducedSize)*5.0;
//		gamma = iGamma;
//		ADMMParam = iADMMParam;
//		graphFile = igraphFile;
//
//	}
//
//	RowMajorMatirx Data;
//	Graph graph = Graph::line;
//	std::string graphFile;
//	fs::path mTrainfilePath;
//	std::string mDataname;
//	RowMajorMatirx Data_lin;
//	LabelVector Label;
//	Row NoDataPoints;
//	int p_test{ -1 };
//	int L_test{-1};
//	int nrNodes{ -1 };
//	double JC_test{ -1 };
//	double EtaTest{ -1 };
//	double gamma{ -1 };
//	int mWorld_rank{ -1 }, mWorld_size{ -1 };
//	bool ADMMParam;
//	bool mShuffle{ false };
//	std::string mScaling{ "no" };
//	WeighVector wj_result;
//	Bias bj_result=-1000;
//	Multipliers initLambda;
//	WeighVector initWj;
//	Bias init_bj = 5;
//	WeighVector average_wi; 
//	Bias average_bi = 0.0;
//	Neighbors neighborIndexes;
//	Multipliers lambda;
//	TypeForDesign Ttype;// = TypeForDesign::identity;
//	DimReductionMat consensusT;
//	TypeForKernel kerneltype;
//	string name;
//	string testTerminationCriteria ="AllThree";
//};
//
//enum class ProblemToTest
//{
//	SimpleLinear2_first6Values,
//	SimpleLinear2,
//	Simple1,
//	fourclass1,
//	ForBroad,
//	splice20GradTest,
//	fourclass100GradTest, 
//	heart,
//	splice,
//	spliceForHitRate,
//	simpleLinear2_16,
//	neighborsCheck
//};
//
//auto createTestProblem(ProblemToTest iVal)//int type)
//{
//	if (ProblemToTest::SimpleLinear2_first6Values == iVal)
//	{
//		std::string name = "simpleLinear2_first6values.dat";
//		double JC_test = 200.0;//Best value for simpleLinear2 is eta=1 and JC=10, 200
//		double EtaTest = 1;
//		double gamma =0.0;
//		bool ADMMPar = false;
//		std::string GraphFolder = "Regular";
//		int Tsize = 5;
//		TestProblem retProblem(Tsize, JC_test, EtaTest, gamma, ADMMPar, GraphFolder, name);
//		retProblem.name = "simpleLinear2_first6values";
//		retProblem.graphFile = GraphFolder;
//		retProblem.p_test = 3;
//		retProblem.L_test = 2;
//		retProblem.lambda = Multipliers::Ones(retProblem.NoDataPoints);
//		retProblem.nrNodes = 1;
//		retProblem.wj_result << -1, -1;
//		retProblem.bj_result = 3;
//		retProblem.mScaling = "no";
//		retProblem.initLambda = Multipliers::Ones(retProblem.NoDataPoints)*5.0;
//		retProblem.initWj = WeighVector::Ones(retProblem.L_test)*5.0;
//		retProblem.init_bj = 5.0;
//		retProblem.kerneltype= TypeForKernel::linear;
//		retProblem.neighborIndexes = Network(Graph::hardCoded, retProblem.nrNodes, retProblem.graphFile,true);
//		//retProblem.neighborIndexes = getNeighbors(retProblem.nrNodes);
//		retProblem.consensusT = ConsensusT(retProblem.L_test, retProblem.Data, TypeForDesign::identity, EtaTest);
//		return retProblem;
//	}
//	//else if (ProblemToTest::SimpleLinear2 == iVal)
//	//{
//
//	//	std::string Dataset = ReadDataSetFile();
//	//	Dataset.append("simpleLinear2.dat");
//
//	//	fs::path filePath = fs::u8path(Dataset);
//	//	const auto prob = FileReader(filePath);
//	//	double JC_test = 200.0;
//	//	double EtaTest = 1.0;
//	//	double gamma = 0.5;
//	//	bool ADMMPar=false;
//	//	std::string GraphFolder = "Regular";
//	//	TestProblem retProblem((int)prob.Data.rows(), (int)prob.Data.cols(),(int)prob.Data.cols(), JC_test, EtaTest, gamma, ADMMPar, GraphFolder);
//	//	retProblem.graphFile = GraphFolder;
//	//	retProblem.Data = prob.Data;
//	//	retProblem.name = "simpleLinear2";
//	//	retProblem.Label = prob.lables;
//	//	retProblem.L_test = 2;
//	//	retProblem.lambda = Multipliers::Ones(retProblem.NoDataPoints);
//	//	retProblem.nrNodes = 1;
//	//	retProblem.wj_result << -1, -1;
//	//	retProblem.bj_result = 3;
//	//	retProblem.initLambda = Multipliers::Ones((int)prob.Data.rows())*5.0;
//	//	retProblem.initWj = WeighVector::Ones(retProblem.L_test)*5.0;
//	//	retProblem.init_bj = 5.0;
//	//	retProblem.kerneltype = TypeForKernel::linear;
//	//	retProblem.neighborIndexes = Network(Graph::hardCoded, retProblem.nrNodes, retProblem.graphFile, true); //getNeighbors(retProblem.nrNodes);
//	//	retProblem.Ttype = TypeForDesign::identity;
//	//	retProblem.consensusT = ConsensusT(retProblem.L_test, retProblem.Data, retProblem.Ttype, EtaTest);
//	//	retProblem.mWorld_size = 4;
//	//	return retProblem;
//	//}
//	//else if (ProblemToTest::simpleLinear2_16 == iVal)
//	//{
//
//	//	std::string Dataset = ReadDataSetFile();
//	//	Dataset.append("simpleLinear2_16.dat");
//
//	//	fs::path filePath = fs::u8path(Dataset);
//	//	const auto prob = FileReader(filePath);
//	//	double JC_test = 200.0;
//	//	double EtaTest = 1.0;
//	//	double gamma = 0.5;
//	//	bool ADMMPar = false;
//	//	std::string GraphFolder = "Regular";
//	//	TestProblem retProblem((int)prob.Data.rows(), (int)prob.Data.cols(), (int)prob.Data.cols(), JC_test, EtaTest, gamma, ADMMPar, GraphFolder);
//	//	retProblem.graphFile = GraphFolder;
//	//	retProblem.graph = Graph::cycle;
//	//	retProblem.Data = prob.Data;
//	//	retProblem.name = "simpleLinear2";
//	//	retProblem.Label = prob.lables;
//	//	retProblem.L_test = 4;
//	//	retProblem.lambda = Multipliers::Ones(retProblem.NoDataPoints);
//	//	retProblem.nrNodes = 4;
//	//	//retProblem.wj_result << -1, -1;
//	//	//retProblem.bj_result = 3;
//	//	retProblem.initLambda = Multipliers::Ones((int)prob.Data.rows())*5.0;
//	//	retProblem.initWj = WeighVector::Ones(retProblem.L_test)*5.0;
//	//	retProblem.init_bj = 5.0;
//	//	retProblem.kerneltype = TypeForKernel::linear;
//	//	retProblem.neighborIndexes = Network(retProblem.graph, retProblem.nrNodes, retProblem.graphFile, true); //getNeighbors(retProblem.nrNodes);
//	//	retProblem.Ttype = TypeForDesign::identity;
//	//	retProblem.consensusT = ConsensusT(retProblem.L_test, retProblem.Data, retProblem.Ttype, EtaTest);
//	//	return retProblem;
//	//}
//	//else if (ProblemToTest::splice20GradTest == iVal)
//	//{
//
//	//	std::string Dataset = ReadDataSetFile();
//	//	Dataset.append("splice20.dat");
//
//	//	fs::path filePath = fs::u8path(Dataset);
//	//	const auto prob = FileReader(filePath);
//	//	double JC_test = 10;
//	//	double EtaTest = 1;
//	//	double gamma = 0.25;
//	//	int reducedSizeofT = 2;
//	//	bool ADMMPar = false;
//	//	std::string GraphFolder = "Regular";
//	//	TestProblem retProblem((int)prob.Data.rows(), (int)prob.Data.cols(), reducedSizeofT, JC_test, EtaTest, gamma, ADMMPar, GraphFolder);
//	//	retProblem.graphFile = GraphFolder;
//	//	retProblem.Data = prob.Data;
//	//	retProblem.name = "splice20";
//	//	retProblem.Label = prob.lables;
//	//	retProblem.L_test = reducedSizeofT;
//	//	retProblem.nrNodes = 1;
//	//	retProblem.initLambda = Multipliers::Ones((int)prob.Data.rows())*5.0;
//	//	retProblem.initWj = WeighVector::Ones(retProblem.L_test)*5.0;
//	//	retProblem.init_bj = 5.0;
//	//	retProblem.kerneltype = TypeForKernel::rbf;
//	//	retProblem.Ttype = TypeForDesign::random;
//	//	retProblem.neighborIndexes = Network(Graph::hardCoded, retProblem.nrNodes, retProblem.graphFile, true); //getNeighbors(retProblem.nrNodes);
//	//	retProblem.consensusT = ConsensusT(retProblem.L_test, retProblem.Data, retProblem.Ttype, EtaTest);
//	//	return retProblem;
//	//}
//	//else if (ProblemToTest::fourclass100GradTest == iVal)
//	//{
//
//	//	std::string Dataset = ReadDataSetFile();
//	//	Dataset.append("fourclass100.dat");
//
//	//	fs::path filePath = fs::u8path(Dataset);
//	//	const auto prob = FileReader(filePath);
//	//	double JC_test = 10;
//	//	double EtaTest = 1;
//	//	double gamma = 0.25;
//	//	int reducedSizeofT = 2;
//	//	bool ADMMPar = false;
//	//	std::string GraphFolder = "Regular";
//	//	TestProblem retProblem((int)prob.Data.rows(), (int)prob.Data.cols(), reducedSizeofT, JC_test, EtaTest,gamma, ADMMPar, GraphFolder);
//	//	retProblem.graphFile = GraphFolder;
//	//	retProblem.Data = prob.Data;
//	//	retProblem.name = "fourclass100";
//	//	retProblem.Label = prob.lables;
//	//	retProblem.L_test = reducedSizeofT;
//	//	retProblem.nrNodes = 1;
//	//	retProblem.initLambda = Multipliers::Ones((int)prob.Data.rows())*5.0;
//	//	retProblem.initWj = WeighVector::Ones(retProblem.L_test)*5.0;
//	//	retProblem.init_bj = 5.0;
//	//	retProblem.kerneltype = TypeForKernel::rbf;
//	//	retProblem.Ttype = TypeForDesign::random;
//	//	retProblem.neighborIndexes = Network(Graph::hardCoded, retProblem.nrNodes, retProblem.graphFile, true); //getNeighbors(retProblem.nrNodes);
//	//	retProblem.consensusT = ConsensusT(retProblem.L_test, retProblem.Data, retProblem.Ttype, EtaTest);
//	//	return retProblem;
//	//}
//	//else if (ProblemToTest::Simple1==iVal)
//	//{
//
//	//	std::string Dataset = ReadDataSetFile();
//	//	Dataset.append("simple1.dat");
//
//	//	fs::path filePath = fs::u8path(Dataset);
//	//	const auto prob = FileReader(filePath);
//	//	double JC_test = 100;
//	//	double EtaTest = 1;
//	//	double gamma = 0.25;
//	//	int reducedSizeofT = 2;
//	//	bool ADMMPar = false;
//	//	std::string GraphFolder = "Regular";
//	//	TestProblem retProblem((int)prob.Data.rows(), (int)prob.Data.cols(), reducedSizeofT, JC_test, EtaTest, gamma, ADMMPar, GraphFolder);
//	//	retProblem.graphFile = GraphFolder;
//	//	retProblem.Data = prob.Data;
//	//	retProblem.name = "simple1";
//	//	retProblem.Label = prob.lables;
//	//	retProblem.L_test = reducedSizeofT;
//	//	retProblem.nrNodes = 1;
//	//	retProblem.wj_result << 0, 0;
//	//	retProblem.bj_result = 0;
//	//	retProblem.initLambda = Multipliers::Ones((int)prob.Data.rows())*5.0;
//	//	retProblem.initWj = WeighVector::Ones(retProblem.L_test)*5.0;
//	//	retProblem.init_bj = 5.0;
//	//	retProblem.kerneltype = TypeForKernel::rbf;
//	//	retProblem.Ttype = TypeForDesign::random;
//	//	retProblem.neighborIndexes = Network(Graph::hardCoded, retProblem.nrNodes, retProblem.graphFile, true); //getNeighbors(retProblem.nrNodes);
//	//	retProblem.consensusT= ConsensusT(retProblem.L_test, retProblem.Data, retProblem.Ttype, EtaTest);
//	//	return retProblem;
//	//}
//	//else if (ProblemToTest::fourclass1 == iVal)
//	//{
//
//	//	std::string Dataset = ReadDataSetFile();
//	//	Dataset.append("fourclass1.dat");
//
//	//	fs::path filePath = fs::u8path(Dataset);
//	//	const auto prob = FileReader(filePath);
//	//	double JC_test = 5;
//	//	double EtaTest = 0.01;
//	//	double gamma = 0.25;
//	//	int reducedSizeofT = 2;
//	//	bool ADMMPar = false;
//	//	std::string GraphFolder = "Regular";
//	//	TestProblem retProblem((int)prob.Data.rows(), (int)prob.Data.cols(), reducedSizeofT, JC_test, EtaTest, gamma, ADMMPar, GraphFolder);
//	//	retProblem.graphFile = GraphFolder;
//	//	retProblem.Data = prob.Data;
//	//	retProblem.name = "fourclass1";
//	//	retProblem.Label = prob.lables;
//	//	retProblem.L_test = reducedSizeofT;
//	//	retProblem.nrNodes = 1;
//	//	retProblem.initLambda = Multipliers::Ones((int)prob.Data.rows())*5.0;
//	//	retProblem.initWj = WeighVector::Ones(retProblem.L_test)*5.0;
//	//	retProblem.init_bj = 5.0;
//	//	retProblem.kerneltype = TypeForKernel::rbf;
//	//	retProblem.Ttype = TypeForDesign::random;
//	//	retProblem.neighborIndexes = Network(Graph::hardCoded, retProblem.nrNodes, retProblem.graphFile, true); //getNeighbors(retProblem.nrNodes);
//	//	retProblem.consensusT = ConsensusT(retProblem.L_test, retProblem.Data, retProblem.Ttype, EtaTest);
//	//	return retProblem;
//	//}
//	//else if (ProblemToTest::heart == iVal)
//	//{
//	//	std::string Dataset = ReadDataSetFile();
//	//	Dataset.append("heart.dat");
//
//	//	fs::path filePath = fs::u8path(Dataset);
//	//	const auto prob = FileReader(filePath);
//	//	double JC_test = 10;
//	//	double EtaTest = 1;
//	//	double gamma = 0.25;
//	//	int reducedSizeofT = 2;
//	//	bool ADMMPar = false;
//	//	std::string GraphFolder = "Regular";
//	//	TestProblem retProblem((int)prob.Data.rows(), (int)prob.Data.cols(), reducedSizeofT, JC_test, EtaTest, gamma, ADMMPar, GraphFolder);
//	//	retProblem.graphFile = GraphFolder;
//	//	retProblem.Data = prob.Data;
//	//	retProblem.name = "heart";
//	//	retProblem.Label = prob.lables;
//	//	retProblem.L_test = reducedSizeofT;
//	//	retProblem.nrNodes = 1;
//	//	retProblem.initLambda = Multipliers::Ones((int)prob.Data.rows())*5.0;
//	//	retProblem.initWj = WeighVector::Ones(retProblem.L_test)*5.0;
//	//	retProblem.init_bj = 5.0;
//	//	retProblem.kerneltype = TypeForKernel::rbf;
//	//	retProblem.Ttype = TypeForDesign::random;
//	//	retProblem.neighborIndexes = Network(Graph::hardCoded, retProblem.nrNodes, retProblem.graphFile, true); //getNeighbors(retProblem.nrNodes);
//	//	retProblem.consensusT = ConsensusT(retProblem.L_test, retProblem.Data, retProblem.Ttype, EtaTest);
//	//	return retProblem;
//	//}
//	//else if (ProblemToTest::splice == iVal)
//	//{
//	//	std::string Dataset = ReadDataSetFile();
//	//	Dataset.append("splice.dat");
//
//	//	fs::path filePath = fs::u8path(Dataset);
//	//	const auto prob = FileReader(filePath);
//	//	double JC_test = 100;
//	//	double EtaTest = 0.1;
//	//	double gamma = 0.025;
//	//	int reducedSizeofT = (int)prob.Data.rows() / 2;
//	//	bool ADMMPar = false;
//	//	std::string GraphFolder = "Regular";
//	//	TestProblem retProblem((int)prob.Data.rows(), (int)prob.Data.cols(), reducedSizeofT, JC_test, EtaTest, gamma, ADMMPar, GraphFolder);
//	//	retProblem.graphFile = GraphFolder;
//	//	retProblem.Data = prob.Data;
//	//	retProblem.graph = Graph::hardCoded;
//	//	retProblem.name = "splice";
//	//	retProblem.Label = prob.lables;
//	//	retProblem.L_test = reducedSizeofT;
//	//	retProblem.gamma = gamma;
//	//	retProblem.nrNodes = 1;
//	//	retProblem.initLambda = Multipliers::Ones((int)prob.Data.rows())*5.0;
//	//	retProblem.initWj = WeighVector::Ones(retProblem.L_test)*5.0;
//	//	retProblem.init_bj = 5.0;
//	//	retProblem.kerneltype = TypeForKernel::rbf;
//	//	retProblem.Ttype = TypeForDesign::random;
//	//	//retProblem.neighborIndexes = Network(Graph::hardCoded, retProblem.nrNodes, retProblem.graphFile, true); //getNeighbors(retProblem.nrNodes);
//	//	retProblem.consensusT = ConsensusT(retProblem.L_test, retProblem.Data, retProblem.Ttype, EtaTest);
//	//	return retProblem;
//	//}
//	//else if (ProblemToTest::spliceForHitRate == iVal)
//	//{
//	//	std::string Dataset = ReadDataSetFile();
//	//	Dataset.append("splice.dat");
//	//	fs::path filePath = fs::u8path(Dataset);
//	//	const auto prob = FileReader(filePath);
//	//	double JC_test = 100;
//	//	double EtaTest = 0.1;
//	//	double gamma = 0.025;
//	//	int reducedSizeofT = 10;
//	//	bool ADMMPar = false;
//	//	std::string GraphFolder = "Regular";
//	//	TestProblem retProblem((int)prob.Data.rows(), (int)prob.Data.cols(), reducedSizeofT, JC_test, EtaTest, gamma, ADMMPar, GraphFolder);
//	//	retProblem.graphFile = GraphFolder;
//	//	retProblem.Data = prob.Data;
//	//	retProblem.graph = Graph::self;
//	//	retProblem.name = "splice";
//	//	retProblem.Label = prob.lables;
//	//	retProblem.L_test = reducedSizeofT;
//	//	retProblem.gamma = gamma;
//	//	retProblem.nrNodes = 4;
//	//	retProblem.initLambda = Multipliers::Ones((int)prob.Data.rows())*5.0;
//	//	retProblem.initWj = WeighVector::Ones(retProblem.L_test)*5.0;
//	//	retProblem.init_bj = 5.0;
//	//	retProblem.kerneltype = TypeForKernel::rbf;
//	//	retProblem.Ttype = TypeForDesign::random;
//	//	//retProblem.neighborIndexes = Network(Graph::self, retProblem.nrNodes, GraphFolder, true); //getNeighbors(retProblem.nrNodes);
//	//	retProblem.consensusT = ConsensusT(retProblem.L_test, retProblem.Data, retProblem.Ttype, EtaTest);
//	//	return retProblem;
//	//}
//	//else {
//	//	throw ExceptionError("No Problem is chosen!!!");
//	//}
//}
//
///**************** Start of tests *******************/
//bool DistributedDataTest(std::string iDataName, int iWorld_size, int iWorld_rank) { // this is for 4 nodes
//
//	ProblemStatment iProb(0, 0); // this is for reading from file directly
//	ProblemStatment TestProb(0, 0);
//	std::string scaling = "no";
//	std::string Testset = ReadDataSetFile();
//	Testset.append(iDataName);
//	fs::path TestfilePath = fs::u8path(Testset);
//	bool shuffle = false, finishTest=false;
//	int LocalTest = 0;
//	iProb = FileReader(Testset);
//	auto DataPart = DistributedData(iWorld_size, iWorld_rank, TestProb, shuffle, scaling, TestfilePath);
//	int partSize = (int)(iProb.Data.rows() / iWorld_size);
//	std::vector<RowMajorMatirx> A(iWorld_size);
//	//if (world_rank==0) {
//	//	int mat = 0;
//	//	for (int i = 0; i < (int)iProb.Data.rows(); i += partSize) {
//	//		 A[mat] = iProb.Data.block(i, 0, partSize, (int)iProb.Data.cols());
//	//		/*cout << "A[" << mat << "]\n" << A[mat] << endl;*/
//	//		mat++;
//	//	}
//	//}
//	A[iWorld_rank] = iProb.Data.block(iWorld_rank*partSize, 0, partSize, (int)iProb.Data.cols());
//	//if (iWorld_rank ==1) {
//	//	A[1] = RowMajorMatirx::Identity(partSize, (int)iProb.Data.cols())*2.0;
//	//}
//	if (!A[iWorld_rank].isApprox(std::get<0>(DataPart))) {
//		cout << "A[" << iWorld_rank << "] and Data[" << iWorld_rank << "] is not the same!!!!" << endl;
//		LocalTest = 0;
//	}
//	else {
//		LocalTest = 1;
//	}
//	int GlobalTest = 0;
//
//	MPI_Allreduce(&LocalTest, &GlobalTest, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
//	if (GlobalTest != iWorld_size) {
//		if (iWorld_rank == 0) {
//			cout << "Distributing data for " << iDataName << " is NOT passed!!!!" << endl;
//			std::cout << "------------------------------------------------" << "\n";
//		}
//		finishTest= false;
//	}
//	else {
//		finishTest = true;
//		if (iWorld_rank == 0) {
//			cout << "Distributing data for " << iDataName << " is PASSED!!!!" << endl;
//			std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++" << "\n";
//		}
//	}
//
//
//	return finishTest;
//}
//
//
///***************************************************/
//void testOfCode() //Compares a linear and nonlinear(with linear kernel) with each other
//{
//	int world_size{ -1 }, world_rank{ -1 };
//	MPI_Init(NULL, NULL);
//	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
//	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
//
//	DistributedDataTest("simpleLinear2.dat", world_size, world_rank);
//	DistributedDataTest("svmguid.dat", world_size, world_rank);
//	DistributedDataTest("higgs10k.dat", world_size, world_rank);
//
//	MPI_Finalize();
//
//
//	//std::cout << endl;
//	//MPI_Finalize();
//}
