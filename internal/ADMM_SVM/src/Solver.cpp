#include "Solver.h"
#include "NodeSolver.h"
#include "utils.h"
#include <math.h>
#include <iostream>
#include <list>
#include "PrintData.h"
#include <chrono>
#include <ctime>
#include "TestPhase.h"
#include "MatrixVectorUtils.h"
#include "FindSVs.h"
#define NOMINMAX
#include <algorithm>
#include <array>
#include "ADMM_penalty.h"
#include "Printer.h"
#include <fmt/format.h>
#include "Printer.h"
#include "TestingPhase.h"
#ifndef _MSC_VER
#ifndef MPIRUN
#include "pstl/execution"
#include "pstl/algorithm"
namespace my_paral = pstl;
#endif
#else
#include <execution>

#include "stdio.h"
#include "stdlib.h"
#ifdef _WIN32
#include "C:\Program Files (x86)\Microsoft SDKs\MPI\Include\mpi.h"
#else
#include "mpi.h"
#endif



namespace my_paral = std;
#endif

using namespace std;
//using concurrency::parallel_for_each;

//void printToMatlab(const std::vector < std::tuple < DataMatrix, LabelVector, IndexOfNeighbors >>& partRes);//forwards function
void CreateMydirectory(std::string directoryName) {
	if (!fs::is_directory(directoryName) || !fs::exists(directoryName)) { // Check if src folder exists
		fs::create_directory(directoryName); // create src folder
	}
}



std::vector<int> PartitionSize(int group, int MatRow) {
	std::vector<int> GroupSize(group);
	int Perrow;
	if (MatRow % group != 0) {
		Perrow = (int)std::floor(MatRow / group);
		for (auto&val : GroupSize) val = Perrow;
		int m = 0;
		for (int i = 0; i < MatRow % group; i++) {
			GroupSize[m++] = GroupSize[m] + 1;
		};
	}
	else {
		Perrow = MatRow / group;
		for (auto&val : GroupSize) val = Perrow;
	}
	return GroupSize;
}

void MPI_Await(std::vector<MPI_Request>& iReq)
{
	int total = 0;
	while (total != iReq.size())
	{
		total = 0;
		for (int i = 0; i < iReq.size(); i++)
		{
			int test;
			MPI_Status stat;
			MPI_Test(&iReq[i], &test, &stat);
			total += test;
		}
	}
}

void Broadcast_wj_bj_MPI(std::vector<int> iLocalNeighVec, int iRank, std::pair<std::unique_ptr<nodesolver>, NodeData>& iNode)
{
	//MPI_Barrier(MPI_COMM_WORLD);
	Eigen::VectorXd tmpWj = iNode.first->getWTilde();
	double tmp_bj = iNode.first->getbj();
	int WiSize = (int)tmpWj.size();
	int LocalNeigSize = (int)iLocalNeighVec.size();
	std::vector<WeighVector> wi(LocalNeigSize, Eigen::VectorXd::Zero(WiSize)); //puting the w of neighbors in vector of vectors
	std::vector<double> bi(LocalNeigSize);

	WeighVector SumOfWi = Eigen::VectorXd::Zero(WiSize);
	double SumOfbi = 0;


	std::vector<MPI_Request> recReq1(LocalNeigSize), recReq2(LocalNeigSize);
	std::vector<MPI_Request> sendReq1(LocalNeigSize), sendReq2(LocalNeigSize);

	for (int i = 0; i < LocalNeigSize; i++)
	{
		MPI_Isend(&tmpWj[0], WiSize, MPI_DOUBLE, iLocalNeighVec[i], 0, MPI_COMM_WORLD, &sendReq1[i]);
		MPI_Isend(&tmp_bj, 1, MPI_DOUBLE, iLocalNeighVec[i], 0, MPI_COMM_WORLD, &sendReq2[i]);
	}
	MPI_Await(sendReq1);
	MPI_Await(sendReq2);
	for (int i = 0; i < LocalNeigSize; i++)
	{
		MPI_Irecv(&wi[i][0], WiSize, MPI_DOUBLE, iLocalNeighVec[i], 0, MPI_COMM_WORLD, &recReq1[i]);
		MPI_Irecv(&bi[i], 1, MPI_DOUBLE, iLocalNeighVec[i], 0, MPI_COMM_WORLD, &recReq2[i]);
	}
	MPI_Await(recReq1);
	MPI_Await(recReq2);

	for (int i = 0; i < LocalNeigSize; i++) {
		SumOfWi = SumOfWi + wi[i];
		SumOfbi = SumOfbi + bi[i];
	}
	iNode.second.SumOfWi = SumOfWi;
	iNode.second.SumOfbi = SumOfbi;
}

std::pair<IndexOfNeighbors, int> distributeNeighbors(svmTypes::Neighbors const& iNeighborsVectors, int world_size)
{
	IndexOfNeighbors localNeigVec = iNeighborsVectors[0]; //this is only for the master local neighbors
	int localNeigSize = (int)localNeigVec.size();
	int slaveSize = world_size - 1;
	std::vector<MPI_Request> sendReq1(slaveSize);
	for (int i = 0; i < slaveSize; i++) {
		int proc = i + 1;
		int sizeVec = (int)iNeighborsVectors[proc].size();
		MPI_Isend(&sizeVec, 1, MPI_INT, proc, 0, MPI_COMM_WORLD, &sendReq1[i]);
	}
	MPI_Await(sendReq1);
	for (int i = 0; i < slaveSize; i++) {
		int proc = i + 1;
		int sizeLocalNeigh = (int)iNeighborsVectors[proc].size();
		MPI_Isend(&iNeighborsVectors[proc][0], sizeLocalNeigh, MPI_INT, proc, 0, MPI_COMM_WORLD, &sendReq1[i]);
	}
	MPI_Await(sendReq1);
	return { localNeigVec,localNeigSize };
}

IndexOfNeighbors DistrNeighToAllNnodes(int iWorld_rank, Graph igraph, int iWorld_size, std::string iGraphFile, bool iCrossVal, double iJC, double iGamma, double iEta) {
	IndexOfNeighbors localNeigVec;
	int localNeigSize = -1;

	if (iWorld_rank == 0) {
		Neighbors neighborsVectors = Network(igraph, iWorld_size, iGraphFile, iCrossVal);
		l::log("...........................  Before Solving ........................... ");
		l::log("JC = {} , Gamma = {} , Eta = {} ", iJC, iGamma, iEta);
		std::tie(localNeigVec, localNeigSize) = distributeNeighbors(neighborsVectors, iWorld_size);
	}
	else {
		MPI_Status st;
		MPI_Recv(&localNeigSize, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &st);
		localNeigVec.resize(localNeigSize);
		MPI_Recv(&localNeigVec[0], localNeigSize, MPI_INT, 0, 0, MPI_COMM_WORLD, &st);
	}
	return localNeigVec;
}

//std::pair<Eigen::MatrixXd, Eigen::MatrixXd> getTestKernel(const RowMajorMatirx& iTestData, const DataMatrix& PartialDataMat, const  DimReductionMat& iT, KrnlMatrx iNormT_squared, KrnlMatrx iNorm_Xsquared, double iGamma) {
//
//	int EvaluationSize = (int)iTestData.rows();
//	int Datasize = (int)PartialDataMat.rows();
//	int ConsensusRow = (int)iT.rows();
//	Eigen::MatrixXd KrnlTestDatXj(EvaluationSize, Datasize);
//	Eigen::MatrixXd KrnlTestAndT(EvaluationSize, ConsensusRow);
//	Eigen::MatrixXd normEvalData_squared(EvaluationSize, EvaluationSize);
//	Eigen::MatrixXd norm_EvalDataXj(EvaluationSize, Datasize);
//	Eigen::MatrixXd norm_EvalDataT(EvaluationSize, ConsensusRow);
//
//	for (int i = 0; i < EvaluationSize; i++) {
//		/************ norm_squared(iTestData) *****************/
//		for (int j = i; j < EvaluationSize; j++) {
//			normEvalData_squared(i, j) = (double)iTestData.row(i).dot(iTestData.row(j));
//			normEvalData_squared(j, i) = normEvalData_squared(i, j);
//		}
//		/************ kernel(iTestData,Xj) *****************/
//		for (int j = 0; j < Datasize; j++) {
//			norm_EvalDataXj(i, j) = (double)iTestData.row(i).dot(PartialDataMat.row(j));
//			KrnlTestDatXj(i, j) = exp((normEvalData_squared(i, i) - 2 * norm_EvalDataXj(i, j) + iNorm_Xsquared(j, j))*(-iGamma));
//		}
//		/************ kernel(iTestData,iT) *****************/
//		for (int j = 0; j < ConsensusRow; j++) {
//			norm_EvalDataT(i, j) = (double)iTestData.row(i).dot(iT.row(j));
//			KrnlTestAndT(i, j) = exp((normEvalData_squared(i, i) - 2 * norm_EvalDataT(i, j) + iNormT_squared(j, j))*(-iGamma));
//		}
//	}
//
//	return { KrnlTestDatXj,KrnlTestAndT };
//}
std::vector<std::tuple<HyperPlaneScalars, Bias, HyperPlaneScalars, hitRate, TotalSolverTime, InnerSolverTime, CommunTime, Residuals, Residuals, Residuals, Residuals, Residuals, KrnlMatrx, KrnlMatrx, ItrTime>> nonlinear_svm::solve(int &iteration, const ProblemStatment& iData, const ProblemStatment& iEvaluation, const  DimReductionMat& iT, const TypeForKernel iKerneltype, const int TotalNumberofNodes, const double iJC, const double iEta, const double iGamma, bool RandomOn, Graph igraph, bool iADMM_update, std::string iGraphFile, int iworld_size, int iworld_rank, double iepsilonVal, double iABSTOL, double iRELTOL, std::string iTerminationCriteria, bool iCrossVal, string iName)
{
	auto startTotalTimeOfSolver = TimeType::now();
	TotalSolverTime Solve_time = 0.0;
	InnerSolverTime NLOPT_time = 0.0;
	CommunTime Communi = 0.0;
	CommunTime Communi3 = 0.0;
	CommunTime Communi4 = 0.0;
	CommunTime Communi5 = 0.0;

	Residuals primal_w = -1000000;
	Residuals primal_b = -1000000;
	Residuals dual_w = -1000000;
	Residuals dual_b = -1000000;
	Residuals AllNode_vs_iNode = -1000000;
	Residuals Vj_Vi{ -10000.0 };

	int world_rank = iworld_rank;
	int NumberofNodes;
	int Max_iter = -1;
	TimePrec iterationEndTime;
	TimePrec iterationEndTime_perNode;
	double iterationEndTime_perNode_TMP;
	double NLOPT_time_per_iter;

	//double iterationEndError_Vj_Vall;
	//double iterationEndError_Vj_Vi;
	//double iterationEndNLOPT_time;
	//int MaxNLOPT_iteration;
#ifdef MPIRUN
	NumberofNodes = 1;
#else
	my_paral::execution::parallel_policy policy = my_paral::execution::par;
	world_rank = -1;
	NumberofNodes = TotalNumberofNodes;
#endif
	//std::string FileName = fmt::format("Time_Accuracy_Parameters_for_{}_{}.time", igraph);
	//std::fstream TimeInfo;
	//if (iworld_rank == 0) {
	//	TimeInfo.open(FileName, std::fstream::out);
	//}
	std::vector<std::tuple<HyperPlaneScalars, Bias, HyperPlaneScalars, hitRate, TotalSolverTime, InnerSolverTime, CommunTime, Residuals, Residuals, Residuals, Residuals, Residuals, KrnlMatrx, KrnlMatrx, ItrTime>> ret_Coeff(NumberofNodes);
	try {
		std::get<3>(ret_Coeff[0]) = -20000;

		/* ABSTOL and RELTOL are from paper "Scalable Classifiers with ADMM and Transpose Reduction" */
		//double ABSTOL = 5e-3;// 1E-4 recomended by Boyd and , for termination criteria; (paper:Scalable Classifiers with ADMM	and Transpose Reduction)
		//double RELTOL = 1e-2; // 1e-2 recomended, for termination criteria;
		std::vector < std::tuple < DataMatrix, LabelVector>> partRes(NumberofNodes);
#ifdef MPIRUN
		get<0>(partRes[0]) = iData.Data;
		get<1>(partRes[0]) = iData.lables;
		int Datasize = (int)iData.Data.rows();
#else
		partRes = Partitioning(iData.Data, iData.lables, NumberofNodes, RandomOn, igraph); //divide the data and labels among nodes
#endif
		std::vector<std::pair<std::unique_ptr<nodesolver>, NodeData>> NodeInfo(NumberofNodes); // Saving the info of each node
		Neighbors neighborsVectors;
		int consensusRow = (int)iT.rows();
		IndexOfNeighbors localNeigVec;
		int localNeigSize = -1;

#ifdef MPIRUN
		localNeigVec = DistrNeighToAllNnodes(iworld_rank, igraph, iworld_size, iGraphFile, iCrossVal, iJC, iGamma, iEta);
		localNeigSize = (int)localNeigVec.size();

#else
		neighborsVectors = Network(igraph, NumberofNodes, iGraphFile);
#endif

#ifdef MPIRUN
		auto& partMat = get<0>(partRes[0]);
		auto& partLab = get<1>(partRes[0]);

		NodeInfo[0] = std::make_pair(std::make_unique<nodesolver>(partMat, partLab, iT, iKerneltype, localNeigVec, iJC, iEta, iGamma, iworld_rank),
			NodeData(consensusRow, (int)partMat.rows(), (int)partMat.cols())
		);
#else
		std::transform(policy, cbegin(partRes), cend(partRes), cbegin(neighborsVectors), begin(NodeInfo), [&](const auto& data, const auto& neighbors) {
			auto& partMat = get<0>(data);
			auto& partLab = get<1>(data);
			return std::make_pair(std::make_unique<nodesolver>(partMat, partLab, iT, iKerneltype, neighbors, iJC, iEta, iGamma),
				NodeData(consensusRow, (int)partMat.rows(), (int)partMat.cols())
			);
		});
#endif
		/************* Preparing for evaluation and saving some matrices for reusibility ***********************/
		int EvaluationSize = (int)iEvaluation.Data.rows();
		Eigen::MatrixXd normEvalData_squared(EvaluationSize, EvaluationSize), norm_EvalDataXj(EvaluationSize, Datasize), norm_EvalDataT(EvaluationSize, consensusRow);
		Eigen::MatrixXd Krnl_EvalDatXj(EvaluationSize, Datasize);
		Eigen::MatrixXd K_EvalTAndT(EvaluationSize, consensusRow);
#ifdef MPIRUN
		if (iKerneltype == TypeForKernel::rbf) {
			std::tie(Krnl_EvalDatXj, K_EvalTAndT) = getTestKernel(iEvaluation.Data, partMat, iT, NodeInfo[0].first->m_normT_squared, NodeInfo[0].first->m_norm_Xsquared, iGamma);
		}
		else {
			K_EvalTAndT = K(iEvaluation.Data, iT, iKerneltype, iGamma);
			Krnl_EvalDatXj = K(iEvaluation.Data, partMat, iKerneltype, iGamma);
		}
#else
		std::transform(policy, cbegin(partRes), cend(partRes), begin(TestInfo), [&](const auto& data) {
			return K(iTestDataColum, get<0>(data), iKerneltype, iGamma);
		});
#endif

		/**************** returning nomT and normXj from solver*******************/
		std::get<12>(ret_Coeff[0]) = NodeInfo[0].first->m_normT_squared;
		std::get<13>(ret_Coeff[0]) = NodeInfo[0].first->m_norm_Xsquared;

		CreateMydirectory(iName);
		stringstream Result_info_for_each_node;
		stringstream ErrorFile;
		ErrorFile << iName << "/" << igraph << "_" << TotalNumberofNodes << "N";
		CreateMydirectory(ErrorFile.str());

		Result_info_for_each_node << ErrorFile.str() << "/" << "Per_Node_Results";

		CreateMydirectory(Result_info_for_each_node.str());

		//string All_Result_info_string;
		//string Result_info_for_each_node_string;
		//std::vector<tuple <int, double, double, double, double, double, double, double, double, double, TimePrec >> All_Result_info_vector;

		///////////////////////////// Start ---- average of all nodes
		//std::string All_Result_Info = fmt::format("{}/All Max and AVG Erros for {}.error", ErrorFile.str(), igraph);
		//std::fstream All_Result_Info_File;

		//if (world_rank == 0) {
		//	All_Result_Info_File.open(All_Result_Info, std::fstream::out);
		//	All_Result_info_string = fmt::format("Iter Wp Bp Wd Bd max(Vall-Vj) max(Vj-Vi) Avg(Vall-Vj) Avg(Vj-Vi) Avg(time) MaxSolvIter SolverTime\n");
		//	All_Result_Info_File << All_Result_info_string;
		//}
		////////////////////////// End ---- average of all nodes /////////////////

		////////////// for each node
		//std::string Result_info_for_each_node_string;
		//std::string Result_info_for_each_node_Name = fmt::format("{}/ResultsForNode[{}] and {}.result", Result_info_for_each_node.str(), world_rank, igraph);
		//std::fstream Result_info_for_each_node_file;
		//Result_info_for_each_node_file.open(Result_info_for_each_node_Name, std::fstream::out);
		//Result_info_for_each_node_string = fmt::format("Iter Vall-Vj Vj-Vi Time_perNode NLOP_time NLOPT_time_perIter\n");
		//Result_info_for_each_node_file << Result_info_for_each_node_string;

		//////****************************************
		std::string All_Result_info_string2, Result_info_for_each_node_string2;
		std::string All_Result_Info2 = fmt::format("{}/All Max and AVG Erros for {}.error2", ErrorFile.str(), igraph);
		std::fstream All_Result_Info_File2;

		if (world_rank == 0) {
			All_Result_Info_File2.open(All_Result_Info2, std::fstream::out);
			All_Result_Info_File2 << fmt::format("Iter Wp Bp Wd Bd max(Vall-Vj) max(Vj-Vi) Avg(Vall-Vj) Avg(Vj-Vi) Avg(time) MaxSolvIter SolverTime\n");
			All_Result_Info_File2.flush();
		}

		std::string Result_info_for_each_node_Name2 = fmt::format("{}/ResultsForNode[{}] and {}.result2", Result_info_for_each_node.str(), world_rank, igraph);
		std::fstream Result_info_for_each_node_file2;
		Result_info_for_each_node_file2.open(Result_info_for_each_node_Name2, std::fstream::out);
		Result_info_for_each_node_file2 << fmt::format("Iter [Vall-Vj] [Vj-Vi] Time_perNode NLOP_time NLOPT_time_perIter\n");
		Result_info_for_each_node_file2.flush();
		//********************************************

		if (iCrossVal == true) {
			Max_iter = 100;
		}
		else {
			Max_iter = iteration;
		}

		TimePointType const solver_start = TimeType::now();

		for (int t = 0; t < Max_iter; t++) {
			if (world_rank == 0) l::log("********[ Node: {}, iteration: {} ]********", world_rank, t + 1);

#ifdef MPIRUN
			NodeInfo[0].second.Wi_old = NodeInfo[0].second.SumOfWi;
			NodeInfo[0].second.bi_old = NodeInfo[0].second.SumOfbi;
			NeighCount sizeNeigh = (int)NodeInfo[0].first->getNeigbors().size();

			/***************** lambda, w, and b calculation */
			auto start_nlopt_t = TimeType::now();

			MPI_Barrier(MPI_COMM_WORLD);
			int localSolverSuccess{ 0 }, GlbSolverSuccess{ 0 };
			localSolverSuccess = NodeInfo[0].first->compute_lambda_w_b(sizeNeigh, NodeInfo[0].second.SumOfbi, NodeInfo[0].second.SumOfWi, world_rank, t);
			MPI_Allreduce(&localSolverSuccess, &GlbSolverSuccess, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
			if (GlbSolverSuccess < iworld_size - 3 && iCrossVal == true) {
				if (iworld_rank == 0) l::log(">>> {} nodes failed OPT solver!!!!!", iworld_size - GlbSolverSuccess);
				return ret_Coeff;
			}
			int maxSolverIteration{ 0 };
			MPI_Allreduce(&NodeInfo[0].first->mSolverCounter, &maxSolverIteration, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

			NLOPT_time = NLOPT_time + CalcTime(start_nlopt_t, TimeType::now());// std::chrono::duration_cast<TimePrec>(tota_nlopt_time).count();
			NLOPT_time_per_iter = (CalcTime(start_nlopt_t, TimeType::now()));

			int AverageSvs = 0;
			int numSv = (int)getNumSV(NodeInfo[0].first->getLambdaj());
			MPI_Allreduce(&numSv, &AverageSvs, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#else
			for_each(policy, begin(NodeInfo), end(NodeInfo), [](auto& Node) {
				Node.second.Vi_old = Node.second.Vi;
				Node.second.ViWithEtas_old = Node.second.ViWithEtas;
				NeighCount sizeNeigh = (int)Node.first->getNeigbors().size();
				Node.first->compute_lambda_w_b(sizeNeigh, Node.second.bi, Node.second.Wi, -1); //WithEtas, Node.second.wiWithEtas);
			});
#endif
			auto startCommuni = TimeType::now();
			Eigen::VectorXd WjAvrgAllNodes = Eigen::VectorXd::Zero((int)iT.rows());
			double bjAvrgAllNodes = 0;
#ifdef MPIRUN
			MPI_Barrier(MPI_COMM_WORLD);
			///*broad casting the wj and bj to neighboors, line 7 in Algorithm 3 MoM - NDSVM Will be done in parallel version*/
			Broadcast_wj_bj_MPI(localNeigVec, world_rank, NodeInfo[0]); //.second.biWithEtas and second.wiWithEtas are updated here
#endif
			/************* GETTING THE AVERAGE OF ALL WJ and bj OF ALL PROCESSORS*/
#ifdef MPIRUN

			MPI_Allreduce(&(NodeInfo[0].first->getWTilde())[0], &WjAvrgAllNodes[0], consensusRow, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			MPI_Allreduce(&(NodeInfo[0].first->getbj()), &bjAvrgAllNodes, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

			WjAvrgAllNodes = WjAvrgAllNodes / TotalNumberofNodes;
			bjAvrgAllNodes = bjAvrgAllNodes / TotalNumberofNodes;
			Eigen::VectorXd VjAvrgAllNodes(consensusRow + 1);
			VjAvrgAllNodes << WjAvrgAllNodes, bjAvrgAllNodes;
#else
			for (auto& Node : NodeInfo) {
				Broadcast_wj_bj(Node, NodeInfo);
				WjAvrgAllNodes = WjAvrgAllNodes + Node.first->getWTilde();
				bjAvrgAllNodes = bjAvrgAllNodes + Node.first->getbj();
			}
			WjAvrgAllNodes = WjAvrgAllNodes / NumberofNodes;
			bjAvrgAllNodes = bjAvrgAllNodes / NumberofNodes;
#endif

			/****line 10 in Algorithm 3 MoM-NDSVM *********for both mpi and non_mpi version*/
			for (int i = 0; i < NumberofNodes; i++) {
				auto& Node = NodeInfo[i];
				Node.first->compute_alpha_Beta(Node.second.SumOfbi, Node.second.SumOfWi);

				int NoOfNeigbors = (int)Node.first->getNeigbors().size();

				Node.second.w_primal = (Node.first->getWTilde() - Node.second.SumOfWi / NoOfNeigbors).norm() / Node.first->getWTilde().norm();
				Node.second.b_primal = abs(Node.first->getbj() - Node.second.SumOfbi / NoOfNeigbors);
				Node.second.w_dual = ((Node.second.SumOfWi - Node.second.Wi_old) / NoOfNeigbors).norm() / (Node.second.SumOfWi / NoOfNeigbors).norm(); //Node.first->mEta*(Node.second.Wi - Node.second.Wi_old).norm() / Node.second.Wi.norm();
				Node.second.b_dual = abs(Node.first->mEta*((Node.second.SumOfbi - Node.second.bi_old) / NoOfNeigbors)); //abs(Node.first->mEta*(Node.second.bi - Node.second.bi_old));

				WeightsAndBias tmpVj(consensusRow + 1);
				tmpVj << Node.first->getWTilde(), Node.first->getbj();
				Node.second.Vall_vs_Vj = (VjAvrgAllNodes - tmpVj).norm() / tmpVj.norm();


				WeightsAndBias tmpVi(consensusRow + 1);
				tmpVi << Node.second.SumOfWi / NoOfNeigbors, Node.second.SumOfbi / NoOfNeigbors;
				Node.second.Vj_Vj = (tmpVj - tmpVi).norm() / tmpVj.norm();

#ifdef MPIRUN

				MPI_Allreduce(&Node.second.w_primal, &primal_w, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD); //BOYD page 51 https://web.stanford.edu/~boyd/papers/pdf/admm_distr_stats.pdf
				MPI_Allreduce(&Node.second.b_primal, &primal_b, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
				MPI_Allreduce(&Node.second.w_dual, &dual_w, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
				MPI_Allreduce(&Node.second.b_dual, &dual_b, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
				MPI_Allreduce(&Node.second.Vall_vs_Vj, &AllNode_vs_iNode, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
				MPI_Allreduce(&Node.second.Vj_Vj, &Vj_Vi, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif
			}

			Communi = Communi + CalcTime(startCommuni, TimeType::now());
			/***************************EVALUATION in every iteration of ADMM********/	//	/*************Checking the correct and wrong hits*****************/
			double AverageHitRate{ 0.0 };
			auto[a, b, c] = NodeInfo[0].first->compute_a_b_c(NodeInfo[0].second.SumOfbi, NodeInfo[0].second.SumOfWi);
			std::get<0>(ret_Coeff[0]) = a;
			std::get<1>(ret_Coeff[0]) = b;
			std::get<2>(ret_Coeff[0]) = c;
			std::get<3>(ret_Coeff[0]) = Testing(iEvaluation.lables, Krnl_EvalDatXj, K_EvalTAndT, a, b, c, false, world_rank, iCrossVal);
			MPI_Allreduce(&std::get<3>(ret_Coeff[0]), &AverageHitRate, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			AverageHitRate = AverageHitRate / TotalNumberofNodes;
			if (iCrossVal == true && AverageHitRate < 60 && t > 40) {
				if (iworld_rank == 0) {
					l::log("----------------------------");
					l::log("NOT 50% accuracy achieved!!!");
					l::log("----------------------------");
				}
				return ret_Coeff;
			}

			if (world_rank == 0) l::log(" maxNLOPT_iter: {}",  maxSolverIteration);

			double iterTime = CalcTime(solver_start, TimeType::now());

			iterationEndTime_perNode = (CalcDiffTime_TimePrec(solver_start, TimeType::now()));
			iterationEndTime_perNode_TMP = (CalcTime(solver_start, TimeType::now()));

			if (world_rank == 0) iterationEndTime = (CalcDiffTime_TimePrec(solver_start, TimeType::now()));

			for (int i = 0; i < NumberofNodes; i++) {
				if (NodeInfo[i].first->mSolverResult != " SUCCESS") l::log("The NLOPT solver does not succeed to solve the problem: {}", NodeInfo[i].first->mSolverResult);
			}

			Solve_time = Solve_time + CalcTime(startTotalTimeOfSolver, TimeType::now());// std::chrono::duration_cast<TimePrec>(timeTotalSolver).count();
#ifdef MPIRUN

			std::get<4>(ret_Coeff[0]) = Solve_time;
			std::get<5>(ret_Coeff[0]) = NLOPT_time;
			std::get<6>(ret_Coeff[0]) = Communi;// +Communi3 + Communi4 + Communi5;
			std::get<7>(ret_Coeff[0]) = NodeInfo[0].second.w_primal;
			std::get<8>(ret_Coeff[0]) = NodeInfo[0].second.b_primal;
			std::get<9>(ret_Coeff[0]) = NodeInfo[0].second.w_dual;
			std::get<10>(ret_Coeff[0]) = NodeInfo[0].second.b_dual;
			std::get<11>(ret_Coeff[0]) = NodeInfo[0].second.Vall_vs_Vj;
			std::get<14>(ret_Coeff[0]) = iterTime;
#endif
			/***************** for each node information **************/
			 std::string tmpPerNodeData= fmt::format("{} {} {} {} {} {}\n", t + 1, NodeInfo[0].second.Vall_vs_Vj, NodeInfo[0].second.Vj_Vj, iterationEndTime_perNode.count(), NLOPT_time, NLOPT_time_per_iter);
			 //Result_info_for_each_node_string += tmpPerNodeData;
			 /********* Writing in each iteration *****/
			Result_info_for_each_node_file2 << tmpPerNodeData;
			Result_info_for_each_node_file2.flush();

			double iterationEndError_Vj_Vall_allNodes{ 0.0 };
			double iterationEndError_Vj_Vi_allNodes{ 0.0 };
			double iterationEndError_time_allnode{ 0.0 };

			MPI_Allreduce(&NodeInfo[0].second.Vall_vs_Vj, &iterationEndError_Vj_Vall_allNodes, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			MPI_Allreduce(&NodeInfo[0].second.Vj_Vj, &iterationEndError_Vj_Vi_allNodes, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			MPI_Allreduce(&iterationEndTime_perNode_TMP, &iterationEndError_time_allnode, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			//MPI_Barrier(MPI_COMM_WORLD);
			///*********puting the average results in a file in folders****/
			if (world_rank == 0) {
				iterationEndError_Vj_Vall_allNodes=iterationEndError_Vj_Vall_allNodes / TotalNumberofNodes;
				iterationEndError_Vj_Vi_allNodes=iterationEndError_Vj_Vi_allNodes / TotalNumberofNodes;
				iterationEndError_time_allnode=iterationEndError_time_allnode / TotalNumberofNodes;
				std::string TMPdata = fmt::format("{} {} {} {} {} {} {} {} {} {} {} {}\n", t + 1, primal_w, primal_b, dual_w, dual_b, AllNode_vs_iNode, Vj_Vi, iterationEndError_Vj_Vall_allNodes, iterationEndError_Vj_Vi_allNodes, iterationEndError_time_allnode, maxSolverIteration, iterationEndTime.count());
				//All_Result_info_string += TMPdata;
				//All_Result_Info_File << All_Result_info_string;

				/********* Writing in each iteration *****/
				All_Result_Info_File2 << TMPdata;
				All_Result_Info_File2.flush();

			}
		}// End of MaxIteration

		if (world_rank == 0) {
			//std::string All_Result_Info = fmt::format("{}/All Max and AVG Erros for {}.error", ErrorFile.str(), igraph);
			//std::fstream All_Result_Info_File;

			//All_Result_Info_File.open(All_Result_Info, std::fstream::out);
			//All_Result_Info_File << fmt::format("Iter Wp Bp Wd Bd max(Vall-Vj) max(Vj-Vi) Avg(Vall-Vj) Avg(Vj-Vi) Avg(time) MaxSolvIter SolverTime\n");
			//All_Result_Info_File << All_Result_info_string;
			//All_Result_Info_File.close();

			/********* Writing in each iteration *****/
			All_Result_Info_File2.flush();
			All_Result_Info_File2.close();
		}

		//std::string Result_info_for_each_node_Name = fmt::format("{}/ResultsForNode[{}] and {}.result", Result_info_for_each_node.str(), world_rank, igraph);
		//std::fstream Result_info_for_each_node_file;
		//Result_info_for_each_node_file.open(Result_info_for_each_node_Name, std::fstream::out);
		//Result_info_for_each_node_file << fmt::format("Iter Vall-Vj Vj-Vi Time_perNode NLOP_time NLOPT_time_perIter\n");
		//Result_info_for_each_node_file << Result_info_for_each_node_string;
		//Result_info_for_each_node_file.close();

		Result_info_for_each_node_file2.flush();
		Result_info_for_each_node_file2.close();

		std::vector<double> SVtemp;
		for (int i = 0; i < Datasize; i++) {
			auto& tmpLAmbda = NodeInfo[0].first->getLambdaj();
			if (tmpLAmbda(i) > (1e-8)) {
				SVtemp.push_back(tmpLAmbda(i));
			}
		}
	} // End of Try
	catch (std::exception& ie)
	{
		cout << "Ended early in the SOLVER_c++!!!!" << "\n";
		l::err(ie.what());
		l::err("Ended early in the SOLVER_c++!!!!");
	}

	if (world_rank == 0) 	l::log("User-defined max iteratons ({}) is reached!!!!!!!!!!!!!!!!!!!!!!!!!!!!", Max_iter);
	return ret_Coeff;
}
