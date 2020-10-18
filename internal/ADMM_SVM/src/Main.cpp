#include "Solver.h"
#include "FileReader.h"
#include "SolverData.h"
#include "NodeSolver.h"
#include "NodeSolver_linear.h"
#include "DistributedData.h"
#include <Eigen/Dense>
//#include <random>
#include "PrintData.h"
#include <cmath>
#include "ChoosingT.h"
#include "Kernel.h"
#include "Types.h"
#include <string>
#include "TestPhase.h"
#include "SequencesOfTests.h"
#include "Partitioning.h"
#include "ProblemSettings.h"
#include "CrossValidation.h"
#include "NetworkCreation.h"
#include "ReadDataset.h"
#include "Randomizing.h"
#include "Scaling.h"
#include "fmt/format.h"
#include "Printer.h"
#include <iostream>
#include <fstream>
#include "utils.h"
#include "TestingPhase.h"
/* MPI implementation */
#include "stdio.h"
#include "stdlib.h"
#include "FoldData.h"

#ifdef _WIN32
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "C:\Program Files (x86)\Microsoft SDKs\MPI\Include\mpi.h"
#else
#include "mpi.h"
#endif

//#define RUNTESTS
#ifndef RUNTESTS

//std::tuple < DataMatrix, LabelVector> Partitioning_MPI(const RowMajorMatirx& _data, const LabelVector& _label, const int group, bool RandomOn) {
//	std::tuple < DataMatrix, LabelVector> Matrix_labl;
//	const int nrows = (int)_data.rows();
//	const int ncols = (int)_data.cols();
//	std::vector<int>  PerCountMat(group), PerCountLabel(group), disp(group), disp2(group);
//	PerCountLabel = PartSize(group, nrows);
//	PerCountMat = PerCountLabel;
//	std::for_each(PerCountMat.begin(), PerCountMat.end(), [ncols](int &el) {el *= ncols; });
//
//	disp[0] = 0;
//	disp2[0] = 0;
//	for (int i = 1; i < group; i++) {
//		disp[i] = disp[i - 1] + PerCountMat[i - 1];
//		disp2[i] = disp2[i - 1] + PerCountLabel[i - 1];
//	}
//	int maxCount = PerCountMat[0];
//	int maxCount2 = PerCountLabel[0];
//	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>  localsMatrix(maxCount2, ncols); //Row Major stored localmatricx for each node
//	Eigen::VectorXd LocalLabels(maxCount2);
//
//	MPI_Scatterv(_data.data(), &PerCountMat[0], &disp[0], MPI_DOUBLE, &localsMatrix(0, 0), maxCount, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//	MPI_Scatterv(&_label[0], &PerCountLabel[0], &disp2[0], MPI_DOUBLE, &LocalLabels[0], maxCount2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//	get<0>(Matrix_labl) = localsMatrix;
//	get<1>(Matrix_labl) = LocalLabels;
//	return Matrix_labl;
//}
void ScalingT(Eigen::VectorXd& iMean, Eigen::VectorXd& iSD, DimReductionMat& iData, string iType) {
	Eigen::VectorXd meanVec;
	if (iType == "standard") {
		for (int i = 0; i < (int)iData.cols(); i++) {
			meanVec = iMean[i] * Eigen::VectorXd::Ones(iData.rows());
			iData.col(i) = (iData.col(i) - meanVec) / iSD[i]; //(iData.col(i) / iSD[i]) - (meanVec / iSD[i]);
		}
	}
	else if (iType == "normal") {
		for (int i = 0; i < (int)iData.cols(); i++) {
			meanVec = iMean[i] * Eigen::VectorXd::Ones(iData.rows());
			auto MaxMin = (iSD[i] - iMean[i]);
			iData.col(i) = (iData.col(i) - meanVec) / MaxMin;
		}
	}
	else if (iType == "normalize") {
		for (int i = 0; i < (int)iData.cols(); i++) {
			iData.col(i) = iData.col(i) / iMean[i];
		}
	}
	else if (iType == "no") {
		//cout << " No scaling is requested!!" << "\n";
	}
	else
	{
		throw ExceptionError("Cant find scalingTrain type");
	}
}


int main(int argc, char* argv[]) {
	auto start_MasterTime = TimeType::now();
	auto startInitMPI = TimeType::now();
#ifdef MPIRUN
	bool MPIflag = true;
	MPI_Init(NULL, NULL); //initialize MPI
	int world_size{ -1 }, world_rank{ -1 };
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
#endif

	auto endInitMPI = TimeType::now();
	CommunTime InitMPI_global = CalcTime(startInitMPI, endInitMPI);
	try
	{
		auto startTotalTime = TimeType::now(); //Calculating time for the whole code
		auto startMiscTime1 = TimeType::now();
			/************* USER DEFINED PARAMETERS**************/
		std::string DataSet = "miRNAs";
		if (argc > 1) {
			DataSet = std::string(argv[1]);
		}
		DataSet.append(".dat");
		std::string Trainset;
		Trainset = ReadDataSetFile();
		std::string BoxPath = Trainset; //the address of all datasets
		Trainset.append(DataSet);
		fs::path TrainfilePath = fs::u8path(Trainset);

		auto settingPath = TrainfilePath;
		settingPath.replace_extension(".config");
		bool foundSettingFile = false;
		ProblemSetting problem(settingPath, foundSettingFile);

		if (argc > 3) {//overide file settings if number of nodes and grap type is specified
			problem.NameOfOutput = std::string(argv[1]);
			problem.NumberofNodes = atoi(argv[2]);
			problem.GraphType = std::string(argv[3]);
			problem.graphFile = std::string(argv[4]);
			problem.graph = returnGraph(problem.GraphType);//converts the user-defined graphtype string to graph enum
		}

		ProblemStatment prob(0, 0), probEvalutaion(0, 0);

		stringstream name;
		name << problem.graphFile << "_" << problem.NameOfOutput;

		int ProblemColSize{ 0 }, ProblemRowSize{ 0 };
		MinVec Min_Scaling(1), Mins(1);
		MaxVec Max_Scaling(1), Maxs(1);
		size_t seed{ 0 };

		if (problem.linearOn == true) {
			problem.kernelTyp = TypeForKernel::linear;
			problem.ReducedMatrixType = TypeForDesign::identity;
		}
		else {
			problem.kernelTyp = TypeForKernel::rbf;
			problem.ReducedMatrixType = TypeForDesign::random;
		}
#ifdef MPIRUN
		if (world_rank == 0) {
			l::log( "**********************************************************************" );
			l::log("------------------------ DISTRIBUTED MEMORY {} ------------------------   ", problem.NameOfOutput );
			l::log("Termination Criteria is {}", problem.TerminationCriteria );
			if (problem.TerminationCriteria == "allthree") {
				l::log("epsilon= {} , ABSTOL= {} , RELTOL= {}", problem.EpsilonVal , problem.ABSTOL , problem.RELTOL );
			}
			else if (problem.TerminationCriteria == "pri-dual") {
				cout << "ABSTOL= " << problem.ABSTOL << " , RELTOL= " << problem.RELTOL << endl;
			}
			else {
				cout << "epsilon= " << problem.EpsilonVal << endl;
			}
			cout << "Data scaling is : " << problem.Scaling << "\n";
			l::log( "...... Iteration defined {}", problem.iter, "....." );

			//********** creating a seed by master and sending to all threads fr creating matrix T
			auto RandSeed = std::random_device();
			seed = (size_t) RandSeed.entropy();
		}

		MPI_Bcast(&seed, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		DimReductionMat T(problem.ReducedRankL, ProblemColSize);
		double Misc1_Time = CalcTime(startMiscTime1, TimeType::now());

		auto startDstrDataTime = TimeType::now(); /*******  Calculating the Distributing data time of the code *****/
		MinVec MinScaling(0);
		MaxVec MaxScaling(0);
		//auto const& DataPart = DistributedData(world_size, world_rank, problem.Shuffle, problem.Scaling, TrainfilePath, problem.ReducedRankL, problem.preDefinedT, problem.ReducedMatrixType,  seed);
		std::tie(prob, T, probEvalutaion, MinScaling, MaxScaling) = DistributedData(world_size, world_rank, problem.Shuffle, problem.Scaling, TrainfilePath, problem.ReducedRankL, problem.preDefinedT, problem.ReducedMatrixType, seed);

		ProblemRowSize = (int)prob.Data.rows();// (int)std::get<0>(DataPart).Data.rows();
		ProblemColSize = (int)prob.Data.cols();//(int)std::get<0>(DataPart).Data.cols();
#else
		else {
			prob = FileReader(TrainfilePath); //Reading the file and put it in the Data matrix and Label vector
			ProblemColSize = (int)prob.Data.cols();
			ProblemRowSize = (int)prob.Data.rows();
		}
#endif
		/*******  Calculating the distrtibuting data time of the code *****/
		double DstrData_local = CalcTime(startDstrDataTime, TimeType::now());
		auto startMiscTime2 = TimeType::now();
#ifdef MPIRUN
		double DstrData_global{ 0 };
		MPI_Allreduce(&DstrData_local, &DstrData_global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#else
		TimeForAllNodes = Elapsed_time_for_entire_solver;
#endif
		/*****************************************************/
		//if (problem.SpilitData == true) {
		//	int rowPick = (int)(prob_pre.Data.rows() * 0.15);
		//	ShuffleRandomly(prob_pre.Data, prob_pre.lables); //Shuffeling Data and corresponding labels
		//	prob.Data = prob_pre.Data.block(0, 0, prob_pre.Data.rows() - rowPick, prob_pre.Data.cols());
		//	prob.lables = prob_pre.lables.block(0, 0, prob_pre.Data.rows() - rowPick, prob_pre.Data.cols());
		//	probTest.Data = prob_pre.Data.block(prob_pre.Data.rows() - rowPick , 0, rowPick, prob_pre.Data.cols());
		//	probTest.lables = prob_pre.lables.block(prob_pre.Data.rows() - rowPick , 0, rowPick, prob_pre.Data.cols());
		//	Eigen::MatrixXd tmp = Eigen::MatrixXd(rowPick, (int)prob_pre.Data.cols()+1);
		//	tmp.col(0) = probTest.lables;
		//	//cout << "tmp" << "\n" << tmp << "\n";
		//	tmp.block(0, 1, (int)probTest.Data.rows(), (int)probTest.Data.cols()) = probTest.Data;
		//	//cout << "tmp" << "\n"<< tmp << "\n";
		//	auto fileValidation = TrainfilePath;
		//	fileValidation.replace_extension(".validation");
		//	PrintData(tmp, fileValidation.u8string());
		//}
		//else{
		//if (world_rank==0) {
		//	cout << "prob Data" << "\n" << prob.Data.block(0, 0, 20, 18) << "\n";
		//	cout << "prob LAbel" << "\n" << prob.lables.block(0, 0, 20, 1) << "\n";

		//}

//auto TestfilePath = TrainfilePath;
//TestfilePath.replace_extension(".validation");
////ProblemStatment probTest_tmp(0, 0);
//probTest = FileReader(TestfilePath);
////probTest.Data = PCA(probTest_tmp.Data);
////probTest.lables = probTest_tmp.lables;
//if (problem.Scaling != "no") {
//	ScalingTest(Min_Scaling, Max_Scaling, probTest.Data, problem.Scaling);
//}

		if (world_rank == 0) {
			l::log( "Average Training smpls for each node: {}",ProblemRowSize);
			Eigen::FullPivLU<DimReductionMat> lu_decomp(T);
			auto rank = lu_decomp.rank();
			l::log( "Matrix T >>> Size=  {} and Rank: {}" ,problem.ReducedRankL, rank );
			l::log( "Number of Nodes are : {}", world_size );
			l::log("Graph Topology >>> {} from Graph folder: {} ", problem.GraphType , problem.graphFile );
		}


		if (problem.CrossValOn == false) {
			CreateMydirectory(name.str());
		}
		std::string Outputname = fmt::format("{}/{}_{}samples_{}Nodes_{}Graph_JC{}_Gamma{}_Eta{}_{}.output", name.str(), problem.NameOfOutput, ProblemRowSize, problem.NumberofNodes, problem.GraphType, problem.JC, problem.Gamma, problem.Eta, problem.graphFile);//writing the output file name based on the settings
		l::setStyle(l::Style::Console, Outputname);


		if (problem.CrossValOn == true) {
			if (world_rank == 0) {
				l::log( "Cross-Validation is running ..... " );
			}
			//PartData = FoldData(problem.NoFold, prob.Data, prob.lables, problem.RandomOn);//Dividing the whole Data matrix and label into training and testing matrices and labels and save them as tuple


			//auto tmpTuple/*[JC, Gamma, Eta, tmp, HitRate]*/ = Cross_validation(problem.NoFold, prob.Data, prob.lables, T, problem.kernelTyp, problem.NumberofNodes, problem.RandomOn, problem.Eta, problem.Gamma, problem.NameOfOutput, problem.Interval, problem.graph, problem.ADMMParaUpdate, problem.graphFile, world_size, world_rank, MPIflag, problem.EpsilonVal, problem.ABSTOL, problem.RELTOL, problem.TerminationCriteria);
			auto tmpTuple = Cross_validation(problem.NoFold, prob, probEvalutaion, T, problem.kernelTyp, problem.NumberofNodes, world_rank, problem.RandomOn, problem.Eta, problem.Gamma, problem.NameOfOutput, problem.Interval, problem.graph, problem.ADMMParaUpdate, problem.graphFile, world_size, world_rank, MPIflag, problem.EpsilonVal, problem.ABSTOL, problem.RELTOL, problem.TerminationCriteria, problem.CrossValOn, problem.GammaStep, problem.EtaStep, problem.JcStep);

			double AvrgHitRates_crossVal{ 0.0 }, CorrectHitRate{ -1 };
			MPI_Allreduce(&std::get<3>(tmpTuple), &AvrgHitRates_crossVal, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			AvrgHitRates_crossVal = AvrgHitRates_crossVal / world_size;
			problem.JC = std::get<0>(tmpTuple);
			problem.Gamma = std::get<1>(tmpTuple);
			problem.Eta = std::get<2>(tmpTuple);
			CorrectHitRate = std::get<3>(tmpTuple);
			l::log("After Cross-Validation ........................");
			l::log("JC = {} , Gamma = {} , Eta = {} , CorrectHitRate = {}, ID = {} ", problem.JC, problem.Gamma, problem.Eta, CorrectHitRate, world_rank);
		}
		else {
			nonlinear_svm solver;

			if (ProblemRowSize == 0 || ProblemColSize == 0)
				throw ExceptionError("It does not read Training Matrix and labels!!");

			double Misc2_Time = CalcTime(startMiscTime2, TimeType::now());

			auto startTotalSolver = TimeType::now();
			int iter=problem.iter;
			auto solving = solver.solve(iter, prob, probEvalutaion, T, problem.kernelTyp, problem.NumberofNodes, problem.JC, problem.Eta, problem.Gamma, problem.RandomOn, problem.graph, problem.ADMMParaUpdate, problem.graphFile, world_size, world_rank, problem.EpsilonVal, problem.ABSTOL, problem.RELTOL, problem.TerminationCriteria, problem.CrossValOn, problem.NameOfOutput);
			auto endTotalSolver = TimeType::now();
			double Solver_solveTime_local = CalcTime(startTotalSolver, endTotalSolver);
			double Solver_solveTime_global{ 0.0 }, Solver_solveTime_global2{ 0.0 };

			auto startMiscTime3 = TimeType::now();
			if (world_rank == 0) {
				l::log( "------------- Stats AFTER the solver-------------" );
				l::log("Problem >>> {} with {} samples and {} features. ", problem.NameOfOutput , ProblemRowSize ,  ProblemColSize );
				l::log("JC = {}, Gamma = {}, Eta = {},  Reduced Dimension to = {}" ,problem.JC , problem.Gamma , problem.Eta , problem.ReducedRankL );
				l::log("Graph Type = {},  ADMM parameter update is: {}, , Nodes: {}, Nr of iterations: {}.", problem.GraphType, problem.ADMMParaUpdate , problem.NumberofNodes, iter + 1 );
				l::log("Shuffling data= {}, Scaling Type= {}, predefined T : {}",problem.Shuffle , problem.Scaling , problem.preDefinedT );
				l::log("------------- Training SVM Finished --------------------------" );
			}

			double MiscTime_global{ 0.0 };
			double MiscTime_local = Misc1_Time + Misc2_Time;
			/*******  Calculating the total time of the code *****/
			auto endTotalTime = TimeType::now();
			double TrainingTime_global{ 0.0 }, TrainingTime_global2{ 0.0 };
			double TrainingTime_local = CalcTime(startTotalTime, endTotalTime);
			/*****************************************************/

			TotalSolverTime totalSolve_time{ 0.0 }, totalSolve_time2{ 0.0 };
			InnerSolverTime totalNLOPT_time{ 0.0 }, totalNLOPT_time2{ 0.0 };
			CommunTime totalCommun_time{ 0.0 }, totalCommun_time2{ 0.0 };
			hitRate AverageHitRates = 0.0;
			Residuals primal_r = std::numeric_limits<double>::max();
			Residuals primal_s = std::numeric_limits<double>::max();
			Residuals dual_r = std::numeric_limits<double>::max();
			Residuals dual_s = std::numeric_limits<double>::max();
			Residuals AllNode_vs_iNode = std::numeric_limits<double>::max();
#ifdef MPIRUN
			//MPI_Allreduce(&std::get<3>(solving[0]), &AverageHitRates, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); //This is hitRate for evaluation Data ,not Testing data
			MPI_Allreduce(&std::get<4>(solving[0]), &totalSolve_time, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
			MPI_Allreduce(&std::get<5>(solving[0]), &totalNLOPT_time, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
			MPI_Allreduce(&std::get<6>(solving[0]), &totalCommun_time, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

			MPI_Allreduce(&std::get<4>(solving[0]), &totalSolve_time2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			MPI_Allreduce(&std::get<5>(solving[0]), &totalNLOPT_time2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			MPI_Allreduce(&std::get<6>(solving[0]), &totalCommun_time2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

			MPI_Allreduce(&std::get<7>(solving[0]), &primal_r, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
			MPI_Allreduce(&std::get<8>(solving[0]), &primal_s, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
			MPI_Allreduce(&std::get<9>(solving[0]), &dual_r, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
			MPI_Allreduce(&std::get<10>(solving[0]), &dual_s, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
			MPI_Allreduce(&std::get<11>(solving[0]), &AllNode_vs_iNode, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

			MPI_Allreduce(&TrainingTime_local, &TrainingTime_global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
			MPI_Allreduce(&Solver_solveTime_local, &Solver_solveTime_global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

			MPI_Allreduce(&TrainingTime_local, &TrainingTime_global2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			MPI_Allreduce(&Solver_solveTime_local, &Solver_solveTime_global2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

			MPI_Allreduce(&MiscTime_local, &MiscTime_global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

#else
			totalSolve_time = std::get<4>(solving[0]);
			totalNLOPT_time = std::get<5>(solving[0]);
			totalCommun_time = std::get<6>(solving[0]);
			TimeForAllNodes = Elapsed_time_for_entire_solver;
#endif
			/*********************  Writing model parameters in file for each node *******/

			HyperPlaneScalars A_vec(ProblemColSize);
			WeightScalar B_bias;
			HyperPlaneScalars C_vec(problem.ReducedRankL);

			MPI_Allreduce(&std::get<0>(solving[0])[0], &A_vec[0], ProblemColSize, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			MPI_Allreduce(&std::get<1>(solving[0]), &B_bias, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			MPI_Allreduce(&std::get<2>(solving[0])[0], &C_vec[0], problem.ReducedRankL, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			if (world_rank == 0) {
				stringstream ErrorFile;
				ErrorFile << problem.NameOfOutput << "/" << problem.GraphType << "_" << world_size << "N";
				std::string ModelFileName = fmt::format("{}/Average_Model_Parameters_for_{}_{}.model", ErrorFile.str(), problem.NameOfOutput, problem.GraphType);
				std::fstream modelInfo;
				modelInfo.open(ModelFileName, std::fstream::out);
				modelInfo << "a:" << endl << A_vec / problem.NumberofNodes << endl;
				modelInfo << "b:" << endl << B_bias / problem.NumberofNodes << endl;
				modelInfo << "c:" << endl << C_vec / problem.NumberofNodes << endl;
				modelInfo.close();
			}
			/********************** Testing after recieved model *******************/
			double Testing_global{ 0.0 };
			auto startTesting = TimeType::now();
			ProblemStatment TestProblem(0, 0);
			bool PCA_flag = false;
			TestProblem = ReadingTestData(TrainfilePath, problem.Scaling, MinScaling, MaxScaling, PCA_flag, world_rank);

			hitRate hit = TestingPhase(TestProblem, T, prob.Data, std::get<12>(solving[0]), std::get<13>(solving[0]), TrainfilePath, problem.Scaling, Min_Scaling, Max_Scaling, PCA_flag, std::get<0>(solving[0]), std::get<1>(solving[0]), std::get<2>(solving[0]), problem.kernelTyp, world_rank, problem.Gamma);
			MPI_Allreduce(&hit, &AverageHitRates, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

			double Testing_local = CalcTime(startTesting, TimeType::now());

#ifdef MPIRUN
			MPI_Allreduce(&Testing_local, &Testing_global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#else
			Testing_global = Testing_local;
#endif
			/********************** End of testing *****************************/
			/********************** end of writing parameters***************************/

			MPI_Barrier(MPI_COMM_WORLD);
			auto end_MasterTime = TimeType::now();
			double MasterTime = 0.0;

			if (world_rank == 0) {
				MasterTime = CalcTime(start_MasterTime, end_MasterTime);

				AverageHitRates = AverageHitRates / problem.NumberofNodes;
				TrainingTime_global = (TrainingTime_global / 1000.0);
				totalSolve_time = (totalSolve_time / 1000.0);
				totalNLOPT_time = (totalNLOPT_time / 1000.0);
				totalCommun_time = (totalCommun_time + InitMPI_global) / 1000.0;
				//totalCommun_time = (totalCommun_time / 1000.0);
				DstrData_global = (DstrData_global / 1000.0);
				MiscTime_global = (MiscTime_global / 1000.0);
				Testing_global = (Testing_global / 1000.0);
				Solver_solveTime_global = (Solver_solveTime_global / 1000);
				MasterTime = MasterTime/1000.0;

				TrainingTime_global2 = (TrainingTime_global2 / 1000.0);
				Solver_solveTime_global2 = (Solver_solveTime_global2 / 1000.0);
				totalSolve_time2 = (totalSolve_time2 / 1000.0);
				totalNLOPT_time2 = (totalNLOPT_time2 / 1000.0);
				totalCommun_time2 = (totalCommun_time2 + InitMPI_global) / 1000.0;

				l::log( "_________________ Average Time for all Nodes for entire solver: {} _________________", TrainingTime_global );
				std::string FinalName = fmt::format("{}/{}_{}_{}_n {}_{}_MPI_itr {}_MasterT {}_TotTime {}_TotSolverTime {}_SlvT {}_DstrData {}_Test {}_MscT {}_NLOPT {}_CommuT {}.final", name.str(), problem.graphFile, problem.NameOfOutput, ProblemRowSize, world_size, problem.GraphType, iter + 1, MasterTime, TrainingTime_global, Solver_solveTime_global, totalSolve_time, DstrData_global, Testing_global, MiscTime_global, totalNLOPT_time, totalCommun_time);//writing the output file name based on the settings
				ofstream FinalResultName;
				FinalResultName.open(FinalName);
				FinalResultName << iter + 1 << " " << MasterTime << " " << TrainingTime_global << " " << Solver_solveTime_global << " " << totalSolve_time << " " << DstrData_global << " " << Testing_global << " " << MiscTime_global <<  " " << totalNLOPT_time << " " << totalCommun_time <<   " \n";// << primal_r << " " << primal_s << " " << dual_r << " " << dual_s << " " << AllNode_vs_iNode << "\n";
				FinalResultName << iter + 1 << " " << TrainingTime_global2 << " " << Solver_solveTime_global2 << " " << totalSolve_time2 << " " << totalNLOPT_time2 << " " << totalCommun_time2 << " \n";

				FinalResultName.close();
				l::log( "************ Residuals ***********" );
				std::cout.precision(17);
				l::log( "Prim_Res[{},{}] <= {}" , primal_r , primal_s , problem.ABSTOL );
				l::log("Dual_Res[{},{}] <= {}" , dual_r , dual_s , problem.RELTOL );
				l::log("Vall_Vj_Res[{}] <= {}", AllNode_vs_iNode , problem.EpsilonVal );
				l::log("MAX: {} {} {} {} {} {} {} {} {} {} ",iter + 1 , MasterTime, TrainingTime_global , Solver_solveTime_global, totalSolve_time, DstrData_global , Testing_global, MiscTime_global , totalNLOPT_time , totalCommun_time );// << primal_r << " " << primal_s << " " << dual_r << " " << dual_s << " " << AllNode_vs_iNode << "\n";
				l::log("SUM:  {} {} {} {} {} {}", iter + 1, TrainingTime_global2, Solver_solveTime_global2, totalSolve_time2, totalNLOPT_time2, totalCommun_time2);
				l::log("Itr: {},MstrTime:{}, total: {}, slvr_main: {}, Solver: {}, Dstr: {}, test: {}, misc: {}, NLopt: {}, Commu: {}", iter + 1 , MasterTime, TrainingTime_global, Solver_solveTime_global , totalSolve_time, DstrData_global , Testing_global, MiscTime_global ,  totalNLOPT_time , totalCommun_time  );// << primal_r << " " << primal_s << " " << dual_r << " " << dual_s << " " << AllNode_vs_iNode << "\n";

			}
		}
	}//end of Try
	catch (std::exception& ie)
	{
		cout << "Ended early because of exception!!!!" << "\n";
		l::err(ie.what());
		l::err("Ended early because of exception");
	}
#ifdef MPIRUN
	MPI_Finalize();
#endif
	return 0;
}
#else

GTEST_API_ int main(int argc, char **argv) {
	int result = 0;

	testing::InitGoogleTest(&argc, argv);
	//testing::InitGoogleMock(&argc, argv);
	MPI_Init(&argc, &argv);
	int world_size{ -1 }, world_rank{ -1 };
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	if (world_rank == 0)
		printf("Running main() from %s\n", __FILE__);
	testing::TestEventListeners& listeners =
		::testing::UnitTest::GetInstance()->listeners();
	if (world_rank != 0) {
		delete listeners.Release(listeners.default_result_printer());
	}
	result = RUN_ALL_TESTS();

	MPI_Finalize();
	return result;

}
#endif // !RUNTESTS