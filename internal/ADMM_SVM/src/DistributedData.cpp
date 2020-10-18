#include "DistributedData.h"
#include <iostream>
#include <vector>
#include "PrintData.h"
#include "ChoosingT.h"
#include "Scaling.h"
#include "Printer.h"

#ifdef _WIN32
#include "C:\Program Files (x86)\Microsoft SDKs\MPI\Include\mpi.h"
#else
#include "mpi.h"
#endif

Eigen::RowVectorXd normalize(Eigen::RowVectorXd vec, Eigen::RowVectorXd iMaxCoeffs, Eigen::RowVectorXd iScalingFactors) {
	for (int i = 0; i < vec.size(); i++) { // normalize each feature.
		vec[i] = (vec[i] - iMaxCoeffs[i]) / iScalingFactors[i];
	}
	return vec;
}
/******************************************/
std::tuple<RowMajorMatirx, RowMajorMatirx> PCA_Train(RowMajorMatirx traindata, int iReducedCol) {

	RowMajorMatirx transformedData((int)traindata.rows(), iReducedCol);
	// Calculate normalization coefficients (globals of type Eigen::VectorXf). 
	Eigen::RowVectorXd maxCoeffs = traindata.colwise().maxCoeff();
	Eigen::RowVectorXd minCoeffs = traindata.colwise().minCoeff();
	Eigen::RowVectorXd scalingFactors = maxCoeffs - minCoeffs;

	// For each datapoint.
	for (int i = 0; i < traindata.rows(); i++) { // Normalize each datapoint.
		(Eigen::RowVectorXd) traindata.row(i) = normalize((Eigen::RowVectorXd) traindata.row(i), maxCoeffs, scalingFactors);
	}

	// Mean centering data.
	Eigen::RowVectorXd featureMeans = traindata.colwise().mean();
	Eigen::RowVectorXd centered = traindata.rowwise() - featureMeans;

	// Compute the covariance matrix.
	RowMajorMatirx cov = centered.adjoint() * centered;
	cov = cov / (traindata.rows() - 1);

	Eigen::SelfAdjointEigenSolver<RowMajorMatirx> eig(cov);
	// Normalize eigenvalues to make them represent percentages.
	Eigen::RowVectorXd normalizedEigenValues = eig.eigenvalues() / eig.eigenvalues().sum();


	// Get the two major eigenvectors and omit the others.
	RowMajorMatirx evecs = eig.eigenvectors();
	RowMajorMatirx pcaTransform = evecs.rightCols(iReducedCol);


	// Map the dataset in the new two dimensional space.
	transformedData = traindata * pcaTransform;
	return { transformedData , pcaTransform };
}

std::tuple<int, int> PosNegRatio(LabelVector iLabl) {
	int negClass{ 0 }, posClass{ 0 };
	int Rows = (int)iLabl.rows();
	for (int ilabe = 0; ilabe < Rows; ilabe++) {
		if (std::signbit(iLabl[ilabe]) == true)
			negClass++;
		else if (std::signbit(iLabl[ilabe]) == false)
			posClass++;
		else
			throw ExceptionError("WRONG LABELS ARE FOUND!!");
	}
	if (posClass + negClass != Rows)
		std::cout << "WARNING: the pos and negatives count is not correct!!!" << "\n";
	return { posClass, negClass };
}
/****************************************/
//std::tuple<DataMatrix, LabelVector, Row, Column, MinVec, MaxVec, MinVec, MaxVec, DimReductionMat, ProblemStatment> DistributedData(int iWorld_size, int iWorldRank, ProblemStatment iProb, bool iShuffle, std::string iScaling, fs::path iTrainfilePath, int ReducedL, bool preDefT, TypeForDesign iReduceMatType, size_t iSeed) {
//std::tuple<DataMatrix, LabelVector, Row, Column, MinVec, MaxVec, MinVec, MaxVec, DimReductionMat, RowMajorMatirx, ProblemStatment> DistributedData(int iWorld_size, int iWorldRank, ProblemStatment iProb, bool iShuffle, std::string iScaling, fs::path iTrainfilePath, int ReducedL, bool preDefT, TypeForDesign iReduceMatType, size_t iSeed) {
std::tuple<ProblemStatment, DimReductionMat, ProblemStatment, MinVec, MaxVec > DistributedData(int iWorld_size, int iWorldRank, bool iShuffle, std::string iScaling, fs::path iTrainfilePath, int ReducedL, bool preDefT, TypeForDesign iReduceMatType, size_t iSeed) {

	DimReductionMat T_mat;
	ProblemStatment Prob_evalution(0, 0), GlobalProb(0, 0), LocalProb(0,0);
	RowMajorMatirx TransformaMat;
	bool GlobalT = true;
	bool PCA_flag = false;
	int TransformRow{ 0 }, TransformCol{ 0 };

	int ProblemColSize{ 0 }, ProblemRowSize{ 0 }, evalutionRows{ 0 }, evaluationCols{ 0 };
	MinVec MinScaling(0), MinsVec(0);
	MaxVec MaxScaling(0), MaxsVec(0);

	if (iWorldRank == 0) { //masterStart
		int slaveNumber = iWorldRank - 1;
		if (PCA_flag == true) {
			ProblemStatment iProp_temp = FileReader(iTrainfilePath);
			GlobalProb.Data = std::get<0>(PCA_Train(iProp_temp.Data, (int)iProp_temp.Data.cols() / 2));
			GlobalProb.lables = iProp_temp.lables;
			TransformaMat = std::get<1>(PCA_Train(iProp_temp.Data, (int)iProp_temp.Data.cols() / 2));
			TransformRow = (int)TransformaMat.rows();
			TransformCol = (int)TransformaMat.cols();
		}
		else {
			GlobalProb = FileReader(iTrainfilePath);

			/***************** Taking the %10 of training data for evaluation ****************/
			ProblemStatment iProb_eval_tmp = GlobalProb;
			ShuffleRandomly(iProb_eval_tmp.Data, iProb_eval_tmp.lables);
			evalutionRows = (int)round(GlobalProb.Data.rows()*0.1);
			if (evalutionRows > 1000) { //no more than 1000 samples to evaluate since for large problems this gets expensive
				evalutionRows = 1000;
			}
			evaluationCols = (int)GlobalProb.Data.cols();
			Prob_evalution.Data = iProb_eval_tmp.Data.block(0, 0, evalutionRows, evaluationCols);
			Prob_evalution.lables = iProb_eval_tmp.lables.block(0, 0, evalutionRows, 1);
		}

		if (iShuffle == true) {
			ShuffleRandomly(GlobalProb.Data, GlobalProb.lables); //Shuffeling Data and corresponding labels
			/* Shuffle data and put it in a file, I use this to create a new file that i can use for all d-regular graph, this should only run once for all graphs*/
			auto filepathTMP = iTrainfilePath;
			filepathTMP.replace_extension(".New");
			Eigen::MatrixXd tmpData((int)GlobalProb.Data.rows(), (int)GlobalProb.Data.cols() + 1);
			tmpData.col(0) = GlobalProb.lables;
			tmpData.block(0, 1, (int)GlobalProb.Data.rows(), (int)GlobalProb.Data.cols()) = GlobalProb.Data;
			PrintData(tmpData, filepathTMP.u8string());
		}
		ProblemColSize = (int)GlobalProb.Data.cols();
		ProblemRowSize = (int)GlobalProb.Data.rows();

		auto ratio = PosNegRatio(GlobalProb.lables);
		int posClass = std::get<0>(ratio);
		int negClass = std::get<1>(ratio);
		l::log("Training Data >>> sampls: {}, features: {} with {} (%{}) positives and {} (%{}) negatives", ProblemRowSize, ProblemColSize, posClass, (posClass * 100) / ProblemRowSize, negClass, (negClass * 100) / ProblemRowSize);
		MinScaling = Eigen::VectorXd::Zero(ProblemColSize);
		MaxScaling = Eigen::VectorXd::Zero(ProblemColSize);
		MinsVec = Eigen::VectorXd::Zero(ProblemColSize);
		MaxsVec = Eigen::VectorXd::Zero(ProblemColSize);
		for (int j = 0; j < ProblemColSize; j++) {
			MaxsVec[j] = GlobalProb.Data.col(j).maxCoeff();
			MinsVec[j] = GlobalProb.Data.col(j).minCoeff();
		}
		if (iScaling != "no") {
			ScalingTrain(MinScaling, MaxScaling, GlobalProb.Data, iScaling); //scale the datadets and save the scaling range in the Mins and MAxs vectors for scaling the test datasets
		}

		T_mat = DimReductionMat::Zero(ReducedL, ProblemColSize);
		//Eigen::FullPivLU<DimReductionMat> lu_decomp(GlobalProb.Data);
		//auto rank1 = lu_decomp.rank();
		//cout << "Global Data Rank: " << rank1 << endl;
		auto ratio_eval = PosNegRatio(Prob_evalution.lables);
		int posClass_eval = std::get<0>(ratio_eval);
		int negClass_eval = std::get<1>(ratio_eval);
		l::log("Evaluation >>> sampls: {}, features {} with {} (%{}) positives and {} (%{}) negatives ", Prob_evalution.Data.rows(), Prob_evalution.Data.cols(), posClass_eval, (posClass_eval * 100) / Prob_evalution.Data.rows(), negClass_eval, (negClass_eval * 100) / Prob_evalution.Data.rows() );

		if (preDefT == false && GlobalT == true) {
			l::log("Creating T from shuffled input Data .....");
			auto propTemp = GlobalProb;
			ShuffleRandomly(propTemp.Data, propTemp.lables);
			T_mat = propTemp.Data.block(0, 0, ReducedL, ProblemColSize);

			/**** Writting matrix T into the file ****/
			auto filePathout = iTrainfilePath;
			filePathout.replace_extension(".T");
			PrintData(T_mat, filePathout.u8string());
		}
		else if (preDefT == true) {
			l::log("Reading T from .T file .....");
			iTrainfilePath.replace_extension(".T");
			T_mat = ReadT(iTrainfilePath, ReducedL, ProblemColSize);
			//ScalingT(Mins, Maxs, T, problem.Scaling); //scaling T
		}
		else {
			l::log("Creating T from function randomly .....");
			T_mat = ConsensusT(ReducedL, GlobalProb.Data, iReduceMatType, MaxsVec, MinsVec, iSeed);
			/**** Writting matrix T into the file ****/
			auto filePathout = iTrainfilePath;
			filePathout.replace_extension(".T");
			PrintData(T_mat, filePathout.u8string());
		}

	} //end of master
	/************* Broadcasting the evaluation data ************/
	MPI_Bcast(&evalutionRows, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&evaluationCols, 1, MPI_INT, 0, MPI_COMM_WORLD);
	Prob_evalution.Data.resize(evalutionRows, evaluationCols);
	Prob_evalution.lables.resize(evalutionRows);
	MPI_Bcast(&Prob_evalution.Data(0, 0), evalutionRows*evaluationCols, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&Prob_evalution.lables[0], evalutionRows, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	/************* end of Broadcasting Data*****************/
	/************* Sending the transformation matrix from PCA for testing data *********/
	if (PCA_flag == true) {
		MPI_Bcast(&TransformRow, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&TransformCol, 1, MPI_INT, 0, MPI_COMM_WORLD);
		TransformaMat.resize(TransformRow, TransformCol);
		MPI_Bcast(&TransformaMat(0, 0), TransformRow*TransformCol, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	}
	MPI_Bcast(&ProblemColSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&ProblemRowSize, 1, MPI_INT, 0, MPI_COMM_WORLD);

	MinScaling.resize(ProblemColSize);
	MaxScaling.resize(ProblemColSize);

	MinsVec.resize(ProblemColSize);
	MaxsVec.resize(ProblemColSize);

	MPI_Bcast(&MinScaling[0], ProblemColSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&MaxScaling[0], ProblemColSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	MPI_Bcast(&MinsVec[0], ProblemColSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&MaxsVec[0], ProblemColSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	/*************** Partitioning the data into processors***************/
	std::vector<int>  PerCountMat(iWorld_size), PerCountLabel(iWorld_size), disp(iWorld_size), disp2(iWorld_size);
	PerCountLabel = PartSize(iWorld_size, ProblemRowSize);


	//std::copy(PerCountLabel.begin(), PerCountLabel.end(), std::ostream_iterator<int>(std::cout, "PartSize: "));

	PerCountMat = PerCountLabel;
	std::for_each(PerCountMat.begin(), PerCountMat.end(), [ProblemColSize](int &el) {el *= ProblemColSize; });

	disp[0] = 0;
	disp2[0] = 0;
	for (int i = 1; i < iWorld_size; i++) {
		disp[i] = disp[i - 1] + PerCountMat[i - 1];
		disp2[i] = disp2[i - 1] + PerCountLabel[i - 1];
	}
	int localMatRowSize = PerCountLabel[iWorldRank];

	LocalProb.Data.resize(localMatRowSize,ProblemColSize);
	LocalProb.lables.resize(localMatRowSize);
	MPI_Scatterv(GlobalProb.Data.data(), &PerCountMat[0], &disp[0], MPI_DOUBLE, &LocalProb.Data(0, 0), PerCountLabel[iWorldRank] * ProblemColSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Scatterv(GlobalProb.lables.data(), &PerCountLabel[0], &disp2[0], MPI_DOUBLE, &LocalProb.lables[0], PerCountLabel[iWorldRank], MPI_DOUBLE, 0, MPI_COMM_WORLD);


	/*********** T Matrix ******/
	T_mat.resize(ReducedL, ProblemColSize);
	MPI_Bcast(&T_mat(0, 0), ReducedL*ProblemColSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	/**************  Reading Testing data*************/
	//iTrainfilePath.replace_extension(".validation");

	//ProblemStatment TestProblem(0, 0);
	//if (PCA_flag == true) {
	//	ProblemStatment TestProblem_tmp = FileReader(iTrainfilePath);
	//	TestProblem.Data = TestProblem_tmp.Data *TransformaMat;
	//	TestProblem.lables = TestProblem_tmp.lables;
	//}
	//else {
	//	TestProblem = FileReader(iTrainfilePath);
	//}

	//if (iScaling != "no") {
	//	ScalingTest(MinScaling, MaxScaling, TestProblem.Data, iScaling);
	//}

	//auto ratioTest = PosNegRatio(TestProblem.lables);
	//int positiveClass_test = std::get<0>(ratioTest);
	//int negativeClass_test = std::get<1>(ratioTest);
	//if (iWorldRank == 0)
	//	std::cout << "Testing with " << TestProblem.Data.rows() << " smpls, " << TestProblem.Data.cols() << " feature, " << positiveClass_test << " (%" << ceil((positiveClass_test * 100) / (int)TestProblem.lables.size()) << ")" << " possitives, " << negativeClass_test << " (%" << ceil((negativeClass_test * 100) / (int)TestProblem.lables.size()) << ") negatives" << endl;

	////return { localsMatrix ,LocalLabels, ProblemRowSize,ProblemColSize, MinScaling ,MaxScaling , MinsVec ,MaxsVec ,T_mat ,TestProblem };
	////return { localsMatrix ,LocalLabels, ProblemRowSize,ProblemColSize, MinScaling ,MaxScaling , MinsVec ,MaxsVec ,T_mat,TestProblem ,iProb_evalution };
	//return { LocalProb, T_mat, TestProblem ,Prob_evalution };
	return { LocalProb, T_mat, Prob_evalution, MinScaling ,MaxScaling };
};