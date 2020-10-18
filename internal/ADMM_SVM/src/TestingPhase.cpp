#include "TestingPhase.h"
#include "TestPhase.h"
#include "Scaling.h"
#include "DistributedData.h"


std::pair<Eigen::MatrixXd, Eigen::MatrixXd> getTestKernel(const RowMajorMatirx& i_TestData, const DataMatrix& PartialDataMat, const  DimReductionMat& iT, KrnlMatrx iNormT_squared, KrnlMatrx iNorm_Xsquared, double iGamma) {
	DataMatrix iTestData= i_TestData;
	int testSize = (int)iTestData.rows();
	int Datasize = (int)PartialDataMat.rows();
	int ConsensusRow = (int)iT.rows();
	Eigen::MatrixXd KrnlTestDatXj(testSize, Datasize);
	Eigen::MatrixXd KrnlTestAndT(testSize, ConsensusRow);
	Eigen::MatrixXd normTestData_squared(testSize, testSize);
	Eigen::MatrixXd norm_tDataXj(testSize, Datasize);
	Eigen::MatrixXd norm_tDataT(testSize, ConsensusRow);

	for (int i = 0; i < testSize; i++) {
		/************ norm_squared(iTestData) *****************/
		for (int j = i; j < testSize; j++) {
			normTestData_squared(i, j) = (double)iTestData.row(i).dot(iTestData.row(j));
			normTestData_squared(j, i) = normTestData_squared(i, j);
		}
		/************ kernel(iTestData,Xj) *****************/
		for (int j = 0; j < Datasize; j++) {
			norm_tDataXj(i, j) = (double)iTestData.row(i).dot(PartialDataMat.row(j));
			KrnlTestDatXj(i, j) = exp((normTestData_squared(i, i) - 2 * norm_tDataXj(i, j) + iNorm_Xsquared(j, j))*(-iGamma));
		}
		/************ kernel(iTestData,iT) *****************/
		for (int j = 0; j < ConsensusRow; j++) {
			norm_tDataT(i, j) = (double)iTestData.row(i).dot(iT.row(j));
			KrnlTestAndT(i, j) = exp((normTestData_squared(i, i) - 2 * norm_tDataT(i, j) + iNormT_squared(j, j))*(-iGamma));
		}
	}

	return { KrnlTestDatXj,KrnlTestAndT };
}
/**************  Reading Testing data*************/
ProblemStatment ReadingTestData(fs::path iTrainfilePath,std::string iScaling, Eigen::VectorXd& iMinScaling, Eigen::VectorXd& iMaxScaling, bool iPCA_flag, int iRank) {
	ProblemStatment TestProblem(0, 0);
	iTrainfilePath.replace_extension(".test");
	TestProblem = FileReader(iTrainfilePath);

	if (TestProblem.Data.rows() == 0 || TestProblem.Data.cols() == 0) {
		l::log("Test data is wrong!!!!");
		throw ExceptionError("It does not read TESTING data (Matrix and labels)!!");
	}

	if (iScaling != "no") {
		ScalingTest(iMinScaling, iMaxScaling, TestProblem.Data, iScaling);
	}

	//if (iPCA_flag == true) {
	//	ProblemStatment TestProblem_tmp = FileReader(iTrainfilePath);
	//	TestProblem.Data = TestProblem_tmp.Data *iTransformaMat;
	//	TestProblem.lables = TestProblem_tmp.lables;
	//}
	//else {
	//TestProblem = FileReader(iTrainfilePath);
	//}



	auto ratioTest = PosNegRatio(TestProblem.lables);
	int positiveClass_test = std::get<0>(ratioTest);
	int negativeClass_test = std::get<1>(ratioTest);
	if (iRank == 0) {
		l::log("Testing with {} samples and {} features and {} (%{}) positives and {} (%{}) negatives samples.", TestProblem.Data.rows(), TestProblem.Data.cols(), positiveClass_test, ceil((positiveClass_test * 100) / (int)TestProblem.lables.size()), negativeClass_test, ceil((negativeClass_test * 100) / (int)TestProblem.lables.size()));
	}
	
	return TestProblem;
}
hitRate TestingPhase(const ProblemStatment& iTestProblem, const  Eigen::MatrixXd& iT, const DataMatrix& iDataMat,const KrnlMatrx& iT_norm, const KrnlMatrx& iX_norm, fs::path iTrainfilePath, std::string iScaling,Eigen::VectorXd& iMinScaling, Eigen::VectorXd&iMaxScaling, bool iPCA_flag, const  HyperPlaneScalars& ia,WeightScalar& ib, const  HyperPlaneScalars& ic, TypeForKernel iKernelType, int iRank, double iGamma) {
	/**************  Reading Testing data*************/

	//auto& TestProblem = ReadingTestData(iTrainfilePath, iTransformaMat, iScaling, iMinScaling, iMaxScaling, iPCA_flag, iRank);

	int testSize = (int)iTestProblem.Data.rows();
	int consensusRow = (int) iT.rows();
	int Datasize =(int) iDataMat.rows();
	Eigen::MatrixXd normTestData_squared(testSize, testSize), norm_tDataXj(testSize, Datasize), norm_tDataT(testSize, consensusRow);
	Eigen::MatrixXd Krnl_TestDatXj(testSize, Datasize);
	Eigen::MatrixXd  K_iTestAndT(testSize, consensusRow);

	if (iKernelType == TypeForKernel::rbf) {
		tie(Krnl_TestDatXj, K_iTestAndT)= getTestKernel(iTestProblem.Data, iDataMat, iT, iT_norm, iX_norm, iGamma);
	}
	else {
		K_iTestAndT = K(iTestProblem.Data, iT, iKernelType, iGamma);
		Krnl_TestDatXj = K(iTestProblem.Data, iDataMat, iKernelType, iGamma);
	}
	hitRate HitRate = Testing(iTestProblem.lables, Krnl_TestDatXj, K_iTestAndT, ia, ib, ic, false, iRank, false);
	return HitRate;
}