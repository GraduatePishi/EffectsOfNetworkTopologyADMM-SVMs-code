#include "TestPhase.h"
#include "Printer.h"

hitRate Testing(const Eigen::VectorXd& iTestLabel, const Eigen::MatrixXd& iFirst,  const Eigen::MatrixXd& iSecond,const  HyperPlaneScalars& ia, WeightScalar& ib, const HyperPlaneScalars& ic, bool iValidFlag, int iRank, bool iCrossVal) {
	int TestSize =  (int)iTestLabel.size();
	/******** equ (31) for the local discriminant function ****************/
	Eigen::VectorXd discriminantFunc = (iFirst * ia) + (iSecond * ic) + (ib*Eigen::VectorXd::Ones(TestSize));
	///*************Checking the correct and wrong hits*****************/
	int nrCorrect = 0 ;
	hitRate correctHit = -2;
	for (int j = 0; j <(int) discriminantFunc.size(); j++) {
		if (std::signbit(discriminantFunc(j)) == std::signbit(iTestLabel(j))) {
			nrCorrect++;
		}
	}
	correctHit = (100 * nrCorrect) / TestSize;

	/************************************************/
		double positiveClass = 0.0, negativeClass = 0.0, TruePos = 0.0, TrueNeg = 0.0, FalsePos = 0.0, FalseNeg = 0.0, Sensivity = 0.0, Precision = 0.0, Accuracy = 0.0;
		double Specificity=0.0, NegPerdicRate=0.0;

		for (int labelind = 0; labelind < TestSize; labelind++) {
			//std::cout << "element: " << labelind << " ,estim: " << std::signbit(iTestLabel[labelind]) << ", actual label: " << std::signbit(discriminantFunc[labelind])<< std::endl;
			if (std::signbit(iTestLabel[labelind]) == false) {
				positiveClass++;
				if (std::signbit(discriminantFunc[labelind]) == false)
					TruePos++;
				else
					FalseNeg++;
			}
			else if (std::signbit(iTestLabel[labelind]) == true) {
				negativeClass++;
				if (std::signbit(discriminantFunc[labelind]) == true)
					TrueNeg++;
				else
					FalsePos++;
			}
			else
				throw ExceptionError("WRONG LABELS ARE FOUND!!");
		}
		Sensivity = (TruePos / (TruePos + FalseNeg));
		Precision = (TruePos / (TruePos + FalsePos));
		Specificity=TrueNeg/negativeClass;
		NegPerdicRate=TrueNeg/(TrueNeg+FalseNeg);
		Accuracy = ((TruePos + TrueNeg) / (TruePos + FalsePos + FalseNeg + TrueNeg));
		if ((TruePos + FalseNeg) != positiveClass) {
			l::log("[{}]_Pos {}+Neg{} != Tot{}", iRank, TruePos, FalseNeg, positiveClass);
			throw ExceptionError("P !=Tp+FN");
		}
		if ((TrueNeg + FalsePos) != negativeClass) {
			l::log("[{}]_Pos {}+Neg{} != Tot{}", iRank, TrueNeg, FalsePos, negativeClass);
			throw ExceptionError(" N!=TN+FP");
		}
		if (positiveClass + negativeClass != iTestLabel.size()) {
			l::log("[{}]_Pos {}+Neg{} != Tot{}", iRank, positiveClass, negativeClass, TestSize);
			throw ExceptionError("Wrong number of + and - points");
		}
		if (iValidFlag ==true && iCrossVal==false && iRank==0) {
			l::log("******* STATISTICS FOR Node[{}] *******", iRank);
			l::log("{} == {}  ,  {} == {}", positiveClass, TruePos+FalseNeg, negativeClass , TrueNeg+FalsePos);
			l::log("Positive points : {}, Negative Points : {}", positiveClass, negativeClass);
			l::log("positive instances: %{}, negative instances: %{}", (positiveClass * 100) / TestSize, (negativeClass * 100) / TestSize);
			l::log("TP: {}, TN:{} , FP:{}, FN:{}", TruePos, TrueNeg, FalsePos, FalseNeg);
			l::log("Sensivity / ture positive rate:{} ", Sensivity);
			l::log("Precision / positive predictive rate:{}", Precision);
			l::log("Accuracy: {}", Accuracy);
			l::log("Specificity / True Negative Rate: {}", Specificity);
			l::log("Negative Predictive Rate: {}", NegPerdicRate);
			l::log("***************************");
		}
	return correctHit;
}

Eigen::VectorXi Predicting(const Eigen::MatrixXd& iTestData, const DataMatrix& DataMat, const  Eigen::MatrixXd& iT, const TypeForKernel& iKernelType, HyperPlaneScalars& ia, WeightScalar& ib, HyperPlaneScalars& ic, double iGamma) {
	Eigen::VectorXd first = K(iTestData, DataMat, iKernelType,iGamma)*ia;
	Eigen::VectorXd second = K(iTestData, iT, iKernelType,iGamma)*ic;
	Eigen::VectorXd Ones = Eigen::VectorXd::Ones((int)iTestData.rows());
	Eigen::VectorXd third = ib * Ones;
	Eigen::VectorXd discriminantFunc = first + second + third;
	Eigen::VectorXi PredicLbl((int)discriminantFunc.size());

	for (int j = 0; j < (int)discriminantFunc.size(); j++) {
		if (std::signbit(discriminantFunc(j))) {
			PredicLbl(j) = -1;
		}
		else {
			PredicLbl(j) = 1;
		}
	}

	return PredicLbl;
}