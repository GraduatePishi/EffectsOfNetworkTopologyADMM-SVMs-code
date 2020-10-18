#include <math.h>
#include "Scaling.h"
#include <iostream>
#include "utils.h"
#include <exception>
#include "fmt/format.h"



void ScalingTrain(Eigen::VectorXd& iMean_or_Min, Eigen::VectorXd& iSD_or_Max, RowMajorMatirx& iData, string iType) {
	try {
		if (iType == "standard") {
			for (int i = 0; i < (int)iData.cols(); i++) {
				iMean_or_Min[i] = iData.col(i).sum() / iData.rows();
				Eigen::VectorXd meanVec = iMean_or_Min[i] * Eigen::VectorXd::Ones(iData.rows());
				iSD_or_Max[i] = sqrt((1.0 / (double)iData.rows())*((iData.col(i) - meanVec).cwiseProduct(iData.col(i) - meanVec)).sum()); //based on I tool away -1 from (iData.rows()-1) in denominator
				if (isZero(abs(iSD_or_Max[i]))) { /******* isZero is defined in utils.h based on the numerical limit ********/
					cout << "Training: Zero denominator: " << iSD_or_Max[i] << "!!!!!!!!!!!!!!!!!!!!!" << endl;
					//throw ExceptionError("The SD Value is zero!!!!!!!!!!!!!!");
				}
				iData.col(i) = (iData.col(i) - meanVec) / iSD_or_Max[i];
			}
		}
		else if (iType == "normal")
		{ //rescaling an input variable to the range between 0 and 1.
			for (int i = 0; i < (int)iData.cols(); i++) {
				iMean_or_Min[i] = iData.col(i).minCoeff(); //iMean is same as mins
				iSD_or_Max[i] = iData.col(i).maxCoeff();    //iSD is the same as Maxs
				Eigen::VectorXd meanVec = iMean_or_Min[i] * Eigen::VectorXd::Ones(iData.rows());
				double MaxMin = (iSD_or_Max[i] - iMean_or_Min[i]);
				if (isZero(abs(MaxMin))) {/******* isZero is defined in utils.h based on the numerical limit ********/
					cout << "Training: Minus in the denominator!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
					//throw ExceptionError("The MaxMin Value is zero!!!!!!!!!!!!!!");
				}
				iData.col(i) = (iData.col(i) - meanVec) / MaxMin;
			}
		}
		else if (iType == "normalize")
		{
			for (int i = 0; i < (int)iData.cols(); i++) {
				iMean_or_Min[i] = iData.col(i).transpose().norm();
				if (isZero(abs(iMean_or_Min[i]))) {/******* isZero is defined in utils.h based on the numerical limit ********/
					cout << "Traing: Minus in the denominator!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
					//throw ExceptionError("The iMean_or_Min Value is zero!!!!!!!!!!!!!!");
				}
				//iData.col(i).transpose().normalize();
				iData.col(i) = iData.col(i) / iMean_or_Min[i];
			}
		}
		else {
			cout << " No scaling is requested!!" << "\n";
		}
	}
	catch (std::exception& ie) {
		{
			cout << "Cant find scalingTrain type" << ie.what() << endl;		
			throw ExceptionError(fmt::format("Cant find scalingTrain type: {}", ie.what()));

		}
	}
}

void ScalingTest(Eigen::VectorXd& iMean_or_Min, Eigen::VectorXd& iSD_or_Max, RowMajorMatirx& iData, string iType) {
	try {
		Eigen::VectorXd meanVec;
		if (iType == "standard") {
			for (int i = 0; i < (int)iData.cols(); i++) {
				meanVec = iMean_or_Min[i] * Eigen::VectorXd::Ones(iData.rows());
				if (isZero(iSD_or_Max[i])) {
					cout << "Test: Zero denominator[i: " << i << "]: " << iSD_or_Max[i] << "!!!!!!!!!!!!!!!!!!!!!" << endl;
					//throw ExceptionError("The iMean_or_Min Value is zero!!!!!!!!!!!!!!");
				}
				iData.col(i) = (iData.col(i) - meanVec) / iSD_or_Max[i];
			}
		}
		else if (iType == "normal") {
			for (int i = 0; i < (int)iData.cols(); i++) {
				meanVec = iMean_or_Min[i] * Eigen::VectorXd::Ones(iData.rows());
				double MaxMin = (iSD_or_Max[i] - iMean_or_Min[i]);
				if (isZero(abs(MaxMin))) {
					cout << "Testing: Minus in the denominator!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
					//throw ExceptionError("The iMean_or_Min Value is zero!!!!!!!!!!!!!!");
				}
				iData.col(i) = (iData.col(i) - meanVec) / MaxMin;
			}
		}
		else if (iType == "normalize") {
			for (int i = 0; i < (int)iData.cols(); i++) {
				if (isZero(abs(iMean_or_Min[i]))) {
					cout << "Testing: Minus in the denominator!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
					//throw ExceptionError("The iMean_or_Min Value is zero!!!!!!!!!!!!!!");
				}
				iData.col(i) = iData.col(i) / iMean_or_Min[i];
			}
		}
		else {
			cout << " No scaling is requested!!" << "\n";
		}
	}
	catch (std::exception& ie)
	{
		cout << "Cant find scalingTest type" << ie.what() << endl;
		throw ExceptionError(fmt::format("Cant find scalingTest type {}",ie.what()));
	}
}