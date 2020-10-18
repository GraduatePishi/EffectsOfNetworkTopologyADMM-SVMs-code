#pragma once
#include "Types.h"
using namespace std;
void ScalingTrain(Eigen::VectorXd& iMean, Eigen::VectorXd& iSD, RowMajorMatirx& iData, std::string iType);
void ScalingTest(Eigen::VectorXd& iMean, Eigen::VectorXd& iSD, RowMajorMatirx& iData, string iType);