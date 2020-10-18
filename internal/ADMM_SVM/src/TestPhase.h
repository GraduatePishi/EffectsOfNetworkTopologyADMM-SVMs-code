#pragma once
#include "Kernel.h"
#include "Types.h"
#include <iostream>

hitRate Testing(const Eigen::VectorXd& iTestLabel, const Eigen::MatrixXd& iFirst, const Eigen::MatrixXd& iSecond, const  HyperPlaneScalars& ia, WeightScalar& ib, const HyperPlaneScalars& ic, bool iValidFlag, int iRank, bool iCrossVal);
Eigen::VectorXi Predicting(const Eigen::MatrixXd& iTestData, const DataMatrix& DataMat, const  Eigen::MatrixXd& iT, const TypeForKernel& iKernelType, HyperPlaneScalars& ia, WeightScalar& ib, HyperPlaneScalars& ic, double iGamma);
