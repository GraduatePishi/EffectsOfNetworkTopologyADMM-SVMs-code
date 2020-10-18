#pragma once
#include <Eigen/Core>
#include <Eigen/QR>
#include <string>
#include <math.h>
#include "Types.h"

enum class TypeForKernel
{
	linear,
	rbf,
	poly
};
Eigen::MatrixXd K(const Eigen::MatrixXd& ia, const Eigen::MatrixXd& ib, TypeForKernel iKerneltype, double iGamma);
Eigen::MatrixXd K_Tilda(const Eigen::MatrixXd& a, const Eigen::MatrixXd& b, const Eigen::MatrixXd& iT, const QRdecomposition& iUTilda, const int Bj, const double iEta, TypeForKernel iKerneltype, double iGamma);