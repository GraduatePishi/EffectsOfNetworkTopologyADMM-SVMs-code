#pragma once
#ifndef INVERSE_MAT
#define INVERSE_MAT

#include <Eigen/IterativeLinearSolvers>
#include <Eigen/unsupported/IterativeSolvers>

Eigen::MatrixXd InversMat(const Eigen::MatrixXd& a);

#endif