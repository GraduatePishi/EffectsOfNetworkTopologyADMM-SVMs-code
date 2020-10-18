#pragma once
#include <Eigen/Core>
#include <random>
#include <string>
#include "Types.h"

enum class TypeForDesign
{
	random,
	grid,
	identity
};
//Eigen::MatrixXd ConsensusT(int rowT, Eigen::MatrixXd _A, TypeForDesign designType, double iEta);
DimReductionMat ConsensusT(int rowT, RowMajorMatirx _A, TypeForDesign designType, MaxVec iMax, MinVec iMin, size_t iSeed);