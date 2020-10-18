#pragma once
#include <iostream>    // std::cout, cin
#include <fstream>
#include <sstream> // for stringstream
#include <vector>
#include <Eigen/Dense> 
#include "Types.h"

using RowMajorMatirx = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

struct ProblemStatment
{
	ProblemStatment() = delete;
	explicit ProblemStatment(int rows, int cols) :Data(RowMajorMatirx::Zero(rows, cols)), lables(Eigen::VectorXd::Zero(rows))
	{
	}
	RowMajorMatirx Data;
	Eigen::VectorXd lables;
};

ProblemStatment FileReader(const fs::path&);
DimReductionMat ReadT(const fs::path& iFilePath, int iReducedRank, int iCol);
//Eigen::MatrixXd ReadGraph(fstream& iFileStream);
Eigen::MatrixXd GraphReader(const fs::path& iFilePath);