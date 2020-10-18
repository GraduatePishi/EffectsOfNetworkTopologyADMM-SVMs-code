#pragma once
#include <Eigen/Core>
#include <vector>
#include <chrono>
#include <filesystem>
#define MPIRUN 
constexpr bool DATACHECKS = false;
constexpr bool EVALUATION = true;
constexpr bool Printing = true;
constexpr bool GoodPrinting = true;
constexpr bool PrintLastResults = true;
constexpr bool MatlabCheck = false;
constexpr bool DuplicatingData = false;
using ExceptionError = std::runtime_error;
namespace fs = std::filesystem;

namespace svmTypes
{
	using HyperPlaneScalars = Eigen::VectorXd;
	using Multipliers = Eigen::VectorXd;
	using WeighVector= Eigen::VectorXd;
	using MinVec = Eigen::VectorXd;
	using MaxVec = Eigen::VectorXd;
	using indxSV = std::vector<int>;
	using VectorSVs= Eigen::VectorXd;
	using WeightsAndBias = Eigen::VectorXd;
	using LabelVector = Eigen::VectorXd;
	using WeightScalar = double;
	using Bias = double;

	using Row = int;
	using Column = int;
	using RowMajorMatirx = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
	using DataMatrix = Eigen::MatrixXd;
	using DimReductionMat = Eigen::MatrixXd;
	using krnlParam = double;
	using TrainDataLabel = std::tuple<RowMajorMatirx, LabelVector, Row, Column, MinVec, MaxVec>;
	//using DistributedInfo = std::tuple<DataMatrix, LabelVector, Row, Column, MinVec, MaxVec, MinVec, MaxVec , DimReductionMat , ProblemStatment> ;
	using Graph = std::string;

	using IndexOfNeighbors = std::vector<int>;
	using Neighbors = std::vector<IndexOfNeighbors>;
	using KrnlMatrx = Eigen::MatrixXd;
	using QRdecomposition = Eigen::HouseholderQR<Eigen::MatrixXd>;

	using KrnlTildeMatrx = Eigen::MatrixXd;
	using NeighCount = int;
	using ReducedDim = int;
	using hitRate = double;
	using TestTrainData = std::vector<std::tuple<DataMatrix, LabelVector, DataMatrix,  LabelVector>>;// TestingMat,TestLab, TraingMat, Traininglab 
	using GammMin = double;
	using GammMax = double;
	using EtaMin = double;
	using EtaMax = double;
	using JCmax = double;
	using JCmin = double;
	using CrossValIntrval = std::tuple<EtaMin, EtaMax, GammMin, GammMax, JCmax, JCmin>;
	using DurationTime = std::chrono::duration<double>;
	using TotalSolverTime = double;
	using InnerSolverTime = double;
	using ItrTime = double;
	using CommunTime = double;
	using Residuals = double;
	using TimeType= std::chrono::high_resolution_clock;
	using TimePointType = std::chrono::time_point<TimeType>;
	using TimePrec = std::chrono::milliseconds;
}
using namespace::svmTypes;

