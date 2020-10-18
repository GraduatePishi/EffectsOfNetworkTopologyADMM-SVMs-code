#pragma once
#include <Eigen/Core>
#include <vector>
#include <memory>
#include "Kernel.h"
#include "Types.h"
#include "Partitioning.h"
#include "FileReader.h"
#include "stdio.h"
#include "stdlib.h"
#ifdef _WIN32
#include "C:\Program Files (x86)\Microsoft SDKs\MPI\Include\mpi.h"
#else
#include "mpi.h"
#endif

class nodesolver;
struct NodeData;
//using RowMajorMatirx = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

void Compute_Weigh_Margin(HyperPlaneScalars _w, WeightScalar _b, HyperPlaneScalars _lambda, Eigen::MatrixXd _FixedMatrix);
void Broadcast_wj_bj(std::pair<std::unique_ptr<nodesolver>, NodeData>& iNode, const std::vector<std::pair<std::unique_ptr<nodesolver>, NodeData>>& iNeighbors);
void MPI_Await(std::vector<MPI_Request>& iReq);
std::pair<Eigen::MatrixXd, Eigen::MatrixXd> getTestKernel(const RowMajorMatirx& iTestData, const DataMatrix& PartialDataMat, const  DimReductionMat& iT, KrnlMatrx iNormT_squared, KrnlMatrx iNorm_Xsquared, double iGamma);
IndexOfNeighbors DistrNeighToAllNnodes(int iWorld_rank, Graph igraph, int iWorld_size, std::string iGraphFile, bool iCrossVal, double iJC, double iGamma, double iEta);
struct NodeData
{
	NodeData() = default;
	NodeData(int consensus, size_t DataRow, size_t DataCol) {
		SumOfWi = WeighVector::Zero(consensus);
		Wi_old = SumOfWi;
		SumOfbi = 0.0;
		Vi = WeightsAndBias::Zero(consensus + 1);
		LambdaPossitive = Eigen::VectorXd::Zero(DataRow);
		DataSV = DataMatrix::Zero(DataRow, DataCol);
		LabelSV = LabelVector::Zero(DataRow);
	};
	WeighVector SumOfWi,Wi_old;
	Bias SumOfbi, bi_old;
	WeightsAndBias Vi; //for termination
	HyperPlaneScalars Lambda_i, Lambda_i_old;

	Eigen::VectorXd LambdaPossitive;

	DataMatrix DataSV;
	LabelVector LabelSV;
	IndexOfNeighbors SVIndx;
	int NrSVs;
	int terminationFlag = 0;
	int convergenceFlagNei = 0;
	double w_primal= std::numeric_limits<double>::max();
	double b_primal = std::numeric_limits<double>::max();
	double w_dual = std::numeric_limits<double>::max();
	double b_dual = std::numeric_limits<double>::max();
	double Vall_vs_Vj = std::numeric_limits<double>::max();
	double Vj_Vj= std::numeric_limits<double>::max();
	double esp_prim=1e-4;
	double esp_dual=1e-3;
};
class nonlinear_svm
{
public:
//	std::vector<std::tuple<HyperPlaneScalars, Bias, HyperPlaneScalars, hitRate, TotalSolverTime, InnerSolverTime, CommunTime, Residuals, Residuals, Residuals, Residuals, Residuals>> solve(int &iteration, const ProblemStatment& iData, const ProblemStatment& iEvaluation, const  DimReductionMat& iT, const TypeForKernel iKerneltype, const int NumberofNodes,  const double iJC, const double iEta, const double iGamma, bool RandomOn, Graph igraph, bool iADMM_update, std::string iGraphFile, int iworld_size, int iworld_rank, double iepsilonVal, double iABSTol, double iRELTOL, std::string iTerminationCriteria, bool iCrossVal);
	std::vector<std::tuple<HyperPlaneScalars, Bias, HyperPlaneScalars, hitRate, TotalSolverTime, InnerSolverTime, CommunTime, Residuals, Residuals, Residuals, Residuals, Residuals, KrnlMatrx, KrnlMatrx, ItrTime>> solve(int& iteration, const ProblemStatment& iData, const ProblemStatment& iEvaluation, const  DimReductionMat& iT, const TypeForKernel iKerneltype, const int TotalNumberofNodes, const double iJC, const double iEta, const double iGamma, bool RandomOn, Graph igraph, bool iADMM_update, std::string iGraphFile, int iworld_size, int iworld_rank, double iepsilonVal, double iABSTOL, double iRELTOL, std::string iTerminationCriteria, bool iCrossVal, std::string iName);

private:

};

void CreateMydirectory(std::string directoryName);

