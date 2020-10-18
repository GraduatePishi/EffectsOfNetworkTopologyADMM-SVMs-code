#pragma once
#include <Eigen/Core>
#include <Eigen/QR>
#include <vector>
#include "Types.h"
#include "Kernel.h"
#include "InverseOfMat.h"
#pragma warning( push )
#pragma warning( disable : 4267)
#include <nlopt/nlopt.hpp>
#pragma warning( pop )

double objeFunc(HyperPlaneScalars iLambda, void *myData);
double myfunc(const std::vector<double> &x, std::vector<double> &grad, void *myData);
struct FuncVariables
{
	FuncVariables(const Eigen::MatrixXd& iFixedMatrix, const Eigen::MatrixXd& ikernelTXj_kernelTildaTXj, const LabelVector& iYj, const Eigen::VectorXd& ifTilda, double iHj, int iNrNeigh, double iEta)
		: FixedMatrix(iFixedMatrix), fTilda(ifTilda),kernelTXj_kernelTildaTXj(ikernelTXj_kernelTildaTXj), Yj(iYj), hj(iHj),NrNeighbors(iNrNeigh), Eta(iEta)
	{
		Eigen::RowVectorXd FtK = (fTilda.transpose()*kernelTXj_kernelTildaTXj);
		Row rows = (int)kernelTXj_kernelTildaTXj.cols();
		Eigen::RowVectorXd h = (hj / (2.0*Eta*NrNeighbors))* (Eigen::RowVectorXd::Ones(rows));
		FtK_h = (FtK + h).transpose().cwiseProduct(Yj);
		FtK_h_PlusOne = FtK_h + Eigen::VectorXd::Ones(FtK_h.size());
	};
	const Eigen::MatrixXd& FixedMatrix;
	const Eigen::VectorXd& fTilda;
	const Eigen::MatrixXd& kernelTXj_kernelTildaTXj;
	Eigen::VectorXd FtK_h;
	Eigen::VectorXd FtK_h_PlusOne;
	LabelVector Yj;
	double hj;
	NeighCount NrNeighbors;
	double Eta;
	int SolverCounter{0};
};
std::tuple<std::vector<double>, double> SolverFunc(nlopt::vfunc f, FuncVariables iFuncData, double JC, int iRow);
class nodesolver
{
	friend double myfunc(const std::vector<double> &lambda, std::vector<double> &grad, void *myNode);
	friend bool BroadcastingTest();
public:
	double ObjectivFunc{ 0 }, kj{ 0 };
	nodesolver( const DataMatrix& Xj, const LabelVector& iY, const DimReductionMat& iT, const TypeForKernel iKerneltype, const IndexOfNeighbors NrOfNeighbors, const double iJC, const double iEta, const double iGamma, const int iRank);
	const auto& getaj() const { return maj; };
	const auto& getbj() const { return mbj; };
	const auto& getcj() const { return mcj; };
	const auto& getLambdaj() const { return mLambdaj; };
	const auto& getWTilde() const { return mwTildaj; };
	void setWTilde(Eigen::VectorXd iW) { mwTildaj = iW; };

	const IndexOfNeighbors& getNeigbors() const { return mNeighbors; }
	auto getAlphaj() const { return malphaj; };
	double getBeta() const { return mBetaj; };
	int compute_lambda_w_b(NeighCount iNrOfNeigh, Bias iSumOfbi, WeighVector iSumOfWi, int iRank, int iter);
	void compute_alpha_Beta(Bias iSumOfbi, WeighVector iSumOfWi);
	std::tuple<HyperPlaneScalars, Bias, HyperPlaneScalars> compute_a_b_c(Bias iSumOfbi, WeighVector iSumOfWi);
	std::tuple< KrnlMatrx, KrnlMatrx, KrnlMatrx, KrnlMatrx, KrnlMatrx, KrnlMatrx, KrnlMatrx> calcGramMat(const DataMatrix& iX, const LabelVector& iy, const DataMatrix& iT, const TypeForKernel iKerneltype, int iNrOfNeighbors, double iEta, double iGamma);
	HyperPlaneScalars calc_c(const Multipliers& ilambda, const Eigen::VectorXd& ifTilda) const;

	TypeForKernel mKerneltype;
	double conditionNO;
	double mEta{ 0 };
	const double mGamma{ -1 };
	int mSolverCounter{ 0};
	HyperPlaneScalars mfTilde;
	WeightScalar mhj;
	LabelVector mYj;
	DataMatrix mXj;
	//double mElapsed_nlopt_time, mElapsed_NodeSolver_time;
	std::string mSolverResult;
	int mNLOPTiter = -1;
	int mRank=-2;
	Eigen::MatrixXd mFixedMatrix;
	double estCondNum_fixedMat{ 0.0 };
	KrnlMatrx mUinvKrnTXj;
	KrnlMatrx mUinvKrnlTT;
	KrnlMatrx m_normT_squared, m_norm_Xsquared;
	KrnlMatrx mKernel_XjXj;
	DimReductionMat mkernelTXj_mkernelTildaTXj; //temporarly here otherwise should be in private
	Eigen::MatrixXd mkernelTT_mkernelTildaTT;
	nodesolver(const IndexOfNeighbors& iNeig) : mNeighbors(iNeig) {};

protected:


private:
	Bias mbj;
	WeighVector mwTildaj ;

	const double JC{ 0 };
	Bias Calc_b(const Multipliers& lambda, const double hj) const;

	WeighVector Calc_Wtilda(const Multipliers& ilambda, const Eigen::VectorXd& ifTilda) const;
	const Eigen::MatrixXd& mT= Eigen::MatrixXd::Zero((int)mwTildaj.rows(), (int)mwTildaj.rows());
	const IndexOfNeighbors mNeighbors;

	HyperPlaneScalars maj;
	HyperPlaneScalars mcj;
	HyperPlaneScalars malphaj;
	Eigen::VectorXd mLambdaj;
	WeightScalar mBetaj;

};