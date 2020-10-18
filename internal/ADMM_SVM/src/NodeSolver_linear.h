#pragma once
#pragma once
#include <Eigen/Core>
#include <vector>
#include "Types.h"
#include "Kernel.h"
#include "InverseOfMat.h"

struct FuncVariables_lin
{
	FuncVariables_lin(const Eigen::MatrixXd& iFixedMatrix, const DataMatrix iXj, const Eigen::MatrixXd iUinv, const LabelVector& iYj, const Eigen::VectorXd& i_f, int iNrNeigh, double iEta)
		: FixedMatrix_lin(iFixedMatrix), mf(i_f),  Yj_lin(iYj), mXj_lin(iXj), mUinv_lin(iUinv), NrNeighbors(iNrNeigh), mEta(iEta) {};
	const Eigen::MatrixXd& FixedMatrix_lin;
	const Eigen::VectorXd& mf;
	LabelVector Yj_lin;
	DataMatrix mXj_lin;
	Eigen::MatrixXd mUinv_lin;
	NeighCount NrNeighbors;
	double mEta;
};

class nodesolver_linear
{
	friend double myfunc_lin(const std::vector<double> &lambda, std::vector<double> &grad, void *myNode);
	
public:
	nodesolver_linear(const DataMatrix& Xj, const LabelVector& iY,  IndexOfNeighbors NrOfNeighbors,  const double iJC, const double iEta);
	//auto getbj() const { return mbj; };
	auto getLambdaj_lin() const { return lambdaj_lin; };
	auto getV_lin() const { return Vj; };
	const IndexOfNeighbors& getNeigbors() const { return mNeighbors; }
	auto getAlphaj_lin() const { return malphaj_lin; };
	void compute_lambda_v(NeighCount iNrOfNeigh, WeightsAndBias iAverage_vi);
	void compute_alpha_lin(WeightsAndBias iAverage_vi);
	std::tuple<Eigen::MatrixXd, Eigen::MatrixXd> calculateKernal(const DataMatrix & iX, const LabelVector & iy, int iNrOfNeighbors, double iEta);
	std::tuple<double, double, double, double> ObejctiveFunc(DataMatrix iGramMat, Eigen::VectorXd ifj, Multipliers iLambda, LabelVector iY, DataMatrix iX, Eigen::MatrixXd iUinv, Row iDataPoints, double iEta, NeighCount iNrNeigh);
	Eigen::VectorXd f_lin(const Eigen::VectorXd& alpha, const WeightsAndBias& vj, const WeightsAndBias& vi, double iEta, double iNrOfNeighrs);
	const double mEta;
	LabelVector Yj_lin;
	DataMatrix mXj_lin;
	Eigen::MatrixXd mFixedMatrix_lin;

	WeightsAndBias Vj;
	const double JC;
	//Bias Calc_b(const Multipliers& lambda) const;
	WeighVector Calc_V(const Multipliers& lambda, const Eigen::VectorXd& mf) const;
	const IndexOfNeighbors mNeighbors;
	Eigen::MatrixXd mU_lin;
	Eigen::MatrixXd mUinv_lin;
	HyperPlaneScalars malphaj_lin;
	Eigen::VectorXd lambdaj_lin;
	//Eigen::VectorXd mfTilda;
	Eigen::VectorXd mfj;
	IndexOfNeighbors indxSV_lin;
};