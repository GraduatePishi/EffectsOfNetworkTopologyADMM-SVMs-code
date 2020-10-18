#include "NodeSolver_linear.h"
#include<iostream>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/unsupported/IterativeSolvers>
#pragma warning( push )
#pragma warning( disable : 4267)
#include <nlopt/nlopt.hpp>
#pragma warning( pop )
using namespace std;
using ExceptionError = std::runtime_error;
Eigen::VectorXd nodesolver_linear::f_lin(const Eigen::VectorXd& alpha, const WeightsAndBias& vj, const WeightsAndBias& vi, double iEta, double iNrOfNeighrs) {
	return 2.0*alpha - iEta * iNrOfNeighrs*(vj + vi);
}

// formula (35) page 1677
double myfunc_lin(const std::vector<double> &x, std::vector<double> &grad, void *myData) {
	const auto data = static_cast<const FuncVariables_lin*>(myData);
	Eigen::VectorXd _lambda = Eigen::VectorXd::Zero(x.size());
	for (int i = 0; i < x.size(); i++) {
		_lambda[i] = x[i];
	}

	Eigen::VectorXd FixedLambdaMult = data->FixedMatrix_lin*_lambda;
	Eigen::VectorXd FixedThirdPart = ((data->mXj_lin*data->mUinv_lin*data->mf).cwiseProduct(data->Yj_lin)).transpose();
	if (!grad.empty()) {
		Eigen::VectorXd firstPartVec = FixedLambdaMult;
		Eigen::VectorXd SecondpartVec = Eigen::VectorXd::Ones(x.size());
		Eigen::VectorXd thirdPartVec = FixedThirdPart;
		Eigen::VectorXd gradEig = firstPartVec - SecondpartVec - thirdPartVec;
		for (int i = 0; i < x.size(); i++) {
			grad[i] = gradEig(i);
		}
	}
	double first = 0.5*_lambda.dot(FixedLambdaMult);
	double second = _lambda.sum();
	double third = FixedThirdPart.dot(_lambda);
	double total = first - second - third;
	return total;
}

// formula (17) page 1671
WeighVector nodesolver_linear::Calc_V(const Multipliers& lambda, const Eigen::VectorXd& mf) const {
	return (nodesolver_linear::mUinv_lin)*( mXj_lin.transpose()*Yj_lin.cwiseProduct(lambda) - mf);
}

std::tuple<Eigen::MatrixXd, Eigen::MatrixXd>  nodesolver_linear::calculateKernal(const DataMatrix& iX, const LabelVector& iy, int iNrOfNeighbors, double iEta)
{
	const Column Attri = (int)iX.cols(); // dimension of the problem + 1 because of the lable columns
	DataMatrix pi_mat = DataMatrix::Zero(Attri, Attri);
	pi_mat(Attri - 1, Attri - 1) = 1.0;
	//cout << "Pi matrix is " << "\n"<< pi_mat<< "\n" ;
	/*Eigen::MatrixXd*/ mU_lin = (1 + 2.0*iEta*iNrOfNeighbors)*Eigen::MatrixXd::Identity(Attri, Attri) - pi_mat;
	/*Eigen::MatrixXd*/ mUinv_lin = InversMat(mU_lin);
	/*Eigen::MatrixXd*/ mFixedMatrix_lin = (iX*mUinv_lin*iX.transpose());
	for (int i = 0; i < iy.size(); i++)
	{
		mFixedMatrix_lin.row(i) = mFixedMatrix_lin.row(i)*iy(i);
	}
	for (int i = 0; i < iy.size(); i++)
	{
		mFixedMatrix_lin.col(i) = mFixedMatrix_lin.col(i)*iy(i);
	}
	//FixedMatrix_lin = FixedMatrix_lin + (1e-10)*Eigen::MatrixXd::Identity(DataRow, DataRow); //In Steves paper he mentioned "% avoid problems when Hessian is badly conditioned "

	//	/* CHECKING THE CONDITION NUMBER OF THE FIXED BIG ATRIX*/
	Eigen::JacobiSVD<Eigen::MatrixXd> svdU(mU_lin);
#ifdef _MSC_VER
	if constexpr(DATACHECKS)
	{
		Eigen::BDCSVD<Eigen::MatrixXd> svdF(mFixedMatrix_lin);
		double conditionNO_lin = svdF.singularValues()(0) / svdF.singularValues()(svdF.singularValues().size() - 1);

		//should be a check here of condit
	}
#endif
	return { mUinv_lin, mFixedMatrix_lin };
}

std::tuple<double, double, double, double> nodesolver_linear::ObejctiveFunc(DataMatrix iGramMat, Eigen::VectorXd ifj, Multipliers iLambda, LabelVector iY, DataMatrix iX, Eigen::MatrixXd iUinv, Row iDataPoints, double /*iEta*/, NeighCount /*iNrNeigh*/)
{
	double first = 0.5*iLambda.transpose()*iGramMat*iLambda;
	double second = Eigen::VectorXd::Ones(iDataPoints).transpose()*iLambda;
	double third = (iX*iUinv*ifj).transpose()*iY.cwiseProduct(iLambda);
	return { first, second, third, (first - second - third) };
}

nodesolver_linear::nodesolver_linear(const DataMatrix& iX, const LabelVector& iy, IndexOfNeighbors iNeighbors, const double iJC, const double iEta) :  mNeighbors(iNeighbors),  JC(iJC), mEta(iEta)
{
	NeighCount NrOfNeighbors = (int)mNeighbors.size();
	const Row DataRow = (int)iX.rows();
	const Column Attri = (int)iX.cols(); // dimension of the problem + 1 because of the lable columns

	mXj_lin = iX; // should be a part matrix
	Yj_lin = iy; //puting lables in a diagonal matrix Yj
	std::tie(mUinv_lin, mFixedMatrix_lin) = calculateKernal(mXj_lin, iy, (int)mNeighbors.size(),mEta);

//	/* INITIALIZATION */
	malphaj_lin = Multipliers::Zero(Attri);

	/* Arbitrary initialization*/
	lambdaj_lin = Multipliers::Ones(DataRow) * 5.0/*rand()*/;
	Vj = WeighVector::Ones(Attri) * 5.0 /*rand()*/;

	mfj = f_lin(malphaj_lin, Vj, WeightsAndBias::Zero(Vj.size()), mEta, NrOfNeighbors);
}
///******** compute lambdaj via(35) **********/
void nodesolver_linear::compute_lambda_v(NeighCount iNrOfNeigh, WeightsAndBias iAverage_vi)
{

	Row DataRow = (int)mXj_lin.rows();
	//Column DataCol = (int)mXj_lin.cols();

	mfj = f_lin(malphaj_lin, Vj, iAverage_vi, mEta, iNrOfNeigh);

	//if constexpr(Printing) {
	//	cout << "----------------- Linear: Before Lambda Calculation-----------------" << "\n\n";
	//	cout << "Vj : " << "\n" << Vj << "\n";
	//}

	nlopt::opt optSolver_lin;
	//nlopt::opt opt(nlopt::GN_ISRES, DataRow);
	//optSolver_lin = nlopt::opt(nlopt::LN_NELDERMEAD, DataRow);
	optSolver_lin = nlopt::opt(nlopt::LD_TNEWTON_RESTART, DataRow);

	FuncVariables_lin tmpFuncData_lin(mFixedMatrix_lin,mXj_lin,  mUinv_lin, Yj_lin, mfj, iNrOfNeigh, mEta);
	optSolver_lin.set_min_objective(myfunc_lin, &tmpFuncData_lin);
	std::vector<double> lb(DataRow, 0.0); //Lower bound of lambda
	optSolver_lin.set_lower_bounds(lb);
	std::vector<double> ub; //Upper bound of lambda
	ub.assign(DataRow, JC);
	optSolver_lin.set_upper_bounds(ub);
	optSolver_lin.set_xtol_rel(1e-4);
	optSolver_lin.set_xtol_abs(1e-4);
	double minf{  };

	std::vector<double> xSolv_lin(DataRow, 0);

	nlopt::result result = optSolver_lin.optimize(xSolv_lin, minf);

	for (int i = 0; i < DataRow; i++) {
		lambdaj_lin[i] = xSolv_lin[i];
		if (lambdaj_lin[i] >(1e-5)) {
			indxSV_lin.push_back(i);
		}
	}

	//cout << "Number of SVs is : " << indxSV.size() << endl;
	//if constexpr(Printing) {
	//	if (result == -5) std::cout << " The status of the NLOPT solver is : FORCED STOP " << "\n";
	//	else if (result == -4) std::cout << " The status of the NLOPT solver is : ROUNDOFF LIMITED " << "\n";
	//	else if (result == -3) std::cout << " The status of the NLOPT solver is : OUT OF MEMORY " << "\n";
	//	else if (result == -2) std::cout << " The status of the NLOPT solver is : INVALID ARGUMENTS " << "\n";
	//	else if (result == -1) std::cout << " The status of the NLOPT solver is : FAILURE " << "\n";
	//	else if (result == 1) std::cout << " The status of the NLOPT solver is : SUCCESS" << "\n";
	//	else if (result == 2) std::cout << " The status of the NLOPT solver is : STOPVAL REACHED" << "\n";
	//	else if (result == 3) std::cout << " The status of the NLOPT solver is : FTOL REACHED" << "\n";
	//	else if (result == 4) std::cout << " The status of the NLOPT solver is : XTOL REACHED" << "\n";
	//	else if (result == 5) std::cout << " The status of the NLOPT solver is :MAXEVAL REACHED" << "\n";
	//	else if (result == 6) std::cout << " The status of the NLOPT solver is : MAXTIME REACHED" << "\n";
	//}
	if (result != 1)
	{
		std::cout << " The NLOPT solver returns " << result << " which is not finished successfully!!" << "\n\n";
		throw ExceptionError(" The NLOPT solver is not finished successfully!!");
	}
	/********  compute wtildaj via (36) **********/
	Vj = Calc_V(lambdaj_lin, mfj);
	/*******************************************/
}

void nodesolver_linear::compute_alpha_lin(WeightsAndBias iAverage_vi)
{
	auto iNrOfNeigh = (int)mNeighbors.size();

	// compute alphaj via (18), iAverage_vi is the Average of Vi of all neighbours
	malphaj_lin = malphaj_lin + (mEta / 2.0)*iNrOfNeigh*(Vj - iAverage_vi);
}