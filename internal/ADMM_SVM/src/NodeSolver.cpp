#include "NodeSolver.h"
#include<iostream>
#include <chrono>
#include "Printer.h"
#include <Eigen/unsupported/IterativeSolvers>
#pragma warning( push )
#pragma warning( disable : 4267)
#include <nlopt/nlopt.hpp>
#ifdef _WIN32
#include "C:\Program Files (x86)\Microsoft SDKs\MPI\Include\mpi.h"
#else
#include "mpi.h"
#endif
#pragma warning( pop )
using namespace std;

namespace {//called a anonomus namespace, it will save us from naming conflict with other files
	auto prettyPrintSolverStatus(int iResult) {
		std::string SolveriResult;
		if (iResult == -5) SolveriResult = " FORCED STOP ";
		else if (iResult == -4) SolveriResult = " ROUNDOFF LIMITED ";
		else if (iResult == -3) SolveriResult = " OUT OF MEMORY ";
		else if (iResult == -2) SolveriResult = " INVALID ARGUMENTS ";
		else if (iResult == -1) SolveriResult = " FAILURE ";
		else if (iResult == 1) SolveriResult = " SUCCESS";
		else if (iResult == 2) SolveriResult = " STOPVAL REACHED";
		else if (iResult == 3) SolveriResult = " FTOL REACHED";
		else if (iResult == 4) SolveriResult = " XTOL REACHED";
		else if (iResult == 5) SolveriResult = " MAXEVAL REACHED";
		else if (iResult == 6) SolveriResult = " MAXTIME REACHED";
		return SolveriResult;
	}
}

//double objeFunc(HyperPlaneScalars iLambda, FuncVariables iFuncData, ) {
//	const auto data = static_cast<const FuncVariables*>(myData);
//	Eigen::VectorXd FixedLambdaMult = data->FixedMatrix*iLambda;
//	double first = 0.5*iLambda.dot(FixedLambdaMult);
//	double second = iLambda.sum();
//	double third = data->FtK_h.dot(iLambda); //I take away the _lambda from here inorder to use in my gradient too without recalculating the matrix
//	double total = first - second - third;
//	return total;
//}

std::tuple<std::vector<double>, double> SolverFunc(nlopt::vfunc f, FuncVariables iFuncData, double JC, int iRow) {
	std::vector<double> ixSolv(iRow, 0);
	double minf{ 0.0 };
	nlopt::opt ioptSolver;
	ioptSolver = nlopt::opt(nlopt::LD_TNEWTON_RESTART, iRow);
	ioptSolver.set_min_objective(f, &iFuncData);
	std::vector<double> lb(iRow, 0.0); //Lower bound of lambda
	ioptSolver.set_lower_bounds(lb);
	std::vector<double> ub; //Upper bound of lambda
	ub.assign(iRow, JC);
	ioptSolver.set_upper_bounds(ub);
	/*nlopt::result iresult =*/ ioptSolver.optimize(ixSolv, minf);
	return { ixSolv , minf };
}


//double Calc_hj(const WeightScalar Betaj, const Bias bj, const Bias Avergbi, const double iEta, int NrOfNeigh) {
//	return 2.0 * Betaj - iEta*NrOfNeigh*(bj + Avergbi);
//}
double Calc_hj(const WeightScalar iBetaj, const Bias ibj, const Bias iSumOfbi, const double iEta, double iNrOfNeigh) {
	return 2.0 * iBetaj - iEta * (iNrOfNeigh*ibj + iSumOfbi);
}
//Eigen::VectorXd f_Tilda(const Eigen::VectorXd& alpha, const Eigen::VectorXd& wTildaj, const Eigen::VectorXd& wTildai, double iEta, double iNrOfNeighrs) {
//	return 2.0*alpha - iEta*iNrOfNeighrs*(wTildaj+ wTildai);
//}
Eigen::VectorXd f_Tilda(const Eigen::VectorXd& ialpha, const Eigen::VectorXd& iwTildaj, const Eigen::VectorXd& iSumOfWi, const double iEta, double iNrOfNeighrs) {
	return 2.0*ialpha - iEta * (iNrOfNeighrs*iwTildaj + iSumOfWi);
}
double objeFunc(HyperPlaneScalars iLambda, void *myData) {
	const auto data = static_cast<const FuncVariables*>(myData);
	double func{ 0 };
	double first = 0.5*iLambda.dot(data->FixedMatrix*iLambda);
	double second = iLambda.sum();
	double third = data->FtK_h.dot(iLambda); //I take away the _lambda from here inorder to use in my gradient too without recalculating the matrix
	func = first - second - third;

	return func;
}

// formula (35) page 1677
double myfunc(const std::vector<double> &x, std::vector<double> &grad, void *myData) {
	auto data = static_cast< FuncVariables*>(myData); //<const FuncVariables*>(myData);
	const Eigen::Map<const Eigen::VectorXd> _lambda(x.data(), x.size());//only masks x as a eigen vector
	Eigen::Map<Eigen::VectorXd> grad_eig(grad.data(), grad.size());//only masks x as a eigen vector
	//int *ptr; //This is to go arounf the fact that data is const
	//ptr = const_cast<int*> (&data->SolverCounter);
	//*ptr = data->SolverCounter + 1;
	grad_eig.noalias() = -1 * data->FixedMatrix*_lambda;
	auto total = (0.5*grad_eig + data->FtK_h_PlusOne).dot(_lambda);

	data->SolverCounter +=1;

	if (!grad.empty()) {
		grad_eig.noalias() += data->FtK_h_PlusOne;
	}
	/********** WARNING: THERE IS A BUG IN THE FORMULATION, IN THE MAX OPTIMIZATION
	THE THIRD PARANTESIS THE SIGN IS + INSTEAD OF - *************************/
	return total;
}

// formula (36) page 1677
WeighVector nodesolver::Calc_Wtilda(const Multipliers& ilambda, const Eigen::VectorXd& ifTilda) const {
	return (mkernelTXj_mkernelTildaTXj)*mYj.cwiseProduct(ilambda) - (mkernelTT_mkernelTildaTT)*ifTilda;
}

std::tuple< KrnlMatrx, KrnlMatrx, KrnlMatrx, KrnlMatrx, KrnlMatrx, KrnlMatrx, KrnlMatrx> nodesolver::calcGramMat(const DataMatrix& iX, const LabelVector& iy, const DataMatrix& iT, const TypeForKernel iKerneltype, int iNrOfNeighbors, double iEta, double iGamma)
{
	int ConsensusRow = (int)iT.rows();
	int nrDataPoints = (int)iX.rows();
	KrnlMatrx Kernl_TT(ConsensusRow, ConsensusRow), Kernl_TXj(ConsensusRow, nrDataPoints), Kernl_XjXj(nrDataPoints, nrDataPoints);
	KrnlMatrx  norm_Xsquared(nrDataPoints, nrDataPoints), norm_TXj(ConsensusRow, nrDataPoints), normT_squared(ConsensusRow, ConsensusRow);

	if (iKerneltype == TypeForKernel::rbf) {
		/********** K(T,T)**********/
		for (int i = 0; i < ConsensusRow; i++) {
			for (int j = i; j < ConsensusRow; j++) {
				normT_squared(i, j) = (double)iT.row(i).dot(iT.row(j));
				normT_squared(j, i) = normT_squared(i, j);
			}
		}
		for (int i = 0; i < ConsensusRow; i++) {
			for (int j = i; j < ConsensusRow; j++) {
				Kernl_TT(i, j) = exp((normT_squared(i, i) - 2 * normT_squared(i, j) + normT_squared(j, j))*(-iGamma));
				Kernl_TT(j, i) = Kernl_TT(i, j);
			}
		}
		/********** K(Xj,Xj)**********/
		for (int i = 0; i < nrDataPoints; i++) {
			for (int j = i; j < nrDataPoints; j++) {
				norm_Xsquared(i, j) = (double)iX.row(i).dot(iX.row(j));
				norm_Xsquared(j, i) = norm_Xsquared(i, j);
			}
		}
		for (int i = 0; i < nrDataPoints; i++) {
			for (int j = i; j < nrDataPoints; j++) {
				Kernl_XjXj(i, j) = exp((norm_Xsquared(i, i) - 2 * norm_Xsquared(i, j) + norm_Xsquared(j, j))*(-iGamma));
				Kernl_XjXj(j, i) = Kernl_XjXj(i, j);
			}
		}
		/********** K(T,Xj)**********/
		for (int i = 0; i < ConsensusRow; i++) {
			for (int j = 0; j < nrDataPoints; j++) {
				norm_TXj(i, j) = (double)iT.row(i).dot(iX.row(j));
				Kernl_TXj(i, j) = exp((normT_squared(i, i) - 2 * norm_TXj(i, j) + norm_Xsquared(j, j))*(-iGamma));
			}
		}
	}
	else {
		Kernl_TT = K(iT, iT, iKerneltype, iGamma);
		Kernl_XjXj = K(iX, iX, iKerneltype, iGamma);
		Kernl_TXj = K(iT, iX, iKerneltype, iGamma);
	}

	KrnlMatrx UTilda = Eigen::MatrixXd::Identity(ConsensusRow, ConsensusRow) + 2.0*iEta*iNrOfNeighbors*Kernl_TT;
	auto& Uinv = UTilda.householderQr();
	auto& UinvKrnlTT = Uinv.solve(Kernl_TT);
	KrnlMatrx kernelTT_kernelTildaTT = Kernl_TT - 2.0*iEta*(iNrOfNeighbors)*Kernl_TT*UinvKrnlTT;


	auto& UinvKrnTXj = Uinv.solve(Kernl_TXj);
	KrnlMatrx KernlTilda_TXj = 2.0*iEta*(iNrOfNeighbors)*Kernl_TT*UinvKrnTXj;
	KrnlMatrx kernelTXj_kernelTildaTXj = Kernl_TXj - KernlTilda_TXj;

	KrnlMatrx KernlTile_XjXj = 2.0*iEta*(iNrOfNeighbors)*(Kernl_TXj.transpose())*UinvKrnTXj;
	KrnlMatrx FixedMatrix = Kernl_XjXj - KernlTile_XjXj + Eigen::MatrixXd::Constant(nrDataPoints, nrDataPoints, 1.0 / (2.0 * iEta*iNrOfNeighbors));


	for (int i = 0; i < iy.size(); i++)
	{
		FixedMatrix.row(i) = FixedMatrix.row(i)*iy(i);
	}
	for (int i = 0; i < iy.size(); i++)
	{
		FixedMatrix.col(i) = FixedMatrix.col(i)*iy(i);
	}

	/*************** Condition Number of Matrices *********/

	//Eigen::VectorXd D = FixedMatrix.ldlt().vectorD().cwiseAbs();
	//bool is_B_singular = D.maxCoeff() >= (D.minCoeff() * 1e8);
	Eigen::VectorXd D = FixedMatrix.ldlt().vectorD();
	bool is_B_singular = isZero(D.minCoeff()) || isGreater(D.minCoeff()) ;
	if (!is_B_singular ) {
		cout << "negative eigenvalue of fixed matrix in node " << mRank << " is " << D.minCoeff() << endl;
	}
	estCondNum_fixedMat = D.maxCoeff() / D.minCoeff();
	//Eigen::EigenSolver<DataMatrix> es(FixedMatrix);
	//cout << "Eigenvalues of Fixed matrix: " << es.eigenvalues().transpose() << endl;


//double HighCondK_T_Xj{ 0.0 }, LowCondK_T_Xj{ 0.0 }, HighCondFixed{ 0.0 }, LowCondFixed{ 0.0 };
//MPI_Allreduce(&cond_TXj, &HighCondK_T_Xj, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
//MPI_Allreduce(&cond_TXj, &LowCondK_T_Xj, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

//MPI_Allreduce(&cond_Fixed, &HighCondFixed, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
//MPI_Allreduce(&cond_Fixed, &LowCondFixed, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
//cout << "I am node " << mRank << " in Gram 186 solver!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11" << endl;
//if (mRank == 0) {
//	//cout << "iT.norm: " << endl << iT.norm() << endl;
//	Eigen::JacobiSVD<DataMatrix> svd_iT(iT);
//	double cond_iT = svd_iT.singularValues()(0) / svd_iT.singularValues()(svd_iT.singularValues().size() - 1);
//	cout << "cond_T in NodeSolver: " << cond_iT << endl;


//	Eigen::JacobiSVD<DataMatrix> svd_Kernl_TT(Kernl_TT);
//	double cond_Kernl_TT = svd_Kernl_TT.singularValues()(0) / svd_Kernl_TT.singularValues()(svd_Kernl_TT.singularValues().size() - 1);
//	cout << "cond_Kernl_TT: " << cond_Kernl_TT << endl;

//	Eigen::JacobiSVD<DataMatrix> svd_UTilda(UTilda);
//	double cond_UTilda = svd_UTilda.singularValues()(0) / svd_UTilda.singularValues()(svd_UTilda.singularValues().size() - 1);
//	cout << "condNum_UTilda: " << cond_UTilda << endl;

//	//cout << "condNum_Kernl_TXj: " << cond_TXj << endl;
//	cout << "Cond_Kernl_TXj[max: " << HighCondK_T_Xj << ", min: " << LowCondK_T_Xj <<"]"<< endl;
//	cout << "Cond_FixedMatrix[max: " << HighCondFixed << ", min: " << LowCondFixed << "]" << endl;
//}



//Eigen::JacobiSVD<DataMatrix> svd_UinvKrnTXj(UinvKrnTXj);
//double cond_UinvKrnTXj = svd_UinvKrnTXj.singularValues()(0) / svd_UinvKrnTXj.singularValues()(svd_UinvKrnTXj.singularValues().size() - 1);
//cout << "condNum_UinvKrnTXj: " << cond_UinvKrnTXj << endl;

//Eigen::JacobiSVD<DataMatrix> svd_KernlTile_XjXj(KernlTile_XjXj);
//double cond_KernlTile_XjXj = svd_KernlTile_XjXj.singularValues()(0) / svd_KernlTile_XjXj.singularValues()(svd_KernlTile_XjXj.singularValues().size() - 1);
//cout << "condNum_KernlTile_XjXj: " << cond_KernlTile_XjXj << endl;

//Eigen::JacobiSVD<DataMatrix> svd_Kernl_XjXj(Kernl_XjXj);
//double cond_Kernl_XjXj = svd_Kernl_XjXj.singularValues()(0) / svd_Kernl_XjXj.singularValues()(svd_Kernl_XjXj.singularValues().size() - 1);
//cout << "condNum_Kernl_XjXj: " << cond_Kernl_XjXj << endl;


//Eigen::JacobiSVD<DataMatrix> svd_fixedMat(FixedMatrix);
//double cond_fixedMat = svd_fixedMat.singularValues()(0) / svd_fixedMat.singularValues()(svd_fixedMat.singularValues().size() - 1);
//cout << "condNum_fixedMat: " << cond_fixedMat << endl;

/*****************************/
	return{ FixedMatrix,UinvKrnTXj, kernelTXj_kernelTildaTXj, normT_squared, norm_Xsquared, UinvKrnlTT, kernelTT_kernelTildaTT };
}

nodesolver::nodesolver(const DataMatrix& iX, const LabelVector& iy, const DimReductionMat& iT, const TypeForKernel iKerneltype, IndexOfNeighbors iNeighbors, const double iJC, const double iEta, const double iGamma, const int iRank) : mT(iT), mKerneltype(iKerneltype), mNeighbors(iNeighbors), JC(iJC), mEta(iEta), mGamma(iGamma), mRank(iRank)
{
	auto start = TimeType::now();
	const NeighCount NrOfNeighbors = (int)mNeighbors.size();
	const Row DataRow = (int)iX.rows();
	ReducedDim ConsensusRow = (int)iT.rows();

	mXj = iX; // should be a part matrix
	mYj = iy;//performance: do not create matrix, can easaly be calculated by element multiplication //puting lables in a diagonal matrix Yj

	std::tie(mFixedMatrix, mUinvKrnTXj, mkernelTXj_mkernelTildaTXj, m_normT_squared, m_norm_Xsquared, mUinvKrnlTT, mkernelTT_mkernelTildaTT) = calcGramMat(mXj, iy, mT, mKerneltype, (int)mNeighbors.size(), mEta, mGamma);

	/* INITIALIZATION */
	malphaj = Multipliers::Zero(ConsensusRow);
	mBetaj = 0.0;
	maj = HyperPlaneScalars::Zero(DataRow);
	mcj = HyperPlaneScalars::Zero(ConsensusRow);

	/* Arbitrary initialization*/
	mLambdaj = Multipliers::Ones(DataRow);// *5.0;
	mwTildaj = WeighVector::Ones(ConsensusRow)*mRank;// *5.0;
	mbj = mRank;// 5.0;// avoid randomizing to control the code

	//mfTilde = f_Tilda(malphaj, wTildaj, HyperPlaneScalars::Zero(ConsensusRow), mEta, NrOfNeighbors);
	//mhj = Calc_hj(mBetaj, mbj, 0, mEta, NrOfNeighbors);
	auto end = TimeType::now();
	DurationTime time = (end - start);
	//std::cout << "Elapsed Time For Kernel Calculations in nodesolver : " <<  time.count() << "\n";
}


/******** compute lambdaj via(35) **********/
//mfTilde, mhj
int nodesolver::compute_lambda_w_b(NeighCount iNrOfNeigh, Bias iSumOfbi, WeighVector iSumOfWi, int iRank, int iter)
{
	Row DataRow = (int)mXj.rows();
	mfTilde = f_Tilda(malphaj, mwTildaj, iSumOfWi, mEta, (double)iNrOfNeigh);
	mhj = Calc_hj(mBetaj, mbj, iSumOfbi, mEta, (double)iNrOfNeigh);
	nlopt::opt optSolver;
	/* G: Global optimization, L: Local optimization, N: NonDerivative-based algorithm, D: Derivative-based algorithm*/
	if (mKerneltype == TypeForKernel::linear) {
		optSolver = nlopt::opt(nlopt::LD_TNEWTON_RESTART, DataRow);
		//optSolver= nlopt::opt(nlopt::GN_ISRES, DataRow);
		//optSolver =nlopt::opt(nlopt::LN_NELDERMEAD, DataRow);
		//optSolver = nlopt::opt(nlopt::GN_DIRECT, DataRow);
		//optSolver = nlopt::opt(nlopt::GD_STOGO_RAND, DataRow); 
	}
	else {
		optSolver = nlopt::opt(nlopt::LD_TNEWTON_RESTART, DataRow);
		/******** the last one used optSolver = nlopt::opt(nlopt::LD_TNEWTON_RESTART, DataRow);*/
		//optSolver = nlopt::opt(nlopt::GD_MLSL, DataRow);
		//optSolver = nlopt::opt(nlopt::GN_ISRES, DataRow);NLOPT_AUGLAG
	}
	//optSolver.set_stopval(optSolver.get_stopval()*0.1);

	FuncVariables tmpFuncData(mFixedMatrix, mkernelTXj_mkernelTildaTXj, mYj, mfTilde, mhj, iNrOfNeigh, mEta);
	optSolver.set_max_objective(myfunc, &tmpFuncData);

	std::vector<double> lb; /****** Lower bound of lambda ****/
	lb.assign(DataRow, 0.0);
	optSolver.set_lower_bounds(lb);
	std::vector<double> ub; /*****Upper bound of lambda ****/
	ub.assign(DataRow, JC);
	optSolver.set_upper_bounds(ub);
	//optSolver.set_xtol_rel(1e-14);
	//optSolver.set_ftol_rel(1e-7);
	//optSolver.set_xtol_abs(1e-4);
	//optSolver.set_maxtime(500);
	//optSolver.set_maxeval(5);	
	/**** Time calculation for nlopt solver ****/
	//auto start = TimeType::now();
	double maxf{ 0.0 };
	std::vector<double> xSolv(DataRow, 0);
	nlopt::result result;
	result = nlopt::result::FAILURE;
	mSolverCounter = 0;
	try {
		result = optSolver.optimize(xSolv, maxf); //shows the status of the optimization whether it was successfull or not: this is not the solution of the optimization
		mSolverResult = prettyPrintSolverStatus(result);
		mSolverCounter = tmpFuncData.SolverCounter;
		//l::log("NLOPT iter in node[{}] is {}", iRank, mSolverCounter);
	}
	catch (std::exception& ie) {
		l::log( "NLOPT solver in node {}  not finished successfully, result is {}.", iRank , prettyPrintSolverStatus(result) );
		std::copy(xSolv.begin(), xSolv.end(), begin(mLambdaj));
		l::log( (ie.what()) );
	}

	int localSuccess=0 ;
	if (result == 1) {
		localSuccess = 1;
	}

	//mNLOPTiter = optSolver.get_numevals;
	//auto end = TimeType::now();
	//mElapsed_nlopt_time = CalcTime(start, end);

	std::copy(xSolv.begin(), xSolv.end(), begin(mLambdaj));
	if (mLambdaj.maxCoeff() > JC || mLambdaj.minCoeff() < 0.0)
		throw ExceptionError("Lambdas are out of range!!");
	/********  compute wtildaj via (36) **********/
	mwTildaj = Calc_Wtilda(mLambdaj, mfTilde);
	/********  compute bj via (37) **********/
	mbj = Calc_b(mLambdaj, mhj);

	return localSuccess; //returns whether solver has been successfull, 1 mean Success
}

/****** formula (34,37) page 1677 ******/
Bias nodesolver::Calc_b(const Multipliers& ilambda, const double ihj) const {
	auto iNrOfNeigh = mNeighbors.size();
	return (1.0 / (2.0*mEta*iNrOfNeigh))*((mYj.cwiseProduct(ilambda)).sum() - ihj);
}
/****** formula (33) page 1677 *****/
HyperPlaneScalars nodesolver::calc_c(const Multipliers& ilambda, const Eigen::VectorXd& ifTilda)  const {
	int iNrOfNeigh = (int)mNeighbors.size();
	Eigen::VectorXd result(ifTilda.size());
	result = mUinvKrnlTT * ifTilda - mUinvKrnTXj * mYj.cwiseProduct(ilambda);
	return 2.0*mEta*iNrOfNeigh*result - ifTilda;
}

void nodesolver::compute_alpha_Beta(Bias iSumOfbi, WeighVector iSumOfWi)
{
	int NrOfNeigh = (int)mNeighbors.size();

	/************ compute alphaj via (38), iAverage_wi is the Average of Wi of all neighbours***/
	malphaj = malphaj + 0.5*mEta*(NrOfNeigh*mwTildaj - iSumOfWi);
	/************ compute alphaj via (39), iAvarage_bi is the Average of bi of all neighbours***/
	mBetaj = mBetaj + 0.5*mEta*(NrOfNeigh*mbj - iSumOfbi);
}
std::tuple<HyperPlaneScalars, Bias, HyperPlaneScalars> nodesolver::compute_a_b_c(Bias iSumOfbi, WeighVector iSumOfWi)
{
	HyperPlaneScalars aj = HyperPlaneScalars::Zero((int)mXj.rows());
	HyperPlaneScalars cj = HyperPlaneScalars::Zero((int)mT.rows());
	Bias bj = 0.0;
	int NrOfNeigh = (int)mNeighbors.size();

	auto fTilde = f_Tilda(malphaj, mwTildaj, iSumOfWi, mEta, NrOfNeigh);
	auto hj = Calc_hj(mBetaj, mbj, iSumOfbi, mEta, NrOfNeigh);

	aj = mYj.cwiseProduct(mLambdaj);	/******** compute aj via (32) not used and overwriten ****/
	cj = calc_c(mLambdaj, fTilde);		/******** compute cj via (33) not used and overwriten ****/
	bj = Calc_b(mLambdaj, hj);			/******** compute bj via (34) ****/
	return { aj,bj,cj };
}
