#include "Kernel.h"
#include <iostream>
#include <numeric>
#include <assert.h> 
#include <Eigen/Dense>
#include "Types.h"
using ExceptionError = std::runtime_error;
template<typename Derived>
double _Kernel(const Derived&& a, const  Derived&& b, TypeForKernel kerneltype, double iGamma) {

	switch (kerneltype) {
	case TypeForKernel::linear:
		if (a.dot(b) != a.cwiseProduct(b).sum()) {
			std::cout << "asize : " << a.size() << "\n" << "bsize " << b.size() << "\n";
			throw ExceptionError("Stop!!");
		}
		return (double)(a.dot(b));
		break;
	case TypeForKernel::rbf:
		return exp((a - b).squaredNorm()*(-iGamma));
		break;
	default:
		std::cout << "Type of the consensus matrix is not chosen!!" << " \n";
		return std::numeric_limits<double>::max();
		break;
	}
}
Eigen::MatrixXd K(const Eigen::MatrixXd& ia, const Eigen::MatrixXd& ib, TypeForKernel iKerneltype, double iGamma) {

	if (ia.cols() != ib.cols()) {
		std::cout << "The size of two vectors are not the same in the Kernel!!"<<"\n";
		throw std::invalid_argument("The size of two vectors are not the same!!");
	}
	Eigen::MatrixXd krnl(ia.rows(), ib.rows());
	for (int i = 0; i < ia.rows(); i++) {
		for (int j = 0; j < ib.rows(); j++) {
			krnl(i, j) = _Kernel(ia.row(i), ib.row(j), iKerneltype, iGamma);
			//krnl(i,j) = exp((ia.row(i) - ib.row(j)).squaredNorm()*(-iGamma));
		}
	}
	return krnl;
}

class nodesolver;
Eigen::MatrixXd K_Tilda(const Eigen::MatrixXd& a, const Eigen::MatrixXd& b, const Eigen::MatrixXd& iT, const QRdecomposition& iUTilda, const int Bj, const double iEtaAverage, TypeForKernel iKerneltype, double iGamma)  {
	KrnlTildeMatrx result((int)iT.rows(),(int)b.rows());
	result= iUTilda.solve(K(iT, b, iKerneltype, iGamma));

	return 2.0*iEtaAverage*(Bj)*K(a, iT, iKerneltype, iGamma)*result;
	//return 2.0*iEta*(Bj)*K(a, iT, iKerneltype)*iUTildainv*K(iT, b, iKerneltype);

}
