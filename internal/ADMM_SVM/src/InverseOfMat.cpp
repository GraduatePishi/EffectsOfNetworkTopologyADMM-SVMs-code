#include "InverseOfMat.h"

Eigen::MatrixXd InversMat(const Eigen::MatrixXd& a) {
	Eigen::GMRES<Eigen::MatrixXd> gmres;
	gmres.compute(a);
	Eigen::MatrixXd I(a.rows(), a.cols());
	I.setIdentity();
	return gmres.solve(I);
}

