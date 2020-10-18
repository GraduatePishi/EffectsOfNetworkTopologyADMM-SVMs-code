#pragma once
#include <Eigen/Core>

double FindCondiNumber(Eigen::MatrixXd const& iMat) {
	Eigen::JacobiSVD<Eigen::MatrixXd> svd_Mat(iMat);
	double cond_Mat = svd_Mat.singularValues()(0) / svd_Mat.singularValues()(svd_Mat.singularValues().size() - 1);
	return cond_Mat;
}
template <typename T>
constexpr void PrintStdVectors(std::string iName, std::vector<T> iVector) {
	cout << iName << " = [";
	std::copy(iVector.begin(), iVector.end(), std::ostream_iterator<int>(std::cout, " , "));
	cout << "]" << endl;
}