#pragma once
#include <fstream>
#include <Eigen/Core>
template<typename MType>

void PrintData(const Eigen::DenseBase<MType>& iMatrix, std::string Filename) {
	std::fstream dataForMatlab;
	dataForMatlab.open(Filename, std::fstream::out);

	dataForMatlab << iMatrix;

	dataForMatlab.close();
}