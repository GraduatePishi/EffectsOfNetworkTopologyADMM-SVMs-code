#include "ChoosingT.h"
#include <iostream>
#include <cmath>


/* Choosing the consensus matrix from unaugmented data matrix, So the consensus matrix does not have any label*/
//Eigen::MatrixXd ConsensusT(int rowT, Eigen::MatrixXd _A, TypeForDesign designType, double iEta) {
//	//std::cout << "rowT : " << rowT << ", _A.cols()" << _A.cols() << "\n";
//	Eigen::MatrixXd Tmat = Eigen::MatrixXd::Zero(rowT, (int)_A.cols());
//	//std::cout << "T" << "\n" << Tmat << "\n";
//	for (int i = 0; i < _A.cols() ; i++) {
//		auto max_val = _A.col(i).maxCoeff();
//		auto min_val = _A.col(i).minCoeff();
//		
//		std::random_device rd;
//		std::mt19937 eng(rd());
//		std::uniform_real_distribution<double> distr(min_val, max_val);
//
//		switch (designType) {
//		case TypeForDesign::random:
//
//			for (int j = 0; j < rowT; j++) {
//				//if ( (std::isnan(max_val) && std::isnan(min_val)) || (max_val==0 && min_val==0) ) {
//				//	Tmat(j, i) = 0;
//				//}
//				//else {
//					Tmat(j, i) = distr(eng);
//					if (std::isnan(Tmat(j, i))) {
//						std::cout << "Incorrect value of NAN in T!" << "\n";
//						throw ExceptionError("Incorrect value of NAN!");
//					}
//				//}
//			}
//			break;
//		case TypeForDesign::grid:
//			//for (int j = 0; j < rowT; j++) {			}
//			std::cout << "Type of the consensus matrix is not chosen!!" << std::endl;
//			break;
//		case TypeForDesign::identity:
//				Tmat = Eigen::MatrixXd::Identity((int)rowT, (int)_A.cols());
//				break;
//		default:
//			std::cout << "Type of the consensus matrix is not chosen!!" << std::endl;
//			break;
//
//		}
//	}
//	//std::cout << "T" << "\n" << Tmat << "\n";
//
//return Tmat;
//}
DimReductionMat ConsensusT(int rowT, RowMajorMatirx _A, TypeForDesign designType, MaxVec iMax, MinVec iMin, size_t iSeed) {
	Eigen::MatrixXd Tmat = Eigen::MatrixXd::Zero(rowT, (int)_A.cols());
	for (int j = 0; j < _A.cols(); j++) {
		std::mt19937 eng((unsigned int)iSeed);//std::mt19937 eng((unsigned int) iSeed);
		std::uniform_real_distribution<double> distr(iMin[j], iMax[j]);

		switch (designType) {
		case TypeForDesign::random:
			for (int i = 0; i < rowT; i++) {
				Tmat(i, j) = distr(eng);
				if (std::isnan(Tmat(i, j))) {
					std::cout << "Incorrect value of NAN in T!" << "\n";
					throw ExceptionError("Incorrect value of NAN!");
				}
			}
			break;
		case TypeForDesign::grid:
			std::cout << "Type of the consensus matrix is not chosen!!" << std::endl;
			break;
		case TypeForDesign::identity:
			Tmat = Eigen::MatrixXd::Identity((int)rowT, (int)_A.cols());
			break;
		default:
			std::cout << "Type of the consensus matrix is not chosen!!" << std::endl;
			break;
		}
	}

	return Tmat;
}