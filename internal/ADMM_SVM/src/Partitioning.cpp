#include "Solver.h"
#include<iostream>
#include<list>
#include<string>
#include<algorithm>
#include <vector>
#include <cmath>
#include <cstdlib>
#include "Printer.h"
#include "stdio.h"
#include "stdlib.h"
#ifdef _WIN32
#include "C:\Program Files (x86)\Microsoft SDKs\MPI\Include\mpi.h"
#else
#include "mpi.h"
#endif
using namespace std;

std::vector<int> PartSize(int group, int MatRow) {
	std::vector<int> GroupSize(group);
	int Perrow;
	if (MatRow % group != 0) {
		Perrow = (int)std::floor(MatRow / group);
		for (auto&val : GroupSize) val = Perrow;
		int m = 0;
		for (int i = 0; i < MatRow % group; i++) {
			GroupSize[m++] = GroupSize[m] + 1;
		};
	}
	else {
		Perrow = MatRow / group;
		for (auto&val : GroupSize) val = Perrow;
	}
	return GroupSize;
}
// std::tuple < DataMatrix, LabelVector> Partitioning_MPI(const RowMajorMatirx& _data, const LabelVector& _label, const int group, bool RandomOn) {
//	std::tuple < DataMatrix, LabelVector> Matrix_labl;
//	const int nrows = (int)_data.rows();
//	const int ncols = (int)_data.cols();
//	std::vector<int>  PerCountMat(group), PerCountLabel(group), disp(group), disp2(group);
//	PerCountLabel = PartSize(group, nrows);
//	PerCountMat = PerCountLabel;
//	std::for_each(PerCountMat.begin(), PerCountMat.end(), [ncols](int &el) {el *= ncols; });
//
//	disp[0] = 0;
//	disp2[0] = 0;
//	for (int i = 1; i < group; i++) {
//		disp[i] = disp[i - 1] + PerCountMat[i - 1];
//		disp2[i] = disp2[i - 1] + PerCountLabel[i - 1];
//	}
//	int maxCount = PerCountMat[0];
//	int maxCount2 = PerCountLabel[0];
//	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>  localsMatrix(maxCount2, ncols); //Row Major stored localmatricx for each node
//	Eigen::VectorXd LocalLabels(maxCount2);
//
//	MPI_Scatterv(_data.data(), &PerCountMat[0], &disp[0], MPI_DOUBLE, &localsMatrix(0, 0), maxCount, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//	MPI_Scatterv(&_label[0], &PerCountLabel[0], &disp2[0], MPI_DOUBLE, &LocalLabels[0], maxCount2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//	get<0>(Matrix_labl) = localsMatrix;
//	get<1>(Matrix_labl) = LocalLabels;
//	return Matrix_labl;
//}


std::vector < std::tuple < DataMatrix, LabelVector>> Partitioning(const DataMatrix& _data, const LabelVector& _label, const int group, bool RandomOn,Graph igraph) {
	int nrrows = (int)_data.rows();
	std::vector<int> a(nrrows);
		for (int i = 0; i < nrrows; ++i)
		a[i] = i;
		if (RandomOn == true) {
			random_unique(a.begin(), a.end(), nrrows); //shuffles the dataset
		}
	auto GroupSize = PartSize(group, nrrows);
	
	std::vector<std::tuple<DataMatrix, LabelVector>> ret_Parts(group);

	/*******Duplicating the dataset into two identical sets*******/
	if(DuplicatingData) {
		for (int i = 0; i < group; i++) {
			std::get<0>(ret_Parts[i]) = DataMatrix(nrrows, (int)_data.cols());
			std::get<1>(ret_Parts[i]) = LabelVector(nrrows);
			//std::get<2>(ret_Parts[i]) = neighbors[i];
			for (int j = 0; j < nrrows; j++) {
				std::get<0>(ret_Parts[i]).row(j) = _data.row(j);
				std::get<1>(ret_Parts[i]).row(j) = _label.row(j);
			}
		}
	}
	else{
	/*******For divding the dataset into different parts*******/
		int jj = 0;
		for (int i = 0; i < group; i++) {
			std::get<0>(ret_Parts[i]) = DataMatrix(GroupSize[i], (int)_data.cols());
			std::get<1>(ret_Parts[i]) = LabelVector(GroupSize[i]);
			for (int j = 0; j < GroupSize[i]; j++) {
				std::get<0>(ret_Parts[i]).row(j) = _data.row(a[jj]);
				std::get<1>(ret_Parts[i]).row(j) = _label.row(a[jj]);
				jj++;
			}
		}
	}
	return ret_Parts;
}

void Compute_Weigh_Margin(HyperPlaneScalars _w, WeightScalar /*_b*/, HyperPlaneScalars _lambda, Eigen::MatrixXd _FixedMatrix) {
	double weight{ 0 }, margin{ 0 };
	//weight = _w.transpose()*_w;
	weight = _lambda.transpose()*_FixedMatrix*_lambda;
	margin = 2.0 / sqrt(weight);
	
	if(Printing) {
		double ww = _w.transpose()*_w;
		l::log("-----------------------------------------" );
		l::log("||w||^2 = sum(lambda) =                  {}", _lambda.sum() ) ;
		l::log("||w||^2 = lambda'*kernel*lambda =        {}" , weight) ;
		l::log("||w||^2 = w*w =                          {}" ,  ww);
		l::log("margin = 2/sqrt(lambda'*kernel*lambda) = {}", margin );
		l::log("-----------------------------------------" );
	}
}