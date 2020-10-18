#include "Types.h"
#include "Partitioning.h"
#include <iostream>
#include "FoldData.h"

struct FoldInfo
{
	DataMatrix TrainMat;
	LabelVector TrainLab;

	DataMatrix TestMat;
	LabelVector TestLab;

};
TestTrainData FoldData(int iNoFold, const DataMatrix& iWholeMat, const LabelVector& iWholeLab, bool iRandomOn) {
	TestTrainData PartData(iNoFold);
	auto DividedMat = Partitioning(iWholeMat, iWholeLab, iNoFold, iRandomOn,"cycle");
	std::vector<int> crossNumbers(iNoFold);
	std::vector<int> crossTmp(iNoFold);
	Column col = (Column)iWholeMat.cols();

	for (int i = 0; i < iNoFold; i++) {
		crossTmp[i] = i;
		crossNumbers[i] = i;
	}
	for (int i = 0; i < iNoFold; i++)
	{
		std::get<0>(PartData[i]) = std::get<0>(DividedMat[i]);//TestingMat
		std::get<1>(PartData[i]) = std::get<1>(DividedMat[i]);//Testing Labels
		crossNumbers.erase(crossNumbers.begin() + i);
		//DataMatrix trnData(0, iDataMat.cols());
		DataMatrix trnData(0, col);
		LabelVector trnLab(0);
		Row block = 0;//row

		for (auto &cross : crossNumbers) {
			Row row = (Row)std::get<0>(DividedMat[cross]).rows();
			trnData.conservativeResize(trnData.rows() + row, Eigen::NoChange);
			trnData.block(block, 0, row, col) = std::get<0>(DividedMat[cross]);
			trnLab.conservativeResize(trnLab.size() + row);
			trnLab.tail(row) = std::get<1>(DividedMat[cross]);
			block = block + row; 
		}
		std::get<2>(PartData[i]) = trnData;
		std::get<3>(PartData[i]) = trnLab;

		//std::cout << "----------Fold " << i << "------------" << "\n";
		//std::cout << "Test Data" << "\n" << std::get<0>(PartData[i]) << "\n";
		//std::cout << "Test LAbel" << "\n" << std::get<1>(PartData[i]) << "\n";

		//std::cout << "trnData[" << trnData.rows() << " , " << trnData.cols() << "]\n" << trnData << "\n";
		//std::cout << "trnLab[" << trnLab.size() << "]\n" << trnLab << "\n";


		crossNumbers = crossTmp;
	}
	return PartData;
}