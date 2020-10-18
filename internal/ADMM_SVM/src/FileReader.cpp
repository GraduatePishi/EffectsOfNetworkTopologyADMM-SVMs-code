#include "FileReader.h"
#include <string>
#include <list>
#include <fmt/format.h>
#include <iostream>
#include <sstream>
#include <iterator>
#include <numeric>
using namespace std;
using ExceptionError = std::runtime_error;
ProblemStatment datReader(fstream& iFileStream)
{
	std::vector<std::tuple<int, int, double>> readData; //<rowCount, column, value>
	std::vector<int> readLable;
	int nrRows = 0, nrCols = 0;
	int rowCount = 1;
	string r;
	getline(iFileStream, r);
	stringstream header(r);
	double headerlable = std::numeric_limits<double>::max();
	double headerval;
	int headercol = -1;
	char headerdumpColon = '.';
	bool withColon = false;
	header >> headerlable >> headercol >> headerdumpColon >> headerval;

	if (headerdumpColon == ':') {
		withColon = true;
	}

	int it = 0;
	while (!iFileStream.eof())
	{
		if(it != 0)
			getline(iFileStream, r);
		stringstream tmp1(r);
		double lable = std::numeric_limits<double>::max();
		double val;
		int col = -1;
		char dumpColon = '.';
		tmp1 >> lable;
		if (withColon) {
			while (tmp1.peek() != EOF) {
				col = -1;
				tmp1 >> col >> dumpColon >> val;
				if (col != -1) { //chacks if the last row exists
					readData.push_back(std::make_tuple(rowCount, col, val));
					nrCols = std::max(col, nrCols);
				}
			}
		}
		else {
			col = 1;
			while (tmp1.peek() != EOF) {
					double value = std::numeric_limits<double>::max();
					tmp1 >> value;
					if (value != std::numeric_limits<double>::max()) {
						readData.push_back(std::make_tuple(rowCount, col++, value));
						nrCols = col - 1;
					}
			}
			
		}

		if (lable != std::numeric_limits<double>::max())
		{
			rowCount++;
			readLable.push_back((int)lable);
		}
		it++;
	}
	nrRows = rowCount - 1;
	ProblemStatment ret(nrRows, nrCols);

	//std::copy(readLable.begin(), readLable.end(), std::ostream_iterator<int>(std::cout, " "));
	//cout << "\n";
	std::list<double> nrClass(readLable.begin(), readLable.end());
	nrClass.sort();
	nrClass.unique();
	
	for (const auto& tmpval : readData)
	{ //because file is of index 1 we use -1
		ret.Data((int)std::get<0>(tmpval) - 1, std::get<1>(tmpval) - 1) = std::get<2>(tmpval);
	}
	int rows = 0;

	for (auto& tmpval : readLable)
	{
		auto list_front = nrClass.begin();

		if (tmpval == *list_front) {
			if (tmpval==1) {
				tmpval = 1;
			}
			else {
				tmpval = -1;
			}
			ret.lables[rows++] = tmpval;
		}
		else {
			std::advance(list_front, 1);
			if (tmpval == 1) {
				tmpval = 1;
			}
			else {
				tmpval = -1;
			}
			ret.lables[rows++] = tmpval;
		}

		//ret.lables[rows++] = val;
	}

	//for (auto val : nrClass) cout << "new list : " << val << endl;
	bool ZeroRowCleaning = false;
	int RowCount;
	std::vector<int> IndexVec;
	if (ZeroRowCleaning == true) {
		Eigen::MatrixXd ZeroVec = Eigen::MatrixXd::Zero(1, ret.Data.cols());
		for (int i = 0; i < ret.Data.rows(); i++) {
			if (!ret.Data.row(i).isApprox(ZeroVec))
				IndexVec.push_back(i);
		}
		RowCount = (int)IndexVec.size();
	}
	else {
		RowCount = (int)ret.Data.rows();
		for (int i = 0; i < ret.Data.rows(); i++) {
			IndexVec.push_back(i);
		}
	}
	//cout << "(int)IndexVec.size()" << (int)IndexVec.size() << "\n" <<"(int)ret.Data.cols()" << (int)ret.Data.cols() <<"\n";

	//std::copy(IndexVec.begin(), IndexVec.end(), std::ostream_iterator<int>(std::cout, " "));
	ProblemStatment returnData(RowCount, (int)ret.Data.cols());
	int rowCounter = 0;

		for (auto& ind : IndexVec) {
			returnData.Data.row(rowCounter) = ret.Data.row(ind);
			//cout << "ret.lables["<<ind<<"] : " <<  ret.lables[ind] << "\n";
			if (ret.lables[ind] == 1) {
				returnData.lables[rowCounter] = 1;
			}
			else {
				returnData.lables[rowCounter] = -1;
			}

			rowCounter++;
		}

	//std::cout << "Data size before cleaning: " << ret.Data.rows() << " x " << ret.Data.cols() << "\n";
	//std::cout << "Data size after cleaning: " << returnData.Data.rows() << " x " << returnData.Data.cols() << "\n";
	//std::cout << "Data matrix" << "\n" << returnData.Data.block(0,0,10, 10 )<< "\n";
	//std::cout << "Data label" << "\n" << returnData.lables.block(0,0,10, 1)<< "\n";

	return returnData;
}

ProblemStatment FileReader(const fs::path& iFilePath) {
	fstream dataFile;
	
	dataFile.open(iFilePath);
	ProblemStatment ret(0,0);
	if (!dataFile) {
		cout << "Path File Not Found!!!" << iFilePath.u8string() << endl;
		throw ExceptionError(fmt::format("Path File Not Found!!! {}", iFilePath.u8string()).c_str());
	}
	else {
		/***counting rows and columns of the input matrix from the file**/
		if (iFilePath.extension() == fs::u8path(".dat") || iFilePath.extension() == fs::u8path(".test"))
		{
			ret = datReader(dataFile);
		}
		if (iFilePath.extension() == fs::u8path(".txt"))
		{
			cout << "not .dat or .validation extension ofr the file!!!!!!!!1" << "\n";
			//problem = txtReader;
		}
		dataFile.close();
	};

	return ret;
}

DimReductionMat ReadT(const fs::path& iFilePath, int iReducedRank, int iCol) {
	fstream dataFile;
	dataFile.open(iFilePath);
	DimReductionMat T(iReducedRank, iCol);
	if (!dataFile) {
		cout << "File T Not Found!!!" << iFilePath << endl;
		throw ExceptionError(fmt::format("File Not Found!!! {}",iFilePath.u8string()).c_str());
	}
	else {
		for (int i = 0; i < iReducedRank; i++) {
			for (int j = 0; j < iCol; j++) {
				dataFile >> T(i, j);
			}
		}
	}
	dataFile.close();
	return T;
}
Eigen::MatrixXd ReadGraph(fstream& iFileStream)
{
	std::vector<std::tuple<int, int, double>> readData; //<rowCount, column, value>
	std::vector<int> readLable;
	int nrRows = 0, nrCols = 0;
	int rowCount = 1;
	string r;


	while (!iFileStream.eof())
	{
		getline(iFileStream, r);
		stringstream tmp1(r);
		int col = 0;
		while (tmp1.peek() != EOF) {
			double value = std::numeric_limits<double>::max();
			tmp1 >> value;
			if (value != std::numeric_limits<double>::max()) {
				readData.push_back(std::make_tuple(rowCount, ++col, value));
				nrCols = col;
			}
		}
		if (col > 0)
		{
			rowCount++;
		}
	}
	nrRows = rowCount - 1;
	Eigen::MatrixXd ret(nrRows, nrCols);

	for (const auto& tmpval : readData)
	{ //because file is of index 1 we use -1
		ret((int)std::get<0>(tmpval) - 1, std::get<1>(tmpval) - 1) = std::get<2>(tmpval);
	}

	return ret;
}
Eigen::MatrixXd GraphReader(const fs::path& iFilePath) {
	fstream dataFile;

	dataFile.open(iFilePath);
	Eigen::MatrixXd ret(0, 0);
	if (!dataFile) {
		cout << "Graph File Not Found!!!" << iFilePath.u8string() << endl;
		throw ExceptionError(fmt::format("Graph File Not Found!!! {}", iFilePath.u8string()).c_str());
	}
	else {
		ret = ReadGraph(dataFile);
		dataFile.close();
	}

	return ret;
}