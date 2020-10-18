#include "ReadDataset.h"
#include "JSon/json.hpp"
#include <fstream>
#include <string>
#include "Types.h"

using nlohmann::json;
using namespace std;

std::string ReadDataSetFile() {
	std::string DataAdrress = "DatasetAddress.dat"; //This make the datasets files to be readable from BOX using different machines
	fs::path DatasetPath = fs::u8path(DataAdrress);
	std::ifstream DataFile(DatasetPath);

	if (!DataFile) {
		cout << "Dataset File Not Found!!!" << endl;
		throw ExceptionError("DataSet File Not Found From BOX!!!");
	}
	json root;
	DataFile >> root;
	std::string Trainset ;
	if (root.count("DatasetAddress") != 0) {
		std::string Tmpset = root["DatasetAddress"];
		Trainset = Tmpset;
	}
	else {
		throw ExceptionError("DataSet Address is incorrect!!!");
	}
	return Trainset;
}