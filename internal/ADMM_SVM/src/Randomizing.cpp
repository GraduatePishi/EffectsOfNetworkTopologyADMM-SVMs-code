#include "Randomizing.h"
#include "Partitioning.h"
#include <iostream>
#include <random>
#include<iterator>
using namespace std;

void ShuffleRandomly(RowMajorMatirx& iData, LabelVector& iLabel) {
	int NrRows = (int)iData.rows();
	auto CopyData = iData;
	auto CopyLabel = iLabel;
	std::vector<int> a(NrRows);
	for (int i = 0; i < NrRows; ++i) {
		a[i] = i;
	}
	random_unique(a.begin(), a.end(), NrRows); //shuffles the dataset

	for (int i = 0; i < NrRows; ++i) {
		iData.row(i) = CopyData.row((int) a[i]);
		iLabel.row(i) = CopyLabel.row((int) a[i]);
	}

}