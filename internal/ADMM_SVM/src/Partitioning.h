#pragma once
#include "NetworkCreation.h"
//std::vector<IndexOfNeighbors> getNeighbors(const int NumberofNodes);
std::vector<int> PartSize(int group, int MatRow);
//std::tuple < DataMatrix, LabelVector> Partitioning_MPI(const RowMajorMatirx& _data, const LabelVector& _label, const int group, bool RandomOn);
std::vector<std::tuple<DataMatrix, LabelVector>> Partitioning(const DataMatrix& A, const LabelVector& _label, const int N, bool iRandomOn, Graph igraph);

template<class bidiiter>
bidiiter random_unique(bidiiter begin, bidiiter end, int num_random) {
	int left = (int)std::distance(begin, end);
	while (num_random--) {
		bidiiter r = begin;
		std::advance(r, rand() % left);
		std::swap(*begin, *r);
		++begin;
		--left;
	}
	return begin;
}