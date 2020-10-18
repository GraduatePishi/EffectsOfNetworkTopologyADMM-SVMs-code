#pragma once
#include "Types.h"
#include "FileReader.h"
#include "Printer.h"
#include "Partitioning.h"
#include "Randomizing.h"
#include "ChoosingT.h"
#include <string>
#ifndef GTEST_SAMPLES_SAMPLE1_H_
#define GTEST_SAMPLES_SAMPLE1_H_

//std::tuple<DataMatrix, LabelVector, Row, Column, MinVec, MaxVec, MinVec, MaxVec, DimReductionMat, ProblemStatment> DistributedData(int iWorld_size, int iWorldRank, ProblemStatment iProb, bool iShuffle, std::string iScaling, fs::path iTrainfilePath , int ReducedL, bool preDefT, TypeForDesign designType, size_t iSeed);
std::tuple<ProblemStatment, DimReductionMat, ProblemStatment, MinVec, MaxVec > DistributedData(int iWorld_size, int iWorldRank, bool iShuffle, std::string iScaling, fs::path iTrainfilePath, int ReducedL, bool preDefT, TypeForDesign designType, size_t iSeed);
std::tuple<RowMajorMatirx, RowMajorMatirx> PCA_Train(RowMajorMatirx traindata, int iReducedCo);
std::tuple<int, int> PosNegRatio(LabelVector iLabl);
#endif  // GTEST_SAMPLES_SAMPLE1_H_