#pragma once
#include <vector>
#include "Types.h"
#include "Kernel.h"
#include "NetworkCreation.h"
#include "FileReader.h"
std::tuple<double, double, double, double> Cross_validation(int iNoFold, ProblemStatment& iData, ProblemStatment& iTest, DimReductionMat iT, TypeForKernel iKrnlTyp, int iNoNodes, int iNodeID, bool iRandomOn, double iEta, double iGamma, std::string iDataSetName, CrossValIntrval Interval, Graph igraph, bool iADMM_update, std::string iGraphFile, int iworld_size, int iworld_rank, bool iMPIflag, double iEpsilonVal, double iABSTOL, double iRELTOL, std::string iTerminationCriteria, bool iCross_val, double iGammaStep, double iEtaStep, double iJcStep);