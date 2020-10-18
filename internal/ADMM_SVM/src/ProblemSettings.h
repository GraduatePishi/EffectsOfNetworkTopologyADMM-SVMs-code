#pragma once
#include "Types.h"
#include "Kernel.h"
#include "ChoosingT.h"
#include "NetworkCreation.h"

using namespace std;
Graph returnGraph(std::string igrType);

struct ProblemSetting
{
	ProblemSetting() =default;
	ProblemSetting(const fs::path& iPath, bool& success);
	//void save2File(const fs::path& iPath);

	Graph graph;
	string GraphType;
	bool ADMMParaUpdate;
	bool SpilitData{ false };
	TypeForKernel kernelTyp;		//Type of the Kenel
	double JC;					//C is the penality parameter and J is the number of nodes
	double Eta;			
	string graphFile;//
	TypeForDesign ReducedMatrixType;			// Common matrix for reducing the feature space rank
	ReducedDim ReducedRankL;				// The reduced dimension value
	double Gamma;			//rbf kernel parameter
	double GammaMin, GammaMax;
	double EtaMin, EtaMax;
	double JCMin, JCMax;
	double GammaStep, EtaStep, JcStep;
	bool linear;				//If we choose linear or nonlinear kernels
	std::string TstSetAddrs;			//Testing dataset
	int NumberofNodes{10 };
	bool linearOn;
	int NoFold;
	bool RandomOn;
	bool preDefinedT;
	bool CrossValOn;
	bool Shuffle;
	double EpsilonVal, ABSTOL, RELTOL;
	std::string NameOfOutput;
	std::string Scaling;
	std::string TerminationCriteria = "allthree"; //3 termination criteria, 1. "AllThree", 2."Pri-Dual", and 3. "onlyVj_VjAll";
	CrossValIntrval Interval;
	//Interval= std::make_tuple(EtaMin, EtaMax, GammaMin, GammaMax, JCMin, JCMax);
	double elapsed_nlopt_time;
	int iter{ 200 };
};

