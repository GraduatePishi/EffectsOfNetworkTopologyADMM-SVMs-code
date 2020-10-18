#include "ProblemSettings.h"
#include "Kernel.h"
#include <fstream>
#include <iostream>
#include <JSon/json.hpp>

using nlohmann::json;
Graph returnGraph(std::string igrType) {
	Graph graph;
	std::string key = igrType;
	//std::transform(key.begin(), key.end(), key.begin(), ::tolower);
	graph = key;
	//if (key == "hardcoded")
	//	graph = Graph::hardCoded;
	//else if (key == "alone")
	//	graph = Graph::Alone;
	//else if (key == "line")
	//	graph = Graph::line;
	//else if (key == "star")
	//	graph = Graph::Star;
	//else if (key == "cycle")
	//	graph = Graph::cycle;
	//else if (key == "self")
	//	graph = Graph::self;
	//else if (key == "shirin_threereg")
	//	graph = Graph::Shirin_ThreeReg;
	//else if (key == "shirin_fivereg")
	//	graph = Graph::Shirin_FiveReg;
	//else if (key == "comb")
	//	graph = Graph::Comb;
	//else if (key == "island")
	//	graph = Graph::Island;
	//else if (key == "unconnected")
	//	graph = Graph::unconnected;
	//else if (key == "unconnislands")
	//	graph = Graph::unConnIslands;
	//else if (key == "tworeg")
	//	graph = Graph::TwoReg;
	//else if (key == "threereg")
	//	graph = Graph::ThreeReg;
	//else if (key == "fivereg")
	//	graph = Graph::FiveReg;
	//else if (key == "sevenreg")
	//	graph = Graph::SevenReg;
	//else if (key == "ninereg")
	//	graph = Graph::NineReg;
	//else if (key == "elevenreg")
	//	graph = Graph::ElevenReg;
	//else if (key == "thirteenreg")
	//	graph = Graph::ThirteenReg;
	//else if (key == "fifteenreg")
	//	graph = Graph::FifteenReg;
	//else if (key == "nineteenreg")
	//	graph = Graph::NineteenReg;
	//else if (key == "twentythreereg")
	//	graph = Graph::TwentyThreeReg;
	//else if (key == "twentysevenreg")
	//	graph = Graph::TwentySevenReg;
	//else if (key == "thirtyonereg")
	//	graph = Graph::ThirtyOneReg;
	//else if (key == "thirtythreereg")
	//	graph = Graph::ThirtyThreeReg;
	//else if (key == "thirtyfivereg")
	//	graph = Graph::ThirtyFiveReg;
	//else if (key == "thirtyninereg")
	//	graph = Graph::ThirtyNineReg;
	///*****Network for big data******/
	//else if (key == "tenreg")
	//	graph = Graph::TenReg;
	//else if (key == "twentyreg")
	//	graph = Graph::TwentyReg;
	//else if (key == "twentyfivereg")
	//	graph = Graph::TwentyFiveReg;
	//else if (key == "thirtyreg")
	//	graph = Graph::ThirtyReg;
	//else if (key == "fourtyreg")
	//	graph = Graph::FourtyReg;
	//else if (key == "fourtyfivereg")
	//	graph = Graph::FourtyFiveReg;
	//else if (key == "fiftyreg")
	//	graph = Graph::FiftyReg;
	//else if (key == "sixtyreg")
	//	graph = Graph::SixtyReg;
	//else if (key == "seventyreg")
	//	graph = Graph::SeventyReg;
	//else if (key == "seventyninereg")
	//	graph = Graph::SeventyNineReg;
	//else if (key == "eightyreg")
	//	graph = Graph::EightyReg;
	//else if (key == "nintyreg")
	//	graph = Graph::NintyReg;
	//else if (key == "hundredreg")
	//	graph = Graph::HundredReg;
	//else if (key == "hundredonereg")
	//	graph = Graph::HundredOneReg;
	//else if (key == "hundredtenreg")
	//	graph = Graph::HundredTenReg;
	//else if (key == "twohundredreg")
	//	graph = Graph::TwoHundredReg;
	//else if (key == "hundrednineteenreg")
	//	graph = Graph::HundredNineteenReg;
	//else if (key == "random")
	//	graph = Graph::Random;
	//else if (key == "thirtyfiverandom")
	//	graph = Graph::ThirtyFiveRandom;
	//else if (key == "thirtyonerandom")
	//	graph = Graph::ThirtyOneRandom;
	//else if (key == "twentysevenrandom")
	//	graph = Graph::TwentySevenRandom;
	//else if (key == "twentythreerandom")
	//	graph = Graph::TwentyThreeRandom;
	//else if (key == "nineteenrandom")
	//	graph = Graph::NineteenRandom;
	//else if (key == "fifteenrandom")
	//	graph = Graph::FifteenRandom;
	//else if (key == "elevenrandom")
	//	graph = Graph::ElevenRandom;
	//else if (key == "sevenrandom")
	//	graph = Graph::SevenRandom;
	//else if (key == "fourtysevenreg")
	//	graph = Graph::FourtySevenReg;
	//else
	//	throw ExceptionError("failed to read node structure " + key);

	return graph;
}

ProblemSetting::ProblemSetting(const fs::path & iPath, bool& success)
{
	std::ifstream fileStream(iPath);
	if (!fileStream) {
		cout << "Json File Not Found!!!" << endl;
		success = false;
		return;
	}
	json root;
	fileStream >> root;
	if (root.count("graphFile") != 0) {
		std::string grpgFl = root["graphFile"];
		std::transform(grpgFl.begin(), grpgFl.end(), grpgFl.begin(), ::tolower);
		graphFile = grpgFl;
		cout << "graph file {}" << graphFile << "\n";
	}

	if (root.count("RedcuedToSize") != 0)
		ReducedRankL = root["RedcuedToSize"];

	if (root.count("EpsilonValue") != 0)
		EpsilonVal = root["EpsilonValue"];

	if (root.count("ABSTOLValue") != 0)
		ABSTOL = root["ABSTOLValue"];

	if (root.count("RELTOLValue") != 0)
		RELTOL = root["RELTOLValue"];

	if (root.count("NameOfOutput") != 0)
	{
		std::string OutputName = root["NameOfOutput"];
		std::transform(OutputName.begin(), OutputName.end(), OutputName.begin(), ::tolower);
		NameOfOutput = OutputName;
	}
	if (root.count("LinearOn") != 0)
	{
		std::string key = root["LinearOn"];
		std::transform(key.begin(), key.end(), key.begin(), ::tolower);
		if (key == "false")
			linearOn = false;
		else if (key == "true")
			linearOn = true;
		else
			throw ExceptionError("failed to read LinearOn");
	}
	if (root.count("Pre-DefinedT") != 0)
	{
		std::string key = root["Pre-DefinedT"];
		std::transform(key.begin(), key.end(), key.begin(), ::tolower);
		if (key == "false")
			preDefinedT = false;
		else if (key == "true")
			preDefinedT = true;
		else
			throw ExceptionError("failed to read preDefinedT");
	}

	if (root.count("RandomOn") != 0)
	{
		std::string key = root["RandomOn"];
		std::transform(key.begin(), key.end(), key.begin(), ::tolower);
		if (key == "false")
			RandomOn = false;
		else if (key == "true")
			RandomOn = true;
		else
			throw ExceptionError("failed to read RandomOn");
	}

	if (root.count("SpilitingData") != 0)
	{
		std::string key = root["SpilitingData"];
		std::transform(key.begin(), key.end(), key.begin(), ::tolower);
		if (key == "false")
			SpilitData = false;
		else if (key == "true")
			SpilitData = true;
		else
			throw ExceptionError("failed to read Spiliting Data");
	}

	if (root.count("CrossValOn") != 0)
	{
		std::string key = root["CrossValOn"];
		std::transform(key.begin(), key.end(), key.begin(), ::tolower);
		if (key == "false")
			CrossValOn = false;
		else if (key == "true")
			CrossValOn = true;
		else
			throw ExceptionError("failed to read CrossOn");
	}
	if (root.count("ShuffleData") != 0)
	{
		std::string key = root["ShuffleData"];
		std::transform(key.begin(), key.end(), key.begin(), ::tolower);
		if (key == "false")
			Shuffle = false;
		if (key == "true")
			Shuffle = true;
	}
	if (root.count("NoNodes") != 0)
		NumberofNodes = root["NoNodes"];

	if (root.count("NoFolds") != 0)
		NoFold = root["NoFolds"];

	if (root.count("EtaMin") != 0)
		std::get<0>(Interval) = root["EtaMin"];
	if (root.count("EtaMax") != 0)
		std::get<1>(Interval) = root["EtaMax"];

	if (root.count("GammaMin") != 0)
		std::get<2>(Interval) = root["GammaMin"];
	if (root.count("GammaMax") != 0)
		std::get<3>(Interval) = root["GammaMax"];
	if (root.count("JCMin") != 0)
		std::get<4>(Interval) = root["JCMin"];
	if (root.count("JCMax") != 0)
		std::get<5>(Interval) = root["JCMax"];


	if (root.count("Eta") != 0)
		Eta = root["Eta"];

	if (root.count("JC") != 0)
		JC = root["JC"];

	if (root.count("Gamma") != 0)
		Gamma = root["Gamma"];

	if (root.count("GammaStep") != 0)
		GammaStep = root["GammaStep"];
	if (root.count("EtaStep") != 0)
		EtaStep = root["EtaStep"];
	if (root.count("JcStep") != 0)
		JcStep = root["JcStep"];

	if (root.count("Iteration") != 0)
		iter = root["Iteration"];

	if (root.count("KernalType") != 0)
	{
		std::string kernelName = root["KernalType"];
		std::transform(kernelName.begin(), kernelName.end(), kernelName.begin(), ::tolower);

		if (kernelName == "rbf")
			kernelTyp = TypeForKernel::rbf;
		else if (kernelName == "linear")
			kernelTyp = TypeForKernel::linear;
		else if (kernelName == "poly")
			kernelTyp = TypeForKernel::poly;
		else
			throw ExceptionError("failed to read KernalType");
	}
	if (root.count("Scaling") != 0) {
		std::string ScaleType = root["Scaling"];
		std::transform(ScaleType.begin(), ScaleType.end(), ScaleType.begin(), ::tolower);
		Scaling = ScaleType;
	}
	if (root.count("TerminationCriteria") != 0) {//I have designed 3 termination criteria, 1. "AllThree", 2."Pri-Dual", and 3. "onlyVj_VjAll"
		std::string Stopping = root["TerminationCriteria"];
		std::transform(Stopping.begin(), Stopping.end(), Stopping.begin(), ::tolower);
		TerminationCriteria = Stopping;
	}

	if (root.count("Graph") != 0)
	{
		std::string grType = root["Graph"];
		std::transform(grType.begin(), grType.end(), grType.begin(), ::tolower);
		GraphType = grType;
		graph = returnGraph(GraphType);
		//if (grType == "hardcoded")
		//	graph = Graph::hardCoded;
		//else if (grType == "line")
		//	graph = Graph::line;
		//else if (grType == "cycle")
		//	graph = Graph::cycle;
		//else if (grType == "self")
		//	graph = Graph::self;
		//else if(grType == "threeregular")
		//	graph = Graph::ThreeRegular;
		//else if (grType == "fiveregular")
		//	graph = Graph::FiveRegular;
		//else if (grType == "comb")
		//	graph = Graph::Comb;
		//else if (grType == "island")
		//	graph = Graph::Island;
		//else if (grType == "unconnected")
		//	graph = Graph::unconnected;
		//else if (grType == "unconnislands")
		//	graph = Graph::unConnIslands;
		//else if (grType == "matlab_threereg")
		//	graph = Graph::Matlab_ThreeReg;
		//else if (grType == "matlab_fivereg")
		//	graph = Graph::Matlab_FiveReg;
		//else if (grType == "matlab_sevenreg")
		//	graph = Graph::Matlab_SevenReg;
		//else if (grType == "matlab_ninereg")
		//	graph = Graph::Matlab_NineReg;
		//else
		//	throw ExceptionError("failed to read node structure");
	}
	if (root.count("ADMMParamUpdate") != 0) {
		std::string key = root["ADMMParamUpdate"];
		std::transform(key.begin(), key.end(), key.begin(), ::tolower);
		if (key == "false")
			ADMMParaUpdate = false;
		else if (key == "true")
			ADMMParaUpdate = true;
		else
			throw ExceptionError("failed to read ADMMParamUpdate");
	}

	success = true;
}

//void ProblemSetting::save2File(const fs::path & iPath)
//{
//
//}
