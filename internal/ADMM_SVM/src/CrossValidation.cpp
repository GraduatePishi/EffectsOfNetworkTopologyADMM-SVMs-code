#include "CrossValidation.h"
#include <iostream>
#include "Solver.h"
#include <fstream>
#include "utils.h"
#include "FoldData.h"

std::tuple<double, double, double, double> Cross_validation(int iNoFold, ProblemStatment& iData, ProblemStatment& iTest, DimReductionMat iT, TypeForKernel iKrnlTyp, int iNoNodes, int iNodeID, bool iRandomOn, double iEta, double /*iGamma*/, std::string iDataSetName, CrossValIntrval iInterval, Graph igraph, bool iADMM_update, std::string iGraphFile, int iworld_size, int iworld_rank, bool iMPIflag, double iEpsilonVal, double iABSTOL, double iRELTOL, std::string iTerminationCriteria, bool iCross_val, double iGammaStep,double iEtaStep, double iJcStep) {
	double Gamma{ -1 };
	double JC{ -1 };
	double Eta{ -1 };
	double  HighestHitRate{ -1 };
	size_t Elapsed_solver_time = 0;

	//TestTrainData DividedMat = FoldData(iNoFold, iDataMat, iLabel, iRandomOn);//Dividing the whole Data matrix and label into training and testing matrices and labels and save them as tuple

	nonlinear_svm solver;
	double MaxGamma = std::get<3>(iInterval);
	double MinGamma = std::get<2>(iInterval);
	double MaxEta = std::get<1>(iInterval);
	double MinEta = std::get<0>(iInterval);
	double _gamma = MinGamma;
	double _eta = MinEta;
	double MaxJC = std::get<5>(iInterval);
	double MinJC = std::get<4>(iInterval);
	double _jc = MinJC;
	int GammStepSize = 0;
	int JcStepSize = 0;
	int EtaStepSize = 0;
	bool manyFile = false;
	//double GammaStep{ 0.0001 }, EtaStep{ 10 }, JcStep{ 1 };
	std::vector<double> GammaVec, JCVec, EtaVec;

	while (_eta <= MaxEta)
	{
		EtaStepSize = EtaStepSize + 1;
		EtaVec.push_back(_eta); //puting each step gamma in a vector
		_eta = _eta * iEtaStep; //change of variable in each step
	}
	while (_gamma <= MaxGamma)
	{
		GammStepSize = GammStepSize + 1;
		GammaVec.push_back(_gamma); //puting each step gamma in a vector
		_gamma = _gamma * iGammaStep; //change of variable in each tep
	}
	while (_jc <= MaxJC)
	{
		JcStepSize = JcStepSize + 1;
		JCVec.push_back(_jc); //puting each step JC in a vector
		_jc = _jc * iJcStep;
	}

	//std::cout << "im node " << iworld_rank << " with JCVec: ";
	//std::copy(JCVec.begin(), JCVec.end(), std::ostream_iterator<double>(std::cout, " "));
	//std::cout << std::endl;
	//std::cout << "im node " << iworld_rank << " with GammaVec: ";
	//std::copy(GammaVec.begin(), GammaVec.end(), std::ostream_iterator<double>(std::cout, " "));
	//std::cout << std::endl;
	//std::cout << "im node " << iworld_rank << " with EtaVec: ";
	//std::copy(EtaVec.begin(), EtaVec.end(), std::ostream_iterator<double>(std::cout, " "));
	//std::cout << std::endl;

	DataMatrix AverageHitMat = DataMatrix::Zero(JcStepSize + 1, GammStepSize + 1), HitRateMat(JcStepSize + 1, GammStepSize + 1);//Vector of matrices, saves hitrates in each fold
	std::stringstream name; //saving the output in a file
	std::fstream outfile;
	if (manyFile == true) {
		name << "HitRateMatrices_" << iDataSetName << "_" << iNodeID << "Node(s).output";
		outfile.open(name.str(), std::fstream::out);
	}
	else {
		if (iworld_rank == 0) {
			name << "AverageHitRate_" << iDataSetName << "_AllNodes" << ".output";
			outfile.open(name.str(), std::fstream::out);
		}
	}
	
	int iter = -1;
	Column col = (Column)iData.Data.cols();
	//TestTrainData PartData = FoldData(iNoFold, iDataMat, iLabel, iRandomOn);//Dividing the whole Data matrix and label into training and testing matrices and labels and save them as tuple
	int NoOfFolds = 1;
	//for (int i = 0; i < iNoFold; i++)
	for (int i = 0; i < NoOfFolds; i++)
	{
		int countCol{ 0 }, countRow{ 0 };
		for (double etaTmp : EtaVec)
		{
			AverageHitMat = DataMatrix::Zero(JcStepSize + 1, GammStepSize + 1); //reset average hitrate for new Eta and Fold
			HitRateMat = DataMatrix::Ones(JcStepSize + 1, GammStepSize + 1)*(-1); //reset the hitrate matrix for all eta and folds

			for (int u = 0; u < GammStepSize; u++)
				AverageHitMat(0, u + 1) = GammaVec[u];
			for (int y = 0; y < JcStepSize; y++)
				AverageHitMat(y + 1, 0) = JCVec[y];

			countCol = 0;

			for (double GammaTmp : GammaVec)
			{
				countCol = countCol + 1;
				countRow = 0;
				/*for (JCtmp =std::get<4>(iInterval); JCtmp<= std::get<5>(iInterval); JCtmp *= JcStep)*/
				for (double JCtmp : JCVec)
				{

					countRow = countRow + 1;
					HitRateMat(0, countCol) = GammaTmp;
					HitRateMat(countRow, 0) = JCtmp;

					try {
						auto start = TimeType::now();
						//auto&& solving = solver.solve(iter, std::get<2>(PartData[i]), std::get<3>(PartData[i]), std::get<0>(PartData[i]), std::get<1>(PartData[i]), iT, iKrnlTyp, iNoNodes, JCtmp, etaTmp, GammaTmp, iRandomOn, igraph, iADMM_update, iGraphFile, iworld_size, iworld_rank, iEpsilonVal, iABSTOL, iRELTOL, iTerminationCriteria);
						auto&& solving = solver.solve(iter, iData, iTest, iT, iKrnlTyp, iNoNodes, JCtmp, etaTmp, GammaTmp, iRandomOn, igraph, iADMM_update, iGraphFile, iworld_size, iworld_rank, iEpsilonVal, iABSTOL, iRELTOL, iTerminationCriteria, iCross_val, "NotForCrossVal");

						auto end = TimeType::now();
						DurationTime time = (end - start);
						Elapsed_solver_time = std::chrono::duration_cast<TimePrec>(time).count();

						if (std::get<3>(solving[0]) > HighestHitRate) {
							HighestHitRate = std::get<3>(solving[0]);
							Gamma = GammaTmp;
							JC = JCtmp;
							Eta = etaTmp;
						}
						HitRateMat(countRow, countCol) = std::get<3>(solving[0]);
						/* getting data from all nodes*/
						double AverageHitRate =0.0;
						MPI_Allreduce(&std::get<3>(solving[0]), &AverageHitRate, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
						AverageHitRate = AverageHitRate / iworld_size ;
						AverageHitMat(countRow, countCol) = AverageHitRate;

						if (manyFile==false){ 
							if (iworld_rank == 0) { //the reason that we write results here is that we can see the reesults before the hitrate matrix is completely full
								outfile << "++++++++++++++++++++++++++" << "\n";
								outfile << "---Eta: " << etaTmp << " ---GammaTmp: " << GammaTmp << "---JCtmp " << JCtmp << "\n";// " ---HighestHitRate " << HighestHitRate << "\n";
								std::cout << "---Eta: " << etaTmp << " ---GammaTmp: " << GammaTmp << "---JCtmp " << JCtmp << "\n";//" ---HighestHitRate " << HighestHitRate << "\n";

								outfile << "AverageHitMat : " <<  "\n" << AverageHitMat << "\n";
							}
						}
						else {
							if (iworld_rank == 0)
								std::cout << " --Eta: " << etaTmp << " ---GammaTmp: " << GammaTmp << "---JCtmp " << JCtmp << " ---HighestHitRate " << HighestHitRate << "\n";
							outfile << "++++++++++++++++++++++++++" << "\n";
							outfile << "HitRate for node " << iNodeID << " is : " << std::get<3>(solving[0]) << "\n";
							outfile << " --Eta: " << etaTmp << " ---GammaTmp: " << GammaTmp << "---JCtmp " << JCtmp << " ---rank " << iNodeID << "\n";
							outfile << "Highest HitRate for node " << iNodeID << " is : " << HighestHitRate << "\n";
							outfile << "HitRate for node " << iNodeID << "\n" << HitRateMat << "\n";
						}
					}
					catch (...) {
						std::cout << "************** An exception in Cross-validation for GammaTmp= " << GammaTmp << ",  JCtmp= " << JCtmp << ",  etaTmp= " << etaTmp << "\n";
						//outfile << "************** An exception is occured for GammaTmp= " << GammaTmp << ",  JCtmp= " << JCtmp << ",  etaTmp= " << etaTmp << "\n";
					}
				}  //End of JC loop!
			}  //End of Gamma loop!
		}  //End of Eta loop!
	}  //End of fold loop
	outfile.close();
	return { JC, Gamma,Eta, HighestHitRate };
}

//std::tuple<double, double, double, TestTrainData, double> Cross_validation_MPI(int iNoFold, DataMatrix iDataMat, LabelVector iLabel, DimReductionMat iT, TypeForKernel iKrnlTyp, int iNoNodes, bool iRandomOn, double iEta, double /*iGamma*/, std::string iDataSetName, CrossValIntrval iInterval, Graph igraph, bool iADMM_update, std::string iGraphFile, int iworld_size, int iworld_rank, bool iMPIflag) {
//	double Gamma{ -1 };
//	double JC{ -1 };
//	double Eta{ -1 };
//	double  HitRate{ -1 };
//	double AverageHitRate{ 0 };
//	double Elapsed_solver_time = 0.0;
//
//	nonlinear_svm solver;
//	double MaxGamma = std::get<3>(iInterval);
//	double MinGamma = std::get<2>(iInterval);
//	double _gamma = MinGamma;
//	double MaxJC = std::get<5>(iInterval);
//	double MinJC = std::get<4>(iInterval);
//	double _jc = MinJC;
//	int GammStepSize = 0;
//	int JcStepSize = 0;
//	double GammaStep{ 5e-6 }, EtaStep{ 10 }, JcStep{ 5 };
//	std::vector<double> GammaVec, JCVec;
//
//	//std::cout << "Gamma  " ;
//	while (_gamma <= MaxGamma)
//	{
//		//std::cout <<  _gamma << "  ";
//		GammStepSize = GammStepSize + 1;
//		GammaVec.push_back(_gamma); //puting each step gamma in a vector
//		_gamma = _gamma + GammaStep;
//	}
//	//std::cout << "\n" << "JC   " << "\n";
//	while (_jc <= MaxJC)
//	{
//		//std::cout << _jc << "\n";
//		JcStepSize = JcStepSize + 1;
//		JCVec.push_back(_jc); //puting each step JC in a vector
//		_jc = _jc + JcStep;
//	}
//
//	//std::cout << "\n" << "Cross-Validation Matrix of size (" << JcStepSize << ", " << GammStepSize << ") is created." << "\n\n";
//	DataMatrix AverageHitMat = DataMatrix::Zero(JcStepSize + 1, GammStepSize + 1);
//
//
//	std::stringstream name; //saving the output in a file
//	name << "HitRateMatrices_" << iDataSetName << "_" << iNoNodes << "Node(s).output";
//	std::fstream outfile;
//	outfile.open(name.str(), std::fstream::out);
//
//	std::vector<DataMatrix> HitRateMat(iNoNodes);//Vector of matrices, saves hitrates in each fold
//
//	int iter = -1;
//	Column col = (Column)iDataMat.cols();
//	//if (world_rank==0) {
//		TestTrainData PartData = FoldData(iNoFold, iDataMat, iLabel, iRandomOn);//Dividing the whole Data matrix and label into training and testing matrices and labels and save them as tuple
//	//}
//	for (int i = 0; i < iNoFold; i++)
//	{
//		double JCtmp{ -1 }, GammaTmp{ -1 }, etaTmp{ -1 };
//		int countCol{ 0 }, countRow{ 0 };
//
//		for (etaTmp = std::get<0>(iInterval); etaTmp <= std::get<1>(iInterval); etaTmp += EtaStep)
//		{
//			AverageHitMat = DataMatrix::Zero(JcStepSize + 1, GammStepSize + 1);
//			for (int u = 0; u < GammStepSize; u++)
//				AverageHitMat(0, u + 1) = GammaVec[u];
//			for (int y = 0; y < JcStepSize; y++)
//				AverageHitMat(y + 1, 0) = JCVec[y];
//
//			for (int d = 0; d < iNoNodes; d++)
//				HitRateMat[d] = DataMatrix::Ones(JcStepSize + 1, GammStepSize + 1)*(-1); //reset the hitrate matrix for all eta and folds
//
//			countCol = 0;
//			for (GammaTmp = std::get<2>(iInterval); GammaTmp <= std::get<3>(iInterval); GammaTmp += GammaStep)
//			{
//				countCol = countCol + 1;
//				countRow = 0;
//				for (JCtmp = std::get<4>(iInterval); JCtmp <= std::get<5>(iInterval); JCtmp += JcStep)
//				{
//					countRow = countRow + 1;
//					for (int z = 0; z < iNoNodes; z++) {
//						HitRateMat[z](0, countCol) = GammaTmp;
//						HitRateMat[z](countRow, 0) = JCtmp;
//					}
//
//					std::cout << "------------- Fold: " << i << " -------- Eta: " << etaTmp << " ----------" << "\n";
//
//					try {
//						auto start = TimeType::now();
//						auto&& solving = solver.solve(iter, std::get<2>(PartData[i]), std::get<3>(PartData[i]), std::get<0>(PartData[i]), std::get<1>(PartData[i]), iT, iKrnlTyp, iNoNodes, JCtmp, etaTmp, GammaTmp, iRandomOn, igraph, iADMM_update, iGraphFile, iworld_size, iworld_rank, iMPIflag);
//						auto end = TimeType::now();
//						DurationTime time = (end - start);
//						Elapsed_solver_time = (double)time.count();
//						int node = 0; //Nr of Nodes
//						AverageHitRate = 0.0;
//						/* getting data from all nodes*/
//						for (auto& slv : solving) {
//							AverageHitRate = AverageHitRate + std::get<3>(slv);
//							if (std::get<3>(slv) > HitRate) {
//								HitRate = std::get<3>(slv);
//								Gamma = GammaTmp;
//								JC = JCtmp;
//								Eta = etaTmp;
//							}
//							HitRateMat[node](countRow, countCol) = std::get<3>(slv);
//							//outfile << "\n" << "----------Fold " << i + 1 << " ----- Eta " << etaTmp << "  ---- Node " << node + 1 << "\n";
//							//outfile << HitRateMat[node] << "\n";
//							node = node + 1;
//						}
//						AverageHitRate = AverageHitRate / iNoNodes;
//						AverageHitMat(countRow, countCol) = AverageHitRate;
//						outfile << "\n" << "----------Fold " << i + 1 << " ----- Eta " << etaTmp << ",GammaTmp = " << GammaTmp << ", JCtmp = " << JCtmp << "\n";
//						std::cout << "Average HitRate for all nodes" << "\n" << AverageHitMat << "\n";
//						outfile << "Average HitRate for all nodes" << "\n" << AverageHitMat << "\n";
//					}
//					catch (...) {
//						std::cout << "************** An exception is occured for GammaTmp= " << GammaTmp << ",  JCtmp= " << JCtmp << ",  etaTmp= " << etaTmp << "\n";
//						outfile << "************** An exception is occured for GammaTmp= " << GammaTmp << ",  JCtmp= " << JCtmp << ",  etaTmp= " << etaTmp << "\n";
//					}
//				}  //End of JC loop!
//			}  //End of Gamma loop!
//			//for (int node = 0; node < iNoNodes; node++) {
//			//	outfile << "\n" << "----------Fold " << i + 1 << " ----- Eta " << etaTmp << "  ---- Node " << node + 1 << "\n";
//			//	outfile << HitRateMat[node] << "\n\n";
//			//}
//		}  //End of Eta loop!
//	}  //End of fold loop
//	outfile.close();
//	return { JC, Gamma,Eta, PartData, HitRate };
//}