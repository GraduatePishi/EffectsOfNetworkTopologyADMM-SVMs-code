#include "ADMM_penalty.h"
#include <iostream>
#include "NodeSolver.h"
//
//using namespace std;
using ExceptionError = std::runtime_error;
//void ADMMMpenaltyUpdate(LabelVector& iY, LabelVector& iF, double PrimalRes, double DoualRes,  ADMMpenalty iMethod, int it, int itMax, std::pair<std::unique_ptr<nodesolver>, NodeData>& iNode, std::vector<std::pair<std::unique_ptr<nodesolver>, NodeData>>& iNodeInfo) {
//
//	double Mu{ 10.0 }, Tu{ 1.0 };
//	if (ADMMpenalty::Formula4 == iMethod) {
//		if (PrimalRes > Mu*DoualRes) {
//			for (int i = 0; i < (int)iNode.first->getNeigbors().size(); i++) {
//				iNode.first->EtaVal[i] = iNode.first->EtaVal[i] * (1.0 + Tu);
//			}
//		}
//		else if (DoualRes > Mu*PrimalRes) {
//			for (int j = 0; j < (int)iNode.first->getNeigbors().size(); j++) {
//				iNode.first->EtaVal[j] = iNode.first->EtaVal[j] / (1.0 + Tu);
//			}
//		}
//
//		//int EtaIndx = 0;
//		//for (auto const& NeighIndex : iNode.first->getNeigbors()) {
//
//		//	if (PrimalRes > Mu*DoualRes) {
//		//		iNode.first->EtaVal[EtaIndx] = iNode.first->EtaVal[EtaIndx] * (1.0 + Tu);
//		//		EtaIndx = EtaIndx + 1;
//		//	}
//		//	else if (DoualRes > Mu*PrimalRes) {
//		//		iNode.first->EtaVal[EtaIndx] = iNode.first->EtaVal[EtaIndx] / (1.0 + Tu);
//		//		EtaIndx = EtaIndx + 1;
//		//	}
//		//	else {
//		//		EtaIndx= EtaIndx +1;
//		//	}
//		//}
//		//if (iNode.first->getNeigbors().size() != EtaIndx+1) 
//		//	throw ExceptionError("The Eta value vector size does not match with neighbor numbers!!!");
//
//		if (GoodPrinting) 
//			cout << " Eta in Eta update: " << iNode.first->EtaVal.transpose() << "\n";
//	}
//
//	if (ADMMpenalty::ADMMAP == iMethod) {
//
//		double fiMax{ -1.0* std::numeric_limits<double>::infinity() }, fiMin{ std::numeric_limits<double>::infinity() };
//
//		const WeighVector& wj = iNode.first->getWTilde();
//		double fi = wj.dot(wj); //Objective function of node itself
//		int EtaIndx = 0;
//		if (it < itMax) {
//			//finding the max and min of one neighbor's objective func
//			for (auto const& NeighIndex : iNode.first->getNeigbors()) {
//				const WeighVector& wi = iNodeInfo[NeighIndex].first->getWTilde();
//				double fij = wi.dot(wi);
//				//iNodeInfo[NeighIndex].first->ObjectivFunc = wi.dot(wi);
//				fiMax = std::max(fi, fij);
//				fiMin = std::min(fi, fij);
//				//iNodeInfo[NeighIndex].first->kj = ((iNodeInfo[NeighIndex].first->ObjectivFunc - fiMin) / (fiMax - fiMin)) + 1;
//				double denominator = (fiMax - fiMin);
//				double kii = ((fi - fiMin) / denominator) + 1.0;
//				double kij = ((fij - fiMin) / denominator) + 1.0;
//
//				double Tu_ij = (kii / kij) - 1.0;
//				//std::cout << "iNodeInfo[iNeighIndex].first->getEta()" <<"\n"<< iNodeInfo[iNeighIndex].first->getEta()[iNeighIndex] << "\n";
//				//iNodeInfo[NeighIndex].first->EtaVal[NeighIndex] = iNode.first->mEta* (1 + Tu_ij);
//				//iNode.first->EtaVal[EtaIndx] = iNode.first->EtaVal[EtaIndx] * (1.0 + Tu_ij);
//				iNode.first->EtaVal[EtaIndx] = iNode.first->mEta * (1.0 + Tu_ij);
//				EtaIndx = EtaIndx + 1;
//			}
//			iNode.first->EtaAverage = iNode.first->EtaVal.sum() / (int)iNode.first->getNeigbors().size();
//		}
//		else
//		{
//			iNode.first->EtaVal = Eigen::VectorXd::Ones((int)iNode.first->getNeigbors().size())*iNode.first->mEta; // change to the initial Eta
//		}
//	}
//
//	if (ADMMpenalty::ADMMCombined == iMethod) { //With primal objective function
//		double fiMax{ -1.0* std::numeric_limits<double>::infinity() }, fiMin{ std::numeric_limits<double>::infinity() };
//
//		const WeighVector& wj = iNode.first->getWTilde();
//		double fi = wj.dot(wj); //Objective function of node itself
//		int EtaIndx = 0;
//		double fij;
//
//		if (it < itMax) {
//			double kii, kij, Tu_ij;
//
//			for (auto const& NeighIndex : iNode.first->getNeigbors()) {
//				const WeighVector& wi = iNodeInfo[NeighIndex].first->getWTilde();
//				fij = wi.dot(wi);
//				fiMax = std::max(fi, fij);
//				fiMin = std::min(fi, fij);
//				kii = ((fi - fiMin) / (fiMax - fiMin)) + 1.0;
//				kij = ((fij - fiMin) / (fiMax - fiMin)) + 1.0;
//				Tu_ij = (kii / kij) - 1.0;
//				if (PrimalRes > Mu*DoualRes) {
//					iNode.first->EtaVal[EtaIndx] = iNode.first->EtaVal[EtaIndx] * (1.0 + Tu_ij)*2.0;
//				}
//				else if ((DoualRes > Mu*PrimalRes)) {
//					iNode.first->EtaVal[EtaIndx] = iNode.first->EtaVal[EtaIndx] * (1.0 + Tu_ij)*0.5;
//				}
//
//				EtaIndx = EtaIndx + 1;
//			}
//			iNode.first->EtaAverage = iNode.first->EtaVal.sum() / (int)iNode.first->getNeigbors().size();
//		}
//	}
	//The below methd is based on the dual objective that seems wrong
	//if (ADMMpenalty::ADMMCombined == iMethod) {
	//	double fiMax{ -1.0* std::numeric_limits<double>::infinity() }, fiMin{ std::numeric_limits<double>::infinity() };

	//	WeighVector& wj = iNode.first->getWTilde();
	//	double fi = wj.dot(wj); //Objective function of node itself
	//	int EtaIndx = 0;
	//	double fii; //The Objective function of the nodej
	//	double fij;
	//	
	//	if (it < itMax) {
	//		double kii, kij, Tu_ij;
	//		FuncVariables tmpFuncData(iNode.first->FixedMatrix, iNode.first->mkernelTXj_mkernelTildaTXj, iNode.first->Yj, iNode.first->mfTilde, iNode.first->mhj, iNode.first->getNeigbors().size(), iNode.first->EtaAverage);
	//		//fii = -1*objeFunc(iNode.first->getLambdaj(), &tmpFuncData);
	//		fii = iNode.first->ObjectivFunc;
	//		cout << "fii : " << fii << "\n";
	//		//cout << "ObjectivFunc : " << iNode.first->ObjectivFunc << "\n";
	//		
	//		if (PrimalRes > Mu*DoualRes) {
	//			cout << "PrimalRes > Mu*DoualRes" << "\n";
	//			for (auto const& NeighIndex : iNode.first->getNeigbors()) {
	//		
	//				//auto& wi = iNodeInfo[NeighIndex].first->getWTilde();
	//				//iNodeInfo[NeighIndex].first->ObjectivFunc = wi.dot(wi);
	//				//fiMax = std::max(fj, iNodeInfo[NeighIndex].first->ObjectivFunc);
	//				//fiMin = std::min(fj, iNodeInfo[NeighIndex].first->ObjectivFunc);
	//				fij =-1* objeFunc(iNodeInfo[NeighIndex].first->getLambdaj(), &tmpFuncData);
	//				cout << "fij : " << fij << "\n";
	//				fiMax = std::max(fii, fij);
	//				cout << "fiMax : " << fiMax << "\n";
	//				fiMin = std::min(fii, fij);
	//				cout << "fiMin : " << fiMin << "\n";
	//				//iNodeInfo[NeighIndex].first->kj = ((iNodeInfo[NeighIndex].first->ObjectivFunc - fiMin) / (fiMax - fiMin)) + 1.0;
	//				kii = ((fii - fiMin) / (fiMax - fiMin)) + 1.0;
	//				kij = ((fij - fiMin) / (fiMax - fiMin)) + 1.0;
	//				Tu_ij = (kii / kij) - 1.0;
	//				cout << "kii : " << kii << "\n";
	//				cout << "kij : " << kij << "\n";
	//				cout << "Tu_ij : " << Tu_ij << "\n";
	//				iNode.first->EtaVal[EtaIndx] = iNode.first->EtaVal[EtaIndx] * (1.0 + Tu_ij)*2.0;
	//				if (GoodPrinting)
	//					cout << "mEta[ " << EtaIndx << " ] " <<"\n"<< iNode.first->EtaVal[EtaIndx] << "\n";
	//				EtaIndx++;

	//			}
	//			iNode.first->EtaAverage = iNode.first->EtaVal.sum() / (int)iNode.first->getNeigbors().size();
	//		}
	//		else if (DoualRes > Mu*PrimalRes)
	//		{
	//			cout << "DoualRes > Mu*PrimalRes" << "\n";
	//			for (auto const& NeighIndex : iNode.first->getNeigbors()) {
	//				fij = -1*objeFunc(iNodeInfo[NeighIndex].first->getLambdaj(), &tmpFuncData);
	//				cout << "fij : " << fij << "\n";
	//				fiMax = std::max(fii, fij);
	//				cout << "fiMax : " << fiMax << "\n";
	//				fiMin = std::min(fii, fij);
	//				cout << "fiMin : " << fiMin << "\n";
	//				double kii = ((fii - fiMin) / (fiMax - fiMin)) + 1.0;
	//				double kij = ((fij - fiMin) / (fiMax - fiMin)) + 1.0;
	//				double Tu_ij = (kii / kij) - 1.0;
	//				cout << "Tu_ij : " << Tu_ij << "\n";
	//				//auto& wi = iNodeInfo[NeighIndex].first->getWTilde();
	//				//iNodeInfo[NeighIndex].first->ObjectivFunc = wi.dot(wi);
	//				//fiMax = std::max(fi, iNodeInfo[NeighIndex].first->ObjectivFunc);
	//				//fiMin = std::min(fi, iNodeInfo[NeighIndex].first->ObjectivFunc);
	//				//iNodeInfo[NeighIndex].first->kj = ((iNodeInfo[NeighIndex].first->ObjectivFunc - fiMin) / (fiMax - fiMin)) + 1;
	//				//double ki = ((fi - fiMin) / (fiMax - fiMin)) + 1.0;
	//				//double Tu_ij = (ki / iNodeInfo[NeighIndex].first->kj) - 1.0;
	//				iNode.first->EtaVal[EtaIndx] = iNode.first->EtaVal[EtaIndx] * (1.0 + Tu_ij) * 0.5;
	//				if (GoodPrinting)
	//					cout << "mEta[ " << EtaIndx << " ] " <<"\n" << iNode.first->EtaVal[EtaIndx] << "\n";
	//				EtaIndx++;
	//			}
	//			iNode.first->EtaAverage = iNode.first->EtaVal.sum() / (int)iNode.first->getNeigbors().size();
	//		}
	//	}
	//	else
	//	{
	//		cout << "Primal = dual" << "\n";
	//		for (auto const& NeighIndex : iNode.first->getNeigbors()) {

	//			iNode.first->EtaVal[EtaIndx] = iNode.first->mEta; // change to the initial Eta
	//			if (GoodPrinting)
	//				cout << "mEta[ " << EtaIndx << " ]" << "\n" << iNode.first->EtaVal[EtaIndx] << "\n";
	//			EtaIndx++;

	//		}
	//		iNode.first->EtaAverage = iNode.first->EtaVal.sum() / (int)iNode.first->getNeigbors().size();
	//	}
	//}

//}
