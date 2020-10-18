#pragma once
#include "Solver.h"
#include "NodeSolver.h"
#include "Types.h"

enum class ADMMpenalty {
	Formula4,
	ADMMAP, //ADMM-AP : Adaptive PEnalty: Fast ADMM Algorithm paper
	ADMMCombined
};
class nodesolver;
struct NodeData;
void ADMMMpenaltyUpdate(LabelVector& iY, LabelVector& iF, double PrimalRes, double DoualRes, ADMMpenalty iMethod, int it, int itMax, std::pair<std::unique_ptr<nodesolver>, NodeData>& iNode, std::vector<std::pair<std::unique_ptr<nodesolver>, NodeData>>& iNodeInfo);
//double ADMMMpenaltyUpdate_Dual(HyperPlaneScalars& iLambdaj, HyperPlaneScalars& iLambdai, HyperPlaneScalars& iLambdaiOld, double Eta, ADMMpenalty iMethod, const int iNeighIndex, int it, int itMax, std::pair<std::unique_ptr<nodesolver>, NodeData>& iNode, std::vector<std::pair<std::unique_ptr<nodesolver>, NodeData>>& iNodeInfo);
