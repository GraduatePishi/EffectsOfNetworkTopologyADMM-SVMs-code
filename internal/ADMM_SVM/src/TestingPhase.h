#include "Kernel.h"
#include "FileReader.h"
#include "Types.h"

ProblemStatment ReadingTestData(fs::path iTrainfilePath, std::string iScaling, Eigen::VectorXd& iMinScaling, Eigen::VectorXd& iMaxScaling, bool iPCA_flag, int iRank);
hitRate TestingPhase(const ProblemStatment& iTestProblem, const  Eigen::MatrixXd& iT, const DataMatrix& iDataMat, const KrnlMatrx& iT_norm, const KrnlMatrx& iX_norm, fs::path iTrainfilePath, std::string iScaling, Eigen::VectorXd& iMinScaling, Eigen::VectorXd&iMaxScaling, bool iPCA_flag, const  HyperPlaneScalars& ia, WeightScalar& ib, const  HyperPlaneScalars& ic, TypeForKernel iKernelType, int iRank, double iGamma);
std::pair<Eigen::MatrixXd, Eigen::MatrixXd> getTestKernel(const RowMajorMatirx& iTestData, const DataMatrix& PartialDataMat, const  DimReductionMat& iT, KrnlMatrx iNormT_squared, KrnlMatrx iNorm_Xsquared, double iGamma);