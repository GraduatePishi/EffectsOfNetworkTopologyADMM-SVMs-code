#pragma once
#include <iterator>
#include <Eigen/Core>
#include "Types.h"

template <typename T>
auto begin(Eigen::MatrixBase<T>& iVec)
{
	EIGEN_STATIC_ASSERT_VECTOR_ONLY(T);
	return iVec.derived().data();
}

template <typename T>
auto end(Eigen::MatrixBase<T>& iVec)
{
	EIGEN_STATIC_ASSERT_VECTOR_ONLY(T);
	return iVec.derived().data() + iVec.derived().size();
}

template <typename T>
auto cbegin(const Eigen::MatrixBase<T>& iVec)
{
	EIGEN_STATIC_ASSERT_VECTOR_ONLY(T);
	return iVec.derived().data();
}

template <typename T>
auto cend(const Eigen::MatrixBase<T>& iVec)
{
	EIGEN_STATIC_ASSERT_VECTOR_ONLY(T);
	return iVec.derived().data() + iVec.derived().size();
}

template <typename type>
constexpr bool isZero(type iVal) {
	return iVal < std::numeric_limits<type>::epsilon() * 10;
}
template <typename type>
constexpr bool isGreater(type iVal) {
	return iVal > std::numeric_limits<type>::epsilon() * 10;
}
//int getNumSV(Eigen::VectorXd const& iLambda);
static inline size_t getNumSV(Eigen::VectorXd const& iLambda)
{
	size_t mNumSV = std::count_if(cbegin(iLambda), cend(iLambda), [](auto& iV) {return iV > 1e-8; });
	return (mNumSV * 100) / iLambda.size();
}

template <typename T>
double CalcTime(const T& StartTime, const T& EndTime)
{
	DurationTime TimeDuration = (EndTime - StartTime);
	double Elapsed_time_for_entire_solver = (double)std::chrono::duration_cast<TimePrec>(TimeDuration).count();

	return Elapsed_time_for_entire_solver;
}

constexpr TimePrec CalcDiffTime_TimePrec(TimePointType const& StartTime, TimePointType const& EndTime) {
	return std::chrono::duration_cast<TimePrec>(EndTime - StartTime);
}

double MachineEpsilon();