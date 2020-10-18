#pragma once
#include <string>
#include <fmt/format.h>
#include <sstream>
#include "utils.h"
namespace l
{
	enum class Style
	{
		Console,
		File,
		Both
	};
	enum class LogLevel
	{
		All,
		WarnAndErr,
		Err,
		None
	};

	class Printer
	{
	public:
		


		virtual void setStyle(Style iStyle, std::string iName) = 0;
		virtual void log(std::string iString) = 0;
		virtual void warn(std::string iString) = 0;
		virtual void err(std::string iString) = 0;
	};

	Printer& getPrint();
	inline void setStyle(Style iStyle, std::string iName)
	{
		getPrint().setStyle(iStyle, iName);
	}

	template <typename... Args>
	void log(const char *format, const Args & ... args) {
		getPrint().log(fmt::format(format, args...));
	}
	template <typename... Args>
	void warn(const char *format, const Args & ... args) {
		getPrint().warn(fmt::format(format, args...));
	}
	template <typename... Args>
	void err(const char *format, const Args & ... args) {
		getPrint().err(fmt::format(format, args...));
	}


	template<class _InIt> inline
	std::string PrintVec(_InIt& start, _InIt& end)
	{
		std::stringstream ss;
		std::copy(start, end, std::ostream_iterator<typename _InIt::value_type>(ss, ", "));
		return ss.str();
	}
	template<class VecT> inline
	std::string PrintVec(const VecT& start)
	{
		std::stringstream ss;
		std::copy(cbegin(start), cend(start), std::ostream_iterator<typename VecT::value_type>(ss, ", "));
		return ss.str();
	}
}