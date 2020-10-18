#include "Printer.h"
#include <iostream>
#include <fstream>
#ifdef _MSC_VER
#include <windows.h>
#endif
namespace l
{
	void SetColor(int value=7) {
		/*1: Blue
			2 : Green
			3 : Cyan
			4 : Red
			5 : Purple
			6 : Yellow(Dark)
			7 : Default white
			8 : Gray / Grey
			9 : Bright blue
			10 : Brigth green
			11 : Bright cyan
			12 : Bright red
			13 : Pink / Magenta
			14 : Yellow
			15 : Bright white*/
#ifdef _MSC_VER
		SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), (WORD)value);
#endif
	}

	class PrivPrinter : public Printer
	{
	public:
		virtual void setStyle(Style iStyle, std::string iName) override final
		{
			mStyle = iStyle;
			if (mStyle == Style::File)
			{
				//iName.append(".Output");
				/*mFile.open("svm.output_global", std::fstream::out);*/
				mFile.open(iName, std::fstream::out);
			}
		}

		virtual void log(std::string iString) override final
		{
			if (mStyle == Style::Console)
				std::cout << iString << "\n";
			if (mStyle == Style::File)
				mFile << "Log  << " << iString << "\n";
		}
		virtual void warn(std::string iString) override final
		{
			
			if (mStyle == Style::Console)
			{
				SetColor(14);
				std::cout << iString << "\n";
				SetColor();
			}
			if (mStyle == Style::File)
				mFile << "Warn << " << iString << "\n";
		}
		virtual void err(std::string iString)override final
		{
			if (mStyle == Style::Console)
			{
				SetColor(12);
				std::cout << iString << "\n";
				SetColor();
			}
			if (mStyle == Style::File)
				mFile << "Err  << " << iString << "\n";
		}

		~PrivPrinter()
		{
			mFile.close();
		}
	private:
		Style mStyle;
		std::fstream mFile;
	};


	//void Printer::operator<<(std::ostream & stream)
	//{
	//	std::cout << stream;
	//}

	PrivPrinter P;

	Printer & getPrint()
	{
		return P;
	}
}
