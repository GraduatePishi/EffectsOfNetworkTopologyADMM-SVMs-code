#include "FindSVs.h"
indxSV FindSVs(Multipliers iLambdas)
{
	indxSV Svs;
	for (int i=0; i < iLambdas.size();i++)
	{
		if (iLambdas[i]> 1e-7)
		{
			Svs.push_back(i);
			//Svs.push_back((int)iLambdas[i]);
		}
	}
	return Svs;
}