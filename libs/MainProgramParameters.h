#pragma once
#include "DebugTools.h"
#include <tchar.h>
#include ".\thecommandline.h"

class CMainProgramParameters
{	
public:
	// The Number of Points in the One Step Diffusion / Normal Quadrature Approximation 
	unsigned int iNoGaussianCubaturePoints;

	// Cubature sets for Gaussian act correctly on polynomials of higher degree
	// this fraction is used to determine the effective multiplier 1< M <=2
	//struct  
	//{
	int iNumerator;
	int iDenominator;
	//};

	CMainProgramParameters(void);
	~CMainProgramParameters(void);
	

	// conceptually these three variables belong to contract parameters
	// rather than pde solver
	double dLocation;
	double dStrike;
	double dBarrier;
	double dTime;	// time 1 < time
	double dRate;
	
	
	double dVol;
	bool bVerbose;
	std::string qtstrLogFile;

	// Cubature sets for Gaussian act correctly on polynomials of higher degree
	// this fraction is used to determine the effective multiplier 1< M <=2
	int INumerator() const
	{
		return iNumerator;
	}
	void INumerator(int val)
	{
		iNumerator = val;
	}
	int IDenominator() const
	{
		return iDenominator;
	}
	void IDenominator(int val)
	{
		iDenominator = val;
	}
	
	// Get the Number of Points in the One Step Diffusion / Normal Quadrature Approximation 
	unsigned int NoGaussianCubaturePoints() const
	{
		return iNoGaussianCubaturePoints;
	}

	// Set the Number of Points in the One Step Diffusion / Normal Quadrature Approximation 
	void NoGaussianCubaturePoints(unsigned int val)
	{
		iNoGaussianCubaturePoints = val;
	}
	
	// Should be incremented each time the integration process evaluates a function
	unsigned int iCountBoundaryFunctionEvaluations;	

	// The maximum depth the solution tree is willing to go before going to the boundary
	unsigned int iMaxEvaluationTreeDepth;
	
	// When forming Clusters the basic diameter is given by the diameter of the cubature used to get to that level - ths modifies it slightly
	double dAdjustClusterDiameter;
	// A measure of the errors used to decide what happens
	double tolerance;
	// The next jump in time as a proportion of the time remaining
	double dProportionOfRemainingTimeConsumedThisStep;
	// the multiple of the minimum error in cubature that will be accepted over the two step test to jump to boundary
	double dErrorAcceptanceFactor;
	
	int MapCommandLineToProgramParameters(int argc, const _TCHAR** argv)
	{
		{
			// Parse the command line
			CTheCommandLine cl;
			cl.Parse(argc, argv);

			// Will show usage or version if there's an error, /?|h|help or /v|version
			if (!cl.Continue())
				return 56;

			//update program parameters
			iNoGaussianCubaturePoints = cl.qiNoGaussianCubaturePoints;
			dAdjustClusterDiameter = cl.qdAdjustClusterDiameter;
			dErrorAcceptanceFactor = cl.qdErrorAcceptanceFactor;
			dProportionOfRemainingTimeConsumedThisStep = cl.qdProportionOfRemainingTimeConsumedThisStep;
			iDenominator = cl.qIDenominator;
			iNumerator = cl.qINumerator;
			iMaxEvaluationTreeDepth = cl.qiMaxEvaluationTreeDepth;
			bVerbose = cl.qbVerbose;
			dLocation = cl.qdLocation;
			dStrike = cl.qdStrike;
			dRate = cl.qdRate;
			dTime = cl.qdTime;
			dVol = cl.qdVol;
			dBarrier = cl.qdBarrier;
			std::string strTemp("log.txt");
			qtstrLogFile = strTemp;//cl.qtstrLogFile;
		
			// extra initialization required
			//if(bVerbose)
			//	dbgLogFile.OpenFile(qtstrLogFile);
		}
		return 0;
	}
	// A log file that can be used for messages
	// CDebugTools dbgLogFile;
};
