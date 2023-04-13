#include "stdafx.h"
#include ".\MainProgramParameters.h"
//#include ".\thecommandline.h"

CMainProgramParameters::CMainProgramParameters(void)
	: iNoGaussianCubaturePoints(8),
	  iCountBoundaryFunctionEvaluations(0),
	  iMaxEvaluationTreeDepth(30),
	  dAdjustClusterDiameter(1.),
	  iNumerator(2),
	  iDenominator(1),
	  tolerance(0),
	  dProportionOfRemainingTimeConsumedThisStep(.5),
	  dErrorAcceptanceFactor(0) {}

CMainProgramParameters::~CMainProgramParameters(void) {}
/*
int CMainProgramParameters::MapCommandLineToProgramParameters(int argc, const _TCHAR** argv)
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
		dTime1 = cl.qdTime1;
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
*/
CTheCommandLine::CTheCommandLine()
	: qiNoGaussianCubaturePoints(__T("NCPts"), __T("No of Gaussian Cubature Points")),
	  qiMaxEvaluationTreeDepth(__T("qiMaxEvaluationTreeDepth"), __T("qiMaxEvaluationTreeDepth value")),
	  qINumerator(__T("qINumerator"), __T("qINumerator value")),
	  qIDenominator(__T("qIDenominator"), __T("qIDenominator value")),
	  qdProportionOfRemainingTimeConsumedThisStep(__T("qdProportionOfRemainingTimeConsumedThisStep"),
												  __T("qdProportionOfRemainingTimeConsumedThisStep value")),
	  qdAdjustClusterDiameter(__T("qdAdjustClusterDiameter"), __T("qdAdjustClusterDiameter value")),
	  qdErrorAcceptanceFactor(__T("qdErrorAcceptanceFactor"), __T("qdErrorAcceptanceFactor value")),
	  qdLocation(__T("qdLocation"), __T("qdLocation value")),
	  qdStrike(__T("qdStrike"), __T("qdStrike value")),
	  qdTime(__T("qdTime"), __T("qdTime value")),
	  qdTime1(__T("qdTime1"), __T("qdTime1 value")),
	  qdRate(__T("qdRate"), __T("qdRate value")),
	  qdVol(__T("qdVol"), __T("qdVol value")),
	  qdBarrier(__T("qdBarrier"), __T("qdBarrier value")),
	  qbVerbose(__T("qbVerbose"), __T("qbVerbose value"))
	//, qtstrLogFile(__T("qtstrLogFile"), __T("qtstrLogFile value"))
	//, zzz11(__T("zzz11"), __T("zzz11 value"))
	//, zzz12(__T("zzz12"), __T("zzz12 value"))
{
	// Add flags matched by name
	// comment out the default value if a value MUST be given

	AddFlag(qiNoGaussianCubaturePoints);
	qiNoGaussianCubaturePoints.SetDefaultValue(16);

// does this need changing? maybe not
	AddFlag(qiMaxEvaluationTreeDepth);
	qiMaxEvaluationTreeDepth.SetDefaultValue(58);	// for AtExpOpt it's 58 // should be 28

	AddFlag(qINumerator);
	qINumerator.SetDefaultValue(8);	

	AddFlag(qIDenominator);
	qIDenominator.SetDefaultValue(16);	

// Proportion Of Remaining Time Consumed this step
	AddFlag(qdProportionOfRemainingTimeConsumedThisStep);
	qdProportionOfRemainingTimeConsumedThisStep.SetDefaultValue(.4);

// adjust cluster diameter
	AddFlag(qdAdjustClusterDiameter);
	qdAdjustClusterDiameter.SetDefaultValue(0.075);	// Terry 0.075		ATExpCall 0.014 is optimal

	AddFlag(qdErrorAcceptanceFactor);
	qdErrorAcceptanceFactor.SetDefaultValue(20.);

	/// redundant begin --- moved to modeldata and optiondata
	AddFlag(qdLocation);
	qdLocation.SetDefaultValue(0.);

	AddFlag(qdStrike);
	qdStrike.SetDefaultValue(1.);

	AddFlag(qdBarrier);
	qdBarrier.SetDefaultValue(2.);

	AddFlag(qdTime);
	qdTime.SetDefaultValue(1.0);

	AddFlag(qdRate);
	qdRate.SetDefaultValue(0.0);
	/// redundant end//////////////////////////

	AddFlag(qdVol);
	qdVol.SetDefaultValue(1.0);

	AddFlag(qbVerbose);

	/*AddFlag(qtstrLogFile);
	std::string strTemp("log.txt");
	qtstrLogFile.SetDefaultValue(strTemp);*/
	//AddFlag(zzz11);
	//zzz11.SetDefaultValue(17);

	//AddFlag(zzz12);
	//zzz12.SetDefaultValue(17);

	//// Add params matched by position, e.g. foo
	//AddParam(files);
}
