#ifndef TheCommandLine_h__
#define TheCommandLine_h__

#include <tchar.h>
#include "DebugTools.h"
#include "CommandLineParser.h"

struct CTheCommandLine : public StandardCommandLineParser
{
public:
	// Create flags and params
	//	MultiValueArg<tstring>          i;
	//	MultiValueArg<FileNameValue>    files;
	ValueArg<unsigned int> qiNoGaussianCubaturePoints;
	ValueArg<unsigned int> qiMaxEvaluationTreeDepth;
	ValueArg<unsigned int> qINumerator;
	ValueArg<unsigned int> qIDenominator;
	ValueArg<double> qdProportionOfRemainingTimeConsumedThisStep;
	ValueArg<double> qdAdjustClusterDiameter;
	ValueArg<double> qdErrorAcceptanceFactor;
	ValueArg<double> qdLocation;
	ValueArg<double> qdStrike;
	ValueArg<double> qdTime;
	ValueArg<double> qdTime1;
	ValueArg<double> qdRate;
	ValueArg<double> qdVol;
	ValueArg<double> qdBarrier;
	FlagArg qbVerbose;
	//ValueArg<std::string> qtstrLogFile;
	//ValueArg<T> zzz11;
	//ValueArg<T> zzz12;

public:
	CTheCommandLine();
	// Set names and descriptions for usage;
};

#endif // TheCommandLine_h__
