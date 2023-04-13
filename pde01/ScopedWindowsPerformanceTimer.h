#pragma once

class CScopedWindowsPerformanceTimer
{
	__int64 lpFrequency;
	__int64 lpPerformanceCountStart;
	__int64 lpPerformanceCountFinish;
	double& secsTimer;
public:
	CScopedWindowsPerformanceTimer(double& dTimeElapsed);
	~CScopedWindowsPerformanceTimer(void);
};
