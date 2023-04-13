#include "StdAfx.h"
#include ".\scopedwindowsperformancetimer.h"
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#pragma comment(lib, "Kernel32.lib")

CScopedWindowsPerformanceTimer::CScopedWindowsPerformanceTimer(double& dTimeElapsed)
	: secsTimer(dTimeElapsed)
{
	QueryPerformanceFrequency(reinterpret_cast< LARGE_INTEGER*> (&lpFrequency));
	QueryPerformanceCounter(reinterpret_cast< LARGE_INTEGER*> (&lpPerformanceCountStart));
}

CScopedWindowsPerformanceTimer::~CScopedWindowsPerformanceTimer(void)
{
	QueryPerformanceCounter(reinterpret_cast< LARGE_INTEGER*> (&lpPerformanceCountFinish));
	secsTimer = double(lpPerformanceCountFinish - lpPerformanceCountStart) / lpFrequency;
}
