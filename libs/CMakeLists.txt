cmake_minimum_required(VERSION 3.8)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

project(libs VERSION 1.0 DESCRIPTION "")


# Tell cmake we want it to automate generating an export stub for the dll
#SET(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)

add_library(libs STATIC 
	CubaturePDESolver.cpp
	CubaturePDESolver2D.cpp
	DebugTools.cpp
	DoublePrecisionGaussianQuadrature.cpp
	ErrorInCubature.cpp
	ExactPDESolution.cpp
	MainProgramParameters.cpp
	StdAfx.cpp
	BasketSpreadQuantLib.h
	CubaturePDESolver.h
	CubaturePDESolver2D.h
	DebugTools.h
	DoublePrecisionGaussianQuadrature.h
	ExactPDESolution.h
	GaussQuadratureSet.h
	MainProgramParameters.h
	SharedFunctionsAndTypes.h
	stdafx.h
	TemplateCubaturePDESolver.h
	TheCommandLine.h
)

target_include_directories(libs PRIVATE
	${CMAKE_SOURCE_DIR}/PointReductionLib_use_dll
	${CMAKE_SOURCE_DIR}/CommandLine
	${CMAKE_SOURCE_DIR}/Morton
	${CMAKE_SOURCE_DIR}/Approximation
	${QUANTLIB_INCLUDE_DIRS}
	${BOOST_INCLUDE_DIRS}
)


