cmake_minimum_required(VERSION 3.8)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

project(PointReductionLib_use_dll VERSION 1.0 DESCRIPTION "")


# Tell cmake we want it to automate generating an export stub for the dll
#SET(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)

add_library(PointReductionLib_use_dll STATIC 
	ArrayOfpDoublesToVectorOfPowers.cpp
	CloudWeightLocation.cpp
	CloudWeightLocation2D.cpp
	DoAbstractPruning.cpp
	DoAbstractPruning2D.cpp
	DoAbstractPruningMultithreaded.cpp
	RdToPowers2.cpp
	stdafx.cpp
	TreeBufferHelper.cpp
    ArrayOfpDoublesToVectorOfPowers.h
	BufferConstructor.h
	CloudWeightLocation.h
	CloudWeightLocation2D.h
	EvaluateAllMonomials.h
	OstreamContainerOverloads.h
	RdToPowers.h
	RdToPowers2.h
	recombine.h
	recombine_helper_fn.h
	stdafx.h
	ToPowersData.h
	TreeBufferHelper.h
)

target_link_libraries(PointReductionLib_use_dll 	
	recombine
	Morton
	libs
)

target_include_directories(
PointReductionLib_use_dll 	PRIVATE
	${CMAKE_SOURCE_DIR}/recombine
	${CMAKE_SOURCE_DIR}/Morton
	${CMAKE_SOURCE_DIR}/libs
	${CMAKE_SOURCE_DIR}/CommandLine
	${BOOST_INCLUDE_DIRS}
)