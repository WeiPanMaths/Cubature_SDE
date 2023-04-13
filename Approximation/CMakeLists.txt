cmake_minimum_required(VERSION 3.8)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

project(Approximation VERSION 1.0 DESCRIPTION "")


# Tell cmake we want it to automate generating an export stub for the dll
#SET(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)

add_library(Approximation STATIC 
	LeastSquaresInterpolation.cpp
	LocalApproximation.cpp
	AdaptiveApproximation.h
	LeastSquaresInterpolation.h
	LocalApproximation.h
	SHOW.h
	Utils.h
)

check_exists(${MKL_LIBRARIES})
check_exists(${MKL_INCLUDE_DIRS})

target_include_directories(Approximation PRIVATE
	${MKL_INCLUDE_DIRS} 
	${BOOST_INCLUDE_DIRS}
	${CMAKE_SOURCE_DIR}/Morton)

target_link_directories(Approximation PRIVATE
	${MKL_LIBRARIES}
	${BOOST_LIBRARIES}
)


#target_compile_definitions(recombine
#	PRIVATE
#	NOSIMPLEX
#	REDUCTION_ALGO=svd
#	MKL
#	RECOMBINE_EXPORTS
#	NDEBUG
#)
