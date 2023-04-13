cmake_minimum_required(VERSION 3.8)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(recombine VERSION 2.0.2 DESCRIPTION "Performs a dynamic Caratheodory process and takes a weighted collection of vectors and identifies by pointers, a subset of minimal cardinality among the vectors and new weights so both empirical measures have the same mean. Software written by Terry Lyons, based on algorithms developed jointly with Christian Litter and then with Maria Tchernychova 2008-2020. Any compiled library built by this cmake script may be used freely without any guarantee or liability. Otherwise, all rights reserved. For further informaton contact Terry Lyons, software@lyonstech.net")

add_library(recombine SHARED
    Cmove1.cpp
    recombine.cpp
    recombine1.cpp
    stdafx.cpp
    TreeBufferHelper.cpp
    LinearAlgebraReductionTool.cpp
	aligned_allocator.h
	BufferConstructor.h
	Cmove1.h
	LinearAlgebraReductionTool.h
	recombine.h
	recombine1.h 
	stdafx.h
	TreeBufferHelper.h
)

target_link_libraries(recombine
	${MKL_LIBRARIES}/mkl_intel_lp64.lib
	${MKL_LIBRARIES}/mkl_intel_thread.lib
	${MKL_LIBRARIES}/mkl_core.lib
	${INTEL_OMP_LIB}
)

target_compile_definitions(recombine PRIVATE
	MKL64
	RECOMBINE_EXPORTS
)
