cmake_minimum_required(VERSION 3.8)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

project(Morton VERSION 1.0 DESCRIPTION "")


# Tell cmake we want it to automate generating an export stub for the dll
#SET(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)

add_library(Morton STATIC 
	bitutilities.cpp
	#Learn.cpp
	splice_char.cpp
	stdafx.cpp
	#test.cpp
	#WeiTest.cpp
    bitutilities.h
	block.h
	log_base2.h
	lsplice.h
	morton.h
	smallest_int_type.h
	splice_lib.h
	stdafx.h
	targetver.h
	unblock.h
	wpUtilities.h
)

