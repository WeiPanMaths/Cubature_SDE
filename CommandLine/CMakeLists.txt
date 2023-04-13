cmake_minimum_required(VERSION 3.8)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

project(CommandLine VERSION 1.0 DESCRIPTION "")


# Tell cmake we want it to automate generating an export stub for the dll
#SET(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)

add_library(CommandLine STATIC 
	CommandLineParser.cpp
	stdafx.cpp
	versioninfo.cpp
	versioninfo.h
	stdafx.h
	CommandLineParser.h
)

