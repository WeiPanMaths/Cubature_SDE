////////////////////////////////////////////////////////
// ALL rights reserved to the author of this software
// Terry Lyons 
// 2011 - 2012
///////////////////////////////////////////////////////
// The following ifdef block is the standard way of creating macros which make exporting 
// from a DLL simpler. All files within this DLL are compiled with the SPLICE_LIB_EXPORTS
// symbol defined on the command line. This symbol should not be defined on any project
// that uses this DLL. This way any other project whose source files include this file see 
// SPLICE_LIB_API functions as being imported from a DLL, whereas this DLL sees symbols
// defined with this macro as being exported.
#ifdef SPLICE_LIB_EXPORTS
#define SPLICE_LIB_API __declspec(dllexport)
#else
#define SPLICE_LIB_API __declspec(dllimport)
#endif

// This class is exported from the splice_lib.dll
class SPLICE_LIB_API Csplice_lib {
public:
	Csplice_lib(void);
	// TODO: add your methods here.
};

extern SPLICE_LIB_API int nsplice_lib;

SPLICE_LIB_API int fnsplice_lib(void);
