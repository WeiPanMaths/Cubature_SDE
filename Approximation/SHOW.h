#ifndef SHOW_h__
#define SHOW_h__
#include <stdlib.h>
#include <iosfwd>
#include <iostream>
#pragma warning( push )
#pragma warning( disable : 4100 )
template <class T>
void SHO_(char* message, const T& value, char* message1, unsigned line, std::ostream& os = std::cout)
{
	os << message << " " << value;	
#ifdef SHO_FILENAME
#pragma warning( push )
#pragma warning( disable : 4996 )
	char drive[_MAX_DRIVE];
	char dir[_MAX_DIR];
	char fname[_MAX_FNAME];
	char ext[_MAX_EXT];

	_splitpath( message1, drive, dir, fname, ext ); // C4996
	// Note: _splitpath is deprecated; consider using _splitpath_s instead
	os << std::endl << "\t" << fname << ext << " Line:"<< line;
#pragma warning( pop )
#endif
	os << std::endl << std::endl;
}
#pragma warning( pop )
// Macro to give trace
#ifndef SHOW
#ifndef SHO_IOS
#define SHO_IOS std::cout
#endif
#define SHOW(arg) SHO_( #arg ,(arg), __FILE__ , __LINE__ , SHO_IOS)
#endif
//#define SHOW(arg)

#endif // SHOW_h__
