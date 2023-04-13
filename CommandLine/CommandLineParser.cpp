#include "stdafx.h"
#include "./CommandLineParser.h"  
// Unicode console output
#ifdef _UNICODE
std::wostream& tcout = std::wcout;
#else
std::ostream& tcout = std::cout;
#endif

bool CommandLineParser::ReadFile( const tstring& sFullPath, tstring& sContents )
{
	// Clean up
	ClearError();
	sContents = _T("");

	FileHandle hFile = CreateFile(sFullPath.c_str(), GENERIC_READ, 0, 0, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, 0);

	if( !hFile )
	{
		_ASSERTE(!"Failed to open parameter file");
		m_sError = _T("Failed to open file ");
		m_sError += sFullPath;

		return false;
	}

	BY_HANDLE_FILE_INFORMATION bhfi;
	::ZeroMemory(&bhfi, sizeof(bhfi));

	if( GetFileInformationByHandle(hFile, &bhfi) == 0 )
	{
		_ASSERTE(!"Failed to get file information");
		m_sError = _T("Failed to get file information for ");
		m_sError += sFullPath;

		return false;
	}

	// A few special cases
	if( bhfi.nFileSizeHigh != 0 )
	{
		_ASSERTE(!"Ridiculously large parameter file");
		m_sError = _T("File ");
		m_sError += sFullPath;
		m_sError += _T(" is too large");

		return false;
	}
	else if( bhfi.nFileSizeLow == 0 )
	{
		// Empty file (but should be accepted anyway)
		return true;
	}

	size_t nDataSize = bhfi.nFileSizeLow;

	_ASSERT( size_t((DWORD)nDataSize)==nDataSize);
	Handle hFileMapping = CreateFileMapping(hFile, 0, PAGE_WRITECOPY, 0, (DWORD) nDataSize, 0);
	if( !hFileMapping )
	{
		_ASSERTE(!"Failed to create file mapping");
		m_sError = _T("Failed to create file mapping on file ");
		m_sError += sFullPath;

		return false;
	}

	LPVOID lpView = MapViewOfFile(hFileMapping, FILE_MAP_COPY, 0, 0, 0);

	if( !lpView )
	{
		_ASSERTE(!"Failed to map view of file");
		m_sError = _T("Failed to map view of file ");
		m_sError += sFullPath;

		return false;
	}

	class MapViewCleanup
	{
	public:
		explicit MapViewCleanup(LPVOID lpv)
			: m_lpView(lpv)
		{}

		~MapViewCleanup()
		{ ::UnmapViewOfFile(m_lpView); }

	private:
		LPVOID m_lpView;
	} mapCleanup(lpView);

	// Try to determine if unicode
	LPVOID  lpData = lpView;
	bool    bUnicode = false;

	if( nDataSize > sizeof(WORD) )
	{
		WORD wMaybeUnicodeMarker = *(reinterpret_cast<LPWORD>(lpView));

		if( wMaybeUnicodeMarker == 0xFEFF )
		{
			_ASSERTE(0 == (nDataSize%sizeof(wchar_t)));
			lpData = &(reinterpret_cast<LPWORD>(lpView)[1]);
			bUnicode = true;
		}
		else if( wMaybeUnicodeMarker == 0xFFFE )
		{
			_ASSERTE(!"Big-endian unicode not supported");
			m_sError = _T("File ");
			m_sError += sFullPath;
			m_sError += _T(" is in big-endian unicode format");

			return false;
		}
	}

	// Convert as necessary
	try
	{
		if( bUnicode ) ConvertStringData(reinterpret_cast<LPCWSTR>(lpData), (nDataSize/sizeof(wchar_t))-1, sContents);
		else ConvertStringData(reinterpret_cast<LPCSTR>(lpData), nDataSize/sizeof(char), sContents);
	}
	catch( std::exception& e )
	{
		_ASSERTE(e.what());
		e;  // Avoid "unreferenced local variable"

		m_sError = _T("Failed to convert character data from file ");
		m_sError += sFullPath;
		m_sError += _T(".");

		return false;
	}

	// Replace CRs & LFs
	std::replace_if(sContents.begin(),
		sContents.end(),
		std::bind2nd(std::equal_to<TCHAR>(), _T('\r')),
		_T(' '));
	std::replace_if(sContents.begin(),
		sContents.end(),
		std::bind2nd(std::equal_to<TCHAR>(), _T('\n')),
		_T(' '));

	return true;
}

bool CommandLineParser::ParseFromFile( const tstring& sFileName )
{
	_ASSERT(sFileName.size());

	// Check if file exists
	DWORD   dwAttrib = GetFileAttributes(sFileName.c_str());
	bool    bFileExists = (dwAttrib != -1) && !(dwAttrib & FILE_ATTRIBUTE_DIRECTORY);

	if( !bFileExists )
	{
		m_sError = __T("");
		m_sError += __T("File not found: ");
		m_sError += sFileName;
		return false;
	}

	// Point current directory at the input file
	TCHAR   szFileDir[MAX_PATH+1];
	LPTSTR  pszFilePart = 0;
	GetFullPathName(sFileName.c_str(), MAX_PATH, szFileDir, &pszFilePart);
	tstring sFullPath = szFileDir;
	_ASSERTE(pszFilePart);
	*pszFilePart = 0;

	// Set (and auto-reset) current working directory
	// We do this so that any file names read out of the file
	// can be relative to the file, not where we're running
	// the app from. We're using the system managed working directory
	// as a "context" for the individual values, i.e. FileNameValue,
	// to be able to compute a complete file name. This is potentially
	// dangerous as the cwd is set per process, not per thread, but
	// since command lines are typically processed before threads are
	// fired off, we should be safe. It saves us from having to pass
	// a virtual cwd to all values as they're parsed.
	class CurrentDir
	{
	public:
		CurrentDir(LPCTSTR pszNewDir)
		{
			GetCurrentDirectory(MAX_PATH, m_szOldDir);
			_ASSERT(pszNewDir);
			BOOL    b = SetCurrentDirectory(pszNewDir);
			_ASSERTE(b);
		}

		~CurrentDir()
		{
			SetCurrentDirectory(m_szOldDir);
		}

	private:
		TCHAR   m_szOldDir[MAX_PATH+1];
	};

	CurrentDir  cd(szFileDir);

	// Read in the contents of the file
	tstring tContents;
	if (!ReadFile(sFullPath, tContents)) return false;

	// Parse the contents of the file like a string
	return Parse(tContents, false, false);
}

void CommandLineParser::ClearError()
{
	m_sError.erase();
}

tstring CommandLineParser::Logo() const
{
	CVersionInfo        vi;
	if( FAILED(vi.LoadVersion()) ) return __T("");

	LPCTSTR pszDesc = vi.QueryValueString(__T("FileDescription"));
	LPCTSTR pszVersion = vi.QueryValueString(__T("FileVersion"));
	LPCTSTR pszCopyright = vi.QueryValueString(__T("LegalCopyright"));

	tstring         sDesc = (pszDesc ? pszDesc : __T(""));
	tstring         sVersion = (pszVersion ? pszVersion : __T(""));
	tstring         sCopyright = (pszCopyright ? pszCopyright : __T(""));
	tstringstream   ss;

	sCopyright = SearchAndReplace(sCopyright, __T("\x0a9"), __T("(c)"));
	sCopyright = SearchAndReplace(sCopyright, __T("\x0ae"), __T("(r)"));

	if( sCopyright[sCopyright.size() - 1] != __T('.') ) sCopyright += __T('.');

	// e.g.
	// Microsoft MIDL Compiler Version 5.01.0164
	// Copyright (c) 1991-1997, Microsoft, Inc. All rights reserved.
	ss  << sDesc << __T(" version ") << sVersion << std::endl
		<< sCopyright << __T(" All rights reserved.") << std::endl;

	return ss.str();
}

bool CommandLineParser::FileExists( const tstring& sFileName )
{
	DWORD dwAttrib = GetFileAttributes(sFileName.c_str());
	return (dwAttrib != -1) && !(dwAttrib & FILE_ATTRIBUTE_DIRECTORY);
}

tstring CommandLineParser::ModuleName()
{
	TCHAR   szPath[MAX_PATH+1];
	GetModuleFileName(0, szPath, MAX_PATH);
	TCHAR   szFile[_MAX_FNAME+1];
	_tsplitpath(szPath, 0, 0, szFile, 0);
	return szFile;
}

bool CommandLineParser::IsWhiteSpace( TCHAR ch )
{
	return ch == __T(' ') || ch == __T('\t') || ch == __T('\n') || ch == __T('\r');
}

tstring CommandLineParser::SearchAndReplace( tstring s, const tstring& sFind, const tstring& sReplace, bool bGlobal /*= true*/ )
{
	_ASSERTE(sFind.size());
	const size_t    cchFind = sFind.size();

	size_t pos = 0;
	while( (pos = s.find(sFind, pos)) != tstring::npos )
	{
		s.replace(pos, cchFind, sReplace);
		if( !bGlobal ) break;
		pos += cchFind;
	}

	return s;
}

bool CommandLineParser::Parse( const std::vector<tstring>& rgsArgs, bool bValidate /*= true*/ )
{
	ClearError();

	for(std::vector<tstring>::const_iterator psArg = rgsArgs.begin(); psArg != rgsArgs.end(); ++psArg)
	{
#ifdef _DEBUG
		OutputDebugString(__T("Processing arg: "));
		OutputDebugString(psArg->c_str());
		OutputDebugString(__T("\n"));
#endif
		Argument*   pArg = 0;
		bool        bIsFlag = false;

		// It's a flag
		if( (psArg->size() > 1) && (psArg->at(0) == __T('/') || psArg->at(0) == __T('-')) )
		{
			// Find the argument by name
			pArg = m_rgArgInfos.FindFlag(psArg->substr(1));
			bIsFlag = true;
		}
		// It's a file name to process parameters from
		else if( psArg->size() > 1 && psArg->at(0) == __T('@') && m_bAllowArgFile )
		{
			if( !ParseFromFile(psArg->substr(1)) ) return false;
			continue;
		}
		// It's a parameter
		else
		{
			// Find the argument by offset
			pArg = m_rgArgInfos.GetNextParam();
		}

		if( !pArg )
		{
			m_sError = __T("Unrecognized argument '") + *psArg + __T("'.");
			return false;
		}

		_ASSERTE(pArg);

		// Argument with a value, e.g. /foo bar
		if( pArg->ExpectsValues() )
		{
			if( bIsFlag && (++psArg == rgsArgs.end()) )
			{
				m_sError = __T("Argument '") + pArg->Name() + __T("' expects a parameter.");
				return false;
			}

			if( !pArg->WantsAValue() )
			{
				m_sError = __T("Argument '") + pArg->Name() + __T("' cannot accept another parameter: '") + *psArg + __T("'.");
				return false;
			}

			if( !pArg->ConsumeValue(*psArg) )
			{
				m_sError = __T("Argument '") + pArg->Name() + __T("' cannot accept parameter: '") + *psArg + __T("'.");
				return false;
			}
		}
		// Argument w/o a value, e.g. /foo
		else
		{
			if( pArg->Found() )
			{
				m_sError = __T("Argument '") + pArg->Name() + __T("' already present.");
				return false;
			}
		}

		pArg->SetFound(true);
	}

	// Check for missing required arguments
	for( ArgInfoContainer::iterator it = m_rgArgInfos.begin(); it != m_rgArgInfos.end(); ++it )
	{
		_ASSERTE(it->pArg);

		if( !it->pArg->IsValid() )
		{
			m_sError = __T("Expected argument '") + it->pArg->Name() + __T("'.");
			return false;
		}
	}

	return true;
}

bool CommandLineParser::Parse( int argc, const TCHAR** argv, bool bValidate /*= true*/, bool bLeadingExe /*= true*/ )
{
	std::vector<tstring> Params;
	for(int i = (bLeadingExe ? 1 : 0); i < argc; ++i) Params.push_back(argv[i]);
	return Parse(Params, bValidate);
}

bool CommandLineParser::Parse( const tstring & CommandLine, bool bValidate /*= true*/, bool bLeadingExe /*= false*/ )
{
	std::vector<tstring> Params;
		tstring CurParam;

		bool bInQuotes = false;
		for(tstring::const_iterator i = CommandLine.begin(); i != CommandLine.end(); ++i)
		{
			if (*i == __T('"'))
			{
				bInQuotes = bInQuotes ? false : true;
				continue;
			}

			if (IsWhiteSpace(*i) && !bInQuotes)
			{
				if (!CurParam.empty())
				{
					if( bLeadingExe ) bLeadingExe = false;
					else Params.push_back(CurParam);
					CurParam = __T("");
				}
				continue;
			}

			CurParam += *i;
		}

		// Check for bLeadingExe to avoid pushing module file
		// name if no parameters were supplied.
		if( !CurParam.empty() && !bLeadingExe )
		{
			Params.push_back(CurParam);
		}

		return Parse(Params, bValidate);
}
Argument* CommandLineParser::ArgInfoContainer::GetNextParam() const
{
	// Give each argument all the values it wants
	// before moving onto the next one
	for( const_iterator it = begin(); it != end(); ++it )
	{
		if( it->type == ARG_PARAM )
		{
			_ASSERTE(it->pArg);
			Argument*   pArg = it->pArg;
			if( pArg->ExpectsValues() && pArg->WantsAValue() )
			{
				return pArg;
			}
		}
	}

	return 0;
}

StandardCommandLineParser::StandardCommandLineParser( bool bAllowArgFile /*= true*/ ) :   CommandLineParser(bAllowArgFile),
help(__T("?"), __T("Show usage.")),
version(__T("v"), __T("Show version."))
{
	AddFlag(help);
	help.AddAlternate(__T("h"));
	help.AddAlternate(__T("help"));

	AddFlag(version);
	version.AddAlternate(__T("version"));
}

bool StandardCommandLineParser::IsConsole()
{
	TCHAR       szPath[MAX_PATH+1];
	GetModuleFileName(0, szPath, MAX_PATH);

	
	SHFILEINFO  sfi = { 0 };
	_ASSERT(size_t((int)&sfi)==size_t(&sfi));
	//tjl DWORD 
	DWORD_PTR dw = SHGetFileInfo(szPath, 0, &sfi, sizeof(sfi), SHGFI_EXETYPE);

	return (LOWORD(dw) == IMAGE_NT_SIGNATURE) && (HIWORD(dw) == 0);
}

void StandardCommandLineParser::Show( const tstring& s ) const
{
	if( IsConsole() )
	{
		// Always send usage to stdout so it's easy to capture the output
		tcout << s.c_str() << std::endl;
	}
	else
	{
		MessageBox(0, s.c_str(), ModuleName().c_str(), MB_SETFOREGROUND);
	}
}

bool StandardCommandLineParser::Continue()
{
	if( version )
	{
		Show(Logo());
		return false;
	}

	if( m_sError.size() || help )
	{
		if( help ) ClearError();
		Show(Usage());
		return false;
	}

	return true;
}

bool StandardCommandLineParser::ParseAndContinue( int argc, const TCHAR* argv[] )
{
	Parse(argc, argv);
	return Continue();
}

bool StandardCommandLineParser::ParseAndContinue( int argc, TCHAR* argv[] )
{
	// To avoid C2664
	Parse(argc, (const TCHAR**)argv);
	return Continue();
}

bool StandardCommandLineParser::ParseAndContinue( LPCTSTR psz )
{
	Parse(psz);
	return Continue();
}

bool FileNameValue::Exists() const
{
	DWORD dwAttrib = GetFileAttributes(c_str());
	return (dwAttrib != -1) && !(dwAttrib & FILE_ATTRIBUTE_DIRECTORY);
}