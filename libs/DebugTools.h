/* Copyright -  All Rights Reserved - Terry Lyons 2008 */
#pragma once
#include <string>
#include <fstream>
#include <vector>
#include <iostream>

class CDebugTools
{
	//	bool m_bFileOpen;
	//
	//	std::string m_filename;
	//
	//public:
	//	std::ofstream m_out;
	//
	//	inline CDebugTools( void ) : 
	//	m_bFileOpen(false)
	//	{
	//
	//	}
	//
	//	inline bool OpenFile( std::string &filename )
	//	{
	//		m_bFileOpen = true;
	//		m_filename = filename;
	//		m_out.open(filename.c_str());
	//		if(!m_out) 
	//		{
	//			m_bFileOpen = false;
	//			std::cerr << "Cannot open the file: " << filename.c_str() <<" Output aborted.\n";  
	//		}
	//		return m_bFileOpen;
	//	}
	//
	//	inline ~CDebugTools( void )
	//	{
	//		if (m_bFileOpen)
	//		{
	//			m_out.close(); 
	//			m_bFileOpen = false;
	//			std::string cmd("start "); 
	//			system( (cmd + m_filename).c_str() );
	//		}
	//	}
	//
	//
	//template < class U >
	//void ShowArithmeticContainer(const U &weights, const std::string &message )
	//{
	//	//double dIntegralBefore(0);
	//	//m_out << std::endl;
	//	//for (U::const_iterator it = weights.begin();it != weights.end() ; ++ it)
	//	//{
	//	//	m_out << "," << *it ;
	//	//	dIntegralBefore += (*it);
	//	//}
	//	//m_out << std::endl << dIntegralBefore << "," << message.c_str() << std::endl;
	//}

};
