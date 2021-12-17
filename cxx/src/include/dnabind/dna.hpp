/*****************************************************************************/
/***************** Copyright (C) 2020-2021, Richard Spinney. *****************/
/*****************************************************************************/
//                                                                           //
//    This program is free software: you can redistribute it and/or modify   //
//    it under the terms of the GNU General Public License as published by   //
//    the Free Software Foundation, either version 3 of the License, or      //
//    (at your option) any later version.                                    //
//                                                                           //
//    This program is distributed in the hope that it will be useful,        //
//    but WITHOUT ANY WARRANTY; without even the implied warranty of         //
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          //
//    GNU General Public License for more details.                           //
//                                                                           //
//    You should have received a copy of the GNU General Public License      //
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.  //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#ifndef DNA_H
#define DNA_H

#include <sstream>

#include "definitions.hpp"
#include "constructs.hpp"
#include "python_interface.hpp"

class Dna{
private:
	mutable PythonInterface m_py;
	std::string reverseString(const std::string&) const;
	std::string baseComplements(const std::string&) const;
public:
	Dna(const std::string a_str):m_py(a_str){}
	std::string getComplement(const std::string&) const;
	bool queryComplementaryBond(const std::string&,const std::string&,uint32_t,uint32_t,const Parameters&) const;  
	std::vector<std::string> bindingPatternAgnostic(const std::string&,const std::string&,const uint32_t,const uint32_t,const Parameters&) const;
	double getFreeEnergy(const Parameters&,const std::string&,const std::string&,const uint32_t&,const uint32_t&,bool&) const;
};

#endif /*DNA_H*/
