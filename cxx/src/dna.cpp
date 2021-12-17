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

#include "dna.hpp"

std::string Dna::baseComplements(const std::string &a_input) const{
	std::string strOut;
	for (const auto &base : a_input){
		if (base=='A') strOut+='T';
		else if (base=='T') strOut+='A';
		else if (base=='C') strOut+='G';
		else if (base=='G') strOut+='C';
	}
	return strOut;
}

std::string Dna::reverseString(const std::string &a_input) const{
	std::string strOut;
	for (auto i = a_input.crbegin();i!=a_input.crend();++i) strOut+=*i;
	return strOut;
}

std::string Dna::getComplement(const std::string &a_input) const{
	return baseComplements(reverseString(a_input));
}

bool Dna::queryComplementaryBond(const std::string &a_sequence,const std::string &a_sequenceC,uint32_t a_index,uint32_t a_indexC,const Parameters &a_params) const{
	std::string subseq = a_sequence.substr(a_index,a_params.bindLength());
	std::string subseqC = a_sequenceC.substr(a_indexC,a_params.bindLength());
	std::string subseqComp = getComplement(subseqC);
	if (subseq.compare(subseqComp)!=0) return false;
	else return true;
}

//agnostic to whether bases are complementary
std::vector<std::string> Dna::bindingPatternAgnostic(const std::string &a_sequence,const std::string &a_sequenceC,const uint32_t a_index,const uint32_t a_indexC,const Parameters &a_params) const{
	std::string str1;
	for (uint32_t i=0;i<a_index;++i) str1+='.';
	for (uint32_t i=0;i<a_params.bindLength();++i) str1+='(';
	for (size_t i=0;i<a_sequence.size()-a_params.bindLength()-a_index;++i) str1+='.';
	std::string str2;
	for (uint32_t i=0;i<a_indexC;++i) str2+='.';
	for (uint32_t i=0;i<a_params.bindLength();++i) str2+=')';
	for (size_t i=0;i<a_sequenceC.size()-a_params.bindLength()-a_indexC;++i) str2+='.';
	return {std::move(str1),std::move(str2)};
}

double Dna::getFreeEnergy(const Parameters &a_params,const std::string &a_seq,const std::string &a_seqC,const uint32_t &a_index,const uint32_t &a_indexC,bool &a_error) const{
	bool match = queryComplementaryBond(a_seq,a_seqC,a_index,a_indexC,a_params);
	if (match){
		std::vector<std::string> sequences = {a_seq,a_seqC};
		std::vector<std::string> bindings = bindingPatternAgnostic(a_seq,a_seqC,a_index,a_indexC,a_params);
		try{//signature is module, function name, then list of arguments expected by the python function
			std::any feAny = m_py.callFunction(std::string("nupack_functions"),std::string("getFreeEnergy"),
			                                	sequences,bindings,a_params.getMaterialString(),
			                                	a_params.getStackingString(),a_params.temperatureC(),a_params.Naconcentration(),a_params.Mgconcentration());
			double fe = std::any_cast<double> (feAny);
			if (a_params.debug()==2){
				std::cout<<sequences[0]<<"+"<<sequences[1]<<std::endl;
				std::cout<<bindings[0]<<"+"<<bindings[1]<<std::endl;
				std::cout<<"free energy cached as: "<<fe<<std::endl<<std::endl;
			}
			return fe;
		}
		catch(const std::exception& ex) {
        	std::cout<<"Exception occurred: \""<<ex.what()<<"\""<<std::endl;
			a_error = true;
			if constexpr(std::numeric_limits<double>::has_quiet_NaN)
				return std::numeric_limits<double>::quiet_NaN();
			else
				return std::nan("");
		}
		catch(...){
			a_error = true;
			if constexpr(std::numeric_limits<double>::has_quiet_NaN)
				return std::numeric_limits<double>::quiet_NaN();
			else
				return std::nan("");	
		}
	}
	else{
		if constexpr(std::numeric_limits<double>::has_quiet_NaN)
			return std::numeric_limits<double>::quiet_NaN();
		else
			return std::nan("");
	}
}

