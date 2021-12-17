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

// global comment:
// when using the CRTP idiom the identity of the base class is dependent, 
// i.e. it depends on <T>, and so is not resolved on the first compiler pass, 
// so we have to use this-> or scope resolve for the compiler to know what we are 
// talking about when accessing parent class members and methods.

#include "derivedmodels.hpp"

template<class T>
void CombinatoricsModel<T>::cacheComplement(const Parameters &a_params){
	m_complementCache.clear();
	m_complementCache.reserve(this->m_sequences.size());
	for (size_t seqIndex=0;seqIndex<this->m_sequences.size();++seqIndex){
		const uint32_t numBindsSeq = this->m_sequences[seqIndex].size()-a_params.bindLength()+1;
		const uint32_t numBindsSeqC = this->m_sequencesC[seqIndex].size()-a_params.bindLength()+1;
		std::vector<std::vector<bool> > outer;
		outer.reserve(numBindsSeq);
		for (uint32_t index=0; index<numBindsSeq;++index){
			std::vector<bool> inner;
			inner.reserve(numBindsSeqC);
			for (uint32_t indexC=0; indexC<numBindsSeqC;++indexC){
				inner.push_back(this->m_dnaMethods.queryComplementaryBond(this->m_sequences[seqIndex],this->m_sequencesC[seqIndex],index,indexC,a_params));
			}
			outer.push_back(std::move(inner));
		}
		m_complementCache.push_back(std::move(outer));
	}
}

template<class T>
std::vector<std::vector<double> > StabilityModel<T>::freeEnergyLoop(uint32_t a_seqIndex,bool &a_error,const Parameters &a_params){
	const uint32_t numBindsSeq = this->m_sequences[a_seqIndex].size()-a_params.bindLength()+1;
	const uint32_t numBindsSeqC = this->m_sequencesC[a_seqIndex].size()-a_params.bindLength()+1;
	std::vector<std::vector<double> > FE;
	FE.reserve(numBindsSeq);
	for (uint32_t index=0;index<numBindsSeq;++index){
		std::vector<double> inner;
		inner.reserve(numBindsSeqC);
		for (uint32_t indexC=0;indexC<numBindsSeqC;++indexC){
			double feVal = this->m_dnaMethods.getFreeEnergy(a_params,this->m_sequences[a_seqIndex],this->m_sequencesC[a_seqIndex],index,indexC,a_error);
			inner.push_back(feVal);
			if (a_error) return FE;	
		}
		FE.push_back(std::move(inner));
	}
	return FE;
}

template<class T>
void StabilityModel<T>::cacheFreeEnergies(const Parameters &a_params){
	std::cout<<"starting free energy cache..."<<std::endl;
	this->m_freeEnergies.clear();
	bool error = false;
	for (size_t seqIndex=0;seqIndex<this->m_sequences.size();++seqIndex) 
		this->m_freeEnergies.push_back(freeEnergyLoop(seqIndex,error,a_params));
	if (error){ 
		std::cerr<<"Fatal error encountered whilst caching free energies - further computation stopped."<<std::endl;
		this->m_fatalError = true;
	}
	else{
		std::cout<<"finished free energy cache"<<std::endl<<std::endl;
	}
}

template<class T>
inline double StabilityModel<T>::logit(const double a_fe,const Solution& a_data) const{
	if (std::isnan(a_fe)) 
		return 0.0;
	else 
		return 1.0/(1.0 + exp(a_data.params().gamma() + a_data.params().beta() * a_fe));
}

template<class T>
void StabilityModel<T>::printFreeEnergies() const{
	for (size_t seqIndex=0;seqIndex<this->m_sequences.size();++seqIndex){
		const uint32_t numBindsSeq = this->m_sequences[seqIndex].size()-this->m_paramsRef.bindLength()+1;
		const uint32_t numBindsSeqC = this->m_sequencesC[seqIndex].size()-this->m_paramsRef.bindLength()+1;
		for (uint32_t index=0;index<numBindsSeq;++index){
			for (uint32_t indexC=0;indexC<numBindsSeqC;++indexC){
				std::vector<std::string> bindings = this->m_dnaMethods.bindingPatternAgnostic(this->m_sequences[seqIndex],this->m_sequencesC[seqIndex],index,indexC,this->m_paramsRef);
				std::cout<<this->m_sequences[seqIndex]<<"+"<<this->m_sequencesC[seqIndex]<<std::endl;
				std::cout<<bindings[0]<<"+"<<bindings[1]<<std::endl;
				std::cout<<"free energy: "<<this->m_freeEnergies[seqIndex][index][indexC]<<std::endl<<std::endl;
			}
		}
	}
}

template class CombinatoricsModel<CombinatoricsModelLength>;
template class CombinatoricsModel<CombinatoricsModelNoLength>;
template class StabilityModel<StabilityModelLength>;
template class StabilityModel<StabilityModelNoLength>;
