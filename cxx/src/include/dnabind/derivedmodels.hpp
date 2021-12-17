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


/********************** CLASS STRUCTURE ***************************************

BindModel         	(abstract [not strictly]) base class containing nupack and structure getting methods + loop functions etc.
	^
	|
	|
	|<-- CombinatoricsModel         (abstract [not strictly]) class that caches and uses whether binding sites are matched but not free energies
	|		|
	|		|<-- CombinatoricsModelNoLength       scores 1 vs 0 for matches - no optimised variables 1 - optimised variable (kappa)
	|		|
	|		|<-- CombinatoricsModelLength     scores 1 vs 0 for matches with length dependence - 2 optimised variables (kappa,alpha)
	|
	|<-- StabilityModel 	(abstract [not strictly]) class containing free energy caching + getting methods
	|		|
	|		|<-- StabilityModelNoLength  class that implements the logit stability model - 2 optimised variables (kappa,gamma)
	|		|
	|		|<-- StabilityModelLength  class that implements the logit stability model with length dependence  - 3 optimised variables (kappa,alpha,gamma)


All subclassses must define a score_function and length_function which sits inside the innermost main loop to determine the model

Efficient calling of the score_function is achieved with static polymorphism through the CRTP idiom

All classes use a shared parameters class - if these need to take specific values set them in the subclass's constructor.

*****************************************************************************/

#ifndef DERIVEDMODELS_H
#define DERIVEDMODELS_H

#include <math.h>
#include <iostream>
#include <string>
#include <iomanip>
#include <vector>
#include <chrono>
#include <thread>
#include <algorithm>
#include <iterator>

#include "threads/ThreadPool.h"
#include "constructs.hpp"
#include "dna.hpp"
#include "bindmodel.hpp"
#include "definitions.hpp"
#include "python_interface.hpp"

//generic templated classes for CRTP idiom

/******************* GENERIC CLASSES COMBINATORIC / STABILITY **********************/

// generic class using bind site matching score function
template<class T>
class CombinatoricsModel:public BindModel<T>{
protected:
	std::vector<std::vector<std::vector<bool > > > m_complementCache;
	void cacheComplement(const Parameters &a_params);
	void debugOut(const uint32_t &a_seqIndex,const uint32_t &a_index,const uint32_t &a_indexC,const Solution &a_modelData) const{
		if ((this->m_paramsRef.debug()==2)&&(this->m_final)){
			std::vector<std::string> bindings = this->m_dnaMethods.bindingPatternAgnostic(this->m_sequences[a_seqIndex],this->m_sequencesC[a_seqIndex],a_index,a_indexC,a_modelData.params());
			std::cout<<this->m_sequences[a_seqIndex]<<" "<<this->m_sequencesC[a_seqIndex]<<std::endl<<bindings[0]<<"+"<<bindings[1]<<std::endl;
			std::cout<<"score: "<<this->m_complementCache[a_seqIndex][a_index][a_indexC]<<std::endl<<std::endl;
		}
	}
	void init(){
		std::this_thread::sleep_for(std::chrono::milliseconds(1000));
		this->m_paramsRef.setStability(false);	
		cacheComplement(this->m_paramsRef);
	}
	//CTOR protected to prevent instantiation of base class
	CombinatoricsModel(const Parameters a_params,const std::vector<std::tuple<std::string,double,double> > a_input,const std::string a_pythonPath,const std::string a_name)
		:BindModel<T>(a_params,a_input,a_pythonPath,a_name){
		if (this->m_fatalError) return;
		init();
	}		
};

// generic class using NUPACK+logit score function
template<class T>
class StabilityModel:public BindModel<T>{
protected:
	std::vector<std::vector<std::vector<double> > > m_freeEnergies;
	std::vector<std::vector<double> > freeEnergyLoop(uint32_t,bool&,const Parameters&);
	void cacheFreeEnergies(const Parameters&);
	double logit(const double,const Solution&) const;
	void debugOut(const uint32_t &a_seqIndex,const uint32_t &a_index,const uint32_t &a_indexC,const Solution &a_modelData) const{
		if ((this->m_paramsRef.debug()==2)&&(this->m_final)){
			std::vector<std::string> bindings = this->m_dnaMethods.bindingPatternAgnostic(this->m_sequences[a_seqIndex],this->m_sequencesC[a_seqIndex],a_index,a_indexC,a_modelData.params());
			std::cout<<this->m_sequences[a_seqIndex]<<" "<<this->m_sequencesC[a_seqIndex]<<std::endl<<bindings[0]<<"+"<<bindings[1]<<std::endl;
			std::cout<<"free energy: "<<this->m_freeEnergies[a_seqIndex][a_index][a_indexC]<<" score: "<<logit(this->m_freeEnergies[a_seqIndex][a_index][a_indexC],a_modelData)<<std::endl<<std::endl;
		}
	}
	void init(){
		std::this_thread::sleep_for(std::chrono::milliseconds(1000));
		this->m_paramsRef.setStability(true);
		cacheFreeEnergies(this->m_paramsRef);
	}
	//CTOR protected to prevent instantiation of base class
	StabilityModel(const Parameters a_params,const std::vector<std::tuple<std::string,double,double> > a_input,const std::string a_pythonPath,const std::string a_name)
		:BindModel<T>(a_params,a_input,a_pythonPath,a_name){
		if (this->m_fatalError) return;
		init();
	}
public:
	void printFreeEnergies() const;
};

/******************* MODEL CLASSES **********************/

class CombinatoricsModelNoLength:public CombinatoricsModel<CombinatoricsModelNoLength>{
private:
	friend class BindModel;
	double scoreFunction(const uint32_t &a_seqIndex,const uint32_t &a_index,const uint32_t &a_indexC,const Solution& a_modelData) const{
		debugOut(a_seqIndex,a_index,a_indexC,a_modelData);
		if (m_complementCache[a_seqIndex][a_index][a_indexC]) return 1.0;
		return 0.0;
	}
	double lengthFunction(const int &a_len,const Solution&) const{
		return std::pow(static_cast<double>(a_len),-2.0);
	}
	void init(){
		m_paramsRef.setLength(false);	
		m_numOptimParams = 1;
		m_optLower = {MIN_KAPPA};
		m_optUpper = {MAX_KAPPA};
		m_optStart = {m_paramsRef.kappa()};
		m_paramNames = {"kappa"};
	}
	void armaToParams(Parameters &a_params,const arma::vec &a_optVars){
		a_params.setKappa(a_optVars(0));
	}
public:
	CombinatoricsModelNoLength(const Parameters a_params,const std::vector<std::tuple<std::string,double,double> > a_input,const std::string a_pythonPath)
		:CombinatoricsModel(a_params,a_input,a_pythonPath,"combinatorics_model_no_length"){
		init();
		printName();
		printParameters();
	}
};

class CombinatoricsModelLength:public CombinatoricsModel<CombinatoricsModelLength>{
private:	
	friend class BindModel;
	double scoreFunction(const uint32_t &a_seqIndex,const uint32_t &a_index,const uint32_t &a_indexC,const Solution& a_modelData) const{
		debugOut(a_seqIndex,a_index,a_indexC,a_modelData);
		if (m_complementCache[a_seqIndex][a_index][a_indexC]) return 1.0;
		return 0.0;
	}
	double lengthFunction(const uint32_t &a_len,const Solution& a_modelData) const{
		return std::pow(static_cast<double>(a_len),a_modelData.params().alpha());
	}
	void init(){	
		m_paramsRef.setLength(true);
		m_numOptimParams = 2;
		m_optLower = {MIN_ALPHA,MIN_KAPPA};
		m_optUpper = {MAX_ALPHA,MAX_KAPPA};
		m_optStart = {m_paramsRef.alpha(),m_paramsRef.kappa()};
		m_paramNames = {"alpha","kappa"};
	}
	void armaToParams(Parameters &a_params,const arma::vec &a_optVars){
		a_params.setAlpha(a_optVars(0));
		a_params.setKappa(a_optVars(1));
	}
public:
	CombinatoricsModelLength(const Parameters a_params,const std::vector<std::tuple<std::string,double,double> > a_input,const std::string a_pythonPath)
		:CombinatoricsModel(a_params,a_input,a_pythonPath,"combinatorics_model_length"){
		init();
		printName();
		printParameters();
	}
};

class StabilityModelNoLength:public StabilityModel<StabilityModelNoLength>{
private:	
	friend class BindModel;
	double scoreFunction(const uint32_t &a_seqIndex,const uint32_t &a_index,const uint32_t &a_indexC,const Solution& a_modelData) const{
		debugOut(a_seqIndex,a_index,a_indexC,a_modelData);
		return logit(m_freeEnergies[a_seqIndex][a_index][a_indexC],a_modelData);
	}
	double lengthFunction(const uint32_t &a_len,const Solution&) const{
		return std::pow(static_cast<double>(a_len),-2.0);
	}
	void init(){	
		m_paramsRef.setLength(false);	
		m_numOptimParams = 2;
		m_optLower = {MIN_GAMMA,MIN_KAPPA};
		m_optUpper = {MAX_GAMMA,MAX_KAPPA};
		m_optStart = {m_paramsRef.gamma(),m_paramsRef.kappa()};
		m_paramNames = {"gamma","kappa"};
	}
	void armaToParams(Parameters &a_params,const arma::vec &a_optVars){			
		a_params.setGamma(a_optVars(0));
		a_params.setKappa(a_optVars(1));
	}
public:
	StabilityModelNoLength(const Parameters a_params,const std::vector<std::tuple<std::string,double,double> > a_input,const std::string a_pythonPath)
		:StabilityModel(a_params,a_input,a_pythonPath,"stability_model_no_length"){
		init();
		printName();
		printParameters();
	}
};

class StabilityModelLength:public StabilityModel<StabilityModelLength>{
private:	
	friend class BindModel;
	double scoreFunction(const uint32_t &a_seqIndex,const uint32_t &a_index,const uint32_t &a_indexC,const Solution& a_modelData) const{				
		debugOut(a_seqIndex,a_index,a_indexC,a_modelData);
		return logit(m_freeEnergies[a_seqIndex][a_index][a_indexC],a_modelData);
	}
	double lengthFunction(const uint32_t &a_len,const Solution& a_modelData) const{
		return std::pow(static_cast<double>(a_len),a_modelData.params().alpha());
	}
	void init(){	
		m_paramsRef.setLength(true);
		m_numOptimParams = 3;
		m_optLower = {MIN_ALPHA,MIN_GAMMA,MIN_KAPPA};
		m_optUpper = {MAX_ALPHA,MAX_GAMMA,MAX_KAPPA};
		m_optStart = {m_paramsRef.alpha(),m_paramsRef.gamma(),m_paramsRef.kappa()};
		m_paramNames = {"alpha","gamma","kappa"};
	}
	void armaToParams(Parameters &a_params,const arma::vec &a_optVars){
		a_params.setAlpha(a_optVars(0));
		a_params.setGamma(a_optVars(1));
		a_params.setKappa(a_optVars(2));
	}
public:
	StabilityModelLength(const Parameters a_params,const std::vector<std::tuple<std::string,double,double> > a_input,const std::string a_pythonPath)
		:StabilityModel(a_params,a_input,a_pythonPath,"stability_model_length"){
		init();
		printName();
		printParameters();
	}	
};

#endif /*DERIVEDMODELS_H*/
