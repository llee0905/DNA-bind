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

#ifndef BINDMODEL_H
#define BINDMODEL_H

#include "boost/math/distributions/binomial.hpp"
#include "boost/random.hpp"
#include "boost/sort/sort.hpp"

#include "threads/ThreadPool.h"
#include "constructs.hpp"
#include "dna.hpp"
#include "definitions.hpp"
#include "random.hpp"

#define OPTIM_ENABLE_ARMA_WRAPPERS
#include "OptimLib/optim.hpp"

//using CRTP idiom to avoid run-time polymorphism with virtual functions
//we template the base class to create a version for each derived class that we might use
template <class T> 
class BindModel{
protected:
	//flags
	bool m_final,m_fatalError,m_distFlag,m_optFlag;
	//name
	std::string m_modelName;
	void printName() const;
	//parameters
	Parameters m_paramsRef;
	void printParameters() const;
	//data
	std::vector<std::tuple<std::string,double,double> > m_inputTuple;
	std::vector<std::string> m_sequences,m_sequencesC;
	std::vector<double> m_ratesRef,m_errorsRef;
	//dna specific methods
	Dna m_dnaMethods;
	//results/optimisation structure
	Solution m_optimisationData;
	std::vector<double> m_corrDist,m_R2Dist;
	std::vector<std::string> m_paramNames;
	//input functions
	void generateComplements();
	void processInput();
	//statistics functions
	std::vector<double> generateBetaBias(const int32_t,const uint64_t,const uint64_t,const double) const;
	void printDescriptiveStatistics(const std::vector<double>&,std::string,std::string) const;
	void permutationAnalysisStatistic(const std::vector<double>&,const std::vector<double> &,double,std::string,uint64_t) const;
	void pValuePoint(const std::vector<double>&,const double,const std::string) const;
	void pValueDist(const std::vector<double>&,const std::vector<double>&,const std::string,const uint64_t) const;
	std::vector<double> getPValueEstimates(const std::vector<double>&,const std::vector<double>&,uint64_t,uint64_t)const;
	//output functions
	void printScores() const;
	void printVarScores(const std::vector<std::vector<double> >&,const std::vector<std::vector<double> >&,const std::vector<std::vector<double> >&) const;
	//optimisation data
	int m_numOptimParams; //needs to be set in relevant constructor
	arma::vec m_optStart,m_optLower,m_optUpper; //needs to be set in relevant constructor
	//optimisation functions
	void initialiseOptimisationData();
	bool optimiseInternal(Solution&);
	void optimisationUpdate(const arma::vec&,Solution*);
	double optimisationFunction(const arma::vec,arma::vec*,void*);
	std::function<double(const std::vector<double>&,const std::vector<double>&)> m_objectiveFunc;
	//threading objects/data 
	const int m_numThreads = (std::thread::hardware_concurrency()==0) ? DEFAULT_THREADS : std::thread::hardware_concurrency();
	mutable ThreadPool m_pool;
	//internal core + utility functions
	void computeScores(Solution&) const;
	//CTOR protected to prevent instantiation of base class
	BindModel(Parameters a_params,std::vector<std::tuple<std::string,double,double> > a_input,std::string a_pythonPath,std::string a_name)
		:m_final(false),m_fatalError(false),m_distFlag(false),m_optFlag(false),
		m_modelName(a_name),m_paramsRef(a_params),m_inputTuple(a_input),m_dnaMethods(a_pythonPath),m_pool(m_numThreads)
	{
		std::cout<<std::endl<<"************************************************************************************************************"<<std::endl;
		std::cout<<"********************************************* INITIALISING MODEL *******************************************"<<std::endl;
		std::cout<<"************************************************************************************************************"<<std::endl<<std::endl;
		processInput(); //complements generated in this method
		if (m_fatalError) std::cout<<"Fatal error encountered. No further computation will proceed."<<std::endl; 			
	}
public:
	//public interface
	void run(); // run once for given parameters;
	void optimise(); // optimise paramters starting from given parameters 
	void errorAnalysis(const uint64_t);
	void permutationAnalysis(const uint64_t,const uint64_t);
};

#endif /*BINDMODEL_H*/
