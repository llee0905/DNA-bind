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

#ifndef CONSTRUCTS_H
#define CONSTRUCTS_H

#include<iostream>
#include<vector>
#include<string>

#define OPTIM_ENABLE_ARMA_WRAPPERS
#include "OptimLib/optim.hpp"

#include "definitions.hpp"

enum OptimisationAlgorithmChoice { NM = 0, DE = 1, PSO = 2};// Nelder-Mead, Differential Evolution, Particle Swarm Optimisation
enum ObjectiveFunction { R_SQUARED = 0, SQUARED_RESIDUALS = 1, RELATIVE_ERROR = 2};
enum StackingOption { NONE = 0, SOME = 1, ALL = 2, STACKING = 3, NOSTACKING = 4}; //0-2 are nupack 3 options, 3-4 are nupack 4 options.
enum StrandMaterial { DNA = 0, RNA = 1, RNA99 = 2};

class Parameters{
private:
    //constants
    constexpr static const double m_KcalToJ = 4185.8;                   //1 kcal in Joules
    constexpr static const double m_NAvTimesKB = 6.0221409*1.38064852; //Avogadros number multiplied by Boltzmann constant
    //debug level
    int m_debug;
    //parameters
    uint16_t m_bindLength;
    double m_temperatureC,m_NaConcentration,m_MgConcentration;
    double m_gamma,m_alpha,m_kappa;
    bool m_stability,m_length;
    StrandMaterial m_material;
    StackingOption m_stacking;
    ObjectiveFunction m_objectiveFunc;
    OptimisationAlgorithmChoice m_optimisationAlgorithm;
    //compound/conditional constants
    double m_beta;
    // vector of parameters for loop access
    arma::vec m_paramList;
public:
    //setters
    void setDebug(const int a_debug){m_debug = a_debug;}
    void setTemperatureC(const double a_temp){
        m_temperatureC = a_temp;
        m_beta = m_KcalToJ / (m_NAvTimesKB*(m_temperatureC +  273.15) );
    }
    void setNaconcentration(const double a_Na){m_NaConcentration = a_Na;}
    void setMgconcentration(const double a_Mg){m_MgConcentration = a_Mg;}
    void setStability(const bool a_stability){m_stability = a_stability;}
    void setLength(const bool a_length){m_length = a_length;}
    void setGamma(const double a_gamma){m_gamma = a_gamma;}
    void setKappa(const double a_kappa){m_kappa = a_kappa;}
    void setAlpha(const double a_alpha){m_alpha = a_alpha;}
    void setOptimisationAlgorithm(const OptimisationAlgorithmChoice a_optimisationAlgorithm){m_optimisationAlgorithm = a_optimisationAlgorithm;}
    void setStacking(const StackingOption a_stacking){m_stacking = a_stacking;}
    void setObjectiveFunc(const ObjectiveFunction a_objectiveFunc){m_objectiveFunc = a_objectiveFunc;}
    void setBindLength(const uint16_t a_bindLength){m_bindLength = a_bindLength;}
    void setMaterial(const StrandMaterial a_material){m_material = a_material;}
    //getters
    const double& beta() const{return m_beta;}
    const int& debug() const{return m_debug;}
    const double& gamma() const{return m_gamma;}
    const double& kappa() const{return m_kappa;}
    const double& alpha() const{return m_alpha;}
    const double& temperatureC() const{return m_temperatureC;}
    const double& Naconcentration() const{return m_NaConcentration;}
    const double& Mgconcentration() const{return m_MgConcentration;}
    const uint16_t& bindLength() const{return m_bindLength;}
    const OptimisationAlgorithmChoice& optimisationAlgorithm() const{return m_optimisationAlgorithm;}
    const ObjectiveFunction& objectiveFunc() const{return m_objectiveFunc;}
    const StrandMaterial& material() const{return m_material;}
    const StackingOption& stacking() const{return m_stacking;}
    const bool& length() const{return m_length;}
    const bool& stability() const{return m_stability;}
    //by reference
    arma::vec& paramList(){return m_paramList;}
    template<typename T> double& paramList(T a_type){static_assert(std::is_integral<T>::value,"Access by integer type only"); return m_paramList(a_type);}
    Parameters()//some defaults
        :m_debug(0),
        m_bindLength(2),m_temperatureC(25.0),m_NaConcentration(0.15),m_MgConcentration(0.0),
        m_gamma(9.5),m_alpha(-2.0),m_kappa(1e7),
        m_stability(false),m_length(false),
        m_material(DNA),m_stacking(STACKING),
        m_objectiveFunc(R_SQUARED),m_optimisationAlgorithm(NM)    
    {}
    std::string getObjectiveFunctionString() const{
        if (m_objectiveFunc==SQUARED_RESIDUALS) return std::string("Squared residuals");
        else if (m_objectiveFunc==R_SQUARED) return std::string("(negative) R-Squared");
        else if (m_objectiveFunc==RELATIVE_ERROR) return std::string("Mean relative error");
        else return std::string("N/A");
    }
    std::string getStackingString() const{
        if (m_stacking==NONE) return std::string("none-nupack3");   
        else if (m_stacking==SOME) return std::string("some-nupack3");   
        else if (m_stacking==ALL) return std::string("all-nupack3");   
        else if (m_stacking==STACKING) return std::string("stacking");
        else if (m_stacking==NOSTACKING) return std::string("nostacking");
        else return std::string("N/A");
    }
    std::string getMaterialString() const{
        if (m_material==DNA) return std::string("dna");
        else if (m_material==RNA) return std::string("rna");
        else if (m_material==RNA99) return std::string("rna1999");  
        else return std::string("N/A");
    }
    std::string getMaterialStringVerbose() const{
        if (m_material==DNA) return std::string("dna2004");
        else if (m_material==RNA) return std::string("rna1995");
        else if (m_material==RNA99) return std::string("rna1999");  
        else return std::string("N/A");
    }
    std::string getOptimisationAlgoString() const{
        if (m_optimisationAlgorithm==NM) return std::string("Nelder-Mead");
        else if(m_optimisationAlgorithm==DE) return std::string("Differential Evolution");
        else if(m_optimisationAlgorithm==PSO) return std::string("Particle Swarm Optimisation");
        else return std::string("N/A");
    }
    void print() const{
        std::cout<<std::endl;
        std::cout<<"bind_length = "<<bindLength()<<std::endl;
        std::cout<<"optimisation algorithm = "<<getOptimisationAlgoString()<<std::endl;
        std::cout<<"objective function = "<<getObjectiveFunctionString()<<std::endl;
        if (stability()){
            std::cout<<"material = "<< getMaterialStringVerbose()<<std::endl;
            std::cout<<"temperatureC = "<<temperatureC()<<std::endl;
            std::cout<<"Naconcentration = "<<Naconcentration()<<std::endl;
            std::cout<<"Mgconcentration = "<<Mgconcentration()<<std::endl;
            std::cout<<"stacking option = "<<getStackingString()<<std::endl;
        }
        std::cout<<"kappa = "<<kappa()<<std::endl;
        if (stability()) std::cout<<"gamma = "<<gamma()<<std::endl; 
        if (length()) std::cout<<"alpha = "<<alpha()<<std::endl;
        std::cout<<std::endl;
    }
    void printVar() const{
        std::cout<<std::endl;
        std::cout<<"kappa = "<<kappa()<<std::endl;
        if (stability()) std::cout<<"gamma = "<<gamma()<<std::endl; 
        if (length()) std::cout<<"alpha = "<<alpha()<<std::endl;
        std::cout<<std::endl;
    }
};

class Solution{
private:
    int m_optimCount;
    std::vector<double> m_scores,m_rates,m_probabilities,m_probabilitiesCond;
    double m_corr,m_R2,m_sqRes,m_meanRelError;
    Parameters m_params;
public:
    //ctors
    Solution():m_optimCount(0){};
    Solution(const std::vector<double> a_rates,const Parameters a_params):m_optimCount(0),m_rates(a_rates),m_params(a_params){};
    //access functions
    int& optimCount(){return m_optimCount;}
    std::vector<double>& scores(){return m_scores;}
    std::vector<double>& rates(){return m_rates;}
    std::vector<double>& probabilities(){return m_probabilities;}
    std::vector<double>& probabilitiesCond(){return m_probabilitiesCond;}
    const std::vector<double>& scores() const{return m_scores;}
    const std::vector<double>& rates() const{return m_rates;}
    const std::vector<double>& probabilities() const{return m_probabilities;}
    const std::vector<double>& probabilitiesCond() const{return m_probabilitiesCond;}
    template<typename T> double& scores(T a_type){static_assert(std::is_integral<T>::value,"Access by integer type only"); return m_scores[a_type];}
    template<typename T> double& rates(T a_type){static_assert(std::is_integral<T>::value,"Access by integer type only"); return m_rates[a_type];}
    template<typename T> double& probabilities(T a_type){static_assert(std::is_integral<T>::value,"Access by integer type only"); return m_probabilities[a_type];}
    template<typename T> double& probabilitiesCond(T a_type){static_assert(std::is_integral<T>::value,"Access by integer type only"); return m_probabilitiesCond[a_type];}
    double& corr(){return m_corr;}
    double& R2(){return m_R2;}
    double& sqRes(){return m_sqRes;}
    double& meanRelError(){return m_meanRelError;}
    Parameters& params(){return m_params;}
    const Parameters& params() const{return m_params;}
};

#endif /*CONSTRUCTS_H*/
