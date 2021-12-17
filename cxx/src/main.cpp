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

// usage:
//
// Specify the parameters of the model using the parameters class/struct. (see lines 88-98 + 117-120)
//
// Specify/identify input, expected as a vector of 3-tuples of types <string,double,double>
// equal to the strand sequence in the 5'-3' direction, the experimental binding rate, and
// standard error in the binding rate. (see line 101 and file "include/sequences.hpp").
//
// Declare an object of a given model type (see lines 125, 130, 135, and 140). Types provided are:
//	1.	CombinatoricsModelNoLength  <--- prefactor kappa only
//	2.	CombinatoricsModelLength    <--- prefactor kappa + length dependence with alpha
//	3.	StabilityModelNoLength      <--- prefactor kappa + stability dependence with gamma
//	4.	StabilityModelLength        <--- prefactor kappa + length dependence with alpha + stability dependence with gamma
//
// Constructor arguments are 1. parameters object 2. input data as vector of 3-tuples 3. string containing relative path to python functions.
//
//	User defined models can be defined by deriving from classes bindmodel, combinatorics model, or stability model
//	all must define a length dependence, score function and function relating optimisation vector to parameter values (see derivedmodels.hpp).
//	Combinatorics_model implements a cache of matched binding locations of length bindingLength (params) on the strand and its complement.
//	Stability_model implements a cache of free energies of matched binding locations of length bindingLength (params),
//	as well as the logit stability function, on the strand and its complement.
//
// There are four public functions for all models which can be run:
//	1. model.run();
//	    This runs the model using the specified parameters and prints detailed results to the console.
//	2. model.optimise();
//	    This optimises the parameters of the model using the specified parameters as initial values,
//	    and the given experimental rates as point estimates (ignoring errors), and prints detailed results to the console.
//	3. model.errorAnalysis(int n);
//	    This optimises the parameters of the model using the specified parameters as initial values
//	    for n copies of the supplied rates with added Gaussian noise with standard deviations equal to the
//	    supplied error values. Correlation and R-squared data (mean, median, std-error, 95% CIs) are printed
//	    to the console.
//	4. model.permutationAnalysis(int n,int m);
//	    This optimises the parameters of the model using the specified parameters as initial values
//	    for n randomly permuted copies of the supplied rates with added Gaussian noise with standard deviation
//	    equal to the supplied error value. These randomly permuted data are treated as the null distribution.
//	    Using any previously calculated correlation and R-squared data from optimise() and error_analysis() calls
//	    inferred p-values and their confidence intervals are then estimated and printed to the console. The p-values
//	    are based on point-estimates and bootstrapped distributions in conjuction with the permutation analysis, respectively.
//	    The value m is used to generate ensembles of beta-distributed p-values to simulate the inevitably incomplete permutations
//	    that are used, in the case of the bootstrapped distribution where off-the-shelf confidence intervals are not available.
//	    It is therefore recommended that permutation_analysis be called after calling optimise() and/or error_analysis().
//
// In addition, stability models have the public function:
//	1. model.printFreeEnergies():
//	    This prints all possible nucleation states, for all sequences and their complements, alongside the free energy of binding
//	    reported by Nupack 4.0 for the parameter values specified in the parameters object, with binding site length specified in the
//	    parameters object

#include <iostream>
#include <chrono>
#include <thread>
#include "derivedmodels.hpp"
#include "sequences.hpp"

int main(int argc, char* argv[]) {

	std::string pythonPath = "./nupack4";
	if (argc>1) pythonPath = std::string(argv[1]);

	std::cout<<std::endl<<std::endl;
	std::cout<<"                 **********************************************************************"<<std::endl;
	std::cout<<"                 *************** Copyright (C) 2020-2021 Richard Spinney **************"<<std::endl;
	std::cout<<"                 **********************************************************************"<<std::endl<<std::endl;
	std::cout<<"                 **** DNA-bind: A simple model for fitting DNA hybridisation rates ****"<<std::endl<<std::endl;
	std::this_thread::sleep_for(std::chrono::milliseconds(1000));

	Parameters params;
	//experimental parameters
	params.setMaterial(DNA);
	params.setNaconcentration(0.15);
	params.setMgconcentration(0.0);
	params.setTemperatureC(25.0);
	params.setStacking(STACKING);
	// model / optimisation parameters
	params.setDebug(0); // 0: no debugging, 1: optimisation details, sorted input, 2: detailed scores for all binding states
	params.setObjectiveFunc(SQUARED_RESIDUALS);//SQUARED_RESIDUALS, RELATIVE_ERROR, R_SQUARED
	params.setOptimisationAlgorithm(NM); // NM, DE, PSO

	//select input out of options in "include/sequences.hpp"
	const auto &input = experimentData;//experimentDataNoOutliers;
	
	//choose to optimise data and perform variance/confidence analysis, and/or print free energies
	bool optimise = true;
	bool printFreeEnergies = false;

	if (optimise){

		uint64_t numResamples       = 200000;        // number of bootstrapped datasets to get error estimates on estimated quantities
		uint64_t numPermutations    = 200000;        // number of optimised permuted datasets to estimate p-values
		uint64_t numBetaResamples   = 2000;          // number of generated beta-distributed p-value samples to replicate uncertainty of l
		                                        // limited numbers of permutations. Only relevant when running permutation_analysis() 
		                                        // after error_analysis()

		for (auto bindLength : {1,2,3,4}){	//loop over bind lengths (for example)

			params.setBindLength(bindLength);
			params.setKappa(1e7*bindLength*bindLength); //initial guess for parameters
			params.setGamma(6.0+2.0*bindLength);
			params.setAlpha(-2.0);			

			// note: an ensemble of size numResamples * numBetaResamples will be generated for the permutation_analysis(),
			// following a call to error_analysis(). Recommended to keep num_binomial_resamples modest.

			CombinatoricsModelNoLength modelA(params,input,pythonPath);
			modelA.optimise();
			modelA.errorAnalysis(numResamples);
			modelA.permutationAnalysis(numPermutations,numBetaResamples);

			CombinatoricsModelLength modelB(params,input,pythonPath);
			modelB.optimise();
			modelB.errorAnalysis(numResamples);
			modelB.permutationAnalysis(numPermutations,numBetaResamples);

			StabilityModelNoLength modelC(params,input,pythonPath);
			modelC.optimise();
			modelC.errorAnalysis(numResamples);
			modelC.permutationAnalysis(numPermutations,numBetaResamples);

			StabilityModelLength modelD(params,input,pythonPath);
			modelD.optimise();
			modelD.errorAnalysis(numResamples);
			modelD.permutationAnalysis(numPermutations,numBetaResamples);
		}
	}
	
	if (printFreeEnergies){
		for (auto bindLength : {1,2,3,4}){
			params.setBindLength(bindLength);
			StabilityModelNoLength model(params,input,pythonPath);
			model.printFreeEnergies();
		}
	}

	return 0;

}

