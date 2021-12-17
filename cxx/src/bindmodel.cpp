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

#include "derivedmodels.hpp"
#include "statistics.hpp"

template<typename T>
void BindModel<T>::printName() const{
	std::cout<<"************************************************************************************************************"<<std::endl;
	std::cout<<"   Model:"<<std::endl;
	std::cout<<"************************************************************************************************************"<<std::endl<<std::endl;
	std::cout<<"  "<<m_modelName<<std::endl<<std::endl;
}

template<typename T>
void BindModel<T>::printParameters() const{
	std::cout<<"************************************************************************************************************"<<std::endl;
	std::cout<<"   Parameters (Initial):"<<std::endl;
	std::cout<<"************************************************************************************************************"<<std::endl;
	m_paramsRef.print();
}

template<typename T>
void BindModel<T>::computeScores(Solution &a_modelData) const{
	const uint32_t numSeq = m_sequences.size();
	std::vector<double> scoresLoc,probabilitiesLoc,probabilitiesCondLoc;
	scoresLoc.reserve(numSeq);
	probabilitiesLoc.reserve(numSeq);
	probabilitiesCondLoc.reserve(numSeq);
	//core model loop
	for (uint32_t seqIndex=0;seqIndex<numSeq;++seqIndex){ //for each sequence
		double totSum = 0.0;
		int64_t count = 0, countCond = 0;
		const uint32_t numBindsSeq = m_sequences[seqIndex].size()-a_modelData.params().bindLength()+1; 
		const uint32_t numBindsSeqC = m_sequencesC[seqIndex].size()-a_modelData.params().bindLength()+1; 	
		for (uint32_t index=0;index<numBindsSeq;++index){ //register on strand
			for (uint32_t indexC=0; indexC<numBindsSeqC;++indexC){ //register on complement
				double score = static_cast<T const*>(this)->scoreFunction(seqIndex,index,indexC,a_modelData); //static polymorphism through CRTP idiom
				totSum += score; 
				++count;
				if (score>0) ++countCond;
			}
		}
		scoresLoc.emplace_back(a_modelData.params().kappa() * totSum * static_cast<T const*>(this)->lengthFunction(numBindsSeq,a_modelData));
		probabilitiesLoc.emplace_back(totSum/static_cast<double>(count));
		probabilitiesCondLoc.emplace_back(totSum/static_cast<double>(countCond));
	}
	//store in optimisation object
	a_modelData.scores() = std::move(scoresLoc);
	a_modelData.probabilities() = std::move(probabilitiesLoc);
	a_modelData.probabilitiesCond() = std::move(probabilitiesCondLoc);
}

template<typename T>
void BindModel<T>::printVarScores(const std::vector<std::vector<double> >&a_scores,const std::vector<std::vector<double> > &a_probabilities,const std::vector<std::vector<double> > &a_probabilitiesCond) const{	
	std::cout.precision(10);
	const std::vector<double> &rates = m_optimisationData.rates();
	const std::vector<double> &errors = m_errorsRef;
	const std::vector<std::string> &sequences = m_sequences;
	int32_t len = a_scores.front().size();
	std::cout<<"************************************************************************************************************"<<std::endl;
	std::cout<<"   Results for normal error distribution:"<<std::endl;
	std::cout<<"************************************************************************************************************"<<std::endl<<std::endl;
	std::cout<<std::left<<std::setw(16)<<"Seqs:"<<",   "<<std::setw(13)<<"Expt:"<<",   "<<std::setw(13)<<"Expt-errors:"<<",   "<<std::setw(20)<<"Model (val,std-dev,95%CI):"<<",   "<<std::setw(13)<<"P-hybridise (val,std-dev,95%CI):"<<",   "<<std::setw(15)<<"P-hybridise | match (val,std-dev,95%CI):"<<std::endl;
	for (size_t seqIndex=0; seqIndex <rates.size();++seqIndex)
		std::cout<<std::setw(16)<< sequences[seqIndex]<<",   "<<std::setw(13)<<rates[seqIndex]<<",   "<<std::setw(13)<<errors[seqIndex]
		<<",   "<<std::setw(10)<<stats::mean(a_scores[seqIndex])<<", "<<std::setw(10)<<stats::standardDeviationSample(a_scores[seqIndex])<<", ["<<a_scores[seqIndex][static_cast<int>(0.5*(1.0-0.95)*len)]<<","<<a_scores[seqIndex][static_cast<int>(0.5*(1.0+0.95)*len-1)]<<"]"
		<<",   "<<std::setw(10)<<stats::mean(a_probabilities[seqIndex])<<", "<<std::setw(10)<<stats::standardDeviationSample(a_probabilities[seqIndex])<<", ["<<a_probabilities[seqIndex][static_cast<int>(0.5*(1.0-0.95)*len)]<<","<<a_probabilities[seqIndex][static_cast<int>(0.5*(1.0+0.95)*len-1)]<<"]"
		<<",   "<<std::setw(10)<<stats::mean(a_probabilitiesCond[seqIndex])<<", "<<std::setw(10)<<stats::standardDeviationSample(a_probabilitiesCond[seqIndex])<<", ["<<a_probabilitiesCond[seqIndex][static_cast<int>(0.5*(1.0-0.95)*len)]<<","<<a_probabilitiesCond[seqIndex][static_cast<int>(0.5*(1.0+0.95)*len-1)]<<"]"
		<<std::endl;
	std::cout<<std::endl;
}


template<typename T>
void BindModel<T>::printScores() const{	
	std::cout.precision(10);
	const std::vector<double> &rates = m_optimisationData.rates();
	const std::vector<double> &scores = m_optimisationData.scores();
	const std::vector<double> &probabilities = m_optimisationData.probabilities();
	const std::vector<double> &probabilities_cond = m_optimisationData.probabilitiesCond();
	std::cout<<"************************************************************************************************************"<<std::endl;
	std::cout<<"   Results for optimised parameters:"<<std::endl;
	std::cout<<"************************************************************************************************************"<<std::endl<<std::endl;
	std::cout<<std::left<<std::setw(18)<<"Seqs:"<<",   "<<std::setw(15)<<"Expt:"<<",   "<<std::setw(15)<<"Expt-errors:"<<",   "<<std::setw(15)<<"Model:"<<",   "<<std::setw(15)<<"P-hybridise:"<<",   "<<std::setw(15)<<"P-hybridise | match:"<<std::endl;
	if (scores.size()==rates.size())
		for (size_t seqIndex=0; seqIndex <rates.size();++seqIndex)
			std::cout<<std::setw(18)<< m_sequences[seqIndex]<<",   "<<std::setw(15)<<rates[seqIndex]<<",   "<<std::setw(15)<<m_errorsRef[seqIndex]<<",   "<<std::setw(15)<<scores[seqIndex]<<",   "<<std::setw(15)<<probabilities[seqIndex]<<",   "<<std::setw(15)<<probabilities_cond[seqIndex]<<std::endl;
	std::cout<<std::endl;
}

template<typename T>
void BindModel<T>::run(){
	std::cout.precision(6);
	m_final=true;
	static_cast<T*>(this)->computeScores(m_optimisationData);
	m_final=false;
	printScores();
	const Solution& optimisationData = m_optimisationData;
	std::cout<<"************************************************************************************************************"<<std::endl;
	std::cout<<"   Summary statistics for point-estimate optimised results against data:"<<std::endl;
	std::cout<<"************************************************************************************************************"<<std::endl<<std::endl;
	std::cout<<" correlation = "<<stats::pearsonR(optimisationData.scores(),optimisationData.rates())<<" R-squared = "<<stats::rSquared(optimisationData.scores(),optimisationData.rates())<<std::endl;
	std::cout<<" mean relative error = "<<stats::relativeError(optimisationData.scores(),optimisationData.rates())<<" squared residuals = "<<stats::squaredResiduals(optimisationData.scores(),optimisationData.rates())<<std::endl<<std::endl;
	m_optimisationData.corr() = stats::pearsonR(optimisationData.scores(),optimisationData.rates());
	m_optimisationData.R2() = stats::rSquared(optimisationData.scores(),optimisationData.rates());
	m_optimisationData.sqRes() = stats::squaredResiduals(optimisationData.scores(),optimisationData.rates());
	m_optimisationData.meanRelError() = stats::relativeError(optimisationData.scores(),optimisationData.rates());
	m_optFlag=true;
}

/// optimisation functions

template<typename T>
void BindModel<T>::optimisationUpdate(const arma::vec &a_optVars,Solution *a_modelData){
	++a_modelData->optimCount();
	static_cast<T*>(this)->armaToParams(a_modelData->params(),a_optVars); //CRTP static poplymorphism 
	a_modelData->params().paramList() = a_optVars;
	static_cast<T*>(this)->computeScores(*a_modelData);
	if (a_modelData->params().debug()){
		std::cout<<"on count "<<a_modelData->optimCount()<<" ";
		for (const auto &var : a_optVars) std::cout<<var<<", ";
		std::cout<<"  R-squared: "<<stats::rSquared(a_modelData->scores(),a_modelData->rates())<<" corr: "<<stats::pearsonR(a_modelData->scores(),a_modelData->rates())<<"  sq-res: "<<stats::squaredResiduals(a_modelData->scores(),a_modelData->rates())<<std::endl;
	}
}

template<typename T>
double BindModel<T>::optimisationFunction(const arma::vec a_optVars,arma::vec* /*not used*/,void* a_modelData){
	Solution *modelData = static_cast<Solution*>(a_modelData);
	optimisationUpdate(a_optVars,modelData);
	return m_objectiveFunc(modelData->scores(),modelData->rates()); //minimise this
}

template<typename T>
bool BindModel<T>::optimiseInternal(Solution &a_modelData){
	bool success = false;
	if (m_numOptimParams){
		a_modelData.optimCount() = 0;
		////////// OPTIMISATION OPTIONS ///////////
		optim::algo_settings_t optimParams;
		optimParams.iter_max = MAXIMUM_ITERATIONS;
		optimParams.rel_objfn_change_tol = ERR_TOL;
		optimParams.rel_sol_change_tol = ERR_TOL;
		optimParams.vals_bound = true;
		optimParams.lower_bounds = m_optLower;
		optimParams.upper_bounds = m_optUpper;	
		///////// INITIAL VECTOR /////////////
		arma::vec optVars = arma::ones(m_numOptimParams,1);
		const arma::vec& optStart = m_optStart;
		for (int paramIndex=0;paramIndex<m_numOptimParams;++paramIndex) optVars(paramIndex) = optStart[paramIndex];
		///////// BIND OBJECTIVE FUNCTION /////////////
		if(a_modelData.params().objectiveFunc() == R_SQUARED) m_objectiveFunc = [](const std::vector<double> a,const std::vector<double> b) -> double {return 1.0-stats::rSquared(a,b);};
		else if(a_modelData.params().objectiveFunc() == SQUARED_RESIDUALS) m_objectiveFunc = stats::squaredResiduals<double,double>;
		else if(a_modelData.params().objectiveFunc() == RELATIVE_ERROR) m_objectiveFunc = stats::relativeError<double,double>;
		else return success;
		///////// RUN CHOSEN OPTIMISATION ALGO /////////////
		auto fn = std::bind(&BindModel::optimisationFunction,this,std::placeholders::_1,std::placeholders::_2,std::placeholders::_3); //make std::function func obj
		if (a_modelData.params().optimisationAlgorithm()== NM) success = optim::nm(optVars,fn,&a_modelData,optimParams);
		else if (a_modelData.params().optimisationAlgorithm() == DE) success = optim::de(optVars,fn,&a_modelData,optimParams);
		else if (a_modelData.params().optimisationAlgorithm() == PSO) success = optim::pso(optVars,fn,&a_modelData,optimParams);
		else return success;
		/////////////////////////////////////////////////
		if (success) static_cast<T*>(this)->armaToParams(a_modelData.params(),optVars); //adjust the parameters to the optimised values
	}
	return success;
}

template<typename T>
void BindModel<T>::initialiseOptimisationData(){
	m_optimisationData.rates() = m_ratesRef;
	m_optimisationData.params() = m_paramsRef;
}

template<typename T>
void BindModel<T>::optimise(){
	if (m_fatalError){
		std::cout<<"Fatal error, not proceeding."<<std::endl;
		return;
	}
	initialiseOptimisationData();
	std::cout<<"************************************************************************************************************"<<std::endl;
	std::cout<<"   Point estimate optimisation:"<<std::endl;
	std::cout<<"************************************************************************************************************"<<std::endl<<std::endl;
	std::cout<<"Optimising using point estimates of hybridisation rates..."<<std::endl;
	if (optimiseInternal(m_optimisationData)){
		std::cout << "Optimisation complete, optimised parameters:"<<std::endl;
		m_optimisationData.params().printVar();
	}
	else {
		std::cerr<<"Something went wrong with the optimisation."<<std::endl<<std::endl;
		return;
	}
	run(); //run for the optmised values
}

template<typename T>
void BindModel<T>::printDescriptiveStatistics(const std::vector<double> &a_sortedStatData,std::string a_statStr,std::string a_distStr) const{
	int len = a_sortedStatData.size();
	double stdDev = stats::standardDeviationSample(a_sortedStatData);
	std::cout<<"Statistics of "<<a_statStr<<" for "<<a_distStr<<" distribution:"<<std::endl<<std::endl;
	std::cout<<"mean "<<a_statStr<<" = "<<stats::mean(a_sortedStatData)<< " +/- "<<stdDev/std::sqrt(static_cast<double>(len))<<" (std-err)"<<std::endl;
	std::cout<<a_statStr<<" std-dev "<< stdDev <<std::endl;
	std::cout<<"median "<<a_statStr<<" = "<<a_sortedStatData[static_cast<int>(0.5*len)]<<std::endl;
	for (const auto CI : {0.5,0.75,0.95,0.99,0.999,0.9999})
		std::cout<<100*CI<<"% CI = ["<<a_sortedStatData[static_cast<int>(0.5*(1.0-CI)*len)]<<","<<a_sortedStatData[static_cast<int>(0.5*(1.0+CI)*len-1)]<<"],"<<std::endl;
	std::cout<<"100% CI = ["<<a_sortedStatData.front()<<","<<a_sortedStatData.back()<<"] (observed)"<<std::endl<<std::endl;
}

// error analysis

template<typename T>
void BindModel<T>::errorAnalysis(const uint64_t a_numSamples){
	if (m_fatalError){
		std::cout<<"Fatal error, not proceeding."<<std::endl;
		return;
	}
	std::cout.precision(6);
	std::cout<<"************************************************************************************************************"<<std::endl;
	std::cout<<"   Error Analysis:"<<std::endl;
	std::cout<<"************************************************************************************************************"<<std::endl<<std::endl;
	std::cout<<"Running error analysis for "<<a_numSamples<<" simulated data-sets"<<std::endl<<std::endl;
	std::cout<<"Achieved by adding zero mean Gaussian noise with std-dev equal to supplied error estimates to the supplied rates."<<std::endl<<std::endl;
	const std::vector<double> &rates = m_ratesRef;
	const std::vector<double> &errors = m_errorsRef;
	std::vector<double> ratesRand = rates;
	Parameters paramsInit = m_paramsRef;
	std::vector<double> corr,Rsq;
	corr.reserve(a_numSamples);
	Rsq.reserve(a_numSamples);
	const uint32_t numSequences = m_sequences.size();
	const uint32_t numParams = m_numOptimParams;
	std::vector<std::vector<double> > paramVals(numParams), modelRates(numSequences), modelProbs(numSequences), modelProbsCond(numSequences);
	RndGauss gaussRand = RndGauss();
	std::vector<std::future<bool> > results;
	std::vector<Solution> dataList;
	results.reserve(a_numSamples);
	dataList.reserve(a_numSamples); // essential to keep references in threads stable whilst extending dataList below
	for (uint64_t sampleIndex=0;sampleIndex<a_numSamples;++sampleIndex){
		for (size_t seqIndex=0; seqIndex<rates.size();++seqIndex) //normal distributed noise added to data
			ratesRand[seqIndex] = rates[seqIndex] + gaussRand(0,errors[seqIndex]);
		static_cast<T*>(this)->init();
		dataList.emplace_back(Solution(ratesRand,paramsInit)); //construct optimisation object
		results.push_back(m_pool.enqueue(&BindModel::optimiseInternal,this,std::ref(dataList.back()))); //queue the optimisation
	}
	for (uint32_t seqIndex=0;seqIndex<numSequences;++seqIndex){
		modelRates[seqIndex].reserve(a_numSamples);
		modelProbs[seqIndex].reserve(a_numSamples);
		modelProbsCond[seqIndex].reserve(a_numSamples);
	}
	for (uint64_t sampleIndex=0;sampleIndex<a_numSamples;++sampleIndex){
		bool success = results[sampleIndex].get();//get result of threaded computation
		if ((m_optimisationData.params().debug()) && (!success)) 
			std::cout<<"There was a problem with one of the optimisations."<<std::endl;
		//package the results
		corr.push_back(stats::pearsonR(dataList[sampleIndex].scores(),dataList[sampleIndex].rates()));
		Rsq.push_back(stats::rSquared(dataList[sampleIndex].scores(),dataList[sampleIndex].rates()));
		for (uint32_t seqIndex=0;seqIndex<numSequences;++seqIndex){
			modelRates[seqIndex].push_back(dataList[sampleIndex].scores(seqIndex));
			modelProbs[seqIndex].push_back(dataList[sampleIndex].probabilities(seqIndex));
			modelProbsCond[seqIndex].push_back(dataList[sampleIndex].probabilitiesCond(seqIndex));
		}
		for (uint32_t paramIndex=0;paramIndex<numParams;++paramIndex) 
			paramVals[paramIndex].push_back(dataList[sampleIndex].params().paramList(paramIndex));
	}
	//sort the results
	boost::sort::block_indirect_sort(corr.begin(),corr.end(),m_numThreads);
	boost::sort::block_indirect_sort(Rsq.begin(),Rsq.end(),m_numThreads);
	for (auto &x : paramVals) boost::sort::block_indirect_sort(x.begin(),x.end(),m_numThreads);
	for (auto &x : modelRates) boost::sort::block_indirect_sort(x.begin(),x.end(),m_numThreads);
	for (auto &x : modelProbs) boost::sort::block_indirect_sort(x.begin(),x.end(),m_numThreads);
	for (auto &x : modelProbsCond) boost::sort::block_indirect_sort(x.begin(),x.end(),m_numThreads);
	if (m_optimisationData.params().debug()>1){
		std::cout<<"normal distributed error ensemble correlation and R-squared values:"<<std::endl;
		for (size_t seqIndex=0;seqIndex<corr.size();++seqIndex) std::cout<<corr[seqIndex]<<", "<<Rsq[seqIndex]<<std::endl;
	}
	//diaplay the results
	printDescriptiveStatistics(corr,std::string("correlation"),std::string("normal error"));
	printDescriptiveStatistics(Rsq,std::string("R-squared"),std::string("normal error"));
	for (size_t paramIndex=0;paramIndex<paramVals.size();++paramIndex) 
		printDescriptiveStatistics(paramVals[paramIndex],m_paramNames[paramIndex],std::string("normal error"));
	printVarScores(modelRates,modelProbs,modelProbsCond);
	//store the goodness of fit distributions
	m_distFlag = true;
	m_corrDist = std::move(corr);
	m_R2Dist = std::move(Rsq);
}

//permutation analysis functions

template<typename T>
std::vector<double> BindModel<T>::generateBetaBias(const int32_t a_seed,const uint64_t a_num,const uint64_t a_n,const double a_p) const{
    boost::random::mt19937 rng(a_seed);
    boost::random::beta_distribution<> dist(static_cast<int>(a_n*a_p)+1,a_n-static_cast<int>(a_n*a_p)+1);
	std::vector<double> ret;
	ret.reserve(a_num);
	for (uint64_t count=0;count<a_num;++count) 
		ret.push_back((static_cast<double>(dist(rng)*a_n)+1)/static_cast<double>(a_n+1)); 
	// convert to conservative biased estimate detailed here: https://www.degruyter.com/document/doi/10.2202/1544-6115.1585/html and https://arxiv.org/abs/1603.05766
	return ret;
}

template<typename T>
std::vector<double> BindModel<T>::getPValueEstimates(const std::vector<double> &a_sortedNullData,const std::vector<double> &a_statData,uint64_t a_startIndex,uint64_t a_endIndex) const{
	uint64_t dataLen = a_sortedNullData.size();
	std::vector<double> pValues;
	pValues.reserve(a_endIndex-a_startIndex);
	for (uint64_t indexData=a_startIndex;indexData<a_endIndex;++indexData){//for each element of distribution
		const double &value = a_statData[indexData];
		int foundIndex = dataLen; // find smallest value of null greater than &value
		for (uint64_t indexCompare=0;indexCompare<dataLen;++indexCompare){
			if (a_sortedNullData[indexCompare]>value){
				foundIndex=indexCompare;
				break;
			}
		}
		pValues.push_back(static_cast<double>(dataLen-foundIndex)/static_cast<double>(dataLen)); //unbiased estimate of p-values here
	}
	return pValues;
}

template<typename T>
void BindModel<T>::pValuePoint(const std::vector<double> &a_sortedNullData,const double a_statVal,const std::string a_statStr) const{
	if (m_optFlag){ //optimisation data exists
		const int64_t dataLen = a_sortedNullData.size(); //num_permutations
		uint64_t foundIndex = dataLen; // find smallest value of null greater than a_statVal
		for (int index=0;index<dataLen;++index){
			if (a_sortedNullData[index]>a_statVal){
				foundIndex=index;
				break;
			}
		}
		double pValue = static_cast<double>(dataLen-foundIndex+1)/static_cast<double>(dataLen+1); // using conservative biased estimate detailed here: https://www.degruyter.com/document/doi/10.2202/1544-6115.1585/html and https://arxiv.org/abs/1603.05766
		double pValStdErr = std::sqrt(pValue*(1-pValue)/static_cast<double>(dataLen));
		std::cout<<"Point estimate "<<a_statStr<<" value of "<<a_statVal<<" greater than "<<100.0*(static_cast<double>(foundIndex)/static_cast<double>(a_sortedNullData.size()))
		<<"% of permuted values based on "<<dataLen<<" permutations."<<std::endl
		<<"This gives an implied p-value of "<<pValue<<" +/- "<<pValStdErr<<" using a conservative estimator"<<std::endl;
		for (const auto CI : {0.95,0.99,0.999}){
			double lower = boost::math::binomial_distribution<>::find_lower_bound_on_p(dataLen,dataLen-foundIndex,(1.0-CI)/2.0,boost::math::binomial_distribution<>::clopper_pearson_exact_interval);
			double upper = boost::math::binomial_distribution<>::find_upper_bound_on_p(dataLen,dataLen-foundIndex,(1.0-CI)/2.0,boost::math::binomial_distribution<>::clopper_pearson_exact_interval);
			std::cout<<"Clopper-Pearson interval at "<<100*CI<<"% = ["<<(lower*dataLen+1.0)/(dataLen+1.0)<<","<<(upper*dataLen+1.0)/(dataLen+1.0)<<"]"<<std::endl; //convert to conservative estimate
		}
		std::cout<<std::endl;
	}
}

template<typename T>
void BindModel<T>::pValueDist(const std::vector<double> &a_sortedNullData,const std::vector<double> &a_statData,const std::string a_statStr,const uint64_t a_numBeta) const{
	if (m_distFlag){ //distribution of resampled optimisation data exists
		uint64_t numPermutations = a_sortedNullData.size();
		uint64_t numData = a_statData.size();
		std::vector<double> pValues;
		pValues.reserve(numData);
		std::vector<std::future<std::vector<double> > > resultsP;
		for (int indexThread=0;indexThread<m_numThreads;++indexThread){ //partition pValue calculation off to threads
			uint64_t indexStartThread = indexThread * numData/m_numThreads;
			uint64_t indexEndThread = (indexThread+1) * numData/m_numThreads;
			if (indexThread==m_numThreads-1) indexEndThread = numData;
			resultsP.push_back(m_pool.enqueue(&BindModel::getPValueEstimates,this,std::ref(a_sortedNullData),std::ref(a_statData),indexStartThread,indexEndThread));
		}
		for (int indexThread=0;indexThread<m_numThreads;++indexThread){
			std::vector<double> temp = resultsP[indexThread].get(); //copy elision should happen
			pValues.insert(pValues.end(),std::make_move_iterator(temp.begin()),std::make_move_iterator(temp.end()));
		}
		std::vector<double> pValuesSample;
		pValuesSample.reserve(numData*a_numBeta);	
		std::vector<std::future<std::vector<double> > > resultsBeta;
		resultsBeta.reserve(numData);
		for (uint64_t indexData=0;indexData<numData;++indexData) // generate sub-ensemble of beta distributed random numbers
			resultsBeta.push_back(m_pool.enqueue(&BindModel::generateBetaBias,this,SEED_BASE+indexData,a_numBeta,numPermutations,pValues[indexData]));
		for (uint64_t indexData=0;indexData<numData;++indexData){
			std::vector<double> temp = resultsBeta[indexData].get(); //copy elision should happen
			pValuesSample.insert(pValuesSample.end(), std::make_move_iterator(temp.begin()),std::make_move_iterator(temp.end()));//store in total ensemble
		}
		boost::sort::block_indirect_sort(pValuesSample.begin(),pValuesSample.end(),m_numThreads); //sort all generated p=values
		printDescriptiveStatistics(pValuesSample,std::string("p-value"),std::string("normal error ")+a_statStr); // print properties
	}
}

template<typename T>
void BindModel<T>::permutationAnalysisStatistic(const std::vector<double> &a_sortedNullData,const std::vector<double> &a_statData,double a_statVal,std::string a_statStr,uint64_t a_numBeta) const{
	std::cout<<"************************************************************************************************************"<<std::endl;
	std::cout<<"   Permutation Analysis - "<<a_statStr<<":"<<std::endl;
	std::cout<<"************************************************************************************************************"<<std::endl<<std::endl;
	printDescriptiveStatistics(a_sortedNullData,a_statStr,std::string("null")); 
	pValuePoint(a_sortedNullData,a_statVal,a_statStr);
	pValueDist(a_sortedNullData,a_statData,a_statStr,a_numBeta);
}

template<typename T>
void BindModel<T>::permutationAnalysis(const uint64_t a_numPermutations,const uint64_t a_numBeta){
	if (m_fatalError){
		std::cout<<"Fatal error, not proceeding."<<std::endl;
		return;
	}
	std::cout.precision(6);
	std::cout<<"************************************************************************************************************"<<std::endl;
	std::cout<<"   Permutation Analysis:"<<std::endl;
	std::cout<<"************************************************************************************************************"<<std::endl<<std::endl;
	std::cout<<"Running permutation analysis for "<<a_numPermutations<<" data set permutations and"<<std::endl<<a_numBeta<<" beta-distributed p-value samples (where applicable)."<<std::endl<<std::endl;
	const std::vector<double> &rates = m_ratesRef;
	const std::vector<double> &errors = m_errorsRef;
	std::vector<double> ratesPerm = rates;
	Parameters paramsInit = m_paramsRef;
	std::vector<double> corr,Rsq;
	corr.reserve(a_numPermutations);
	Rsq.reserve(a_numPermutations);
	std::mt19937 mtGen{std::random_device{}()};
	RndGauss gaussRand = RndGauss();
	std::vector<std::future<bool> > results;
	std::vector<Solution> dataList;
	results.reserve(a_numPermutations);
	dataList.reserve(a_numPermutations); // required to keep vector references intact
	for (uint64_t count=0;count<a_numPermutations;++count){ 
		for (size_t seqIndex=0; seqIndex<rates.size();++seqIndex) //add normal noise to each point
			ratesPerm[seqIndex] = rates[seqIndex] + gaussRand(0,errors[seqIndex]);
		std::shuffle(ratesPerm.begin(),ratesPerm.end(),mtGen); //shuffle the data
		static_cast<T*>(this)->init();
		dataList.emplace_back(Solution(ratesPerm,paramsInit)); //initialise/store optimisation object
		results.push_back(m_pool.enqueue(&BindModel::optimiseInternal,this,std::ref(dataList.back()))); //queue job
	}
	for (uint64_t permIndex=0;permIndex<a_numPermutations;++permIndex){
		bool success = results[permIndex].get(); // get result of threaded computation
		if ((m_optimisationData.params().debug()) && (!success)) 
			std::cout<<"There was a problem with one of the optimisations."<<std::endl;
		//package results
		corr.push_back(stats::pearsonR(dataList[permIndex].scores(),dataList[permIndex].rates()));
		Rsq.push_back(stats::rSquared(dataList[permIndex].scores(),dataList[permIndex].rates()));
	}
	if (m_optimisationData.params().debug()>1){
		std::cout<<"permutation correlation and R-squared values:"<<std::endl;
		for (size_t i=0;i<std::min(corr.size(),Rsq.size());++i) std::cout<<corr[i]<<" "<<Rsq[i]<<std::endl;
	}
	//sort results
	boost::sort::block_indirect_sort(corr.begin(),corr.end(),m_numThreads);
	boost::sort::block_indirect_sort(Rsq.begin(),Rsq.end(),m_numThreads);
	//analyse/display results
	permutationAnalysisStatistic(corr,m_corrDist,m_optimisationData.corr(),std::string("correlation"),a_numBeta);
	permutationAnalysisStatistic(Rsq,m_R2Dist,m_optimisationData.R2(),std::string("R-squared"),a_numBeta);
}

//input processing functions

template<typename T>
void BindModel<T>::generateComplements(){
	for (const auto &seq : m_sequences) 
		m_sequencesC.push_back(m_dnaMethods.getComplement(seq));
}

template <class T>
void BindModel<T>::processInput(){
	for (const auto &seq : m_inputTuple){ //process sequences -  check for non base characters
		const std::string &str = std::get<0>(seq);
		auto it = std::find_if_not(str.begin(),str.end(),
			[](char base){
				if ((base=='A')||(base=='C')||(base=='G')||(base=='T')) return true;
				else return false;
			} 
		);
	    if (it != str.end()){
	        std::cout << "ERROR! Sequences contain characters other than 'A', 'C', 'T', or 'G'."<<std::endl;
	        m_fatalError = true;
	    }
	}
	std::sort(std::begin(m_inputTuple),std::end(m_inputTuple),  //sort the input
		[](std::tuple<std::string,double,double> a, std::tuple<std::string,double,double> b){
			if ( std::get<0>(a).size() != std::get<0>(b).size() ) //sort by length - then by rate
				return std::get<0>(a).size() < std::get<0>(b).size();
			else
				return std::get<1>(a) < std::get<1>(b); 
		}
	);
	//populate internal structures
	m_sequences.clear();
	m_ratesRef.clear();
	m_errorsRef.clear();
	for (const auto &seq : m_inputTuple){
		m_sequences.push_back(std::get<0>(seq));
		m_ratesRef.push_back(std::get<1>(seq));
		m_errorsRef.push_back(std::get<2>(seq));
	}		
	if (m_paramsRef.debug()){
		std::cout<<"Processed input:"<<std::endl<<std::endl;
		for (const auto &seq : m_inputTuple)
			std::cout<<std::left<<std::setw(18)<<std::get<0>(seq)<<std::setw(18)<<std::get<1>(seq)<<std::setw(18)<<std::get<2>(seq)<<std::endl;
		std::cout<<std::endl;
	}
	//generate complentary sequences and initialise the optimisation Solution object
	generateComplements();	
	initialiseOptimisationData();
}

//need to tell the compiler which class templates are going to be used...
template class BindModel<CombinatoricsModelLength>;
template class BindModel<CombinatoricsModelNoLength>;
template class BindModel<StabilityModelLength>;
template class BindModel<StabilityModelNoLength>;
