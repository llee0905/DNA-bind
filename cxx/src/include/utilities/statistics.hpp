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

// Basic stats functions as header only implementation

#ifndef STATISTICS_H
#define STATISTICS_H

#include<iostream>
#include<vector>
#include<numeric>
#include<cmath>

namespace stats{

template<typename T>
inline double mean(const std::vector<T> &a_data){
	static_assert(std::is_arithmetic<T>::value,"statistics functions only valid for numerical types");
	return std::accumulate(a_data.begin(),a_data.end(),0.0)/static_cast<double>(a_data.size());
}

//one pass version currenty for reference only - testing suggests two pass equally fast
template<typename T,typename U>
inline double totalCovarianceOnePass(const std::vector<T> &a_uVec,const std::vector<U> &a_vVec){
	static_assert(std::is_arithmetic<T>::value && std::is_arithmetic<U>::value,"statistics functions only valid for numerical types");
	size_t n = std::min(a_uVec.size(),a_vVec.size());
	if (n<2) return 0.0;
	double cov = 0;
	double uMean = a_uVec.front();
	double vMean = a_vVec.front();
	for (size_t i=1;i<n;++i){  // Based on equation III.9 of "Numerically Stable, Single-Pass, Parallel Statistics Algorithms", Bennet et al. (2009)
		double div = 1.0/static_cast<double>(i+1); //adapt Boost version here for performance
		double uTemp = (a_uVec[i]-uMean)*div;
		double vTemp = a_vVec[i]-vMean;
		cov += (i*uTemp*vTemp);
		uMean += uTemp;
		vMean += (vTemp*div);
	}
	return cov;
}

template<typename T,typename U>
inline double totalCovariance(const std::vector<T> &a_data1,const std::vector<U> &a_data2){
	static_assert(std::is_arithmetic<T>::value && std::is_arithmetic<U>::value,"statistics functions only valid for numerical types");
	double mean1 = stats::mean(a_data1),mean2 = stats::mean(a_data2),cov=0.0;
	for (size_t i=0;i<std::min(a_data1.size(),a_data2.size());++i) 
		cov += (static_cast<double>(a_data1[i])-mean1)*(static_cast<double>(a_data2[i])-mean2);
	return cov;
}

template<typename T,typename U>
inline double covariance(const std::vector<T> &a_data1,const std::vector<U> &a_data2){
	static_assert(std::is_arithmetic<T>::value && std::is_arithmetic<U>::value,"statistics functions only valid for numerical types");
	return stats::totalCovariance(a_data1,a_data2)/static_cast<double>(std::min(a_data1.size(),a_data2.size()));
}

template<typename T,typename U>
inline double covarianceSample(const std::vector<T> &a_data1,const std::vector<U> &a_data2){
	static_assert(std::is_arithmetic<T>::value && std::is_arithmetic<U>::value,"statistics functions only valid for numerical types");
	return stats::totalCovariance(a_data1,a_data2)/static_cast<double>(std::min(a_data1.size(),a_data2.size())-1);
}

template<typename T>
inline double standardDeviation(const std::vector<T> &a_data){
	static_assert(std::is_arithmetic<T>::value,"statistics functions only valid for numerical types");
	return std::sqrt(stats::covariance(a_data,a_data));
}

template<typename T>
inline double standardDeviationSample(const std::vector<T> &a_data){
	static_assert(std::is_arithmetic<T>::value,"statistics functions only valid for numerical types");
	return std::sqrt(stats::covarianceSample(a_data,a_data));
}

template<typename T,typename U>
inline double relativeError(const std::vector<T> &a_scores, const std::vector<U> &a_refScores){
	static_assert(std::is_arithmetic<T>::value && std::is_arithmetic<U>::value,"statistics functions only valid for numerical types");
	double sum = 0.0;
	for (size_t i=0;i<a_scores.size();++i) 
		sum += std::fabs(static_cast<double>(a_scores[i])-static_cast<double>(a_refScores[i])/static_cast<double>(a_refScores[i]));
	return sum / static_cast<double>(a_scores.size());
}

template<typename T,typename U>
inline double squaredResiduals(const std::vector<T> &a_scores, const std::vector<U> &a_refScores){
	static_assert(std::is_arithmetic<T>::value && std::is_arithmetic<U>::value,"statistics functions only valid for numerical types");
	double sum = 0.0;
	for (size_t i=0;i<a_scores.size();++i) 
		sum += (static_cast<double>(a_scores[i])-static_cast<double>(a_refScores[i]))*(static_cast<double>(a_scores[i])-static_cast<double>(a_refScores[i]));
	return sum;
}

template<typename T,typename U>
inline double pearsonR(const std::vector<T> &a_scores, const std::vector<U> &a_refScores){
	static_assert(std::is_arithmetic<T>::value && std::is_arithmetic<U>::value,"statistics functions only valid for numerical types");
	double sig2Ref = stats::totalCovariance(a_refScores,a_refScores);
	double sig2Scores = stats::totalCovariance(a_scores,a_scores);
	double cov2 = stats::totalCovariance(a_refScores,a_scores);
	return cov2 / std::sqrt(sig2Ref*sig2Scores);		
}

template<typename T,typename U>
inline double rSquared(const std::vector<T> &a_scores, const std::vector<U> &a_refScores){
	static_assert(std::is_arithmetic<T>::value && std::is_arithmetic<U>::value,"statistics functions only valid for numerical types");
	double cov = stats::totalCovariance(a_refScores,a_refScores);
	double res = stats::squaredResiduals(a_scores,a_refScores);
	return  1.0 - res / cov;
}

}

#endif /*STATISTICS_H*/
