##############################################################################
################## Copyright (C) 2020-2021, Richard Spinney. #################
##############################################################################
#                                                                            #
#    This program is free software: you can redistribute it and/or modify    #
#    it under the terms of the GNU General Public License as published by    #
#    the Free Software Foundation, either version 3 of the License, or       #
#    (at your option) any later version.                                     #
#                                                                            #
#    This program is distributed in the hope that it will be useful,         #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of          #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           #
#    GNU General Public License for more details.                            #
#                                                                            #
#    You should have received a copy of the GNU General Public License       #
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.   #
#                                                                            #
##############################################################################


####################################################################################################################
##### Functions and globals designed to be imported in file 'nupack_rate.py' which should accompany this code ######
####################################################################################################################

##############################
######### IMPORTS ############
##############################

import numpy as np 
from scipy.optimize import leastsq
from scipy.optimize import minimize
from scipy import stats
import multiprocessing
from joblib import Parallel, delayed
from datetime import datetime
from nupack import * #nupack 4.0

#################################
###### GLOBALS & CONSTANTS ######
#################################

### global constants
misMatchFreeEnergy = 9999999.9      # abitrarily unstable free energy used when binding between non-complementary base pairs
kcalToJ = 4185.8  				    # 1 kcal in Joules
NAvTimesKB = 6.0221409 * 1.38064852 # Avogadros number multiplied by Boltzmann constant

### globals for caching and multi_processing
isCached = 0
optimisationCounter = 0
freeEnergyCache = []
numCores = multiprocessing.cpu_count()
config.parallelism = True

###########################################
############ PARAMETER CLASS ##############
###########################################

class parameters:
  def __init__(self):
  	self.bindLength = 3
  	self.temperatureC = 25
  	self.NaConcentration = 0.15
  	self.MgConcentration = 0.0
  	self.stacking = "stacking"
  	self.debug = 0
  	self.material = "dna"
  
#######################
###### FUNCTIONS ######
#######################

# celsisus to kelvin
def celsiusToKelvin(temp):
	return temp + 273.15

# kelvin to celsisus
def kelvinToCelsius(temp):
	return temp - 273.15

#boltzmann factor in units of kcal/mol
def boltzmannFactor(temp):
	return kcalToJ / (NAvTimesKB*celsiusToKelvin(temp)) 

# returns R^2 value
def Rsquared(scoresReference,modelPrediction):
	meanRef = np.mean(scoresReference)
	residuals=0
	error=0
	for index in range(len(scoresReference)):
		residuals+= (scoresReference[index]-modelPrediction[index])**2
		error+= (scoresReference[index]-meanRef)**2
	return 1.0 - (residuals/error)

#returns mean relative error
def relativeError(scoresReference,modelPrediction):
	val = 0
	for index in range(len(scoresReference)):
		val += abs(scoresReference[index]-modelPrediction[index])/scoresReference[index]
	return val / len(scoresReference)

# returns correlation
def pearsonCorr(scoresReference,modelPrediction):
	meanRef = np.mean(scoresReference)
	meanModel = np.mean(modelPrediction)
	varRef  = 0
	varModel  = 0
	cov = 0
	for index in range(len(scoresReference)):
		varModel  += (modelPrediction[index]-meanModel)**2
		varRef  += (scoresReference[index]-meanRef)**2
		cov += (modelPrediction[index]-meanModel)*(scoresReference[index]-meanRef)
	return (cov / ((varRef*varModel)**0.5))

# returns reversed string. I.e. reverseStr("ABCD") returns "DCBA"
def reverseStr(s): 
  str = "" 
  for i in s: 
    str = i + str
  return str

# returns 1 if strings are complements, 0 otherwise. I.e. isComplement("GCTA","CGAT") returns 1, isComplement("GCTA","CGTT") returns 0
def isComplement(sequence1,sequence2):
	isComp = 1
	for i in range(len(sequence1)):
		if sequence1[i]=="A" and sequence2[i]!="T":
			isComp = 0
			break
		if sequence1[i]=="T" and sequence2[i]!="A":
			isComp = 0
			break
		if sequence1[i]=="C" and sequence2[i]!="G":
			isComp = 0
			break
		if sequence1[i]=="G" and sequence2[i]!="C":
			isComp = 0
			break
	return isComp

# returns crick-watson complement of sequence, getComplement("GCTA") returns "CGAT"
def getComplement(sequence):
	comp = ""
	for i in range(len(sequence)):
		if sequence[i] == "A":
			comp += "T"
		elif sequence[i] =="C":
			comp += "G"
		elif sequence[i] =="T":
			comp += "A"
		elif sequence[i] =="G":
			comp += "C"
		else:
			print("invalid sequence")
	return comp

# given a sequence, binding length and two initial offsets, returns nupack style binding pattern with parentheses
# i.e. getBindingPattern("AAAAAA",0,0,params) with params.bindLength=2 returns {"((....","....))"}
def getBindingPattern(sequence,index,indexC,params):
	length = params.bindLength
	debug = params.debug
	sub1 = sequence[index:index+params.bindLength]
	sub2 = getComplement(sequence[indexC:indexC+params.bindLength])
	bindPattern1 = ""
	bindPattern2 = ""
	if isComplement(sub1,sub2):
		for i in range(len(sequence)):
			if (i>=index) and (i<index+length):
				bindPattern1 += "("
			else:
				bindPattern1 += "."
		for i in range(len(sequence)):
			if (i>=indexC) and (i<indexC+length):
				bindPattern2 += ")"
			else:
				bindPattern2 += "."
		bindPattern2 = reverseStr(bindPattern2)
	else:
		for i in range(len(sequence)):
			bindPattern1 += "."
		for i in range(len(sequence)):
			bindPattern2 += "."
	bindPattern = []
	bindPattern.append(bindPattern1)
	bindPattern.append(bindPattern2)
	return bindPattern

def getBindingPatternNoMatch(sequence,index,indexC,params):
	length = params.bindLength
	debug = params.debug
	bindPattern1 = ""
	bindPattern2 = ""
	for i in range(len(sequence)):
		if (i>=index) and (i<index+length):
			bindPattern1 += "("
		else:
			bindPattern1 += "."
	for i in range(len(sequence)):
		if (i>=indexC) and (i<indexC+length):
			bindPattern2 += ")"
		else:
			bindPattern2 += "."
	bindPattern2 = reverseStr(bindPattern2)
	bindPattern = []
	bindPattern.append(bindPattern1)
	bindPattern.append(bindPattern2)
	return bindPattern

# given a sequence, binding length and two initial offsets, returns all sequences, in order, in the 5-3 direction
# plus the nupack style binding pattern with parentheses
def getStructure(sequence,index,indexC,params):
	length = params.bindLength
	debug = params.debug
	bindingModel = []
	strands = []
	bindingSequence=""
	strands.append(sequence)
	strands.append(reverseStr(getComplement(sequence)))
	bindPattern = getBindingPattern(sequence,index,indexC,params)
	bindingSequence+=bindPattern[0]
	bindingSequence+="+"+bindPattern[1]
	bindingModel.append(strands)
	bindingModel.append(bindingSequence)
	return bindingModel

#Nupack 4.0 interface
def getFreeEnergy(sequences,bindings,mat,stack,temp,na,mg):
	model = Model(material=mat,ensemble=stack,celsius=temp,sodium=na,magnesium=mg)
	struct = bindings[0]
	for index in range(1,len(bindings)):
		struct += "+"
		struct += bindings[index]
	dGstruc = structure_energy(strands=sequences,structure=struct,model=model) #nupack function
	return dGstruc

# given a sequence, binding length and two initial offsets, returns reported free energy difference between secondary structure and unbound ligand of the secondary structure using NUPACK
def dfNupack(sequence,index,indexC,params):
	bindingModel = getStructure(sequence,index,indexC,params)
	val = getFreeEnergy(bindingModel[0],bindingModel[1],params.material,params.stacking,params.temperatureC,params.NaConcentration,params.MgConcentration)
	if bindingModel[1].find("(")==-1:
		val = misMatchFreeEnergy;
	if params.debug == 2:
		print(bindingModel[0][0]+"+"+bindingModel[0][1])
		bonds = getBindingPatternNoMatch(sequence,index,indexC,params)
		print(bonds[0]+"+"+bonds[1])
		print ("df = ",val)
	return val

# logit probability model
def logit(gamma,freeEnergy,params):
	temperatureC = params.temperatureC
	preFactor = boltzmannFactor(temperatureC)
	if gamma+preFactor*freeEnergy > 999999.9:  # avoid overflow
		return 0.0
	elif gamma+preFactor*freeEnergy < -999999.9:
		return 1.0
	return 1.0/(1.0 + np.exp(gamma+preFactor*freeEnergy))

# main loop for stability model
def mainLoop(sequences,gamma,params):
	scores = []
	for seqIndex, sequence in enumerate(sequences):  # loop over sequences
		sumVal = 0
		for ligandStart in range(0,len(sequence)-params.bindLength+1): #loop over ligand site
			for receptorStart in range(0,len(sequence)-params.bindLength+1): #loop over receptor site
				df = 0
				if isCached == 1:
					df = freeEnergyCache[seqIndex][ligandStart][receptorStart]
				else:
					df = dfNupack(sequence,ligandStart,receptorStart,params)
				sumVal += logit(gamma,df,params)
				if params.debug == 2:
					bindingModel = getStructure(sequence,ligandStart,receptorStart,params)
					strBind = bindingModel[0][0]+"+"+bindingModel[0][1]
					print(strBind)
					print(bindingModel[1])
					print("score = ",logit(gamma,df,params)," cumulative score = ",sumVal,"\n")
		scores.append(sumVal)
	return scores

# main loop for combinatorics model
def mainLoopCombo(sequences,params):
	scores = []
	for sequence in sequences:  # loop over sequences
		sumVal = 0
		for ligandStart in range(0,len(sequence)-params.bindLength+1): #loop over ligand site
			for receptorStart in range(0,len(sequence)-params.bindLength+1): #loop over receptor site
				val=1;
				pattern = getBindingPattern(sequence,ligandStart,receptorStart,params)
				if pattern[0].find("(")==-1:
					val = 0
				sumVal += val
				if params.debug == 2:
					strBind = sequence+"+"+reverseStr(getComplement(sequence))
					print(strBind)
					print(pattern[0]+"+"+pattern[1])
					print("score = ",val, "cumulative score = ",sumVal,"\n")
		scores.append(sumVal)
	return scores

# models in terms of min loops

def stabilityModelLength(sequences,optVars,params):
	kappa = optVars[0]
	alpha = optVars[1]
	gamma = optVars[2]
	scores = mainLoop(sequences,gamma,params)
	for seqIndex in range(len(sequences)):
		scores[seqIndex] *= kappa * ((len(sequences[seqIndex])-params.bindLength+1)**(alpha))
	return scores

def stabilityModelNoLength(sequences,optVars,params):
	kappa = optVars[0]
	gamma = optVars[1]
	scores = mainLoop(sequences,gamma,params)
	for seqIndex in range(len(sequences)):
		scores[seqIndex] *= kappa * ((len(sequences[seqIndex])-params.bindLength+1)**(-2.0))
	return scores

def combinatoricsModelLength(sequences,optVars,params):
	kappa = optVars[0]
	alpha = optVars[1]
	scores = mainLoopCombo(sequences,params)
	for seqIndex in range(len(sequences)):
		scores[seqIndex] *= kappa * ((len(sequences[seqIndex])-params.bindLength+1)**(alpha))
	return scores

def combinatoricsModelNoLength(sequences,optVars,params):
	kappa = optVars[0]
	scores = mainLoopCombo(sequences,params)
	for seqIndex in range(len(sequences)):
		scores[seqIndex] *= kappa * ((len(sequences[seqIndex])-params.bindLength+1)**(-2.0))
	return scores

# returns 1 - correlation between model and reference for optimisation
def objectivePearsonCorr(optVars,model,reference,sequences,params):
	global optimisationCounter
	optimisationCounter +=1
	results = model(sequences,optVars,params)
	diff = []
	rho = pearsonCorr(results,reference)
	if params.debug==1:
		print("loop: ",optimisationCounter)
		print("trying: ",x)
		print("p = ",rho)
	return 1.0 - rho
	
# returns sum of residuals between model and reference for optimisation
def objectiveSumSqResiduals(optVars,model,reference,sequences,params):
	global optimisationCounter
	optimisationCounter +=1
	results = model(sequences,optVars,params)
	diff = []
	for seqIndex in range(len(results)):
		diff.append(reference[seqIndex]-results[seqIndex])
	if params.debug==1:
		print("loop: ",optimisationCounter)
		print("trying: ",x)
		print("residuals = ",sum(np.square(diff)))
	return sum(np.square(diff))
	
def freeEnergyInner(sequence,params):
	sumVal = 0
	ligandDf = []
	for ligandStart in range(0,len(sequence)-params.bindLength+1): #loop over register
		receptorDf = []
		for receptorStart in range(0,len(sequence)-params.bindLength+1): #loop over register
			receptorDf.append(dfNupack(sequence,ligandStart,receptorStart,params))
		ligandDf.append(receptorDf)
	return ligandDf

# store free energy calculations for optimisation
def cacheFreeEnergies(sequences,params):
	global isCached
	global freeEnergyCache
	freeEnergyCache = []
	freeEnergyCache = Parallel(n_jobs=numCores, verbose=0)(delayed(freeEnergyInner)(sequence,params) for sequence in sequences)
	isCached = 1

# Wraps convex optimisation given an objective function, model, ground truth etc. Returns optimal parameters, the optimal prediction, the objective function and probabilities of binding
def optimiseModel(objectiveFunc,modelFunc,optVarsInit,groundTruth,sequences,params):
	if (modelFunc == stabilityModelLength) or (modelFunc == stabilityModelNoLength):
		now = datetime.now()
		currentTime = now.strftime("%H:%M:%S")
		print("starting free energy cache ",currentTime)
		cacheFreeEnergies(sequences,params)
		now = datetime.now()
		currentTime = now.strftime("%H:%M:%S")
		print("finished free energy cache ",currentTime)
	
	optMethod="nelder-mead"
	global optimisationCounter
	optimisationCounter = 0
	optRes = minimize(objectiveFunc,optVarsInit,args=(modelFunc,groundTruth,sequences,params), method=optMethod,options={'xatol': 1e-2, 'disp': False})
	modelPrediction = modelFunc(sequences,optRes.x,params)
	R2 = Rsquared(groundTruth,modelPrediction)
	rho = pearsonCorr(groundTruth,modelPrediction)
	relError = relativeError(groundTruth,modelPrediction)
	return optRes.x,modelPrediction,optRes.fun,R2,rho,relError
