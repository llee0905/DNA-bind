##############################################################################
################## Copyright (C) 2020-2021, Richard Spinney. #################
##############################################################################
#                                                                            #
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
#                                                                            #
##############################################################################

# Description:

# Optimises and predicts a binding rate model for DNA strands
#
# k_{bind} = \sum_{i\in A} kappa * ((num)^alpha)/( 1+ exp[gamma + beta * df_i])
#
# where:
# 
# num: number of binding sites on a single sequence. i.e. (L - bindLength + 1), all measured in bases
# A: set of all initial binding nucleation states comprising all pairs of contiguous bindLength bases on both the strand and its complement
# kappa,alpha,gamma: parameters of the model
# L: length of DNA strand in bases
# df_i: free energy of binding structure i
# beta: Boltzmann pre-factor (1/ ( k_B * T))
#
# uses convention: -ve df means strongly bound (same as Nupack)

# Imports:

from dna_bind_functions import *

#####################
###### GLOBALS ######
#####################

########### PARAMETERS #################

params = parameters()
params.bindLength = 3         # length of binding state
params.material = "dna"       # NUPACK parameter: material choice - dna
params.stacking = "stacking"  # NUPACK parameter: stacking, nostacking, (backwards compatability versions: none-nupack3, some-nupck3, all-nupack3)
params.temperatureC = 25      # NUPACK parameter: temperature in celsius
params.NaConcentration = 0.15 # NUPACK parameter: Sodium ion concentration (M)
params.MgConcentration = 0.0  # NUPACK parameter: Magnesium ion concentration (M)
params.debug = 0              # 0: run normally, 1: print errors, 2: print binding states and scores

#################### INPUT: SEQUENCES AND RATES #########################

# reference rates/sequences - these must be in the same order.

sequences = ["GCTGTTCGGTCTAT","GTTCGGTCTA","ACCAAACCACCAAC","CCAAAACCAA","CAACACCAAACAAC","AAACCACACA","CAAAACCCCAACAC","AAACCCACCACACA","CAACACCCAAACAC","ACCAAACCAC","ACACCAAACC","CCAAAACCAACAAC","AAAAACCCAC","AAAAACCCACCCAA","CAACACCCAA","AAACCCACCA","ACAACACCAC","CAAAACCCCA","AACCAACACC","CAACACAACC","AAACCACCCAACAC","ACCAACACCA","CCACCAACAA","CCCAAACCCAACCA","ACCAACACCAACCA","CAACAACACCACCA","AACCAACACCACCA","AACCACCACAAACC","CAACACAACCAACC","CACCACAACCACCA","ACACACACAC","AACCACCACA","CAACCAACCA","CCCCACACAACAAC","CCACCAACAACAAC","ACACACACCA","ACACACACCACACA","CACACACACACACA","ACCAACCAACCAAC","CCCCACACAA","ACACCACCAC"]
scoresReference = [1042172.8,1145949.98,1186656.61,1668153.61,1779671.88,1809459.04,1831342.81,1892296.19,2124654.2,2167674.56,2178513.2,2304252.75,2364640.37,2433452.04,2511513.09,2527035.41,2559279.25,2651666.39,2655669.66,2695608.39,2745622.75,2765864.25,2825046.68,2975580.89,3186915.71,3242444.03,3629903.36,3659749.4,3763039.46,3805399.66,3843880.1,3865948.42,3960894.12,3961922.22,4062257.54,4159346.56,4442077.08,4818114.36,5208461.34,5448523.07,5613147.49]

sequencesNoOutliers = ["GCTGTTCGGTCTAT","GTTCGGTCTA","ACCAAACCACCAAC","CCAAAACCAA","CAACACCAAACAAC","AAACCACACA","CAAAACCCCAACAC","AAACCCACCACACA","CAACACCCAAACAC","ACCAAACCAC","ACACCAAACC","CCAAAACCAACAAC","AAAAACCCAC","AAAAACCCACCCAA","CAACACCCAA","AAACCCACCA","ACAACACCAC","CAAAACCCCA","AACCAACACC","CAACACAACC","AAACCACCCAACAC","ACCAACACCA","CCACCAACAA","CCCAAACCCAACCA","ACCAACACCAACCA","CAACAACACCACCA","AACCAACACCACCA","AACCACCACAAACC","CAACACAACCAACC","CACCACAACCACCA","AACCACCACA","CAACCAACCA","CCCCACACAACAAC","CCACCAACAACAAC","ACACACACCA","ACACACACCACACA","ACCAACCAACCAAC","CCCCACACAA","ACACCACCAC"]
scoresReferenceNoOutliers = [1042172.8,1145949.98,1186656.61,1668153.61,1779671.88,1809459.04,1831342.81,1892296.19,2124654.2,2167674.56,2178513.2,2304252.75,2364640.37,2433452.04,2511513.09,2527035.41,2559279.25,2651666.39,2655669.66,2695608.39,2745622.75,2765864.25,2825046.68,2975580.89,3186915.71,3242444.03,3629903.36,3659749.4,3763039.46,3805399.66,3865948.42,3960894.12,3961922.22,4062257.54,4159346.56,4442077.08,5208461.34,5448523.07,5613147.49]

######################## SORT INPUT #####################################

zipped = sorted(zip(scoresReference, sequences),key=lambda tuple: len(tuple[1]))
sequences = [x[1] for x in zipped]
scoresReference = [x[0] for x in zipped]

zipped = sorted(zip(scoresReferenceNoOutliers, sequencesNoOutliers),key=lambda tuple: len(tuple[1]))
sequencesNoOutliers = [x[1] for x in zipped]
scoresReferenceNoOutliers = [x[0] for x in zipped]

######################### END SEQUENCES #################################

######################### CHOOSE INPUT #################################

seqs = sequencesNoOutliers
scores = scoresReferenceNoOutliers

########################## HEADER/PREAMBLE ##############################

print("\n*****************************************************************")
print("************* Copyright Richard Spinney 2020-2021 ***************")
print("*****************************************************************","\n")
print("DNA-bind:  A simple model of DNA hybridisation rates")

######################### FUNCTION SIGNATURES/RETURN VALUES ############################

# run a model:
#
# 	signature:
#
#		modelName(sequences,modelParameters,params)
#
#			sequences: list of bases sequences - defined above
#			params: model/nupack parameters - defined above
#			modelParameters: vaues of alpha, kappa, gamma for the model 
#
#	 modelName + model parameters options:
#
#		1. combinatoricsModelNoLength(sequences,[kappa],params)
#		2. combinatoricsModelLength(sequences,[kappa,alpha],params)
# 		3. stabilityModelNoLength(sequences,[kappa,gamma],params)
# 		4. stabilityModelLength(sequences,[kappa,alpha,gamma],params)
#
#	returns:
#
#		list of predicted/model values
#
#
# optimise a model:
#
# 	signature:
#
#		optimiseModel(objectiveFunction, modelName,startModelParameters,referenceScores,sequences,params)
#
#			sequences: list of DNA strands as base sequences in 5'-3' direction - defined above
#			reference_scores: list of experimental hybridisation rates matched to sequences - defined above
#			params: model/nupack parameters - defined above
#			start_model_parameters: initial vaues of alpha, kappa, gamma for the model for the optimisation
#
#	objective function options:
#
#		1. objectivePearsonCorr - optimisation maximises correlation
#		2. objectiveSumSqResiduals - optimisation minimses square residuals
#
#	model_name / start_model_parameters choices/combinations:
#
#		1. combinatoricsModelNoLength / [startKappa]
#		2. combinatoricsModelLength / [startKappa,startAlpha]
# 		3. stabilityModelNoLength / [startKappa,startGamma]
# 		4. stabilityModelLength / [startKappa,startAlpha,startGamma]
#
#	returns:
#
#		x,modelPrediction,residuals,R2,rho,re
#
#		x: vector of optimum parameters [kappa,alpha,gamma] - where appropriate
#		modelPrediction: list of predicted/model values
#		residuals: sum of square residuals between model and data
#		R2: R-squared of prediction to data - defined as 1 - ((sum square residuals)/(total sum of squares))
#		rho: pearson correlation coefficient between model and data
#		re: mean relative error between model and data
#	

########################### RUN MODEL OPTIMISATIONS #####################

# optimise all models for binding lengths 1,2,3,4 using ojective_sum_sq_residuals objective function

#### LOOP OVER BIND LENGTHS ####
for length in range(1,5):
	
	params.bindLength = length

	#### LOOP OVER MODEL CHOICE ####
	for modelIndex in range(1,5):

		print("\n\n#####################################################################################")

		###################### OPTIMISATION INIT ############################
		# determines where to start the optimisation - try estimates:
		startKappa = 1e7 * length * length
		startGamma = 6.0 + 2.0 * length
		startAlpha = -2.0
		#################### END OPTIMISATION INIT ##########################

		############### SELECT MODEL + STARTING VECTOR ######################

		if modelIndex == 1:
			print("\nModel: combinatoricsModelNoLength")
			model = combinatoricsModelNoLength
			optVarsInit = [startKappa]
		elif modelIndex == 2:
			print("\nModel: combinatoricsModelLength")
			model = combinatoricsModelLength
			optVarsInit = [startKappa,startAlpha]
		elif modelIndex == 3:	
			print("\nModel: stabilityModelNoLength")
			model = stabilityModelNoLength
			optVarsInit = [startKappa,startGamma]
		elif modelIndex == 4:		
			print("\nModel: stabilityModelLength")
			model = stabilityModelLength
			optVarsInit = [startKappa,startAlpha,startGamma]
		else:
			print("\nOOB. ",modelIndex)
			break

		print("Binding length = ",params.bindLength,'\n')

		############# END SELECT MODEL + STARTING VECTOR ####################

		######################## OPTIMISE MODEL #############################

		optVars,modelPrediction,residuals,R2,rho,re = optimiseModel(objectiveSumSqResiduals,model,optVarsInit,scores,seqs,params)

		###################### END OPTIMISE MODEL ###########################

		######################### PRINT RESULT ##############################

		print("OPTIMISATION RESULT:\n")
		
		if modelIndex == 1:
			print("kappa = ",optVars[0])
		elif modelIndex == 2:
			print("kappa = ",optVars[0]," alpha = ",optVars[1])
		elif modelIndex == 3:	
			print("kappa = ",optVars[0]," gamma = ",optVars[1])
		elif modelIndex == 4:		
			print("kappa = ",optVars[0]," alpha = ",optVars[1]," gamma = ",optVars[2])
		else:
			print("\nOOB. ",modelIndex)
			break

		print("mean relative error = ",re," residuals = ",residuals)
		print("correlation = ",rho," R-squared = ",R2,"\n")
		print("Seq:,        Expt:,       Model:\n")
		
		for i in range(len(modelPrediction)):
			print(seqs[i],", ",scores[i],", ",modelPrediction[i])

		####################### END PRINT RESULT ############################

exit(0)
