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

# functions to get binding energies from nupack 4.0

from nupack import *
config.parallelism = True
def getFreeEnergy(sequences,bindings,mat,stack,temp,na,mg):
	model = Model(material=mat,ensemble=stack,celsius=temp,sodium=na,magnesium=mg)
	struct = bindings[0]
	for i in range(1,len(bindings)):
		struct += "+"
		struct += bindings[i]
	dGstruc = structure_energy(strands=sequences,structure=struct,model=model)
	return dGstruc