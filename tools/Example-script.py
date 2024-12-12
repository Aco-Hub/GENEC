####Small example script to show hwo one can use the functionalitys of the netgen class.
#%%

import netgen_class as ng
import matplotlib.pyplot as plt
# Read in the rate file

input_file = 'vit_GENET48.datCNEO'


#Store all the reactions contained in the input file.
all_reactions = ng.load_all_reactions(input_file)

#%%




############################################
############################################
############################################

#Looking at the reaction type.

############################################
############################################
############################################



#Print a reaciton to see the species involved and the Q values.
print(all_reactions[0])

#Look at the rates a function for T for this reaciton.
T_grid = all_reactions[0].T
rates = all_reactions[0].rates


fig,ax = plt.subplots(figsize=(10,8))

ax.plot(T_grid,rates,lw=2,color='k',label='Rate as a function of T')

ax.minorticks_on()
ax.tick_params(labelsize=12,direction='in')
ax.tick_params('both', length=10, width=2, which='major',right=True, top=True,direction='in')
ax.tick_params('both', length=5, width=1, which='minor',right=True, top=True,direction='in')

ax.set_xlabel(r'T [GK]',fontsize=20)
ax.set_ylabel(r'Rate',fontsize=20)


#%%


############################################
############################################
############################################

#Loading a specific reaction.

############################################
############################################
############################################


#Define a reaction

reaction_test= ng.reaction()

# C12 
reaction_test.A = ng.AZ_to_format(12,6)
#HE4 
reaction_test.B = ng.AZ_to_format(4,2)
#Gamma
reaction_test.C = ng.AZ_to_format(0,0)
#O16
reaction_test.D = ng.AZ_to_format(16,8)

#Set the Q value
reaction_test.Qrad = ng.compute_Qrad(reaction_test) #Compute the Q value from mass excess
reaction_test.Qnu = 0.0  #No neutrino loses. 
print(reaction_test)#Here T grid and rates are not defined yet.

# Note that the rates of the type of the reaction should not be set by hand but rather from a netgen
#output. So the functionality is not implemented to do this by hand.
#However if one wants to load a specific reaction and fill in the details defined by netgen one can run:


reaction_test.load_reaction_specific(input_file)
#Doing so reads through all the reactions until finding the one you want.
print(reaction_test) #Here everything is defined.


#%%


############################################
############################################
############################################

#Change the rate of a reaction in input file

############################################
############################################
############################################

#File path where we output the changed rates (for the test made random variations)
output_file = 'vit_GENET48_changed.datCNEO'

#Contains the rate or rates of the reactions to be changed.
rate_to_change_file = 'vit_C12.datCNEO'


#Takes the input file.
#Finds the reaction or reactions to be changed
#Outputs the updates rates in output_file
ng.change_reactions(input_file,rate_to_change_file,output_file)


#%%

#Check to see how it has worked.

reaction_changed = ng.load_all_reactions(rate_to_change_file)


reaction_before =ng.reaction()
reaction_before.copy_reaction(reaction_changed[0])
reaction_before.load_reaction_specific(input_file)

reaction_after =ng.reaction()
reaction_after.copy_reaction(reaction_changed[0])
reaction_after.load_reaction_specific(output_file)

fig,ax = plt.subplots(figsize=(10,8))

ax.plot(reaction_before.T,reaction_before.rates,lw=2,color='k',label='Rate before change')
ax.plot(reaction_after.T,reaction_after.rates,lw=2,color='r',label='Rate after change')

ax.minorticks_on()
ax.tick_params(labelsize=12,direction='in')
ax.tick_params('both', length=10, width=2, which='major',right=True, top=True,direction='in')
ax.tick_params('both', length=5, width=1, which='minor',right=True, top=True,direction='in')

ax.set_yscale('log')


#%%


############################################
############################################
############################################

#Interpolating temperature into a different grid

############################################
############################################
############################################
import numpy as np


all_reactions  =ng.load_all_reactions(input_file)

T_out = np.logspace(8,9,20)/10**9 #Peculiar grid to interpolate into. Put in GK

input_reaction = ng.reaction()
input_reaction.copy_reaction(all_reactions[0])




input_reaction.interpolate_rate(T_out) #Changes reaction in place. Updates the rates and the temperature grid.

rates_after = input_reaction.rates
T_after = input_reaction.T

rates_before = all_reactions[0].rates  
T_before = all_reactions[0].T

fig,ax = plt.subplots(figsize=(10,8))

ax.plot(T_before,rates_before,lw=2,color='k',label='Rate before interpolation')
ax.plot(T_after,rates_after,lw=2,color='r',label='Rate after interpolation')

ax.minorticks_on()
ax.tick_params(labelsize=12,direction='in')
ax.tick_params('both', length=10, width=2, which='major',right=True, top=True,direction='in')
ax.tick_params('both', length=5, width=1, which='minor',right=True, top=True,direction='in')

ax.set_yscale('log')

#Warning do not mix grid sizes of reactions in the final vitCNEO file it will not work.




# %%


#In summary the netgen class has some basic functionality with the following usage in mind.

## 1. Add reactions to the reaction file that weren't already there.
## 2. Change the T grid of a reaction to fit with others.
## 3. Change the rates of a specific set of reactions and then output them into net vit file.
