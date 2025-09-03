#Class to manipulate the netgen format for the rates. Follows is some basic documentation.
#The reactions are written as A + B --> C + D 
#The netgen format is a header for the reaction in Question with speices A,B,C,D and the Q values.
#Qrad and Qnu where one is the neutrino losses and Qrad is then Qtotal - Qnu and the type identifier.
#Then the temperature grid and the rates are written in a table format.

#%%


from typing import Any
import pandas as pd
import numpy as np
import os

#Some defaults
header_size = 15 
rate_size = 60


#Dictonary of species to identfiy Z with the name of the species.
species = {0:'NEUT',1:'PROT',2:'HE',6:'C',7:'N',8:'O',9:'F',10:'NE',11:'NA',12:'MG',
           13:'AL',14:'SI',15:'P',16:'S',17:'CL',18:'AR',19:'K',20:'CA',22:'TI',
           24:'CR',26:'FE',27:'CO',28:'NI'}

#Make a dictonnary associating the species in netgen format to the mass excess.
mass_excess = {'OOOOO':0.0,'NEUT ':8.0713181,'PROT ':7.288971,'HE  4':2.4249156,'C  12':0.0,'C  13':3.12501129,'C  14':3.01989305,'N  14':2.86341704,
             'N  15':0.10143805,'O  16':-4.73700141,'O  17':-0.808813,'O  18':-0.781522,'F  18':0.873701,'F  19':-1.487386,'NE 20':-7.04193131,
             'NE 21':-5.731776,'NE 22':-8.024715,'NA 23':-9.52985358,'MG 24':-13.933567,'MG 25':-13.192826,'MG 26':-16.214582,'AL 26':-12.210309,
             'AL 27':-17.196658,'SI 28':-21.49279678,'SI 30':-24.432928,'P  31':-24.440885,'S  32':-26.015697,'S  34':-29.931788,'CL 35':-29.01354,
             'AR 36':-30.23154,'AR 38':-34.714551,'K  39':-33.807011,'CA 40':-34.846275,'CA 42':-38.547072,'TI 44':-37.548459,'TI 46':-44.123422,
             'CR 48':-42.81918,'CR 50':-50.259499,'CR 56':-55.281245,'FE 52':-48.331615,'FE 53':-50.945323,'FE 54':-56.252456,'FE 55':-57.479368,
             'FE 56':-60.605352,'CO 55':-54.027557,'CO 56':-56.039352,'CO 57':-59.344204,'NI 56':-53.903674}
            
             
             



def reac_in_file(file_name,reaction0):
    """Checks if a reaction is in a file.
    """

    reactions = load_all_reactions(file_name)


    for i in range(len(reactions)):
        if reaction0 == reactions[i]:
            return(i)
        


    return(-1)



def add_reactions(initial_reactions,reactions_to_add,output_file):

    ini_reactions = load_all_reactions(initial_reactions)
    add_reactions = load_all_reactions(reactions_to_add)

    #If output file exists delete it.
    try:
        os.remove(output_file)
    except:
        pass

    for i in range(len(add_reactions)):
        if add_reactions[i] not in ini_reactions:
            ini_reactions.append(add_reactions[i])
        else:
            print('Reaction already in file: ',add_reactions[i])
    
    for i in range(len(ini_reactions)):
        write_reaction(ini_reactions[i],output_file,'a')

        #Write a # for the last line.
    with open(output_file,'a') as f:
        f.write('#')
        f.close()
    
    return()


def readin_grid(filename):
    """Reads in the temperature grid desired for interpolation of the rates.
    """

    data = pd.read_fwf(filename,header=None)

    T_grid = data.to_numpy().flatten()
    return(T_grid)


def AZ_to_format(A,Z,two=False):
    """Give a A number and Z number and produces the string 
    identifier needed by NetGen.
    Input : 
          A : Mass number int
          Z : Atomic number int
    Optional 
         two : If one requires that the string identifier begins 
         with a 2. Bool
    Output: 
           String identifer 
           
    Remarks: 
            The special case of gammas is treated by A = 0.
            Total size of string should be 7                   """
    

    total_length = 7

    if two:
        amount = '2'
    else:
        amount = '1'

    if A > 0:
        species_name = species[Z]
        species_string = amount + ' ' + species_name 
        #Special cases of neutrons and protons.
        if (Z == 0 or Z == 1):
            return(species_string+' ')
        else:
            remaning_length = 7 - len(species_string) - len(str(A))

            species_string = species_string + remaning_length*' '+str(A)
            return(species_string)
    else:
        return('0 OOOOO')

    
def format_to_AZ(reaction):
    """Takes a reaction and outputs tuples of (A,Z) for each species."""

    try:
        Z1 = species[reaction.A[1:-2].strip()]
    except:
        Z1 = -1
    try:
        Z2 = species[reaction.B[1:-2].strip()]
    except:
        Z2 = -1
    try:
        Z3 = species[reaction.C[1:-2].strip()]
    except:
        Z3 = -1
    try: 
        Z4 = species[reaction.D[1:-2].strip()]
    except:
        Z4 = -1



    try:
        A1 = int(reaction.A[-2:])
    except:
        A1 = 0
    try:
        A2 = int(reaction.B[-2:])
    except:
        A2 = 0
    try:
        A3 = int(reaction.C[-2:])
    except:
        A3 = 0 
    try:
        A4 = int(reaction.D[-2:])
    except:
        A4 = 0


    if reaction.A == '1 PROT ' or reaction.A == '1 NEUT ':
        A1 = 1
    if reaction.B == '1 PROT ' or reaction.B == '1 NEUT ':
        A2 = 1
    
    if reaction.C == '1 PROT ' or reaction.C == '1 NEUT ':
        A3 = 1    
    if reaction.D == '1 PROT ' or reaction.D == '1 NEUT ':
        A4 = 1



    return([(A1,Z1),(A2,Z2),(A3,Z3),(A4,Z4)])




def convert_line_to_species(line,left_pos=14):
    """Converts a line of text from the header into a species object.
    If the rate file is in the normal netgen format then the left_pos is 
    correctly set.
    """

    species = line[left_pos:left_pos+7]

    return(species)

def convert_line_to_Q(line,left_pos=15):
    """Converts a line of text from the header into a species object.
    """

    Q = float(line[left_pos:].strip())
    
    return(Q)

def convert_line_to_type(line,left_pos=15):
    """Converts a line of text from the header into a species object.
    """

    type = line[left_pos:].strip()
    
    return(type)


def write_reaction(reaction,filename,writing_mode='w+'):
    """Writes a reaction object into a file.
    """
    #Header object has a specific sturcture to reproduce.
    with open(filename,writing_mode) as f:
        f.write('#\n')
        f.write('#'+13*' '+reaction.A+'\n')
        f.write('#'+17*' '+'+'+'\n')
        f.write('#'+13*' '+reaction.B+'\n')
        f.write('#'+17*' '+'='+'\n')
        f.write('#'+13*' '+reaction.C+'\n')
        f.write('#'+17*' '+'+'+'\n')
        f.write('#'+13*' '+reaction.D+'\n')
        f.write('#\n')
        f.write('#Qrad(Mev)'+(12-len(f'{reaction.Qrad:.3f}'))*' '+f'{reaction.Qrad:.3f}\n')
        f.write('#Qnu (Mev)'+(12-len(f'{reaction.Qnu:.3f}'))*' '+f'{reaction.Qnu:.3f}\n')
        f.write('#Type/Ne(cm-3)'+(8-len(reaction.type))*' '+reaction.type+'\n')
        f.write('#\n')
        f.write('#'+7*' '+'T8'+'\n')
        f.write('#'+'\n')
        #Proceed to write the temperature grid and the rates
        for i in range(len(reaction.T)):
            size_of_temp = len(f'{reaction.T[i]:.4f}')
            number_of_spaces = 11 - size_of_temp
            f.write(number_of_spaces*' '+ f'{reaction.T[i]:.4f} {  reaction.rates[i]:.3E}\n')
        f.close()

    return()



def update_reaction(reaction,filename):
    """Changes the reaction in the file to the new reaction.
    Load in whole file and rewrite it with new reaction updated.
    """

    reactions = load_all_reactions(filename)
    for i in range(len(reactions)):
        if reaction == reactions[i]:
            reactions[i] = reaction
            break
        
    with open(filename,'w+') as f:
        for i in range(len(reactions)):
            write_reaction(reactions[i],filename,'a')
    return()




def load_all_reactions(filename):
    """Loads all of the reactions in a file and returns a list of reaction objects.
    """

    #First find the number of reactions in the file.
    with open(filename,'r') as f:
        lines = f.readlines()
        number_of_reactions = 0
        #Loop through the lines and find the number of reactions.
        #One can count the number of 'T8' in the whole file.
        for i in range(len(lines)):
            if 'T8' in lines[i]:
                number_of_reactions += 1
        print('Number of reactions in file: ',number_of_reactions)


    #Now load all of the reactions
    reaction_list = []
    indiviual_species=[]
    for i in range(number_of_reactions):
        reaction1 = reaction()
        reaction1.load_reaction_from_number(filename,reaction_number=i)

        if reaction1.A[1:] not in indiviual_species:
            indiviual_species.append(reaction1.A[1:])
        if reaction1.B[1:] not in indiviual_species:
            indiviual_species.append(reaction1.B[1:])
        if reaction1.C[1:] not in indiviual_species:
            indiviual_species.append(reaction1.C[1:])
        if reaction1.D[1:] not in indiviual_species:
            indiviual_species.append(reaction1.D[1:])

        reaction_list.append(reaction1)

    print('Number of species in file: ',len(indiviual_species))

    return(reaction_list)

def change_reactions(rates_file,rates_to_be_changed_file,output_file):
    """In rates_file one has the entire set of rates used in the simulation.
     Rates_to_be_changed contains all of the rates that need to be changed in the main file.
    """

    #Store all the inital rates.
    ini_reaction_list = load_all_reactions(rates_file)
    to_be_changed_reaction_list = load_all_reactions(rates_to_be_changed_file)


    #If output file exists delete it.
    try:
        os.remove(output_file)
    except:
        pass

    for i in range(len(ini_reaction_list)):
        for j in range(len(to_be_changed_reaction_list)):
            if ini_reaction_list[i] == to_be_changed_reaction_list[j]:
                ini_reaction_list[i] = to_be_changed_reaction_list[j]
        write_reaction(ini_reaction_list[i],output_file,'a')


    #Write a # for the last line.
    with open(output_file,'a') as f:
        f.write('#')
        f.close()

    return()

def add_reactions(initial_reactions,reactions_to_add,output_file):

    ini_reactions = load_all_reactions(initial_reactions)
    add_reactions = load_all_reactions(reactions_to_add)

    #If output file exists delete it.
    try:
        os.remove(output_file)
    except:
        pass

    for i in range(len(add_reactions)):
        if add_reactions[i] not in ini_reactions:
            ini_reactions.append(add_reactions[i])
        else:
            print('Reaction already in file: ',add_reactions[i])
    
    for i in range(len(ini_reactions)):
        write_reaction(ini_reactions[i],output_file,'a')

        #Write a # for the last line.
    with open(output_file,'a') as f:
        f.write('#')
        f.close()
    
    return()

    


def compute_Qrad(reaction):
    """Using the mass excesss computed the Q rad of a reaction"""

    m1 = mass_excess[reaction.A[2:]] * int(reaction.A[0])
    m2 = mass_excess[reaction.B[2:]] * int(reaction.B[0])
    m3 = mass_excess[reaction.C[2:]] * int(reaction.C[0])
    m4 = mass_excess[reaction.D[2:]] * int(reaction.D[0])

    Q = m1 + m2 - m3 - m4

    Q = round(Q,3)

    return(Q)


class reaction(object):
    """Reaction class which contains the informtation about a given reaction
    A+B -> C+D it will contain A,B,C,D the Qrad and Qnu values and the type identifier.
    It also contains the temperature grid and reaciton rate grid.
    """

    def __init__(self,grid_size=50):
        self.A = None
        self.B = None
        self.C = None
        self.D = None
        self.Qrad = None
        self.Qnu = None
        self.type = None
        self.T = np.zeros(grid_size)
        self.rates = np.zeros(grid_size)

    def __repr__(self):
        return("Reaction: {} + {} -> {} + {} \n Qrad: {} Qnu: {} Type: {}".format(self.A,self.B,self.C,self.D,self.Qrad,self.Qnu,self.type))
    
    def __eq__(self,other):
        #A reaction is equal to the other is all speices are involved are the same. To be written in the right order.
        if self.A == other.A and self.B == other.B and self.C == other.C and self.D == other.D:
            return(True)
        else:
            return(False)

    def update_rates(self,new_rates):
        """Updates the rates of the reaction.
        """
        self.rates = new_rates

    def reaction_to_header(self,sp1,sp2,sp3,sp4,two=[False,False,False,False]):
        """Creates the reaction header as in netgen. 
            Inputs: 
                   (A,Z) of each species involved
            Outputs:
                    Reaction format with correct netgen
                    identifiers.
        """
        
        #Check reaction
        A_left = sp1[0]+sp2[0]
        A_right = sp3[0]+sp4[0]
        Z_left = sp1[1]+sp2[1]
        Z_right =sp3[1]+sp4[1]

        if (A_left != A_right) or (Z_left != Z_right):
            print("Warning reaction may be badly defined.")
            print("A conservation ",A_left,"?=?",A_right)
            print("Z conservation ",Z_left,"?=?",Z_right)
    

        self.A = AZ_to_format(sp1[0],sp1[1],two[0])
        self.B = AZ_to_format(sp2[0],sp2[1],two[1])
        self.C = AZ_to_format(sp3[0],sp3[1],two[2])
        self.D = AZ_to_format(sp4[0],sp4[1],two[3])

        

    def load_reaction_from_number(self,filename,header_size=15,reaction_number=0):
        """Loads a reaction from a file produced by netgen. The only thing needs to be known is the reaction
        number in the file. The reaction does not need to have any information loaded.
        """


        with open(filename,'r') as f:
            header = f.readlines()[reaction_number*(header_size+rate_size):reaction_number*(header_size+rate_size)+header_size]



        data = pd.read_fwf(filename,header=None,delimiter='  ',skiprows=header_size+reaction_number*(rate_size+header_size),nrows=rate_size,infer_nrows=rate_size)
        T = data[0].to_numpy().flatten()
        rates = data[1].to_numpy().flatten()
        self.T = T
        self.rates = rates

        #Line converted from the header calibarated to netgen format
        self.A = convert_line_to_species(header[1]) 
        self.B = convert_line_to_species(header[3])
        self.C = convert_line_to_species(header[5])
        self.D = convert_line_to_species(header[7])
        self.Qrad = convert_line_to_Q(header[9])
        self.Qnu = convert_line_to_Q(header[10])
        self.type = convert_line_to_type(header[11])





    def load_reaction_specific(self,filename):
        """Given a well defined reaction A B C D load from the file the reaction specifically.
        """
        
        idx_of_reaction = reac_in_file(filename,self)
        if idx_of_reaction >=0:
            with open(filename,'r') as f:
                header = f.readlines()[idx_of_reaction*(header_size+rate_size):idx_of_reaction*(header_size+rate_size)+header_size]


        data = pd.read_fwf(filename,header=None,delimiter='  ',skiprows=header_size+idx_of_reaction*(rate_size+header_size),nrows=rate_size,infer_nrows=rate_size)
        T = data[0].to_numpy().flatten()
        rates = data[1].to_numpy().flatten()

        self.T = T
        self.rates = rates
        self.Qrad = convert_line_to_Q(header[9])
        self.Qnu = convert_line_to_Q(header[10])
        self.type = convert_line_to_type(header[11])




    def interpolate_rate(self,T_out):
        """Interpolates the rates to the new temperature grid and changes the grid at the same time.
        Also NaNs are set to 0 if they appear in the rates.
        """
        where_is_nan = np.isnan(self.rates)

        rates_in_corr = np.where(where_is_nan,0.0,self.rates)


        rates_out = np.interp(T_out,self.T,rates_in_corr).flatten()

        self.T = T_out
        self.rates = rates_out

        return(rates_out)

    def apply_cutoff(self,rate_cutoff=1e-50):
        """Applies a cutoff to the rates below a certain value.
        """
        self.rates = np.where(self.rates<rate_cutoff,1e-98,self.rates)

    def clear_NaNs(self):
        """If there is a NaN in the rate set it to 0.
        """
        self.rates = np.where(np.isnan(self.rates),1e-98,self.rates)

    def copy_reaction(self,reaction):
        """Copies the reaction object.
        """
        self.A = reaction.A
        self.B = reaction.B
        self.C = reaction.C
        self.D = reaction.D
        self.Qrad = reaction.Qrad
        self.Qnu = reaction.Qnu
        self.type = reaction.type
        self.T = reaction.T
        self.rates = reaction.rates
        return()



