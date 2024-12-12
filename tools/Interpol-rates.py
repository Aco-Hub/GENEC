####Script that interpolates energy rates into the temperature table that GENEC requires

#%%



import pandas as pd
import numpy as np





def readin_grid(filename):
    """Reads in the temperature grid that GENEC uses from a file directly
    """

    data = pd.read_fwf(filename,header=None)

    T_grid = data.to_numpy().flatten()
    return(T_grid)





def readin_rate(filename,header_size=15,reaction_number=0,rate_size=40):
    """Read in the rates produced by whatever scheme.
    """


    data = pd.read_fwf(filename,header=None,delimiter='  ',skiprows=header_size+reaction_number*(rate_size+header_size),nrows=rate_size,infer_nrows=rate_size)
    T = data[0].to_numpy().flatten()
    rates = data[1].to_numpy().flatten()

    return(T,rates)

def extract_header(filename,header_size=15,reaction_number=0,rate_size=40):
    """Reads in the text contained in the header of the rates file"""

    start = reaction_number*(header_size+rate_size)

    with open(filename,'r') as f:
        header = f.readlines()[start:start+header_size]


    return(header)


def interpolate_rate(T_in,T_out,rates_in):
    """Interpolates the rates to the new temperature grid
    """
    rates_out = np.interp(T_out,T_in,rates_in).flatten()

    return(rates_out)


def produce_output(output_grid_file,input_rates_file,output_name,number_of_reactions):
    """Interpolates the grid into a new one and produces the output needed"""



    output_temp = readin_grid(output_grid_file)

    with open(output_name,'w+') as f:
        for reaction_number in range(number_of_reactions):

            input_temp, input_rate = readin_rate(input_rates_file,reaction_number=reaction_number)

            output_rate = interpolate_rate(input_temp,output_temp,input_rate)
            header = extract_header(input_rates_file,reaction_number=reaction_number)

            #Write the desired output file


            f.writelines(header)
            for i in range(len(output_temp)):
                size_of_temp = len(f'{output_temp[i]:.4f}')
                number_of_spaces = 11 - size_of_temp
                f.write(number_of_spaces*' '+ f'{output_temp[i]:.4f} {  output_rate[i]:.4E}\n')



    return()

produce_output('temp_grid.txt','approx21_latest.txt','rates_interpolated_approx21.txt',number_of_reactions=20)

# extract_header('input_rates.txt')

#%%



d= pd.read_fwf('approx21_latest.txt',header=None,delimiter='  ',skiprows=15,nrows=40,infer_nrows=40)

d[0]
# %%
