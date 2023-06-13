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





def readin_rates(filename,header_size=15):
    """Read in the rates produced by whatever scheme.
    """

    
    data = pd.read_fwf(filename,header=None,skiprows=header_size)

    T = data[0].to_numpy().flatten()
    rates = data[1].to_numpy().flatten()

    return(T,rates)

def extract_header(filename,header_size=15):
    """Reads in the text contained in the header of the rates file"""

    with open(filename,'r') as f:
        header = f.readlines()[0:header_size]


    return(header)


def interpolate_rates(T_in,T_out,rates_in):
    """Interpolates the rates to the new temperature grid
    """

    rates_out = np.interp(T_out,T_in,rates_in).flatten()

    return(rates_out)


def produce_output(output_grid_file,input_rates_file,output_name):
    """Interpolates the grid into a new one and produces the output needed"""

    output_temp = readin_grid(output_grid_file)

    input_temp, input_rates = readin_rates(input_rates_file)

    output_rates = interpolate_rates(input_temp,output_temp,input_rates)

    header = extract_header(input_rates_file)

    #Write the desired output file

    with open(output_name,'w') as f:
        f.writelines(header)
        for i in range(len(output_temp)):
            size_of_temp = len(f'{output_temp[i]:.4f}')
            number_of_spaces = 11 - size_of_temp
            f.write(number_of_spaces*' '+ f'{output_temp[i]:.4f} {  10**output_rates[i]/10**3:.4E}\n')


    return()

produce_output('temp_grid.txt','new_input_rates.txt','rates_interpolated_1000less.txt')

# extract_header('input_rates.txt')

#%%


# %%
