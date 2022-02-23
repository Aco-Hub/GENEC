import re, sys, os.path
from os import path

#How many lines after a match do we look for the next section
henyey_methode_search_length = 10 #henyey-methode
j_vm_p_t_search_length = 5        #j    vm      p       t
i_Z_i_A_i_search_length = 10      #i) Z(i) A(i):
random_data_search_length = 5     # Search for the random data, I hope it always starts with '1 ' otherwise you could search for the first non-empty line starting with a number

line_number = 0
matched_line_number = 0

uvlpm = 0.00
vltm = 0.00
vlrm = 0.00
uvlpm_compare = 0.00
vltm_compare = 0.00
vlrm_compare = 0.00

#Function that takes a line, and checks if it is the 'uvlpm' line
def find_starting_data(line_data):
    #global is used to tell python that these variables are the global ones, and not to create local - function only - versions of them
    global matched_line_number
    global uvlpm
    global vltm
    global vlrm

    search_data = re.search('uvlpm.*', line_data)
    #If the data we want is not in this line we return False to the caller
    if not search_data:
        return False

    matched_line_number = line_number
    print('Found uvlpm on line', matched_line_number)

    #Multi part line here, it takes the found match 'search_data.group(0)', splits it into an array based on spaces
    #then removes any empty entries (where there were many spaces back to back), then turns this back into a list
    found = list(filter(None, search_data.group(0).split(" ")))
    print (found)
    #If all has gone well, found will be a list in the format ['uvlpm', number, 'vltm', number, 'vlrm', number
    if (found[0] != 'uvlpm='):
        raise "uvlpm not in correct place"
    uvlpm = float(found[1])

    if (found[2] != 'vltm='):
        raise "vltm not in correct place"
    vltm = float(found[3])

    if (found[4] != 'vlrm='):
        raise "vlrm not in correct place"
    vlrm = float(found[5])

    print('uvlpm =', uvlpm)
    print('vltm =', vltm)
    print('vlrm =', vlrm)
    return True

def find_second_data(line_data):
    global matched_line_number
    global uvlpm_compare
    global vltm_compare
    global vlrm_compare

    #This is a very vague search, but hopefully, with the way-markers used before, it should be correct
    search_data = re.search('1 .*', line_data)
    if not search_data:
        return False
    print('Found the data \'1 \' on line', line_number)

    #Multi part line here, it takes the found match 'search_data.group(0)', splits it into an array based on spaces
    #then removes any empty entries (where there were many spaces back to back), then turns this back into a list
    found = list(filter(None, search_data.group(0).split(" ")))

    #If all has gone well, found will be a list in the format ['1', something, 'uvlpm', 'vltm', 'vlrm', Lots of something
    uvlpm_compare = float(found[2])
    vltm_compare = float(found[3])
    vlrm_compare = float(found[4])
    print('uvlpm_compare =', uvlpm_compare)
    print('vltm_compare =', vltm_compare)
    print('vlrm_compare =', vlrm_compare)
    return True

def parse_data_file(file_name):
    global matched_line_number
    global line_number

    #Check the file exists
    if not path.exists(file_name):
        print("Could not open", file_name)
        print("This script needs a file name, if there are spaces in the file name, or the path, put in \"\"")
        return

    print("Opening file",file_name)
    #Open the file, storing a handle to the file as 'f'
    with open(file_name, 'r+') as f:
        #For this section we will read the file in line at a time, each 'for line in f:' block searches for one thing, then moves to the next block
        for line in f:
            line_number = line_number + 1
            #Search the line for the first occurence of our string - .* means read the whole line, from the found part until the end
            if (find_starting_data(line) is True):
                break

        #We now want to find the next line that looks memorable - 'henyey-methode', just as a way-marker
        for line in f:
            line_number = line_number + 1
            if (line_number > (matched_line_number + henyey_methode_search_length)):
                raise "Cannot find henyey_methode within henyey_methode_search_length lines"
            search_data = re.search('henyey-methode.*', line)
            if not search_data:
                continue
            matched_line_number = line_number
            break

        #Next way-marker, the line starting with 'j    vm      p       t'
        for line in f:
            line_number = line_number + 1
            if (line_number > (matched_line_number + j_vm_p_t_search_length)):
                raise "Cannot find j    vm      p       t within j_vm_p_t_search_length lines"
            search_data = re.search('j    vm      p       t.*', line)
            if not search_data:
                continue
            matched_line_number = line_number
            break

        #Final way-marker, the line starting with 'i) Z(i) A(i):'
        for line in f:
            line_number = line_number + 1
            if (line_number > (matched_line_number + i_Z_i_A_i_search_length)):
                raise "Cannot find i) Z(i) A(i): within i_Z_i_A_i_search_length lines"
            #Fun note here, '\' is the escape character.  In regular expressions brackets have meaning, so if you want to search for them, you need to add a '\' before it
            search_data = re.search('i\) Z\(i\) A\(i\).*', line)
            if not search_data:
                continue
            matched_line_number = line_number
            break

        #Now to find the final line, which - we hope - has our numbers
        for line in f:
            line_number = line_number + 1
            if (line_number > (matched_line_number + random_data_search_length)):
                raise "Cannot find the data within random_data_search_length lines"
            if find_second_data(line) is True:
                matched_line_number = line_number
                break

        #All going well, by this point you now have the data, and the _compare data
        print("Found the data!")

#Passes the file name (we hope) to parse_data_file
def main(argv):
    parse_data_file(argv[0])
    if (uvlpm > uvlpm_compare):
        print('Higher pressure in the topmost layer, change FITM')
    elif (vltm > vltm_compare):
        print('Higher temperature in the topmost layer, change FITM')
    elif (vlrm < vlrm_compare):
        print('Outer layer radius is smaller than inner layer, change FITM')

    else:
        print('All good in the hood!')
#If this script is being called directly from the command line, then call main
if __name__ == "__main__":
    if (len(sys.argv) < 2):
        print("needs a file name, if there are spaces in the file, or the path, put in \"\"")
        print("Examples:")
        print("    ", os.path.basename(__file__), "M200Z0V0.4.l0058231")
        print("    ", os.path.basename(__file__), "\"C:\data path\M200Z0V0.4.l0058231\"")
        sys.exit(1)
    main(sys.argv[1:])
