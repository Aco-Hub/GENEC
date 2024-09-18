import argparse
import glob

StarName = glob.glob('*.input')[0][:-6]
new_mod_string1 = 'nouveau pas temporel modele'
new_mod_string2 = 'New timestep, model'

parser = argparse.ArgumentParser(description='Arguments for the search of largest correction', \
                                 usage='largest_corr.py #num' \
                                 '\n--------------------------------------------------------' )

parser.add_argument('numL',help='.l file number.',type=int)
parser.add_argument('-m','--mod_num',help='model number.',type=int,default=0)
args = parser.parse_args()

numL = args.numL
nwmd = args.mod_num
specific_model = False
print('----------------------------------------------------------------')
if nwmd != 0:
    specific_model = True
    new_mod_string1 = new_mod_string1+'{0: 6d}'.format(nwmd)
    new_mod_string2 = new_mod_string2+'{0: 6d}'.format(nwmd)

lname = StarName+'.l{0:07d}'.format(numL)
try:
    with open(lname) as lfile:
        file_content = lfile.readlines()
        old_l = False
        for s in file_content[0:100]:
            if new_mod_string1 in s:
                old_l = True
                break
        if old_l:
            new_mod_string = new_mod_string1
        else:
            new_mod_string = new_mod_string2
        if not specific_model:
            i_goingback = [i for i,x in enumerate(file_content) if 'GOING BACK' in x][0]
            i_mod = [i for i,x in enumerate(file_content[:i_goingback]) if new_mod_string in x][-1]
            nwmd = int(file_content[i_mod].replace(new_mod_string,'').rstrip('\n'))
            print('GOING BACK at line {0} (model {1})\n'.format(i_goingback,nwmd))
            i_largestcorr = [i for i,x in enumerate(file_content[:i_goingback]) if 'biggest correction' in x][-1]
        else:
            i_mod = [i for i,x in enumerate(file_content) if new_mod_string in x][0]
            i_goingback = [i for i,x in enumerate(file_content[i_mod:]) if 'GOING BACK' in x][0]
            print('Model {0} starts at line {1}, GOING BACK at line {2}\n'.format(nwmd,i_mod,i_mod+i_goingback))
            i_largestcorr = i_mod+[i for i,x in enumerate(file_content[i_mod:i_mod+i_goingback]) if 'biggest correction' in x][-1]
        print(file_content[i_largestcorr].lstrip(' ').rstrip('\n').replace(':',':\n'))
        print('----------------------------------------------------------------')
except IOError:
    print('File '+lname+' not found, check and try again.')
except IndexError:
    print('No ELEMENT NEGATIF in this set.')
