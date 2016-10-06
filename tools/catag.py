#! /Users/ekstrom/Library/Enthought/Canopy_64bit/User/bin/python
#=======================================================================
import os
import argparse

wgreductor = 'wg_reductor.py'

parser = argparse.ArgumentParser(description='Arguments for the creation of .wg and .wa files', \
                                 usage='catag.py #star_name' \
                                 '\n--------------------------------------------------------' \
                                 '\n\nYou can create automatically a reduced wg file with option -r' \
                                 '\n\nMore details on the options by calling catag.py -h' \
                                 '\n--------------------------------------------------------\n')

parser.add_argument('StarName',help='Star name.',type=str)
parser.add_argument('-r','--red',help='create a reduced file.',action='store_true')
parser.add_argument('-f','--forced',help='replaces pre-existing file.',action='store_true')
args = parser.parse_args()

StarName = args.StarName
reduce = args.red
forced = args.forced
flag = ''

os.system('find . -name '+StarName+'".g0*" -print0 | sort -z | xargs -0 cat > '+StarName+'.wg')

if reduce:
    if forced:
        flag = ' -f'
    os.system(wgreductor+flag+' '+StarName)
