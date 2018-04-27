import os
import argparse

source_dir = os.path.dirname(os.path.abspath(__file__))
wgreductor = os.path.join(source_dir,'wg_reductor.py')

parser = argparse.ArgumentParser(description='Arguments for the creation of .wg and .wa files', \
                                 usage='catag.py #star_name' \
                                 '\n--------------------------------------------------------' \
                                 '\nOptions:' \
                                 '\n\n-r to create automatically a reduced wg file' \
                                 '\n\n-f to force overwriting of the reduced file' \
                                 '\n\n-p for keeping preZAMS lines in the reduced file' \
                                 '\n\n-h to get more details on the options' \
                                 '\n--------------------------------------------------------\n')

parser.add_argument('StarName',help='Star name.',type=str)
parser.add_argument('-r','--red',help='create a reduced file.',action='store_true')
parser.add_argument('-f','--forced',help='replaces pre-existing file.',action='store_true')
parser.add_argument('-p','--preZAMS',help='keeps preZAMS lines.',action='store_true')
args = parser.parse_args()

StarName = args.StarName
reduce = args.red
forced = args.forced
preMS = args.preZAMS
flag = ''

os.system('find . -name '+StarName+'".g0*" -print0 | sort -z | xargs -0 cat > '+StarName+'.wg')
print('wg file done')
if reduce:
    if forced:
        flag = ' -f'
        if preMS:
            flag = ' -fp'
    elif preMS:
        flag = ' -p'
    os.system('python '+wgreductor+' '+flag+' '+StarName)
