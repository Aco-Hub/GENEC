#! /usr/bin/env python
#=======================================================================
import os
import glob
import argparse
import gzip
import numpy as np

parser = argparse.ArgumentParser(description='Arguments for the creation of an .input file from a .s file', \
                                 usage='MakeInput.py #star_name #num_model' \
                                 '\n--------------------------------------------------------' \
                                 '\n\nYou can create file from data in another directory with option -d [#dir_path]' \
                                 '\n\nMore details on the options by calling MakeInput.py -h' \
                                 '\n--------------------------------------------------------\n')

parser.add_argument('StarName',help='Star name.',type=str)
parser.add_argument('num',help='model number.',type=int)
parser.add_argument('-d','--dir',help='create a .input file from data in another directory.',default='')
args = parser.parse_args()

StarName = args.StarName
num = args.num
dir_path = args.dir

if dir_path:
    dir_path = os.path.expanduser(dir_path)

longNum = '{:0>7d}'.format(num)
sfile = StarName+'.s'+longNum
cleaned = False
zipped = False
modanfs = 'modanf='
nwseqs = 'nwseq='
params_end = '============='

try:
    with open(os.path.join(dir_path,sfile),'r') as sref:
        refs = sref.read()
        params = refs.split(params_end)[0]
except IOError:
    try:
        with gzip.open(os.path.join(dir_path,StarName+'.ws.gz'),'r') as ziprefs:
            cleaned = True
            refs = ziprefs.read()
            split_mods = refs.split('&CharacteristicsParams')
            params = ''
            models = {}
            searchmod = 'nwseq='+str(num)+'\n'
            print(searchmod)
            for i in range(1,len(split_mods)):
                if searchmod in split_mods[i]:
                    params = '&CharacteristicsParams'+split_mods[i].split(params_end)[0]
                key = int(split_mods[i].split('modanf=')[1].split('\n')[0])
                models[key] = split_mods[i].split('nwseq=')[1].split('\n')[0]
            if not params:
                print('problem while searching the parameters')
                exit(0)
    except IOError:
        print(os.path.join(dir_path,sfile)+' file not found')
        exit(0)

ibfile = params.rfind(modanfs)+len(modanfs)
ibfile_end = params[ibfile:].find('\n')
modanf = int(params[ibfile:ibfile+ibfile_end])
bfile = StarName+'.b{:0>5d}'.format(modanf)
print('looking for file '+bfile)
if dir_path:
    bfile = os.path.join(dir_path,bfile)
    if cleaned:
        blist = sorted(glob.glob(os.path.join(dir_path,'*.b[0-9]*')))
        blist_num = []
        for i in range(len(blist)):
            blist_num.append(int(blist[i].split('.b')[1].split('.gz')[0]))
        blist_num = np.array(blist_num)
        ind_closest = np.argmin(abs(blist_num-modanf))
        closestb = blist_num[ind_closest]
        if modanf-closestb < 0:
            nextb = blist_num[ind_closest-1]
        else:
            nextb = blist_num[ind_closest+1]
        if modanf-closestb != 0:
            print('\nBEWARE, the directory '+dir_path+' has been cleaned:')
            print(bfile.rsplit('/',1)[1]+' does not exist anymore\n')
            print("you'll have to start with model {0} or {1}\n".format(models[closestb],models[nextb]))
            exit(0)
        else:
            if bfile not in blist:
                if bfile+'.gz' in blist:
                    bfile = bfile+'.gz'
                    zipped = True
                else:
                    print('{0} not found in {1}'.format(bfile.rsplit('/',1)[1],dir_path))
    os.system('cp '+bfile+' .')
    if zipped:
        os.system('gunzip '+bfile[bfile.rfind('/')+1:])
    os.system('cp '+os.path.join(dir_path,'net*')+' .')

inputfile = open(StarName+'.input','w')
inputfile.write(params)
inputfile.close()
if not dir_path:
  try:
    with open('computation.log','a') as logfile:
      logfile.write('BACK TO MODEL '+str(num)+'\n')
  except IOerror:
    pass
print('input files created')
