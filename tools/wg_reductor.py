#! /usr/bin/env python
#=======================================================================
import os
import sys
import argparse
import numpy as np
from six.moves import input

parser = argparse.ArgumentParser(description='Arguments for the reduction of .wg files', \
                                 usage='wg_reductor.py #star_name' \
                                 '\n--------------------------------------------------------' \
                                 '\nOptions:' \
                                 '\n\n-f to overwrite an existing file' \
                                 '\n\n-p to keep prezams lines' \
                                 '\n\n-h to get more details on the options' \
                                 '\n--------------------------------------------------------\n')

parser.add_argument('StarName',help='Star name.',type=str)
parser.add_argument('-f','--forced',help='replaces pre-existing file.',action='store_true')
parser.add_argument('-p','--preZAMS',help='keeps preMS lines.',action='store_true')
args = parser.parse_args()

StarName = args.StarName
forced=args.forced
preMS=args.preZAMS

skipline = 0
delta = 300
colH1s = 5
colH1c = 21

try:
    answer = ''
    with open(StarName+'.dat','r') as datfile:
        if not forced:
            while answer == '':
                answer = input('File '+StarName+'.dat already exists, overwrite? y/n:')
            if answer in 'nN0':
                sys.exit(0)
            else:
                forced = True
except IOError:
    forced = False
with open(StarName+'.wg','r') as f:
    linesarray = np.array(f.readlines())
try:
	wgfile = np.loadtxt(StarName+'.wg',skiprows=skipline)
except ValueError as VE:
	print('!!! Value error in wgfile: {0}'.format(VE))
	sys.exit(0)
h1c = wgfile[:,colH1c]
hini = wgfile[0,colH1s]
try:
    ind_zams = np.where(abs(h1c-hini)>=3.e-3)[0][0]
    print('ZAMS line: {0}'.format(ind_zams))
except:
    ind_zams = 0
if preMS:
    ind_ini=0
else:
    ind_ini=ind_zams

time = wgfile[ind_ini:,1]
lum = wgfile[ind_ini:,3]
teff = wgfile[ind_ini:,4]
rhoc = wgfile[ind_ini:,19]
tc = wgfile[ind_ini:,20]

diff_time = time[-1]-time
diff_time[-1] = time[-1]-time[-2]
lgtime = np.log10(diff_time)

dt = (np.max(time)-np.min(time))/float(delta)
lgdt = (np.max(lgtime)-np.min(lgtime))/float(delta)
dlum = (np.max(lum)-np.min(lum))/float(delta)
dteff = (np.max(teff)-np.min(teff))/float(delta)
drho = (np.max(rhoc)-np.min(rhoc))/float(delta)
dtc = (np.max(tc)-np.min(tc))/float(delta)

list_index = [0]

for i,vars in enumerate(zip(time,lgtime,lum,teff,rhoc,tc)):
    cond0 = abs(vars[0]-time[list_index[-1]]) > dt
    cond1 = abs(vars[1]-lgtime[list_index[-1]]) > lgdt
    cond2 = abs(vars[2]-lum[list_index[-1]]) > dlum
    cond3 = abs(vars[3]-teff[list_index[-1]]) > dteff
    cond4 = abs(vars[4]-rhoc[list_index[-1]]) > drho
    cond5 = abs(vars[5]-tc[list_index[-1]]) > dtc
    if cond0 or cond1 or cond2 or cond3 or cond4 or cond5:
        list_index.append(i)

print('file with {0} lines reduced to {1} lines'.format(len(time),len(list_index)))
list_index = np.array(list_index)+ind_ini

if not forced:
    wmode = 'w'
else:
    wmode = 'w+'
with open(StarName+'.dat',wmode) as datfile:
    datfile.seek(0)
    for i in list_index:
        datfile.write(linesarray[i])
