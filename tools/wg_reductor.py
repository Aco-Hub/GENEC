#! /Users/ekstrom/Library/Enthought/Canopy_64bit/User/bin/python
#=======================================================================
import os
import argparse
import numpy as np

#prog = '/Users/ekstrom/OBS/Programs/UtilsEvol/filesFormat/Release/filesFormat.paf'

parser = argparse.ArgumentParser(description='Arguments for the reduction of .wg files', \
                                 usage='wg_reductor.py #star_name' \
                                 '\n--------------------------------------------------------' \
                                 '\n\nYou can overwrite an existing file with option -f' \
                                 '\n\nMore details on the options by calling catag.py -h' \
                                 '\n--------------------------------------------------------\n')

parser.add_argument('StarName',help='Star name.',type=str)
parser.add_argument('-f','--forced',help='replaces pre-existing file.',action='store_true')
args = parser.parse_args()

StarName = args.StarName
forced=args.forced

delta = 300

try:
    answer = ''
    datfile = open(StarName+'.dat','r')
    if not forced:
        while answer == '':
            answer = raw_input('File '+StarName+'.dat already exists, overwrite? y/n:')
        if answer in 'nN0':
            sys.exit(0)
        else:
            forced = True
    datfile.close()
except IOError:
    forced = False
with open(StarName+'.wg','r') as f:
    linesarray = np.array(f.readlines())
wgfile = np.loadtxt(StarName+'.wg')
h1c = wgfile[:,21]
hini = wgfile[0,5]
ind_zams = np.where(abs(h1c-hini)>=3.e-3)[0][0]
print 'ZAMS line:',ind_zams

time = wgfile[ind_zams:,1]
lum = wgfile[ind_zams:,3]
teff = wgfile[ind_zams:,4]
rhoc = wgfile[ind_zams:,19]
tc = wgfile[ind_zams:,20]

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

print 'file with',len(time),'lines reduced to',len(list_index),'lines'
list_index = np.array(list_index)+ind_zams

if not forced:
    datfile = open(StarName+'.dat','w')
else:
    datfile = open(StarName+'.dat','w+')
    datfile.seek(0)

for i in list_index:
    datfile.write(linesarray[i])
#print linesarray[list_index]
#datfile.write(linesarray[list_index])
datfile.close()
