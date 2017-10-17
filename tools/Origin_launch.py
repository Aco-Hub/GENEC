#! /usr/bin/env python
#=======================================================================
import sys
import os
import argparse
import shutil
import time
import socket

def mymail(email,message):
	tmpfile = open('tmpf','w')
	tmpfile.write(message)
	tmpfile.close()
	current_dir = os.getcwd()
	MyCommand = 'cat - tmpf << EOF | /usr/sbin/sendmail -r'+email+' -t\nto:'+email+'\nsubject:'+current_dir+'\n\nEOF'
	os.system(MyCommand)
	os.system('rm tmpf')

# TO BE CHANGED BY USER ================================================================
email_adress = 'r.hirschi@keele.ac.uk'
default_prog = '/shen/hirschi/data/GENEC/git/GENEC/code/make/evosolL2016'
#=======================================================================================
current_dir = os.getcwd()
print 'Current dir: ',current_dir

parser = argparse.ArgumentParser(description='Arguments for the launch of stellar models computation', \
                                 usage='Origin_launch.py #star_name' \
                                 '\n--------------------------------------------------------' \
                                 '\nThe default program is:\n   '+default_prog+ \
                                 '\nYou can change it with option -e #prog' \
                                 '\n\nYou can specify a phase (or model) to stop at with option -p #phase (-m #model)' \
                                 '\n\nMore details on the options by calling Origin_launch.py -h' \
                                 '\n--------------------------------------------------------\n')

parser.add_argument('StarName',help='Star name.',type=str)
parser.add_argument('-e','--exe',help='Program to be used.',type=str,default=default_prog)
parser.add_argument('-p','--phase',help='Phase to stop.',type=int)
parser.add_argument('-m','--model',help='Model to stop.',type=int)
parser.add_argument('-i','--initial',help='Initial file.',type=str)
parser.add_argument('-z','--zip',help='Activate zipping all files after a series computation',action='store_true')
parser.add_argument('-c','--calcdir',help='Calculation directory.',type=str)
parser.add_argument('-N','--NoMail',help='Activate Nomail mode',action='store_false')

args = parser.parse_args()

ProgEvol = args.exe
StarName = args.StarName
phase_stop = args.phase
model_stop = args.model
initial_file = args.initial
calc_dir = args.calcdir
MailMode = args.NoMail
Zipping = args.zip

if calc_dir == None:
    calc_dir = ''
if initial_file == None:
    initial_file = ''

modanfs = 'modanf='
nwseqs = 'nwseq='
nzmods = 'nzmod='
phases = 'phase='
xcns = 'xcn='
CommandLaunch = ProgEvol+' < '+StarName+'.input'

print 'Prog: ',ProgEvol
print 'StarName: ',StarName
print 'Calc dir: ',calc_dir
if calc_dir != '':
    Zipping = True
if initial_file == '':
    InputFile = open(StarName+'.input','r')
    Inputs = InputFile.read()
    ibfile = Inputs.rfind(modanfs)+len(modanfs)
    ibfile_end = Inputs[ibfile:].find('\n')
    imod = Inputs.rfind(nwseqs)+len(nwseqs)
    imod_end = Inputs[imod:].find('\n')
    iphase = Inputs.rfind(phases)+len(phases)
    iphase_end = Inputs[iphase:].find('\n')
    modanf = int(Inputs[ibfile:ibfile+ibfile_end])
    nwseq = int(Inputs[imod:imod+imod_end])
    phase = int(Inputs[iphase:iphase+iphase_end])
    InputFile.close()

mytime = time.strftime('%A %d %B %Y at %H:%M:%S')
logfile = open('computation.log','a')
logfile.write('Computation on '+socket.gethostname()+' started on '+mytime+'\n')
logfile.write('Directory: '+current_dir+'\n')
if initial_file == '':
  logfile.write('starting model is '+str(nwseq)+', phase: '+str(phase)+'\n')
if calc_dir != '':
    logfile.write('Computation performed in directory: '+calc_dir+'\n')
logfile.write('Program used: '+ProgEvol+'\n')
if initial_file != '':
    logfile.write('Program launched on initial file '+initial_file+'\n')
logfile.write('----------------------------------------\n')
if phase_stop != None:
	logfile.write('Requested stop at phase '+str(phase_stop)+'\n')
if model_stop != None:
	logfile.write('Requested stop at model '+str(model_stop)+'\n')
logfile.close()

time_start = time.time()
answer = ''

answer = ''

if initial_file != '':
    if StarName+'.input' in os.listdir('.'):
        answer = raw_input('This star seems to be already partially computed.\n'+\
                       'Are you sure you want to proceed from file '+initial_file+\
                       '?\n yes(y) or no(n): ')
    if answer.lower()!='':
        sys.exit()
    elif answer.lower() in 'yes':
        os.system(ProgEvol+' < '+initial_file)
        try:
            runlog = open('runfile','r')
        except:
            sys.exit()
        runstat = runlog.read().strip(' \n\t')
        if runstat != 'running':
            if runstat != '':
			    print 'Program stopped with message: ',runstat
            else:
			    print 'Program aborted...'
            sys.exit()
        else:
	        if Zipping:
	            CommandZip = 'gzip -f '+StarName+'.[l,v,x,y]0000001 '+StarName+'.b00000 '+StarName+'.b00001 '+\
	                     StarName+'_StrucData_'+'0000001.dat'
	            os.system(CommandZip)
logfile.close()

relaunch_advection = True

if calc_dir != '':
    try:
        os.makedirs(calc_dir)
        print 'Calculation directory ',calc_dir,' successfully created'
    except OSError:
        print 'Calculation directory ',calc_dir,' already exists'
    try:
        os.remove(calc_dir+'/*')
    except:
        pass

while True:
    InputFile = open(StarName+'.input','r')
    Inputs = InputFile.read()
    ibfile = Inputs.rfind(modanfs)+len(modanfs)
    ibfile_end = Inputs[ibfile:].find('\n')
    imod = Inputs.rfind(nwseqs)+len(nwseqs)
    imod_end = Inputs[imod:].find('\n')
    iphase = Inputs.rfind(phases)+len(phases)
    iphase_end = Inputs[iphase:].find('\n')
    modanf = int(Inputs[ibfile:ibfile+ibfile_end])
    nwseq = int(Inputs[imod:imod+imod_end])
    phase = int(Inputs[iphase:iphase+iphase_end])
    InputFile.close()

    if phase_stop != None and phase == phase_stop:
        stop_message = 'Phase '+str(phase_stop)+' reached.'
        print stop_message
        if MailMode and len(email_adress) != 0:
            mymail(email_adress,stop_message)
        break

    if model_stop != None and nwseq > model_stop:
        stop_message = 'Model '+str(model_stop)+' reached.'
        print stop_message
        if MailMode and len(email_adress) != 0:
            mymail(email_adress,stop_message)
        break

    print 'New model ',nwseq,' with bfile ',modanf
    print 'Phase: ',phase

    needed_net = [i for i in os.listdir('.') if i[0:3] == 'net']
    needed_for_calc = needed_net+['input_changes.log',StarName+'.input']
    if os.path.isfile(StarName+'.b{0:05d}.gz'.format(modanf)):
        needed_for_calc = needed_for_calc+[StarName+'.b{0:05d}.gz'.format(modanf)]
    else:
        needed_for_calc = needed_for_calc+[StarName+'.b{0:05d}'.format(modanf)]

    if calc_dir != '':
        for file in needed_for_calc:
            shutil.copy2(file,calc_dir)
        os.chdir(calc_dir)

    try:
        os.remove('runfile')
    except:
        pass
    if Zipping:
        if os.path.isfile(StarName+'.b{0:05d}.gz'.format(modanf)):
            os.system('gunzip '+StarName+'.b{0:05d}.gz'.format(modanf))

    os.system(CommandLaunch)
    if Zipping:
        CommandZip = 'gzip -f '+StarName+'.[l,v,x,y]{0:07d} '.format(nwseq)+StarName+'.b{0:05d} '.format(modanf)+\
	                 StarName+'.b{0:05d} '.format(modanf+1)+StarName+'_StrucData_{0:07d}.dat'.format(nwseq)
        os.system(CommandZip)

    if calc_dir != '':
        result_files = [i for i in os.listdir('.') if i[0:4] == StarName[0:4]]
        result_files = result_files+['input_changes.log']
        for file in result_files:
#            os.rename(file,current_dir+'/'+file)
            shutil.move(file,current_dir+'/'+file)
        shutil.copy2('runfile',current_dir)
        os.chdir(current_dir)

    try:
        runlog = open('runfile','r')
    except:
        break
    runstat = runlog.read().strip(' \n\t')
    if runstat != 'running':
        if runstat != '':
            if 'ZAMS reached' in runstat:
                input_file = open(StarName+'.input','r')
                input_card = input_file.read()
                input_file.close()
                nwseq = int(input_card[input_card.rfind(nwseqs)+len(nwseqs):input_card.find('\n modanf')])
                modanf = int(input_card[input_card.rfind(modanfs)+len(modanfs):input_card.find('\n nzmod')])
                print 'ZAMS reached: NWSEQ, MODANF :',nwseq,modanf
                if Zipping:
                    if os.path.isfile(StarName+'.b{0:05d}.gz'.format(modanf)):
                        os.system('gunzip '+StarName+'.b{0:05d}.gz'.format(modanf))
                os.system(CommandLaunch)
                if Zipping:
                    CommandZip = 'gzip -f '+StarName+'.[l,v,x,y]{0:07d} '.format(nwseq)+StarName+'.b{0:05d} '.format(modanf)+\
	                             StarName+'_StrucData_{0:07d}.dat'.format(nwseq)
                    os.system(CommandZip)
                input_file = open(StarName+'.input','r')
                input_card = input_file.read()
                input_file.close()
                nzmod = int(input_card[input_card.rfind(nzmods)+len(nzmods):input_card.find('\n&END')])
                input_card = input_card.replace('nzmod='+str(nzmod),'nzmod=10')
                input_card = input_card.replace('gkorm=.300','gkorm=.100')
                input_file = open(StarName+'.input','w')
                input_file.write(input_card)
                input_file.close()
            elif 'Problem during advection' in runstat:
                if relaunch_advection:
                    input_file = open(StarName+'.input','r')
                    input_card = input_file.read()
                    input_file.close()
                    nwseq = int(input_card[input_card.rfind(nwseqs)+len(nwseqs):input_card.find('\n modanf')])
                    to_replace = "xcn="+input_card[input_card.rfind(xcns)+len(xcns):input_card.find('\n&END',input_card.rfind(xcns)+len(xcns))]
                    xcn = float(input_card[input_card.rfind(xcns)+len(xcns):input_card.find('\n&END',input_card.rfind(xcns)+len(xcns))])
                    input_card = input_card.replace(to_replace,'xcn=0.100')
                    input_file = open(StarName+'.input','w')
                    input_file.write(input_card)
                    input_file.close()
                    relaunch_advection = False
                else:
                    stop_message = 'Program stopped with message: '+runstat
                    break
            else:
                stop_message = 'Program stopped with message: '+runstat
                if MailMode and len(email_adress) != 0:
                    mymail(email_adress,runstat)
                break
        else:
            stop_message = 'Program aborted...'
            break
    
    if 'Problem during advection'  not in runstat:
        relaunch_advection = True

time_stop = time.time()
diff_time = int(round(time_stop-time_start))
m, s = divmod(diff_time, 60)
h, m = divmod(m, 60)
d, h = divmod(h, 24)
computation_time = '{0:d} day(s) {1:d} hour(s) {2:d} minutes and {3:d} seconds'.format(d,h,m,s)
logfile = open('computation.log','a')
logfile.write(stop_message+'\n')
logfile.write('Computation lasted '+computation_time+'\n')
logfile.write('----------------------------------------\n\n')
logfile.close()
