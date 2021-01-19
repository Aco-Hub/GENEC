#! /Users/ekstrom/Library/Enthought/Canopy_64bit/User/bin/python
#=======================================================================
import os
import sys
import argparse

parser = argparse.ArgumentParser(description='Arguments for program newstart.py', \
                                 usage='newstart.py #dir_name' \
                                 '\n--------------------------------------------------------' \
                                 '\n\nMore details on the options by calling newstart.py -h' \
                                 '\n--------------------------------------------------------\n')

parser.add_argument('DirName',help='Directory name for the initial files.',type=str)
args = parser.parse_args()
DirName = args.DirName

try:
  os.system('cp '+os.path.join(DirName,'ini*')+' .')
  os.system('cp '+os.path.join(DirName,'net*')+' .')
except IOError:
  print('problem with files ini* or net*')
