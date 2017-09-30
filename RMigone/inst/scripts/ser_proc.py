#!/usr/bin/env python

##############################
#
# process a series of files to generate apd and related files
#
##############################

import os, re

SCRIPT_PATH = '~/lib/R/library/RMigone/scripts/'

def f1(x):
    return re.search(r'.*\.ini.*',x)

def f2(x):
    return re.search(r'~',x) == None

files = os.listdir('.')
files = filter(f1,files)
files = filter(f2,files)
for i in files:
    n = re.sub(r'\.ini.*$','',i)
    print 'Processing ', n, '...'
    os.system(SCRIPT_PATH + 'pre_process.py ' + n)

R_script = """library(RMigone,lib.loc="~/lib/R/library")
papd <- apd.ser(patt=".*\\\\.apd.*$")
save(papd,file="papd.Rdata")"""

os.system("echo '" + R_script + "'| R --vanilla > r.out")
