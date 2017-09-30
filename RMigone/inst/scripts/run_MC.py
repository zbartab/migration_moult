#!/usr/bin/env python

##############################
#
# run the MC simulation of migra1
#
##############################

#------------------------------ out_manip ------------------------------#
def out_manip(f, y_s, w_s, y_e, w_e, tipus, loc, vol):
    start = y_s*52 + w_s
    end = y_e*52 + w_e
    for t in range(start,end+1):
        y = t / 52
        w = t % 52
        f.write('%s %d %d %d %lf\n' % (tip,y,w,loc,vol))
#------------------------------ out_manip ------------------------------#


import sys
import os
import os.path
import rpy
import re

WORK_HOME = '~/Projects/Migration/one_feather/'
AWKSCRIPTS = WORK_HOME + 'RMigone/inst/scripts'
EXEC_DIR = WORK_HOME + 'build/'

file_types = {'.food':False,'.stb':False,'.uvb':False,'.acb':False,
              '.mob':False,'.ini':False}


usage = """usage:
%s filename [ext] [-N int] [-x start_year start_week end_year end_week loc vol]
   where -N: scaling of the simulated population
         -x: treatments,
             where x:
                   -p predation
                   -f food
                   -r reserves
                   -t feather quality
""" % sys.argv[0]                   

if len(sys.argv) > 1 and sys.argv[1][0] != '-':
    f_name = sys.argv[1];
    del sys.argv[1];
else:
    print usage
    sys.exit(1)
if len(sys.argv) >= 2 and sys.argv[1][0] != '-':
    ext_name = sys.argv[1]
    del sys.argv[1]
else:
    ext_name = 'sim'

manip_name = f_name + '-' + ext_name + '.man'
manip_f = open(manip_name, "w")

N = 1

# read variables from the command line, one by one:
while len(sys.argv) >= 8:
    option = sys.argv[1];       del sys.argv[1]
    if option == '-N':
        N = int(sys.argv[1]); del sys.argv[1]
    else:
        if   option == '-p': #predation
            tip = 'preda'
        if   option == '-f':
            tip = 'food'
        if   option == '-r':
            tip = 'reserve'
        if   option == '-t':
            tip = 'feather'
        s_y = int(sys.argv[1]); del sys.argv[1]
        s_w = int(sys.argv[1]); del sys.argv[1]
        e_y = int(sys.argv[1]); del sys.argv[1]
        e_w = int(sys.argv[1]); del sys.argv[1]
        loc  = int(sys.argv[1]); del sys.argv[1]
        vol  = float(sys.argv[1]); del sys.argv[1]
        out_manip(manip_f, s_y, s_w, e_y, e_w, tip, loc, vol)

manip_f.close()

for ft in file_types.keys():
    if os.path.exists(f_name + ft):
        continue
    elif os.path.exists(f_name + ft + '.gz'):
        os.system("gunzip " + f_name + ft + ".gz")
        file_types[ft] = True
    else:
        print "ERROR: necessary file '%s' does not exists." % (f_name + ft)
        sys.exit(1)

print 'Running MCrunner ...'
failure = os.system(EXEC_DIR + 'MCrunner -N ' +str(N) + ' ' + f_name \
                    + ' ' + ext_name)

for k,zipped in file_types.iteritems():
    if zipped:
        os.system("gzip " + f_name + k)

if failure:
    print 'ERROR: MCrunner failed'
    sys.exit(1)

print 'Done.'
print 'Plotting...'
y = []
count = []
f_data = open(f_name + '-' + ext_name + '.sum')
while True:
    line = f_data.readline()
    if not line:
        break
    if re.match(r'^[0-9]',line):
        l = map(int,line.split(','))
        y.append(l[0])
        count.append(l[1])

rpy.r.plot(y,count,type="l",xlab="years",ylab="pop. size")
rpy.r.abline(v=20)
#raw_input('Please press return to continue...\n')
rpy.r.library('RMigone',lib_loc="~/lib/R/library")
f_tup = (f_name,)*8
print f_tup
R_script="""
papd <- apd.ser()
x11()
layout(matrix(1:4,ncol=2,byrow=TRUE))
plot.apd("^p\\.w[0-3]$",papd,"%s",type="o",main="winter",ylab="proportion")
plot.apd("^p\\.s[0-3]$",papd,"%s",type="o",main="summer",ylab="proportion")
plot.apd("^n\\...$",papd,"%s",type="o",ylab="number of actions")
plot.apd("^p\\.br[0-3]$",papd,"%s",type="o",ylab="prop. of breeding")
layout(1)
x11()
layout(matrix(1:4,ncol=2,byrow=TRUE))
plot.apd("^p\\.mo[0-3]$",papd,"%s",type="o",ylab="prop. of moult")
plot.apd("(^w.br0|^w.mo[0-3])",papd,"%s",type="o",ylab="week")
plot.apd("^jo\\..*",papd,"%s",type="o",ylab="week")
plot.apd(".*\\.st$",papd,"%s",type="o",ylab="week")
layout(1)
""" % f_tup
print R_script
rpy.r(R_script)
raw_input('Please press return to continue...\n')
        
