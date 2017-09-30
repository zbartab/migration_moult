#!/usr/bin/env python

import sys, os.path, os
if sys.version.split()[0]<'2.3':
	True = 1
	False = 0

WORK_HOME = '~/Projects/Migration/one_feather'
AWKSCRIPTS = WORK_HOME + '/RMigone/inst/scripts'
EXEC_DIR = WORK_HOME + '/build'

N = 1
l_argv = len(sys.argv)

if l_argv == 1:
    print "ERROR: a file name should be given!"
    sys.exit(1)
elif l_argv > 1 and l_argv <= 3:
    f_name = sys.argv[1]
else:
    option = sys.argv[1]; del sys.argv[1]
    if option == "-N":
        N = int(sys.argv[1]); del sys.argv[1]
        f_name = sys.argv[1]
    else:
        f_name = option
        if sys.argv[1] == "-N":
            N = int(sys.argv[2])

os.system("gunzip " + f_name + "*gz")

i = 0
while True:
    i += 1
    s = EXEC_DIR + "/MCrunner -N " + str(N) +" "+ f_name
    failure = os.system(s)
    if failure:
        print "ERROR: generating sim.dat file is failed"
        sys.exit(1)
    else:
        s_name = f_name + "-sim"
        
    if (os.path.exists(s_name + ".dat")):
        dat_file = s_name + ".dat"
        CAT = "cat"
    elif (os.path.exists(s_name + ".dat.gz")):
        dat_file = s_name + ".dat.gz"
        CAT = "zcat"
    else:
        print "ERROR: ", s_name, " does not exists"
        sys.exit(1);
        
    s = CAT + " " + dat_file + " " + "|" + AWKSCRIPTS + \
        "/pre_process.awk >> " + s_name + ".pre"

    failure = os.system(s)
    if failure:
        print "ERROR: generating pre file is failed"
        sys.exit(1)
    else:
        k = os.popen("wc " + s_name + ".pre")
        wc_res = int(k.read().split()[0])
        k.close()
        if wc_res > 250 or i > 25:
            break
    
print "pre OK, no. of individuals:", wc_res

s = AWKSCRIPTS + "/pre_pre.py -y 0 " + s_name + ".pre > " + s_name + \
    ".apd"

failure = os.system(s)
if failure:
    print "ERROR: generating first apd file is failed"
    sys.exit(1)

print "apd OK"



#os.system("rm -f " + s_name + ".dat")
os.system("gzip -f " + f_name + "*")
