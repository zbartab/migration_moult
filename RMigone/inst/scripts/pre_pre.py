#!/usr/bin/env python

##############################
#
# .pre file processor
#
##############################


# non migratory activity seems to be OK
# migratory activity: record errors (what is an error?)

import sys,re,string

if sys.version.split()[0]<'2.3':
	True = 1
	False = 0

#------------------------- comp_week -------------------------#
def comp_week(w1,w2): 
	"""Return true if w1 is earlier than w2 taking 
	acount the annual circularity
	"""
	d = abs(w1-w2)
	if (d < 26):
		return w1 < w2
	else:
		return w1 > w2
#------------------------- comp_week -------------------------#

mig_fontos = re.compile(r'[WLMeBadNSFRPH]')

#------------------------- valid_mig -------------------------#
def valid_mig(CI,li,acts):
	i = CI
	while True:
		if i > 0:
			i -= 1
		else:
			i = 0
		code_p, w, l_p, state = acts[i].split(",",3)
		if (mig_fontos.match(code_p) or i == 0):
			break
	i = CI
	while True:
		if i < li-1:
			i += 1
		else:
			i = li-1
		code_n, w, l_n, state = acts[i].split(",",3)
		if (mig_fontos.match(code_n) or i == li-1):
			break
	return not (l_n == l_p)
#------------------------- valid_mig -------------------------#



local_code = re.compile(r'[MeBaT]')
mig_code = re.compile(r'[NS]')
death_code = re.compile(r'[FHRP]')
pos_code = re.compile(r'[WUL]')
#------------------------- proc_acts -------------------------#
def proc_acts(acts):
	r = [{} for i in range(ML)]
	for i in range(len(acts)):
		code, w, l, state = acts[i].split(",",3)
		w, l = map(int, [w,l])
		if local_code.match(code):      # local actions e.g. moult
			if r[l].has_key(code):	# action already occured
				r[l][code][0] += 1
				if comp_week(w, r[l][code][1]):
					r[l][code][1:] = [w, state]
			else:			# action new
				r[l][code] = [1, w, state]
		elif pos_code.match(code):      # position info e.g. midwinter
			if r[l].has_key(code):	# action already occured
				r[l][code][0] += 1
				if comp_week(w, r[l][code][1]):
					r[l][code][1:] = [w, string.join([str(l),state],",")]
			else:			# action new
				r[l][code] = [1, w, string.join([str(l),state],",")]
		elif death_code.match(code):    # cause of death
			if r[l].has_key('D'):	# action already occured
				r[l]['D'][0] += 1
				if comp_week(w, r[l]['D'][1]):
					r[l]['D'][1:] = [w, string.join([code,str(l),state],",")]
			else:			# action new
				r[l]['D'] = [1, w, string.join([code,str(l),state],",")]
		elif mig_code.match(code):	#migration
			if valid_mig(i,len(acts),acts):
				if r[l].has_key(code):	# action already occured
					r[l][code][0] += 1
					if comp_week(w, r[l][code][1]):
						r[l][code][1:] = [w, state]
				else:			# action new
					r[l][code] = [1, w, state]
	return r
#------------------------- proc_acts -------------------------#

#------------------------- var_names -------------------------#
def var_names(res,vn):
	IDs = res.keys()
	for i in IDs:
		for j in range(len(res[i])):
			k = res[i][j].keys()
			for s in k:
				if local_code.match(s) or mig_code.match(s):
					vn[string.join([s,str(j)],"")] = []
				else:
					vn[s] = []
#------------------------- var_names -------------------------#

#------------------------- fill_mat  -------------------------#
def fill_mat(res,vn,m_o):
	var_n = vn.keys()
	var_n.sort()
	IDs = res.keys()
	for i in IDs:
		vn['ID'].append(i)
		vn['mo.over'].append(m_o[i])
		for j in range(len(res[i])):
			res_keys = res[i][j].keys()
			for k in res_keys:
				if local_code.match(k) or mig_code.match(k):
					a = string.join([k,str(j)],"")
				else:
					a = k
				vn[a].append(res[i][j][k])
		for k in var_n:
			if len(vn[k]) < len(vn['ID']):
				vn[k].append(None)
#------------------------- fill_mat -------------------------#

#------------------------- print_header -------------------------#
def print_header(v_n):
	helyi_valt = ['n','w','r','f','u','ua','na']
	mig_valt = ['n','w','r','f','u','ar','af','ua','na']
	pos_valt = ['n','w','loc','res','fea','u','ua','na']
	death_valt = ['n','w','act','loc','res','fea','u','ua','na']
	v_recode = {'M':'mo', 'e':'en', 'B':'br', 'a':'su', 'T':'te',
		    'N':'no','S':'so'}
	for k in v_n:
		a = k[0]
		if local_code.match(a):
			print string.join([string.join([h,'.',v_recode[a],k[1]],\
						       "") for h in helyi_valt]),
		elif mig_code.match(a):
			print string.join([string.join([h,'.',v_recode[a],k[1]],\
						       "") for h in mig_valt]),
		elif re.match(r'ID',k):
			print 'ID',
		elif re.match(r'mo.over',k):
			print 'mo.over',
		elif pos_code.match(a):
			print string.join([string.join([a.lower(),'.',\
							h],"") \
					   for h in pos_valt]),
		elif a == 'D':
			print string.join([string.join([a.lower(),'.',\
							h],"") \
					   for h in death_valt]),
	pass
#------------------------- print_header -------------------------#


#------------------------- print_rec -------------------------#
def print_rec(k,rek):
#	print k, rek
	if rek == None:
		if local_code.match(k) :
			print 0, 'NA '*6,
		elif mig_code.match(k):
			print 0, 'NA '*8,
		elif pos_code.match(k) or death_code.match(k):
			print 0, 'NA '*7,
		elif k == 'D':
			print 0, 'NA A ', 'NA '*6,
	elif k == 'ID' or k == 'mo.over':
		print rek,
	else:
		print rek[0], rek[1], re.sub(r',',' ',rek[2]),
	pass
#------------------------- print_rec -------------------------#

#------------------------- print_res -------------------------#
def print_res(vn):
	var_n = vn.keys()
	var_n.sort()
	print_header(var_n)
	print 
	for i in range(len(vn['ID'])):
		for n in var_n:
			print_rec(n, vn[n][i]);
		print
	pass
#------------------------- print_res -------------------------#


#------------------------- main -------------------------#
ML = 4
ev = 0
f_name = ''
if(len(sys.argv) > 1):
	if (sys.argv[1] == '-y'):
		ev = int(sys.argv[2])
		try:
			f_name = sys.argv[3]
		except:
			f_name = ''
	else:
		f_name = sys.argv[1]
		ev = 0
	
try:
	f_pre = open(f_name,"r")
except:
	f_pre = sys.stdin
#	print "ERROR: " + f_pre + " cannot be opened!"
#	sys.exit(1)

y_split1 = re.compile(r'\|W')
y_split2 = re.compile(r'^\*')

mo_NA = re.compile(r'NA')
mo_ex = re.compile(r'[^MeTBadNS]')
mo_T = re.compile(r'MT+e')
mo_B = re.compile(r'M[Bad]+e')
mo_M = re.compile(r'M[NS]+e')
mo_over = {}

result = {}

while True:
	indiv = f_pre.readline()
	if not indiv:
		break
	ID, behav = indiv.split(" ",1)
	behav = y_split1.sub("*W",behav)
	behav = y_split2.sub("",behav)
	years = behav.split("*")
	if(len(years) > ev):
		acts = years[ev].strip().split("|");
		result[ID] = proc_acts(acts)
		ny = mo_NA.sub("",years[ev])
		ny = mo_ex.sub("",ny)
		mo_str = ''
		if mo_T.search(ny):
			mo_str += 'T'
		if mo_B.search(ny):
			mo_str += 'B'
		if mo_M.search(ny):
			mo_str += 'M'
		if mo_str == '':
			mo_str = 'N'
		mo_over[ID] = mo_str
f_pre.close()

vn = {}
vn['ID'] = []
vn['mo.over'] = []
var_names(result,vn)
fill_mat(result,vn,mo_over)

print_res(vn)
#------------------------- main -------------------------#
