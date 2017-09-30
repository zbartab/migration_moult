#!/usr/bin/awk -f

#####
## Script to calculate the average time of starting breeding and moult
## and migration.
#####

function print_res(behav, n) {
  for(i in n) {
    j = i;
    gsub(SUBSEP, " ", j);
    if(n[i]>0) print behav, j, n[i];
  }
}

function med(n, m){
  max_w = max_e = max_l = 0;
  min_w = min_e = min_l = 10000000;
  for(i in n) {
    j = i;
    gsub(SUBSEP, " ", j);
    ns = split(j,A," ");
    if(A[1] < min_w) min_w = A[1];
    else if(A[1] > max_w) max_w = A[1];
    if(A[2] < min_e) min_e = A[2];
    else if(A[2] > max_e) max_e = A[2];
    if(A[3] < min_l) min_l = A[3];
    else if(A[3] > max_l) max_l = A[3];
  }
  for(e=min_e; e<=max_e; e++)
    for(l=min_l; l<=max_l; l++) {
      cum_sum[min_w] = n[min_w,e,l];
      s = n[min_w,e,l];
      for(w=min_w+1; w<=max_w; w++) {
	cum_sum[w] = cum_sum[w-1]+n[w,e,l];
	s += n[w,e,l];
      }
      for(w=0; w<=max_w; w++) {
	if(cum_sum[w] >= s/2.0) break;
      }
      m[e,l] = w-1.0 + ((s/2.0)-cum_sum[w-1])/(cum_sum[w]-cum_sum[w-1]);
    }
}

NF == 12 {
  i = $1*52+$2;
  j = $0;
  sub("^[0-9]+\t+[0-9]+","",j);
  $0 = i j;
}

$5 ~ /^S$/ {# start of breeding
  b_n[$1%52,$8,$10]++;
  
}

$5 ~ /^No/ {# migrate north
  no_n[$1%52,$8,$10]++;
}

$5 ~ /^So/ {# migrate south
  so_n[$1%52,$8,$10]++;
}

$5 ~ /B$/ {# search for territory
  te_n[$1%52,$8,$10]++;
}

$4 ~ /^S$/ {# start of moult
  m_n[$1%52,$8,$10]++;
}

END {
#   print_res("moult",m_n);
#   print_res("breed", b_n);
#   print_res("north", no_n);
#   print_res("south", so_n);
#   print_res("terri", te_n);
  m[0,0] = 0;
  med(b_n,m);
  print_res("breed",m);
  delete m;
  m[0,0] = 0;
  med(m_n,m);
  print_res("moult",m);
  delete m;
  m[0,0] = 0;
  med(te_n,m);
  print_res("terr",m);
  delete m;
  m[0,0] = 0;
  med(no_n,m);
  print_res("north",m);
  delete m;
  m[0,0] = 0;
  med(so_n,m);
  print_res("south",m);
}
