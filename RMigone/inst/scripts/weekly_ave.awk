#!/usr/bin/awk -f

#####
## Script to count the individuals of different exp classes and at
## different location. 
#####

BEGIN {
  if (y == "") {
    y = -1;
  }
  s[0,0,0,0] = 0;
  p[0,0,0,0] = 0;
  sm[0,0,0,0] = 0;
  pm[0,0,0,0] = 0;
}


$1 == y || y == -1 {
  r[$1,$2,$9,$11] += $7;
  f[$1,$2,$9,$11] += $8;
  u[$1,$2,$9,$11] += $12;
  if($4 ~ /^S$/) s[$1,$2,$9,$11]++;
  if($4 ~ /^P$/) p[$1,$2,$9,$11]++;
  if($4 ~ /^SM$/) sm[$1,$2,$9,$11]++;
  if($4 ~ /^PM$/) pm[$1,$2,$9,$11]++;
  c[$1,$2,$9,$11]++;
}

END {
  for(i in c) {
    j = i;
    gsub(SUBSEP, " ", j);
    if(s[i] == "") s[i] = 0;
    if(p[i] == "") p[i] = 0;
    if(sm[i] == "") sm[i] = 0;
    if(pm[i] == "") pm[i] = 0;
    print j, r[i]/c[i], f[i]/c[i], u[i]/c[i], c[i], s[i], p[i], sm[i], pm[i];
  }
}
