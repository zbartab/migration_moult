#! /usr/bin/gawk -f

##################################
#
#  migra1 `.dat' file preprocessor
#
##################################

#########################
# Codes used for output
# 
# B: start breeding
# M: start moult
# N: migrate north
# S: migrate south
# T: search for territory
# U: midsummer
# W: begining of the year (midwinter, week == 0)
# L: end  of the year (midwinter, week == 51)
# a: end of successful breeding (abandon)
# d: end of unsuccessful breeding (desert)
# e: end of moult
# P: died because of predation
# R: died during migration because of predation
# F: died because of starvation
# H: died during migration because of starvation
#
#########################



BEGIN {
  if (y == "") {
    y = -1;
    year = -1;
    n_year = 0;
  }
}


$1 == y || y == -1 {
  if (year != $1) {
    year = $1;
    n_year++;
  }
  u[$3] = u[$3] + $12;
  n[$3] = n[$3] + 1;
  if ($2 == 0 && $9 == 2 && n_year == 1) {
    gyak[$3]++;
  }
  if($3 in gyak) {
    if ($2 == 0) {
      if (n[$3] > 0) ua = u[$3]/n[$3];
      else ua = "NA";
      sched[$3] = sched[$3] "|W" "," $2 "," $11 "," $7 "," $8 "," $12 "," ua "," n[$3];
      u[$3] = n[$3] = 0;
    } else  if ($2 == 26) {
      if (n[$3] > 0) ua = u[$3]/n[$3];
      else ua = "NA";
      sched[$3] = sched[$3] "|U" "," $2 "," $11 "," $7 "," $8 "," $12 "," ua "," n[$3];
      u[$3] = n[$3] = 0;
    } 
    
    if ($4 ~ /^P$/) {
      if (n[$3] > 0) ua = u[$3]/n[$3];
      else ua = "NA";
      sched[$3] = sched[$3] "|P" "," $2 "," $11 "," $7 "," $8 "," $12 "," ua "," n[$3];
      u[$3] = n[$3] = 0;
    } else  if ($4 ~ /^PM$/) {
      if (n[$3] > 0) ua = u[$3]/n[$3];
      else ua = "NA";
      sched[$3] = sched[$3] "|R" "," $2 "," $11 "," $7 "," $8 "," $12 "," ua "," n[$3];
      u[$3] = n[$3] = 0;
    } else  if ($4 ~ /^S$/) {
      if (n[$3] > 0) ua = u[$3]/n[$3];
      else ua = "NA";
      sched[$3] = sched[$3] "|F" "," $2 "," $11 "," $7 "," $8 "," $12 "," ua "," n[$3];
      u[$3] = n[$3] = 0;
    } else  if ($4 ~ /^SM$/) {
      if (n[$3] > 0) ua = u[$3]/n[$3];
      else ua = "NA";
      sched[$3] = sched[$3] "|H" "," $2 "," $11 "," $7 "," $8 "," $12 "," ua "," n[$3];
      u[$3] = n[$3] = 0;
    } else {  ### process further if the bird is alive
    
      if ($5 ~ /^S$/) {
	if (n[$3] > 0) ua = u[$3]/n[$3];
	else ua = "NA";
	sched[$3] = sched[$3] "|M" "," $2 "," $11 "," $7 "," $8 "," $12 "," ua "," n[$3];
	u[$3] = n[$3] = 0;
      } else if($5 ~ /^E$/) {
	if (n[$3] > 0) ua = u[$3]/n[$3];
	else ua = "NA";
	sched[$3] = sched[$3] "|e" "," $2 "," $11 "," $7 "," $8 "," $12 "," ua "," n[$3];
	u[$3] = n[$3] = 0;
      }
      if($6 ~ /^No./) {
	if (n[$3] > 0) ua = u[$3]/n[$3];
	else ua = "NA";
	sched[$3] = sched[$3] "|N" "," $2 "," $11 "," $7 "," $8 "," $12 "," $13 "," $14 "," ua "," n[$3];
	u[$3] = n[$3] = 0;
      } else if($6 ~ /^So./) {
	if (n[$3] > 0) ua = u[$3]/n[$3];
	else ua = "NA";
	sched[$3] = sched[$3] "|S" "," $2 "," $11 "," $7 "," $8 "," $12 "," $13 "," $14 "," ua "," n[$3];
	u[$3] = n[$3] = 0;
      }
      if ($6 ~ /^S$/) {
	if (n[$3] > 0) ua = u[$3]/n[$3];
	else ua = "NA";
	sched[$3] = sched[$3] "|B" "," $2 "," $11 "," $7 "," $8 "," $12 "," ua "," n[$3];
	u[$3] = n[$3] = 0;
      } else if ($6 ~ /^A$/) {
	if (n[$3] > 0) ua = u[$3]/n[$3];
	else ua = "NA";
	sched[$3] = sched[$3] "|a" "," $2 "," $11 "," $7 "," $8 "," $12 "," ua "," n[$3];
	u[$3] = n[$3] = 0;
      } else if ($6 ~ /^B$/) {
	if (n[$3] > 0) ua = u[$3]/n[$3];
	else ua = "NA";
	sched[$3] = sched[$3] "|T" "," $2 "," $11 "," $7 "," $8 "," $12 "," ua "," n[$3];
	u[$3] = n[$3] = 0;
      } else if ($6 ~ /oB$/) {
	if (n[$3] > 0) ua = u[$3]/n[$3];
	else ua = "NA";
	sched[$3] = sched[$3] "|T" "," $2 "," $15 "," $7 "," $8 "," $12 "," ua "," n[$3];
	u[$3] = n[$3] = 0;
      } else if ($6 ~ /^AT$/) {
	if (n[$3] > 0) ua = u[$3]/n[$3];
	else ua = "NA";
	sched[$3] = sched[$3] "|d" "," $2 "," $11 "," $7 "," $8 "," $12 "," ua "," n[$3];
	u[$3] = n[$3] = 0;
      } 
      if ($2 == 51) {
	if (n[$3] > 0) ua = u[$3]/n[$3];
	else ua = "NA";
	sched[$3] = sched[$3] "|L" "," $2 "," $11 "," $7 "," $8 "," $12 "," ua "," n[$3];
	u[$3] = n[$3] = 0;
      }
    }
  }
#  if ($2 == 51 && gyak[$3] == n_year) surv[$3]++;
}


END {
  for (bird in gyak) {
    printf "%s %s\n", bird, sched[bird];
  }
}
