#!/usr/bin/awk -f

#####
## Script to count the individuals of different exp classes and at
## different location. 
#####


/^ *[0-9]/ {
  c[$1,$8,$10]++;
}

END {
  for(i in c) {
    j = i;
    gsub(SUBSEP, " ", j);
    print j, c[i];
  }
}
