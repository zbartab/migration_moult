#!/usr/bin/awk -f

#####
## Script to extract information about the timing of events from dat files
#####

function print_res(behav, n) {
  for(i in n) {
    j = i;
    gsub(SUBSEP, " ", j);
    if(n[i]>0) print behav, j, n[i];
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
  print_res("moult",m_n);
  print_res("breed", b_n);
  print_res("north", no_n);
  print_res("south", so_n);
  print_res("terri", te_n);
}
