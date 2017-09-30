#!/usr/bin/awk -f

BEGIN {
  max_year = 0;
  max_loc = 0;
  max_ex = 0;
}

/^[0-9]/ {
  if($1>max_year) max_year = $1;
  if($9>max_ex) max_ex = $9;
  if($11>max_loc) max_loc = $11;
  if($2 == 51) census[$1,$9,$11]++;
  if($4 ~ "P") pred[$1,$9,$11]++;
  if($4 ~ "S") starv[$1,$9,$11]++;
}

END {
  for(l = 0; l <= max_loc; l++) {
    for(y = 0; y <= max_year; y++) {
      for(e = 0; e <= max_ex; e++) {
	printf "\t%d\t%d\t%d", census[y,e,l], pred[y,e,l], starv[y,e,l];
      }
      printf "\n";
    }
    printf "\n";
  }
}
