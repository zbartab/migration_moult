#!/usr/bin/awk -f

BEGIN {
  if(y=="") y = 3;
  max_loc = -1;
}

/^[0-9]/ {
  if($1 == y && $2 == 26) c26[$11,$3]++;
  else if($1 == y && $2 == 51) {
    c51[$11,$3]++;
    surv[$3]++;
  }
  if($11 > max_loc) max_loc = $11;
}

END {
  for(l = 0; l <= max_loc; l++) {
    for(bird in surv) {
      nyar[l] += c26[l,bird];
      tel[l] += c51[l,bird];
    }
  }
  for(l = 0; l <= max_loc; l++) print l, nyar[l], tel[l];
}
