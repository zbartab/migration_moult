#! /usr/bin/awk -f


#######
## preprocess migra1's `txt' file
#######

BEGIN {
  last = "";
  lambda = -1;
  sz = -1;
  diff = -1;
}

/^[bf][ao]/ {
  gsub(/=/," ");
  gsub(/,/,"");
  if (last != "" && last !~ $1) print substr(last,1,1), sz, lambda, diff;
  lambda = $6;
  sz = $4;
  diff = $8;
  last = $1;
}

END {
  if(diff != "") print substr(last,1,1), sz, lambda, diff;
}
