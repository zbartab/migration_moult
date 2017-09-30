#! /usr/bin/gawk -f

#########
# migra1 `.dat' file preprocessor
#########


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
  if ($2 == 0) {
    sched[$3] = sprintf("%s|(%s,%d,%d,%d,%d,%d,%f)",sched[$3],"W",\
			$2,$11,$7,$8,$9,$12);
  }
  if ($2 == 26) {
    sched[$3] = sprintf("%s(%s,%d,%d,%d,%d,%d,%f)",sched[$3],"U",$2,$11,$7,$8,$9,$12);
  }
  if ($4 ~ /^S$/) {
    sched[$3] = sprintf("%s(%s,%d,%d,%d,%d,%d,%f)",sched[$3],"E",$2,$11,$7,$8,$9,$12);
  } else if($4 ~ /^P$/) {
    sched[$3] = sprintf("%s(%s,%d,%d,%d,%d,%d,%f)",sched[$3],"P",$2,$11,$7,$8,$9,$12);
  } else {
    if ($6 ~ /^S$/) {
      sched[$3] = sprintf("%s(%s,%d,%d,%d,%d,%d,%f)",sched[$3],"B",$2,$11,$7,$8,$9,$12);
    } else if ($6 ~ /^A$/) {
      sched[$3] = sprintf("%s(%s,%d,%d,%d,%d,%d,%f)",sched[$3],"a",$2,$11,$7,$8,$9,$12);
    } else if ($6 ~ /B$/) {
      sched[$3] = sprintf("%s(%s,%d,%d,%d,%d,%d,%f)",sched[$3],"T",$2,$11,$7,$8,$9,$12);
    } else if ($6 ~ /^AT$/) {
      sched[$3] = sprintf("%s(%s,%d,%d,%d,%d,%d,%f)",sched[$3],"d",$2,$11,$7,$8,$9,$12);
    } else if($6 ~ /^No./) {
      sched[$3] = sprintf("%s(%s,%d,%d,%d,%d,%d,%f)",sched[$3],"N",$2,$11,$7,$8,$9,$12);
    } else if($6 ~ /^So./) {
      sched[$3] = sprintf("%s(%s,%d,%d,%d,%d,%d,%f)",sched[$3],"S",$2,$11,$7,$8,$9,$12);
    }
    if ($5 ~ /^S$/) {
      sched[$3] = sprintf("%s(%s,%d,%d,%d,%d,%d,%f)",sched[$3],"M",$2,$11,$7,$8,$9,$12);
    } else if($5 ~ /^E$/) {
      sched[$3] = sprintf("%s(%s,%d,%d,%d,%d,%d,%f)",sched[$3],"e",$2,$11,$7,$8,$9,$12);
    }
  }
}


END {
  for (bird in sched) {
    if(sched[bird] == "") printf "%s\n", "n";
    else printf "%s\n", sched[bird];
  }
}
