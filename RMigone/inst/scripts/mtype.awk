#!/usr/bin/awk -f

BEGIN {
  if (y == "") {
    y = 2;
  }
}

{
  if ($1 == y) {
    if ($2 == 0 && $11 == 2) {
      gyak[$3]++;
    }
    if ($6 ~ "START") {
      sched[$3] = sprintf("%s%s", sched[$3], "B");
    } else if ($6 ~ "ABAND") {
      sched[$3] = sprintf("%s%s", sched[$3], "a");
    } else if ($6 ~ "ABORT") {
      sched[$3] = sprintf("%s%s", sched[$3], "d");
    }
    if ($5 ~ "START") sched[$3] = sprintf("%s%s", sched[$3], "p");
    if ($2 == 51 && gyak[$3] == 1) surv[$3]++;
  } else if ($1>y && FNR > 1) exit;
}

END {
  for (bird in surv) {
    s_type[sched[bird]]++;
  }
  for (sched_t in s_type) {
    if(sched_t == "") printf "%s\t%d\n", "n", s_type[sched_t];
    else printf "%s\t%d\n", sched_t, s_type[sched_t];
  }
}
