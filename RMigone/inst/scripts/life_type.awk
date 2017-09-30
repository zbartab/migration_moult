#!/usr/bin/awk -f

{
  if ($7 ~ "START") {
    sched[$3] = sprintf("%s%s", sched[$3], "B");
  }
  if ($5 ~ "START") sched[$3] = sprintf("%s%s", sched[$3], "p");
  if ($6 ~ "START") sched[$3] = sprintf("%s%s", sched[$3], "s");
  if ($2 == 51) sched[$3] = sprintf("%s%s", sched[$3], "|");
}

END {
  for (bird in sched) {
    s_type[sched[bird]]++;
  }
  for (sched_t in s_type) {
    if(sched_t == "") printf "%s\t%d\n", "n", s_type[sched_t];
    else printf "%s\t%d\n", sched_t, s_type[sched_t];
  }
}
