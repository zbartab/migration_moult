#!/usr/bin/awk -f

BEGIN {
  if (y == "") {
    y = 2;
  }
  n_brood = 0;	# mean number of breedings
  n_abort = 0;	# mean number of desertions
  n_aband = 0;	# mean number of abandons
  n_m1 = 0;	# mean number of moult 1
  l_b = 0;	# mean length of breeding
  l_m1 = 0;	# mean length of moult 1
  l_mbo = 0;	# mean length of moult-breeding overlap
  n_b = 0;	# number of birds breeding
  n_1 = 0;	# number of birds moulting 1
  n_m = 0;	# number of birds moulting 1 or 2
  m_o = 0;	# mean proportion of moulting overlap
  t_fb = 0;	# mean time of the start of the first breeding
  t_fm1 = 0;	# mean time of the start of the first moult 1
  t_b = 0;	# mean time of the start of the breeding
  t_m1 = 0;	# mean time of the start of the moult 1
  n_s = 0;	# number of survivals untill the end of year
}

{
  if ($1 == y) {
    if ($2 == 0 && $11 == 2) {
      gyak[$3]++;
    }
    if ($2 == 0) { 
      tfb[$3] = -1;
      tfm1[$3] = -1;
      lastbreed[$3] = -1;
      firstmoult[$3] = 1000;
    }
    if ($6 ~ "START") { 
      if (tfb[$3] == -1) tfb[$3] = $2;
      nbrood[$3]++;
      tb[$3] += $2;
    }
    if ($6 ~ "^AB") {
      if (lastbreed[$3] < $2) lastbreed[$3] = $2; 
      if ($6 ~ "ABORT") { 
	nabort[$3]++;
      }
      if ($6 ~ "ABAND") { 
	naband[$3]++;
      }
    }
    if ($5 ~ "START") { 
      if (tfm1[$3] == -1) tfm1[$3] = $2;
      nm1[$3]++;
      tm1[$3] += $2;
      if (firstmoult[$3] > $2) firstmoult[$3] = $2;
    }
    if ($6 ~ "START" || $6 ~ "KEEP") lb[$3]++;
    if ($5 !~ "NO") lm1[$3]++;
    if (($5 !~ "NO") && ($6 ~ "START" || $6 ~ "KEEP")) lmbo[$3]++;
    if ($2 == 51 && gyak[$3] == 1) {
      surv[$3]++;
    }
  } else if ($1>y && FNR > 1) exit;
}

END {
  for (bird in surv) {
    n_brood += nbrood[bird];
    n_abort += nabort[bird];
    n_aband += naband[bird];
    if (nbrood[bird] > 0) {
      n_b++;
      tb[bird] /= nbrood[bird];
    }
    n_m1 += nm1[bird];
    if (nm1[bird] > 0) {
      n_1++;
      tm1[bird] /= nm1[bird];
    }
    if (nm1[bird] > 0) n_m++;
    if (nm1[bird] > 0) m_o += (lmo[bird]/lm[bird]);
    if (lastbreed[bird] > -1 && firstmoult[bird] < 1000) {
      n_mb++;
      l_mb += firstmoult[bird] - lastbreed[bird];
    }
    l_b += lb[bird];
    l_m1 += lm1[bird];
    l_mbo += lmbo[bird];
    t_b += tb[bird];
    t_m1 += tm1[bird];
    if (tfb[bird] >= 0) t_fb += tfb[bird];
    if (tfm1[bird] >= 0) t_fm1 += tfm1[bird];
    n_s++;
  }
  print n_brood/n_s, n_aband/n_s, n_abort/n_s, n_m1/n_s, l_b/n_b, \
    l_m1/n_1, l_mbo/n_b, t_b/n_b, \
    t_fb/n_b, t_m1/n_1, t_fm1/n_1, l_mb/n_mb, n_s;
}
