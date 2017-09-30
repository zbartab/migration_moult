#! /usr/bin/gawk -f

#########
# process the results of the `pre_process.awk' the migra1 `.dat' file preprocessor
#########

function ex_week(act, from,    re) 
{
  re = "^.*" act ",([0-9]+),.*$";
  return gensub(re,"\\1","g",from);
}

function ex_loc(act, from,    re) 
{
  re = "^.*" act ",[0-9]+,([0-9]+),.*$";
  return gensub(re,"\\1","g",from);
}

function ex_res(act, from,    re) 
{
  re = "^.*" act ",[0-9]+,[0-9]+,([0-9]+),.*$";
  return gensub(re,"\\1","g",from);
}

function ex_fea(act, from,    re) 
{
  re = "^.*" act ",[0-9]+,[0-9]+,[0-9]+,([0-9]+),.*$";
  return gensub(re,"\\1","g",from);
}

function ex_u(act, from,    re) 
{
  re = "^.*" act ",[0-9]+,[0-9]+,[0-9]+,[0-9]+,([0-9.]+).*$";
  return gensub(re,"\\1","g",from);
}

function ex_state(act, from,    re)
{
  re = "^.*" act ",[0-9]+,[0-9]+,([0-9]+),([0-9]+),([0-9.]+),([0-9.]+|NA),([0-9]+).*$";
  return gensub(re,"\\1 \\2 \\3 \\4 \\5 ","g",from);
}

function ex_mig_state(act, from,    re)
{
  re = "^.*" act ",[0-9]+,[0-9]+,([0-9]+),([0-9]+),([0-9.]+),([0-9]+),([0-9]+),([0-9.]+|NA),([0-9]+).*$";
  return gensub(re,"\\1 \\2 \\3 \\4 \\5 \\6 \\7 ","g",from);
}

function delarray(arr,    i) 
{
  for(i in arr) delete arr[i];
}

function comp_w(w1, w2) # return true if w2 is earlier than w1
{
  d = w1 - w2;
  d = sqrt(d*d);
  if(d < 26) {
    return w1 > w2
  } else {
    return w2 > w1
  }
}


function do_job(act, from, n_act, w_act, s_act,     l, w)
{
  l = ex_loc(act,from);
  w = strtonum(ex_week(act,from));
  n_act[l]++;
  if (w_act[l] == 100 || comp_w(w_act[l],w)) {
    w_act[l] = w;
    s_act[l] = ex_state(act,from);
  }
}

function do_mig_job(act, from, n_act, w_act, s_act,     l, w)
{
  l = ex_loc(act,from);
  w = strtonum(ex_week(act,from));
  n_act[l]++;
#  if(((eltol+w_act[l])%52) > ((eltol+w)%52)) {
  if (w_act[l] == 100 || comp_w(w_act[l],w)) {
    w_act[l] = w;
    s_act[l] = ex_mig_state(act,from);
  }
}

BEGIN {
  if(MOL == "") MOL = 10;
  if(ML == "") ML = 4;
  printf "ID mo.over ";
  for(l=0; l < ML; l++) {
    printf  " n.no" l " n.te" l " n.br" l " n_su" l " n.mo" l " n.en" l " n.so" l\
      " w.no" l " r.no" l " f.no" l " u.no" l " ar.no" l " af.no" l " ua.no" l " na.no" l \
      " w.te" l " r.te" l " f.te" l " u.te" l " ua.te" l " na.te" l \
      " w.br" l " r.br" l " f.br" l " u.br" l " ua.br" l " na.br" l \
      " w.su" l " r.su" l " f.su" l " u.su" l " ua.su" l " na.su" l \
      " w.mo" l " r.mo" l " f.mo" l " u.mo" l " ua.mo" l " na.mo" l \
      " w.en" l " r.en" l " f.en" l " u.en" l " ua.en" l " na.en" l \
      " w.so" l " r.so" l " f.so" l " u.so" l " ar.so" l " af.so" l " ua.so" l " na.so" l ;
  }
  printf " w.loc " "w.res " "w.fea " "w.u " "s.loc " "s.res " "s.fea " "s.u " \
    "l.loc " "l.res " "l.fea " "l.u " "d.act " "d.w " "d.loc "\
    "d.res " "d.fea " "d.u " "d.ua " "d.na ";
  printf "\n";
}


# if you want to process a specific year add th `-v ev=x' switch to
# the command line, where `x' is the year you want

NR == 1 && ev == "" {
  e = $2;
  gsub(/[^W]/,"",e);
  if(length(e) > 1) ev = 2;
  else ev = 1;
}

{
  d = $2;
  gsub(/NA/,"",d);
  gsub(/[^a-zA-Z]/,"",d);
  if(d ~ /[PRFHL]$/) {
    for(l=0; l < ML; l++) {
      n_south[l] = 0;
      n_north[l] = 0;
      n_terr[l] = 0;
      n_breed[l] = 0;
      n_succ[l] = 0;
      n_moult[l] = 0;
      n_end[l] = 0;
      w_terr[l] = 100;
      s_terr[l] = "NA NA NA NA NA ";
      w_north[l] = 100;
      s_north[l] = "NA NA NA NA NA NA NA ";
      w_south[l] = 100;
      s_south[l] = "NA NA NA NA NA NA NA ";
      w_breed[l] = 100;
      s_breed[l] = "NA NA NA NA NA ";
      w_succ[l] = 100;
      s_succ[l] = "NA NA NA NA NA ";
      w_moult[l] = 100;
      s_moult[l] = "NA NA NA NA NA ";
      w_end[l] = 100;
      s_end[l] = "NA NA NA NA NA ";
    }
    nyar = $2;
    sub(/^\|/,"",nyar);
    gsub(/\|W/,"*W",nyar);
    split(nyar,nyev,"*");
    if(ev in nyev) {
      n_nya_d = split(nyev[ev],nya_d,"|");
      for(i in nya_d) {
	if(nya_d[i] ~ /^T/) {
	  do_job("T",nya_d[i],n_terr,w_terr,s_terr);
	}else if(nya_d[i] ~ /^B/) {
	  do_job("B",nya_d[i],n_breed,w_breed,s_breed);
	}else if(nya_d[i] ~ /^a/) {
	  do_job("a",nya_d[i],n_succ,w_succ,s_succ);
	}else if(nya_d[i] ~ /^N/) {
	  do_mig_job("N",nya_d[i],n_north,w_north,s_north);
	}else if(nya_d[i] ~ /^S/) {
	  do_mig_job("S",nya_d[i],n_south,w_south,s_south);
	}else if(nya_d[i] ~ /^M/) {
	  do_job("M",nya_d[i],n_moult,w_moult,s_moult);
	}else if(nya_d[i] ~ /^e/) {
	  do_job("e",nya_d[i],n_end,w_end,s_end);
	}
      }
      
      ny2 = nyev[ev];
      gsub(/[^MeTBadNS]/,"",ny2);
      m_over = "";
      if(ny2 ~ /MT+e/) m_over = sprintf("%s%s",m_over, "T");
      if(ny2 ~ /M[Bad]+e/) m_over = sprintf("%s%s",m_over, "B");
      if(ny2 ~ /M[NS]+e/) m_over = sprintf("%s%s",m_over, "M");
      if(m_over == "") m_over = "N";
      
      printf $1 " " m_over " ";
      for(l=0; l < ML; l++) {
	if(w_north[l] == 100) w_north[l] = "NA";
	if(w_south[l] == 100) w_south[l] = "NA";
	if(w_terr[l] == 100) w_terr[l] = "NA";
	if(w_breed[l] == 100) w_breed[l] = "NA";
	if(w_succ[l] == 100) w_succ[l] = "NA";
	if(w_moult[l] == 100) w_moult[l] = "NA";
	if(w_end[l] == 100) w_end[l] = "NA";
	printf " " n_north[l] "  " n_terr[l] "  " n_breed[l] "  " n_succ[l] "  " \
	  n_moult[l] "  " n_end[l] "  " n_south[l] "  " w_north[l] "  " s_north[l] "  " w_terr[l] "  " s_terr[l] "  " \
	  w_breed[l] "  " s_breed[l] "  " w_succ[l] "  "\
	  s_succ[l] "  " w_moult[l] " " s_moult[l] "  " w_end[l] "  " s_end[l] "  " w_south[l] "  " s_south[l];
      }
      printf " " ex_loc("W",nyev[ev]) " " ex_res("W",nyev[ev]) " " ex_fea("W",nyev[ev]) "  " ex_u("W",nyev[ev]);
      if (nyev[ev] ~ /U/) printf " " ex_loc("U",nyev[ev]) " " ex_res("U",nyev[ev]) " " ex_fea("U",nyev[ev]) " " ex_u("U",nyev[ev]) " ";
      else printf " NA NA NA NA ";
      if (nyev[ev] ~ /L/) printf ex_loc("L",nyev[ev]) " " ex_res("L",nyev[ev]) " " ex_fea("L",nyev[ev]) " " ex_u("L",nyev[ev]) " ";
      else printf " NA NA NA NA ";
      if(nyev[ev] ~ /[PRFH]/) {
	gsub(/,/," ",nya_d[n_nya_d]);
	printf " " nya_d[n_nya_d];
      } else {
	printf " A NA NA NA NA NA NA NA ";
      }
      printf "\n";
    }
  }
}
