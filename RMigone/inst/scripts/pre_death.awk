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


function do_job(act, from, n_act, w_act, s_act,eltol,     l, w)
{
  l = ex_loc(act,from);
  w = strtonum(ex_week(act,from));
  if(l==ML) eltol=26;
  else eltol=0;
  n_act[l]++;
  if(((eltol+w_act[l])%52) > ((eltol+w)%52)) {
    w_act[l] = w;
    s_act[l] = ex_state(act,from);
  }
}

function do_mig_job(act, from, n_act, w_act, s_act,eltol,     l, w)
{
  l = ex_loc(act,from);
  w = strtonum(ex_week(act,from));
  if(l==ML) eltol=26;
  else eltol=0;
  n_act[l]++;
  if(((eltol+w_act[l])%52) > ((eltol+w)%52)) {
    w_act[l] = w;
    s_act[l] = ex_mig_state(act,from);
  }
}

BEGIN {
  if(MOL == "") MOL = 10;
  if(ID == "") ID = "NO";
  else ID = "YES";
  if (ID == "YES") printf "ID";
  print "p.act", "p.week","p.loc","p.res","p.fea", "p.u","a.res","a.fea","p.ua","p.na","death";
}


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
  if(d ~ /[PRFH]$/) {
    nf = split($2,field,"|");
    n_prev = split(field[nf-1],prev,",");
    if(n_prev < 9) {
      hh = "";
      for(i = 1; i < 7; i++) hh = hh prev[i] " ";
      hh = hh " NA NA " prev[n_prev-1] " " prev[n_prev];
    } else {
      hh = field[nf-1];
    }
    hh = hh " " field[nf];
    gsub(/,/," ",hh);
    if (ID == "YES") print $1, hh;
    else print hh;
  }
}
