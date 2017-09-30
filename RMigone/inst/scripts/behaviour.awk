#!/usr/bin/awk -f

#####
## Script to calculate the average time of starting breeding and moult
## and migration.
#####

BEGIN {
  pi = -1;
  b_pn = 0;
  m_pn = 0;
}

$1 != pi {
  pi = $1;
  if(b_n == 0) {
    if (b_pn > 0) {
      print "breed", $1, b_w/b_pn, b_pn;
    }
    b_pn = b_w = 0;
  }
  b_pn += b_n;
  b_n = 0;
  if(m_n == 0) {
    if (m_pn > 0) {
      print "moult", $1, m_w/m_pn, m_pn;
    }
    m_pn = m_w = 0;
  }
  m_pn += m_n;
  m_n = 0;
}

$5 ~ /^S$/ && $8 == 2 {
#  print "Breed: ", $2;
  b_w += $1;
  b_n++;
}

$4 ~ /^S$/ && $8 == 2 {
  m_w += $1;
  m_n++;
}

END {
  b_pn += b_n;
  m_pn += m_n;
  if(b_pn>0) print "breed", $1, b_w/b_pn, b_pn;
  if(m_pn>0) print "moult", $1, m_w/m_pn, m_pn;
}
