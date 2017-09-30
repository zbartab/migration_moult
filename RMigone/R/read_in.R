
O.sim <- function(nev="baseline",y=3, ext=".dat") {
  if (is.null(ext)) { f.nev=nev;}
  else {f.nev <- paste(nev,ext,sep="");}
  zz1 <- pipe(paste("sed -n -e '1p' -e '1q' ",f.nev,sep=""));
  n.ds <- scan(zz1,what=character(0),quiet=TRUE);
  close(zz1);
  zz2 <- pipe(paste("awk '{if ($1==",y,") print $0}' ",f.nev, sep=""));
  if(length(n.ds)==12) {
    ds <- scan(zz2,what=list(year=integer(0),week=integer(0),
                     bird=character(0),
                     event=character(0),moult1=character(0),
                     action=character(0),res=integer(0),fea1=integer(0),
                     exper=integer(0),brood=integer(0),loc=integer(0),
                     uval=double(0)));
    close(zz2);
  } else {
    ds <- scan(zz2,what=list(year=integer(0),week=integer(0),
                     bird=character(0),
                     event=character(0),moult1=character(0),
                     action=character(0),res=integer(0),fea1=integer(0),
                     exper=integer(0),brood=integer(0),loc=integer(0),
                     uval=double(0),a.res=integer(0),a.fea1=integer(0),
                     a.loc=integer(0)));
  }
}
