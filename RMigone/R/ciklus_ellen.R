# ellenorzi, hogy az olvasott  repval ertekek voltak-e elozoleg irva.




irott <- function(i,x,data) {
  x <- x[i,];
  r <- data[data$sorszam<x$sorszam &
       data$mod=="WRITE" &
       data$tt==x$tt &
       data$rr==x$rr &
       data$f1==x$f1 &
       data$f2==x$f2 &
       data$ee==x$ee &
       data$aa==x$aa,];
  rs <- r$sorszam;
  r <- nrow(r);
  names(r) <- x$sorszam;
  cat(i," ",x$sorszam," ",rs," ",r,"\n");
  r;
#  if (nrow(r)==0) {
#    cat(x$sorszam," is NOT written!\n");
#  }
}

olvasott <- function(i,x,data) {
  x <- x[i,];
  r <- data[data$sorszam<x$sorszam &
       data$mod=="read" &
       data$tt==x$tt &
       data$rr==x$rr &
       data$f1==x$f1 &
       data$f2==x$f2 &
       data$ee==x$ee &
       data$aa==x$aa,];
  r <- nrow(r);
  names(r) <- x$sorszam;
  cat(i," ",x$sorszam," ",r,"\n");
  r;
#  if (nrow(r)==0) {
#    cat(x$sorszam," is NOT written!\n");
#  }
}
