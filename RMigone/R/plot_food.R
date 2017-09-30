
food <- function(fp) {
  w <- 0:51
  f <- fp[1] + fp[2]*sin(pi*((w-13)/26));
  ifelse(f<0,0,f);
}


plot.food <- function(A,e=NULL,...) {
  if(is.null(e)) {
    f <- data.frame(A)
  } else {
    f <- data.frame(A=A,e=e);
  }
  fv <- apply(f,1,food);
  matplot(fv,xlab="Weeks",...)
}

