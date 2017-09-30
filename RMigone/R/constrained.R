con.xy <- function(xn,yn, es, odf=one.df,...) {
  x.y <- c(odf[es&odf$con=="yes",xn],odf[odf$ser=="baseline-con",xn]);
  x.n <- c(odf[es&odf$con=="no", xn],odf[odf$ser=="baseline",xn]);
  y.y <- c(odf[es&odf$con=="yes",yn],odf[odf$ser=="baseline-con",yn]);
  y.n <- c(odf[es&odf$con=="no", yn],odf[odf$ser=="baseline",yn]);
  plot(c(x.y,x.n),c(y.y,y.n),type="n",xlab=xn,ylab=yn,...)
  points(x.n,y.n)
  points(x.y,y.y,pch=3);
  abline(v=odf[odf$ser=="baseline",xn],lty=3,col=grey(0.5));
  abline(h=odf[odf$ser=="baseline",yn],lty=3,col=grey(0.5));
}


rajz.con <- function(xn,es,odf=one.df,...) {
  layout(matrix(1:8,ncol=4,byrow=TRUE));
  on.exit({
    layout(1);
  })
  con.xy(xn,"n.broods",es,odf,...);
  con.xy(xn,"n.abort",es,odf,...);
  con.xy(xn,"n.m1",es,odf,...);
  con.xy(xn,"n.s",es,odf,...);
  con.xy(xn,"t.brood",es,odf,...);
  con.xy(xn,"t.m1",es,odf,...);
  con.xy(xn,"t.bm",es,odf,...);
}
