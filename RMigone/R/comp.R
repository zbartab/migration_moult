comp.res <- function(nev1,nev2,y=2) {
                                        #compare the reserves for two runs

  ds1 <- O.sim(nev=nev1,y=y);
  ds2 <- O.sim(nev=nev2,y=y);
  s.bird1 <- unique(ds1$bird[ds1$week==51]);
  s.bird2 <- unique(ds2$bird[ds2$week==51]);
  es1 <- ds1$bird %in% s.bird1;
  es2 <- ds2$bird %in% s.bird2;
  weeks <- sort(unique(ds1$week));
  plot(0:58,0:58,type="n",ylim=c(0,1),xlab="week",
       ylab="reserves");
  lines(weeks,tapply(ds1$res[es1],ds1$week[es1],mean)/12,col="black",lwd=1);
  lines(weeks,tapply(ds2$res[es2],ds2$week[es2],mean)/12,col="red",lwd=1,lty=2);
  lines(c(52,54),rep(0.33,2),col="black",lty=1)
  lines(c(52,54),rep(0.66,2),col="red",lty=2)
  text(54,0.33,nev1,pos=4,cex=0.6)
  text(54,0.66,nev2,pos=4,cex=0.6)
  title(paste(nev1,"vs",nev2));
}


p.mort <- function(nev="baseline",y=2,ds=NULL) {
  if(is.null(ds)) {
    ds <- O.sim(nev=nev,y=y);
  } else {
    nev <- "";
  }
  starv <- tapply(ds$event=="STARV",ds$week,mean);
  preda <- tapply(ds$event=="PREDA",ds$week,mean);
  max.y <- max(c(starv,preda));
  plot(0:60,rep(0,61),type="n",ylim=c(0,max.y),
       xlab="week",ylab="mortality",main=nev);
  lines(0:(length(starv)-1),starv,lty=1);
  lines(0:(length(preda)-1),preda,lty=2);
  points(0:(length(starv)-1),starv,pch=16);
  points(0:(length(preda)-1),preda);
  y1 <- max.y*0.6;
  y2 <- max.y*0.8;
  lines(c(52,56),rep(y1,2),lty=1);
  lines(c(52,56),rep(y2,2),lty=2);
  points(54,y1,pch=16);
  points(54,y2);
  text(56,y1,"starv.",pos=4,cex=0.6)
  text(56,y2,"pred.",pos=4,cex=0.6)
}
