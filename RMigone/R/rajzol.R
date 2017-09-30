

rajz.pol <- function(nev="baseline") {
  sr <- G.sta(nev=nev,typ="rv");
  sa <- G.sta(nev=nev,typ="ac");
  st <- O.bin(nev=nev,typ="st");
  rv <- O.bin(nev=nev,typ="rv");
  st <- st[sr$T!=52];
  rv <- rv[sr$T!=52];
  rm(sr);
  ac <- O.bin(nev=nev,typ="ac");
  uv <- O.bin(nev=nev,typ="uv");
  mo <- O.bin(nev=nev,typ="mo");
  stt <- tapply(st,sa$T,sum);
  pro.t <- list();
  actions <- unique(ac);
  for(i in actions) {
    pro.t[[as.character(i)]] <- tapply((ac==i)*st,sa$T,sum);
  }
  UV <- tapply(uv*st,sa$T,sum);
  MO <- tapply((sa$F<0)*st,sa$T,sum);
  plot(0:51,0:51,type="n",ylim=c(0,1));
  for(i in actions) {
    lines(0:51,pro.t[[as.character(i)]]/stt,col=i,lty=i);
  }
  lines(0:51,UV/stt,col="orange",lwd=2,lty=2)
  lines(0:51,MO/stt,col="brown",lwd=2,lty=2)
  lines(0:51,stt/max(stt),col="red",lwd=3);
  legend(50,1,xjust=1,legend=actions,col=actions,lty=actions);
  sa$uv <- uv;
  sa$ac <- ac;
  sa$mo <- mo;
  sa$rv <- rv;
  sa$st <- st;
  invisible(sa);
}

rajz.behav <- function(s,es=TRUE,rajz=TRUE) {
  ido <- range(s$T[es]);
  ido <- ido[1]:ido[2]-1;
  beh <- list()
  beh[["ss"]] <- tapply(s$st[es],s$T[es],sum);
  beh[["st"]] <- tapply((s$ac[es]==2)*s$st[es],s$T[es],sum);
  beh[["ke"]] <- tapply((s$ac[es]==3)*s$st[es],s$T[es],sum);
  beh[["ab"]] <- tapply((s$ac[es]==4)*s$st[es],s$T[es],sum);
  beh[["at"]] <- tapply((s$ac[es]==5)*s$st[es],s$T[es],sum);
  beh[["no"]] <- tapply((s$ac[es]==6)*s$st[es],s$T[es],sum);
  beh[["so"]] <- tapply((s$ac[es]==7)*s$st[es],s$T[es],sum);
  beh[["bi"]] <- tapply((s$ac[es]==10)*s$st[es],s$T[es],sum);
  m.beh <- max(unlist(beh));
  l.beh <- length(beh);
  if(rajz) {
    plot(ido,ido,type="n",ylim=c(0,m.beh),xlab="time",ylab="prop");
    l.wd = 1;
  } else {l.wd <- 2;}
  for(i in 1:l.beh) {
    lines(ido,beh[[i]],col=i,lty=i,lwd=l.wd);
  }
  if(rajz)
    legend(max(ido),m.beh,legend=names(beh),col=1:l.beh,lty=1:l.beh,xjust=1);
}
