qp.sim <- function(y=3,nev="baseline",ask=FALSE,print=FALSE,kep="all",
                   ext=".dat", MR=12, MF=10, cim=NULL) {
                                        #quick plot for moult data
  if (is.null(ext)) { f.nev=nev;}
  else {f.nev <- paste(nev,ext,sep="");}
  ds <- NA;
  fig.sim <- NA;
  if(file.exists(f.nev)) {
    if(kep=="all") {layout(matrix(1:2,ncol=1,nrow=2));}
    opar <- par(mar=c(4,4,2,1));
    on.exit({
      if(kep=="all") {layout(matrix(1,ncol=1,nrow=1));}
      par(opar);
    });
    zz1 <- pipe(paste("sed -n -e '1p' -e '1q' ",f.nev,sep=""));
    n.ds <- scan(zz1,what=character(0),quiet=TRUE);
    close(zz1);
    zz2 <- pipe(paste("awk '{if ($1==",y,") print $0}' ",f.nev, sep=""));
    ds <- scan(zz2,what=list(year=integer(0),week=integer(0),bird=integer(0),
                     event=character(0),moult1=character(0),
                     action=character(0),res=integer(0),fea1=integer(0),
                     ml1=integer(0),
                     brood=integer(0),exper=integer(0),uval=double(0)));
    close(zz2);
    if (length(ds$bird)!=0) {
      s.bird <- unique(ds$bird[ds$week==51]);
      e.bird <- unique(ds$bird[ds$exper==2&ds$week==0]);
      if(length(s.bird)!=0) {
        if(is.null(cim)) {
          cim <- nev;
        }
        es <- (ds$bird %in% s.bird) & (ds$bird %in% e.bird);
        weeks <- sort(unique(ds$week));
        if (kep=="all" || kep=="behav") {
          plot(0:58,0:58,type="n",ylim=c(0,1),xlab="week",ylab="Proportion");
          lines(weeks,tapply(ds$action[es]=="START"|ds$action[es]=="KEEP",
                             ds$week[es],mean),col="black",lwd=3);
          
          lines(weeks,tapply(ds$action[es]=="ABAND",
                             ds$week[es],mean),col="yellow",lwd=2);
          lines(weeks,tapply(ds$action[es]=="ABORT",
                             ds$week[es],mean),col="red",lwd=2);
          
          lines(weeks,tapply(ds$moult1[es]=="START"|ds$moult1[es]=="MOULT",
                             ds$week[es], mean),col="blue",lty=1)
          lines(weeks,tapply(ds$uval[es],ds$week[es],mean),
                col="green",lty=3,lwd=2);
          lines(c(52,55),rep(0.25,2),col="black",lwd=3)
          lines(c(52,55),rep(0.5,2),col="blue",lty=1)
          lines(c(52,55),rep(0.75,2),col="green",lty=2)
          text(55,0.25,"breed",pos=4,cex=0.6)
          text(55,0.5,"moult 1",pos=4,cex=0.6)
          text(55,0.75,"u val",pos=4,cex=0.6)
          title(cim);
        }
        if (kep=="all" || kep=="state") {
          plot(0:58,0:58,type="n",ylim=c(0,1),xlab="week",
               ylab="state");
          lines(weeks,tapply(ds$fea1[es],ds$week[es],mean)/MF,col="blue",lty=1);
          lines(weeks,tapply(ds$res[es],ds$week[es],mean)/MR,col="brown",lwd=3);
          lines(c(52,55),rep(0.2,2),col="blue",lty=1)
          lines(c(52,55),rep(0.6,2),col="brown",lwd=3)
          text(55,0.2,"fea1",pos=4,cex=0.6)
          text(55,0.6,"res",pos=4,cex=0.6)
        }
##        if (kep=="all" || kep=="mort") {
##          p.mort(ds=ds);
##       }
        fig.sim <- recordPlot();
        if (print) {dev.print();}
        if (ask) {readline("<RETURN>");}
      ##  return();
      } else {
        cat("No survivor until the end of year ",y,"\n",sep="");
        plot(1,1,type="n",xlab="",ylab="");
        text(1,1,"ERROR");
      }
    } else {
      cat("No data for year ",y,"\n",sep="");
      plot(1,1,type="n",xlab="",ylab="");
      text(1,1,"ERROR");
    }
  } else {
    cat("ERROR: File '",f.nev,"' is not exists.\n");
        plot(1,1,type="n",xlab="",ylab="");
        text(1,1,"ERROR");
  }
  invisible(fig.sim);
}


p.moult.types <- function(var="MAXEPS",x,...) {
  y <- x[,var];
  x <- as.matrix(x[order(y),c("mult","post","pren","intr","susp","arre","inco",
                              "nomo")]);
  x <- t(apply(x,1,function(xx) {xx/sum(xx)}));
  y <- y[order(y)];
  m.y <- max(y);
  nevek <- colnames(x);
  n <- ncol(x);
  my.colour <- rainbow(n+2);
  nl <- length(nevek);
  plot(c(y,1.25*m.y),c(x[,1],1),xlab=var,ylab="proportion",
       type="n",ylim=c(0,1),...)
  title(paste("moult types vs",var));
  j <- 0; k <- 1;
  for(i in nevek) {
    if (i == "post") {
      lwd <- 2; pch <- 1; lty <- 1;
    } else if (i == "pren") {
      lwd <- 2; pch <- 2; lty <- 2;
    } else if (i == "arre") {
      lwd <- 1; pch <- 3+j; lty <- 3+j;
      j <- j+1;
    } else if (i == "susp") {
      lwd <- 1; pch <- 3+j; lty <- 3+j;
      j <- j+1;
    } else if (i == "inco") {
      lwd <- 1; pch <- 3+j; lty <- 3+j;
      j <- j+1;
    } else {
      lwd <- 0.5; pch <- 3+j; lty <- 3+j;
      j <- j+1;
    }
    lines(y,x[,i],lty=lty,lwd=lwd,pch=pch,type="b",cex=0.8,
          );
    ly <- 0.05+(k-1)*(0.95/n);
    lines(c(m.y*1.1,m.y*1.15),c(ly,ly),lwd=lwd,lty=lty,cex=0.8,
          );
    points(m.y*1.125,ly,pch=pch,cex=0.8,lty=lty,lwd=lwd,
           #col=my.colour[j]
           );
    text(m.y*1.2,ly,i);
    k <- k+1;
  }
}

p.time <- function(var="MAXEPS",x,...) {
  y <- x[,var];
  x <- x[order(y),c("t.brood","t.fbrood","t.m1","t.fm1","t.m2","t.fm2","n.s")];
  y <- y[order(y)];
  m.y <- max(y);
  nevek <- names(x);
  n <- ncol(x);
  my.colour <- rainbow(n+2);
  nl <- length(nevek);
  plot(c(y,1.25*m.y),c(x[,1],1),xlab=var,ylab="week",
       type="n",ylim=c(0,52),...)
  title(paste("time vs",var));
  j <- 0; k <- 0;
  for(i in nevek[-nl]) {
    if (i == "t.brood") {
      lwd <- 2; pch <- 1; lty <- 1;
    } else if (i == "t.m1") {
      lwd <- 2; pch <- 2; lty <- 2;
    } else if (i == "t.m2") {
      lwd <- 2; pch <- 3; lty <- 3;
    } else {
      lwd <- 1; pch <- 4+j; lty <- 4+j;
      j <- j+1;
    }
    lines(y,x[,i],lty=lty,lwd=lwd,pch=pch,type="b",cex=0.8,
          );
    ly <- 0.05*52+k*((0.95*52)/n);
    lines(c(m.y*1.075,m.y*1.125),c(ly,ly),lwd=lwd,lty=lty,cex=0.8,
          );
    points(m.y*1.1,ly,pch=pch,cex=0.8,lty=lty,lwd=lwd,
           #col=my.colour[j]
           );
    text(m.y*1.2,ly,i,cex=0.8);
    k <- k+1;
  }
  text(y,0.95,x$n.s,cex=0.7,srt=90);
}

p.number <- function(var="MAXEPS",x,...) {
  y <- x[,var];
  x <- x[order(y),c("n.broods","n.m1","n.m2")];
  y <- y[order(y)];
  m.y <- max(y);
  nevek <- names(x);
  n <- ncol(x);
  my.colour <- rainbow(n+2);
  nl <- length(nevek);
  plot(c(y,1.25*m.y),c(x[,1],1),xlab=var,ylab="number",
       type="n",ylim=c(0,5),...)
  title(paste("numbers vs",var));
  j <- 0; k <- 0;
  for(i in nevek) {
    if (i == "n.broods") {
      lwd <- 2; pch <- 1; lty <- 1;
    } else if (i == "n.m1") {
      lwd <- 2; pch <- 2; lty <- 2;
    } else if (i == "n.m2") {
      lwd <- 2; pch <- 3; lty <- 3;
    } 
    lines(y,x[,i],lty=lty,lwd=lwd,pch=pch,type="b",cex=0.8,
          );
    ly <- 0.05*5+k*((0.95*5)/n);
    lines(c(m.y*1.075,m.y*1.125),c(ly,ly),lwd=lwd,lty=lty,cex=0.8,
          );
    points(m.y*1.1,ly,pch=pch,cex=0.8,lty=lty,lwd=lwd,
           #col=my.colour[j]
           );
    text(m.y*1.2,ly,i,cex=0.8);
    k <- k+1;
  }
}

p.overlap <- function(var="MAXEPS",x,...) {
  y <- x[,var];
  x <- x[order(y),];
  x <- data.frame(p.mbo=x$l.mbo/x$l.brood,p.mo=x$l.mo/x$l.moult,
                  mort.max=x$m.max);
#  x <- t(apply(x,1,function(xx) {xx/sum(xx)}));
  y <- y[order(y)];
  m.y <- max(y);
  nevek <- names(x);
  n <- ncol(x);
  my.colour <- rainbow(n+2);
  nl <- length(nevek);
  plot(c(y,1.25*m.y),c(x[,1],1),xlab=var,ylab="proportion",
       type="n",ylim=c(0,1),...)
  title(paste("overlap vs",var));
  j <- 0; k <- 0;
  for(i in nevek) {
    if (i == "p.mbo") {
      lwd <- 2; pch <- 1; lty <- 1;
    } else if (i == "p.mo") {
      lwd <- 2; pch <- 2; lty <- 2;
    } else {
      lwd <- 1; pch <- 3; lty <- 3;
    }
    lines(y,x[,i],lty=lty,lwd=lwd,pch=pch,type="b",cex=0.8,
          );
    ly <- 0.05*1+k*((0.95*1)/n);
    lines(c(m.y*1.075,m.y*1.125),c(ly,ly),lwd=lwd,lty=lty,cex=0.8,
          );
    points(m.y*1.1,ly,pch=pch,cex=0.8,lty=lty,lwd=lwd,
           #col=my.colour[j]
           );
    text(m.y*1.2,ly,i,cex=0.8);
    k <- k+1;
  }
}

p.wrapp <- function(var,nev,df,...){
  opar <- par(mfcol=c(2,2),mar=c(4,4,2,1));
  on.exit({
    par(opar);
  });
  df <- df[c(grep(nev,rownames(df)),grep("baseline",rownames(df))),];
  p.moult.types(var,df,...);
  p.time(var,df,...);
  p.number(var,df,...);
  p.overlap(var,df,...);
}


p.ser.print <- function(y=2,m.path=".",patt="\\.dat$",f.list=NULL) {
  if(is.null(f.list)) {
    f.list <- list.files(path=m.path,pattern=patt,full.names=TRUE);
    f.list <- gsub("\\.dat$","",f.list);
  }
  for(i in f.list){
    cat(i,"\n");
    qp.sim(y=y,nev=i,print=TRUE);
  }
}



