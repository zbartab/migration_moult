plot.sim1 <- function(nev="baseline",y=3, ext=".dat", exper=2,MR=12, MF=10,
                      MOL=10,cim=NULL) {

  ## quick plot for migra1 data
  ## `y' plot this year
  ## `nev' is the name of the dat file WITHOUT extension
  ## `ext' is the extension used
  ## `MR', `MF' and `MOL' are the same as in the ini file
  ## `cim' is the title of the graph
  
  if (is.null(ext)) { f.nev=nev;}
  else {f.nev <- paste(nev,ext,sep="");}
  ds <- NA;
  fig.sim <- NA;
  if(!file.exists(f.nev)){
    cat("ERROR: File '",f.nev,"' is not exists.\n");
    plot(1,1,type="n",xlab="",ylab="");
    text(1,1,"ERROR");
  } else  {
    opar <- par(mar=c(4,4,2,1));
    on.exit({
      layout(matrix(1,ncol=1,nrow=1));
      par(opar);
    });
    if(.Platform$OS.type=="unix") {
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
    } else {
      ds <- read.table(f.nev,header=TRUE,as.is=TRUE);
      ds <- ds[ds$year==y,];
    }
    if (length(ds$bird)==0) {
      cat("No data for year ",y,"\n",sep="");
      plot(1,1,type="n",xlab="",ylab="");
      text(1,1,"ERROR");
    } else {
      s.bird <- unique(ds$bird[ds$week==51]);
      e.bird <- unique(ds$bird[ds$exper==exper&ds$week==0]);
      n.loc <- length(unique(ds$loc));
      if(length(s.bird)==0 || length(e.bird)==0)  {
        cat("No survivor until the end of year ",y,"\n",sep="");
        plot(1,1,type="n",xlab="",ylab="");
        text(1,1,"ERROR");
      } else {
        if(is.null(cim)) {
          cim <- nev;
        }
        layout(matrix(1:(3*n.loc),ncol=n.loc,byrow=FALSE));
        es.s <- (ds$event=="A") & (ds$bird %in% s.bird) & (ds$bird %in% e.bird);
        for(l in 0:(n.loc-1)) {
          es <- es.s & (ds$loc == l);
          if(sum(es)>0) {
            weeks <- sort(unique(ds$week[es]));
            denz <- tapply(ds$bird[es],ds$week[es],length);
            m.denz <- max(denz);
            plot(0:58,0:58,type="n",ylim=c(0,m.denz),
                 xlab="week",ylab="density");
            lines(weeks,denz,lty=3);
            lines(weeks,tapply(ds$action[es]=="S"|ds$action[es]=="K"|
                               ds$action[es]=="A",ds$week[es],sum),
                  col="black",lwd=1,type="o");
            lines(weeks,tapply(ds$moult1[es]!="N",ds$week[es], sum),
                  col="blue",lty=1,type="o",pch=2)
            lines(weeks,tapply(regexpr("B$",ds$action[es])!=-1,
                               ds$week[es],sum),
                  col="green",lty=1,lwd=1,type="o",pch=3);
            lines(weeks,tapply(regexpr("^No",ds$action[es])!=-1,
                               ds$week[es],sum),
                  col="red",lty=1,lwd=1,type="o",pch=4);
            lines(weeks,tapply(regexpr("^So",ds$action[es])!=-1,
                               ds$week[es],sum),
                  col="brown",lty=1,lwd=1,type="o",pch=5);
            legend(60,m.denz,bty="n",lty=c(3,rep(1,5)),cex=0.8,
                   legend=c("total","breed","moult","terr.","north","south"),
                   col=c("black","black","blue","green","red","brown"),
                   xjust=1,pch=c(NA,1:5))
            title(paste(cim,"l=",l));
            plot(0:58,0:58,type="n",ylim=c(0,1),xlab="week",
                 ylab="state");
            lines(weeks,tapply(ds$fea1[es]-MOL,ds$week[es],mean)/MF,
                  col="blue",lty=1);
            lines(weeks,tapply(ds$res[es],ds$week[es],mean)/MR,
                  col="black",lty=2);
            lines(weeks,tapply(ds$uval[es],ds$week[es],mean),
                  col="red",lty=3);
            legend(60,1,legend=c("fea1","res","u"),col=c("blue","black","red"),
                   xjust=1,bty="n",lty=1:3,cex=0.8);
            es <- ds$bird %in% e.bird & ds$loc == l;
            weeks <- sort(unique(ds$week[es]));
            mort <- tapply(ds$event[es]!="A",ds$week[es],sum,na.rm=TRUE);
            plot(0:58,0:58,type="n",ylim=c(0,max(mort,na.rm=TRUE)),xlab="week",
                 ylab="mortality");
            lines(weeks,mort);
          }
        }
        fig.sim <- recordPlot();
      } 
    } 
  } 
##  invisible(fig.sim);
  invisible(ds);
}


p.ser <- function(y=2,m.path=".",patt="\\.dat$",f.list=NULL,print=FALSE,
                  kep="all",descr=TRUE,types=TRUE,census=TRUE,ask=TRUE,
                  MF=10,MR=12,MOL=8) {
  if(is.null(f.list)) {
    f.list <- list.files(path=m.path,pattern=patt,full.names=TRUE);
    f.list <- gsub("\\.dat$","",f.list);
  }
  fig.list <- list();
  for(i in f.list){
    cat("-----------------------------------------\n");
    cim <- paste("File:", i);
    cat(cim,"\n");
    fig.list[[i]] <- plot.sim1(y=y,nev=i,MF=MF,MR=MR, MOL=MOL,cim=cim);
##    if(descr) {print(moult.descr(y=y,nev=i));
##               cat("\n");}
##    if(types) {print(moult.types(y=y,nev=i));
##               cat("\n");}
##    if(census) {print(moult.census(nev=i));
##                cat("\n");}
    if(ask) {readline("Return..");}
  }
  invisible(fig.list);
}



