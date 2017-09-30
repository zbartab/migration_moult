
O.apd <- function(nev) {
  f.nev <- paste(nev,".apd",sep="");
  if(file.exists(f.nev)) {
    zz <- file(f.nev);
    r <- read.table(zz,header=TRUE)
  } else {
    f.nev <- paste(f.nev,".gz",sep="");
    if(file.exists(f.nev)) {
      zz <- gzfile(f.nev)
      r <- read.table(zz,header=TRUE)
    } else {
      stop(paste("File '",f.nev,"' does not exists!",sep=""));
      r <- NA;
    }
  }
  r.n <- apply(r,2,function(x) all(is.na(x)))
  r <- r[,!r.n]
  invisible(r);
}


ex.apd <- function(vars,apd,L=0:3) {
  v <- paste(vars,L,sep="")
  v.i <- v %in% names(apd)
  if(sum(v.i)>0) {
    r <- apd[,v[v.i],drop=FALSE]
  } else {
    r <- NULL
  }
  if(sum(!v.i)>0) {
    if(!is.null(r)) {
      r <- cbind(r,sapply(v[!v.i], function(i) as.numeric(rep(NA,nrow(apd)))))
    } else {
      r <- sapply(v[!v.i], function(i) as.numeric(rep(NA,nrow(apd))))
    }
  }
  r;
}


proc.apd <- function(apd,MOL=10,ML=3) {
  apd <- apd[apd$d.act=="A",]
  if(nrow(apd) <= 1) {
    n.r <- c("n", "mo.B", "mo.M", "mo.T"    ,
             "p.w0", "p.w1", "p.w2", "p.w3", "p.s0", "p.s1"    ,
             "p.s2", "p.s3", "p.te", "n.te", "p.te0", "p.te1"   ,
             "p.te2", "p.te3", "na.te0", "na.te1", "na.te2", "na.te3"  ,
             "p.br", "n.br", "p.br0", "p.br1", "p.br2", "p.br3"   ,
             "na.br0", "na.br1", "na.br2", "na.br3", "p.su", "n.su"    ,
             "p.su0", "p.su1", "p.su2", "p.su3", "na.su0", "na.su1"  ,
             "na.su2", "na.su3", "p.mo", "n.mo", "p.mo0", "p.mo1"   ,
             "p.mo2", "p.mo3", "na.mo0", "na.mo1", "na.mo2", "na.mo3"  ,
             "w.te0", "w.te1", "w.te2", "w.te3", "w.br0", "w.br1"   ,
             "w.br3", "w.br2", "w.mo0", "w.mo1", "w.mo2", "w.mo3"   ,
             "w.no1", "w.no2", "w.no3", "w.no0", "w.so0", "w.so1"   ,
             "w.so2", "w.so3", "sd.w.te0", "sd.w.te1", "sd.w.te2", "sd.w.te3",
             "sd.w.br0", "sd.w.br1", "sd.w.br3", "sd.w.br2", "sd.w.mo0",
             "sd.w.mo1", "sd.w.mo2", "sd.w.mo3", "sd.w.no1", "sd.w.no2",
             "sd.w.no3", "sd.w.no0", "sd.w.so0", "sd.w.so1", "sd.w.so2",
             "sd.w.so3", "w.f", "s.f", "w.r", "s.r", "w.u", "s.u", "sum.st",
             "win.st", "jo.spr", "jo.aut");
    r <- rep(NA,length(n.r))
    names(r) <- n.r;
    r <- data.frame(t(r));
  } else {
    r <- data.frame(n=nrow(apd));
#    ML <- max(unique(as.numeric(gsub("[^0-9]","",names(apd)))),na.rm=TRUE);
    if(!is.null(apd$mo.over)) {
      ss <- table(apd$mo.over);
      r$mo.B <- sum(ss[regexpr("B",names(ss)) > -1]);
      r$mo.M <- sum(ss[regexpr("M",names(ss)) > -1]);
      r$mo.T <- sum(ss[regexpr("T",names(ss)) > -1]);
    }
    wl <- sapply(0:ML,function(i) sum(apd$w.loc==i,na.rm=TRUE));
    names(wl) <- paste("p.w",0:ML,sep="")
    wl <- wl/sum(wl);
    sapply(names(wl), function(s) r[,s] <<- wl[s]);
    wl <- sapply(0:ML,function(i) sum(apd$u.loc==i,na.rm=TRUE));
    names(wl) <- paste("p.s",0:ML,sep="")
    wl <- wl/sum(wl);
    sapply(names(wl), function(s) r[,s] <<- wl[s]);
    
    v.nev <- "te"
    n <- ex.apd(paste("n.",v.nev,sep=""),apd,0:ML)
    h <- rowSums(n,na.rm=TRUE);
    r[,paste("p.",v.nev,sep="")] <- sum(h>0)/r$n;
    r[,paste("n.",v.nev,sep="")] <- mean(h[h>0],na.rm=TRUE)
    sapply(0:ML,function(i) {
      r[,paste("p.",v.nev,i,sep="")] <<-
        sum(n[,paste("n.",v.nev,i,sep="")]>0,na.rm=TRUE)/r$n;
    })
    sapply(0:ML,function(i) {
      r[,paste("na.",v.nev,i,sep="")] <<-
        mean(n[,paste("n.",v.nev,i,sep="")],na.rm=TRUE)
    })
    
    
    v.nev <- "br"
    n <- ex.apd(paste("n.",v.nev,sep=""),apd,0:ML)
    h <- rowSums(n,na.rm=TRUE);
    h1 <- h;
    r[,paste("p.",v.nev,sep="")] <- sum(h>0)/r$n;
    r[,paste("n.",v.nev,sep="")] <- mean(h[h>0],na.rm=TRUE)
    sapply(0:ML,function(i) {
      r[,paste("p.",v.nev,i,sep="")] <<-
        sum(n[,paste("n.",v.nev,i,sep="")]>0,na.rm=TRUE)/r$n;
    })
    sapply(0:ML,function(i) {
      r[,paste("na.",v.nev,i,sep="")] <<-
        mean(n[,paste("n.",v.nev,i,sep="")],na.rm=TRUE)
    })
    
    v.nev <- "su"
    n <- ex.apd(paste("n.",v.nev,sep=""),apd,0:ML)
    h <- rowSums(n,na.rm=TRUE);
    r[,paste("p.",v.nev,sep="")] <- sum(h>0)/r$n;
    r[,paste("n.",v.nev,sep="")] <- mean(h[h1>0],na.rm=TRUE)
    sapply(0:ML,function(i) {
      r[,paste("p.",v.nev,i,sep="")] <<-
        sum(n[,paste("n.",v.nev,i,sep="")]>0,na.rm=TRUE)/r$n;
    })
    sapply(0:ML,function(i) {
      r[,paste("na.",v.nev,i,sep="")] <<-
        mean(n[,paste("n.",v.nev,i,sep="")],na.rm=TRUE)
    })
    
    v.nev <- "mo"
    n <- ex.apd(paste("n.",v.nev,sep=""),apd,0:ML)
    h <- rowSums(n,na.rm=TRUE);
    r[,paste("p.",v.nev,sep="")] <- sum(h>0)/r$n;
    r[,paste("n.",v.nev,sep="")] <- mean(h[h>0],na.rm=TRUE)
    sapply(0:ML,function(i) {
      r[,paste("p.",v.nev,i,sep="")] <<-
        sum(n[,paste("n.",v.nev,i,sep="")]>0,na.rm=TRUE)/r$n;
    })
    sapply(0:ML,function(i) {
      r[,paste("na.",v.nev,i,sep="")] <<-
        mean(n[,paste("n.",v.nev,i,sep="")],na.rm=TRUE)
    })
    
    
    h <- colMeans(ex.apd("w.te",apd,0:ML),na.rm=TRUE);
    sapply(names(h), function(s) r[,s] <<- h[s]);
    h <- colMeans(ex.apd("w.br",apd,0:ML),na.rm=TRUE);
    sapply(names(h), function(s) r[,s] <<- h[s]);
    h <- colMeans(ex.apd("w.mo",apd,0:ML),na.rm=TRUE);
    sapply(names(h), function(s) r[,s] <<- h[s]);
    h <- colMeans(ex.apd("w.no",apd,0:ML),na.rm=TRUE);
    sapply(names(h), function(s) r[,s] <<- h[s]);
    h <- colMeans(ex.apd("w.so",apd,0:ML),na.rm=TRUE);
    sapply(names(h), function(s) r[,s] <<- h[s]);
    
    h <- apply(ex.apd("w.te",apd,0:ML),2,sd,na.rm=TRUE);
    sapply(names(h), function(s) r[,paste("sd.",s,sep="")] <<- h[s]);
    h <- apply(ex.apd("w.br",apd,0:ML),2,sd,na.rm=TRUE);
    sapply(names(h), function(s) r[,paste("sd.",s,sep="")] <<- h[s]);
    h <- apply(ex.apd("w.mo",apd,0:ML),2,sd,na.rm=TRUE);
    sapply(names(h), function(s) r[,paste("sd.",s,sep="")] <<- h[s]);
    h <- apply(ex.apd("w.no",apd,0:ML),2,sd,na.rm=TRUE);
    sapply(names(h), function(s) r[,paste("sd.",s,sep="")] <<- h[s]);
    h <- apply(ex.apd("w.so",apd,0:ML),2,sd,na.rm=TRUE);
    sapply(names(h), function(s) r[,paste("sd.",s,sep="")] <<- h[s]);
    
    r$w.f <- mean(apd$w.fea-MOL,na.rm=TRUE)
    r$s.f <- mean(apd$u.fea-MOL,na.rm=TRUE)
    r$w.r <- mean(apd$w.res,na.rm=TRUE)
    r$s.r <- mean(apd$u.res,na.rm=TRUE)
    r$w.u <- mean(apd$w.u,na.rm=TRUE)
    r$s.u <- mean(apd$u.u,na.rm=TRUE)
    su.loc <- unique(apd$u.loc[!is.na(apd$u.loc)]);
    x <- numeric(0)
    for(i in 1:length(su.loc)) {
      l <- su.loc[i]
      l1 <- if(l<ML) {l+1} else {l}
      es <- apd$u.loc == l
      n.l <- paste("w.so",l,sep="")
      n.l1 <- paste("w.no",l1,sep="")
      if(sum(c(n.l,n.l1) %in% names(apd),na.rm=TRUE)==2) {
        y <- (apd[es,n.l] - apd[es,n.l1])
      } else {
        y <- rep(52,sum(es,na.rm=TRUE))
      }
      x <- c(x,y)
    }
    r$sum.st <- mean(x,na.rm=TRUE)
    su.loc <- unique(apd$w.loc[!is.na(apd$w.loc)]);
    x <- numeric(0)
    for(i in 1:length(su.loc)) {
      l <- su.loc[i]
      l1 <- if(l>0) {l-1} else {l}
      es <- apd$w.loc == l
      n.l <- paste("w.no",l,sep="")
      n.l1 <- paste("w.so",l1,sep="")
      if(sum(c(n.l,n.l1) %in% names(apd),na.rm=TRUE)==2) {
        y <- (((26+apd[es,n.l])%%52) - ((26+apd[es,n.l1])%%52))
      } else {
        y <- rep(52,sum(es,na.rm=TRUE))
      }
      x <- c(x,y)
    }
    r$win.st <- mean(x,na.rm=TRUE)
    x <- ex.apd("w.no",apd,0:ML)
    r$jo.spr <- mean(apply(x,1,function(xx) {
      a <- abs(outer(xx,xx,"-"));
      if(sum(!is.na(a)) > 0) {
        max(a,na.rm=TRUE)
      } else {
        NA;
      }
    }),na.rm=TRUE);
    x <- ex.apd("w.so",apd,0:ML)
    r$jo.aut <- mean(apply(x,1,function(xx) {
      a <- abs(outer(xx,xx,"-"));
      if(sum(!is.na(a)) > 0) {
        max(a,na.rm=TRUE)
      } else {
        NA;
      }
    }),na.rm=TRUE);
  }
  r;
}


apd.ser <- function(patt="^.*\\.apd.*$",m.path=".",f.list=NULL) {
  if(is.null(f.list)) {
    f.list <- list.files(path=m.path,pattern=patt,full.names=TRUE);
    f.list <- gsub("\\.apd.*$","",f.list);
  }
  START <- TRUE;
  for(i in f.list) {
    print(i);
    if(START) {
      r <- data.frame(bf.converged(i),proc.apd(O.apd(i)));
      START <- FALSE;
    } else {
      r <- rbind(r,cbind(bf.converged(i),proc.apd(O.apd(i))));
    }
  }
  r$file <- f.list;
  invisible(r);
}


split.case <- function(v) {
  v <- sub("[^-]*-([0-9]+\\(_|e-\\)*[0-9]*)(-sim)*$","\\1",v)
  v <- as.numeric(gsub("_",".",v))
  v
}

ex.papd <- function(vars,papd,case) {
  if(length(vars) == 1) {
    n.papd <- names(papd)
    vars <- n.papd[regexpr(vars,n.papd) > -1]
  }
  if(length(vars)==0) stop("given variable(s) not in data.frame");
  es <- regexpr(case,papd$file) > -1;
  v <- split.case(papd$file[es]);
  r <- papd[es,vars,drop=FALSE];
  list(v=v[order(v)],r=r[order(v),]);
}

ex.papd.no.sort <- function(vars,papd,case) {
  if(length(vars) == 1) {
    n.papd <- names(papd)
    vars <- n.papd[regexpr(vars,n.papd) > -1]
  }
  if(length(vars)==0) stop("given variable(s) not in data.frame");
  es <- regexpr(case,papd$file) > -1;
  v <- split.case(papd$file[es]);
  r <- papd[es,vars,drop=FALSE];
  list(v=v,r=r);
}


plot.apd <- function(vars,papd,case,...) {
  r <- ex.papd(vars,papd,case);
  if(sum(!is.na(r$r)) > 0) {
    matplot(r$v,r$r,xlab=case,...)
  } else {
    plot(1,1,type="n");
    text(1,1,"ERROR");
  }
  invisible(r);
}

plot.both <- function(vars,case,p.w=win.df,p.s=sum.df,...) {
  layout(matrix(1:2,ncol=2));
  plot.apd(vars,p.w,case,sub="winter moult",...);
  plot.apd(vars,p.s,case,sub="summer moult",...);
  layout(1);
}


circ.week <- function(w) {
  w.rad <- ((w%%52)/51)*2*pi
  si <- sum(sin(w.rad),na.rm=TRUE)
  co <- sum(cos(w.rad),na.rm=TRUE)
  print(si)
  print(co)
  a.w <- (atan2(si,co)/(2*pi))*51
  a.w
}

comp.week <- function(w1,w2) {
  d <- abs(w1-w2)
  ifelse((d<26),(w1 < w2),(w1 > w2))
}

week.diff <- function(w1,w2) {
  ## return the difference between the two dates (in weeks)
  d <- (w2%%52)-(w1%%52)
  ad <- abs(d)
  ifelse(comp.week(w1,w2),ifelse(ad<26,-ad,-(52-ad)),ifelse(ad<26,ad,(52-ad)))
}
