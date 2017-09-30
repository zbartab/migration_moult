create.ini <- function(path=".", pattern="\\.ini$",y=2) {
  f.inik <- list.files(path=path,pattern=pattern,full.names=TRUE);
  zz1 <- pipe(paste("sed -e 's/		*/	/' ",f.inik[1],sep=""));
  ini1 <- read.table(file=zz1,header=FALSE,sep="\t");
##  close(zz1);
  m.types <- list();
  dat.nev <- sub("\\.ini$",".dat.gz",f.inik[1]);
  if (file.exists(dat.nev)) {
    system(paste("gunzip -f ",dat.nev));
  } 
  dat.nev <- sub("\\.gz$","",dat.nev);
  print(dat.nev);
##  qp.sim(y=y,nev=dat.nev,print=FALSE,kep="both",ext=NULL);
  m.desc <- moult.descr(y=y,nev=dat.nev,ext=NULL);
  m.census <- moult.census(y=y,nev=dat.nev,ext=NULL);
  m.types[[f.inik[1]]] <- moult.types(y=y,nev=dat.nev,
                                      ext=NULL);
  mort <- mortality(y=y,nev=sub("\\.ini$",".C.pop",f.inik[1]),ext=NULL);
  yrv <- young.rv(y=y,nev=sub("\\.ini$",".yrv",f.inik[1]),ext=NULL);
  m.t.desc <- moult.typ.descr(m.types[[f.inik[1]]]);
  m.class <- moult.classify(m.types[[f.inik[1]]]);
  inik.df <- data.frame(matrix(
                               c(as.vector(ini1[,2]),m.desc,m.t.desc,mort,yrv,
                                 m.class,m.census),
                               nrow=1));
  n <- sumo(f.inik[1]);
  nevek <- c(gsub(" +","",gsub("_",".",as.character(ini1[,1]))),
             names(m.desc),names(m.t.desc),names(mort),names(yrv),
             names(m.class),names(m.census));
  names(inik.df) <- make.names(nevek);
##  if (file.exists(dat.nev)) {
##    print("gunzip 2");
##    system(paste("gzip -f ",dat.nev));
##  }
  for(f.nev in f.inik[-1]) {
    print(f.nev);
    dat.nev <- sub("\\.ini$",".dat.gz",f.nev);
    if (file.exists(dat.nev)) {
      print("gunzip 3");
      system(paste("gunzip -f ",dat.nev));
    }
    dat.nev <- sub("\\.gz$","",dat.nev);
    zz1 <- pipe(paste("sed -e 's/		*/	/' ",f.nev,sep=""));
    ini1 <- read.table(file=zz1,header=FALSE,sep="\t");
##    qp.sim(y=y,nev=dat.nev,print=FALSE,kep="both",ext=NULL);
    m.desc <- moult.descr(y=y,nev=dat.nev,ext=NULL);
    m.census <- moult.census(y=y,nev=dat.nev,ext=NULL);
    m.types[[f.nev]] <- moult.types(y=y,nev=dat.nev,
                                        ext=NULL);
    mort <- mortality(y=y,nev=sub("\\.ini$",".C.pop",f.nev),ext=NULL);
    yrv <- young.rv(y=y,nev=sub("\\.ini$",".yrv",f.nev),ext=NULL);
    m.t.desc <- moult.typ.descr(m.types[[f.nev]]);
    m.class <- moult.classify(m.types[[f.nev]]);
    inik.df <- rbind(inik.df,c(ini1[,2],m.desc,m.t.desc,mort,yrv,m.class,
                               m.census));
    n <- c(n,sumo(f.nev));
##    if (file.exists(dat.nev)) {
##      system(paste("gzip -f ",dat.nev));
##    }
  }
  f.inik <- gsub("\\.ini$","",f.inik);
  f.inik <- gsub("^.*/","",f.inik);
  rownames(inik.df) <- f.inik;
  inik.df <- data.frame(inik.df,n.y=as.numeric(n));
  invisible(list(database=inik.df,type.list=m.types))
}

sumo <- function(i.nev) {
  f.nev <- sub("\\.ini$",".sum",i.nev);
  if(file.exists(f.nev)) {
    sum1 <- scan(file=f.nev,
                 what=character(0),sep="\n",quiet=TRUE);
    sum1 <- grep("^bac",sum1,value=TRUE);
    if (length(sum1)!=0) {
      gsub("^.*= ","",sum1);
    } else {
      NA;
    }
  } else {
    NA;
  }
}

moult.types2 <- function(ds) {
### return the moult-breeding type of a bird where letters mark the
### event (B: breeding; P: moulting feather 1; S: moulting feather 2),
### and there order is corresponde the order of events.

  breeding <- ds$week[ds$action=="START"];
  if(length(breeding)==0) breeding <- -1;
  moultP <- ds$week[ds$moult1=="START"];
  if(length(moultP)==0) moultP <- -1;
  moultS <- ds$week[ds$moult2=="START"];
  if(length(moultS)==0) moultS <- -1;
  r <- rbind(cbind(breeding,"B"), cbind(moultP,"p"), cbind(moultS,"s"));
  r <- r[as.numeric(r[,1])!= -1,,drop=FALSE];
  r <- paste(r[order(as.numeric(r[,1])),2],collapse="");
  r;
}

moult.types <- function(y=2,nev="moult",ext=".dat") {
  if (is.null(ext)) { f.nev=nev;}
  else {f.nev <- paste(nev,ext,sep="");}
  if(file.exists(f.nev)) {
    command <- paste(.path.package("Rone"),"/scripts/mtype.awk -v y=",y," ",
                     f.nev,sep="");
    zz1 <- pipe(command);
    r <- scan(zz1,what=list(character(0),integer(0)),quiet=TRUE);
    close(zz1);
    m.type <- r[[2]];
    names(m.type) <- r[[1]];
    m.type <- sort(m.type,decreasing=TRUE);
  } else {
    cat("ERROR: File '",f.nev,"' is not exists.\n");
    m.type <- numeric(0);
#    names(m.type) <- "n";
  }
  m.type;
}

ind.moult.types <- function(y=2,nev="moult",ext=".dat") {
  if (is.null(ext)) { f.nev=nev;}
  else {f.nev <- paste(nev,ext,sep="");}
  if(file.exists(f.nev)) {
    command <- paste(.path.package("Rone"),"/scripts/ind_mtype.awk -v y=",
                     y," ", f.nev,sep="");
    zz1 <- pipe(command);
    r <- scan(zz1,what=list(numeric(0),character(0)),quiet=TRUE);
    close(zz1);
    m.type <- r; ##cbind(r[[1]],r[[2]]);
    names(m.type) <- c("bird","routine");
  } else {
    cat("ERROR: File '",f.nev,"' is not exists.\n");
    m.type <- numeric(0);
#    names(m.type) <- "n";
  }
  m.type;
}

life.types <- function(nev="moult",ext=".dat") {
  if (is.null(ext)) { f.nev=nev;}
  else {f.nev <- paste(nev,ext,sep="");}
  if(file.exists(f.nev)) {
    command <- paste(.path.package("Rone"),"/scripts/life_type.awk ",
                     f.nev,sep="");
    zz1 <- pipe(command);
    r <- scan(zz1,what=list(character(0),integer(0)),quiet=TRUE);
    close(zz1);
    m.type <- r[[2]];
    names(m.type) <- r[[1]];
##    m.nevek <- gsub("B+","B",gsub("(ps|sp)","m",names(m.type)));
    m.nevek <- gsub("B+","B",names(m.type));
    names(m.type) <- m.nevek;
    m.nevek <- unique(grep(".*\\|$",m.nevek,value=TRUE));
    gyak <- numeric(length=length(m.nevek));
    names(gyak) <- m.nevek;
    for(i in 1:length(m.nevek)) {
      gyak[m.nevek[i]] <- gyak[m.nevek[i]] + m.type[m.nevek[i]];
    }
    gyak <- sort(gyak,decreasing=TRUE);
    gyak;
  } else {
    cat("ERROR: File '",f.nev,"' is not exists.\n");
  }
}


moult.typ.descr <- function(m.types){
  n <- length(m.types);
  p <- m.types/sum(m.types);
  H <- -1*sum(log(p)*p);
  E <- H/log(n);
  r <- c(n,H,E,max(p));
  names(r) <- c("S","H","E","pmax");
  r;
}

mortality <- function(y=3,nev="moult",ext=".C.pop",ll=1) {
  if (is.null(ext)) { f.nev=nev;}
  else {f.nev <- paste(nev,ext,sep="");}
  r <- rep(NA,4);
  names(r) <- c("m.min","m.max","t.min","t.max");
  if (file.exists(f.nev)) {
    pop <- read.table(f.nev,header=FALSE);
    v.nevek <- c("alive","breed","start","aban",
                 "moultP","moultS");
    szamok <- 1:ll
    nevek <- paste(v.nevek,rep(szamok,rep(length(v.nevek),length(szamok))),
                   sep="");
    colnames(pop) <- c("week",nevek);
    n <- nrow(pop);
    n.years <- n/52;
    year <- rep(0:(n.years-1),rep(52,n.years));
    pop <- cbind(year=year,pop)
    es <- pop$year==y;
    mort <- sapply(2:length(pop$alive1),
                   function(i) (pop$alive1[i-1]-pop$alive1[i])/pop$alive1[i-1]);
    mort <- c(NA,mort);
    if(sum(pop$alive1[es])!=0) {
      mort <- mort[es];
      m.max <- max(mort,na.rm=TRUE);
      m.min <- min(mort,na.rm=TRUE);
      t.max <- which.max(mort);
      t.min <- which.min(mort);
      r <- c(m.min,m.max,t.min,t.max);
      names(r) <- c("m.min","m.max","t.min","t.max");
    } else {
      cat("ERROR: No survival in year ",y,"\n");
    }
  } else {
    cat("ERROR: File '",f.nev,"' is not exists.\n");
  }
  r;
}  


young.rv <- function(y=3,nev="moult",ext=".yrv") {
  if (is.null(ext)) { f.nev=nev;}
  else {f.nev <- paste(nev,ext,sep="");}
  r <- rep(NA,4);
  names(r) <- c("rv.min","rv.max","rt.min","rt.max");
  if (file.exists(f.nev)) {
    yrv <- scan(file=f.nev,quiet=TRUE);
    m.max <- max(yrv,na.rm=TRUE);
    m.min <- min(yrv,na.rm=TRUE);
    t.max <- which.max(yrv);
    t.min <- which.min(yrv);
    r <- c(m.min,m.max,t.min,t.max);
    names(r) <- c("rv.min","rv.max","rt.min","rt.max");
  } else {
    cat("ERROR: File '",f.nev,"' is not exists.\n");
  }
  r;
}  
   
moult.descr <- function(y=3,nev="moult",ext=".dat") {

  if (is.null(ext)) { f.nev=nev;}
  else {f.nev <- paste(nev,ext,".gz",sep="");}
##  if(file.exists(f.nev)) {
##    print("gunzip 1.0");
##    system(paste("gunzip -f ",f.nev));
##  }
  f.nev <- sub("\\.gz$","",f.nev);
  nevek <-
    c("n.broods","n.aband","n.abort","n.m1","l.brood",
      "l.m1","l.mbo", "t.brood","t.fbrood","t.m1",
      "t.fm1","t.bm","n.s");
  r <- rep(NA,length(nevek));
  names(r) <- nevek;
  if(file.exists(f.nev)) {
    command <- paste(.path.package("Rone"),"/scripts/mdescr.awk -v y=",y," ",
                     f.nev,sep="");
    zz1 <- pipe(command);
    r <- scan(zz1,quiet=TRUE);
    close(zz1);
    names(r) <- nevek;
  } else {
    cat("ERROR: File '",f.nev,"' is not exists.\n");
  }
  r;
}

moult.osztaly <- function(mt) {
  m.t <- unlist(strsplit(mt,""));
  if ((ms <- sum(m.t=="m"))>1) {return("mult");}
  else if (ms == 1) {
    if (regexpr("[ps]",mt) > -1) {return("arre");}
    else if (m.t[length(m.t)] == "m") {return("post");}
    else if (m.t[1] == "m") {return("pren");}
    else {return("intr");}
  } else if (ms == 0) {
    p <- regexpr("p",mt);
    s <- regexpr("s",mt);
    if ((p > -1) && (s > -1)) {
      return("susp");
    } else if ((p > -1) || (s > -1)) {return("inco");}
    else {return("nomo");}
  }
}

moult.classify <- function(mt) {
  n.gyak <- c("mult","post","pren","intr",
              "susp","arre","inco","nomo");
  if(length(mt) > 0) {
    m.types <- names(mt);
    m.types <- gsub("B+","B",gsub("p","m",names(mt)));
    for(i in 1:length(m.types)) {
      m.types[i] <- moult.osztaly(m.types[i]);
    }
    gyak <- numeric(length=length(n.gyak));
    names(gyak) <- n.gyak;
    for(i in 1:length(mt)) {
      gyak[m.types[i]] <- gyak[m.types[i]] + mt[i];
    }
  } else {
    cat("ERROR: variable(mt) has zero length\n");
    gyak <- rep(NA,length(n.gyak));
    names(gyak) <- n.gyak;
  }
  gyak;
}

moult.census <- function(y=NULL,nev="moult",ext=".dat") {

  if (is.null(ext)) { f.nev=nev;}
  else {f.nev <- paste(nev,ext,sep="");}
  r <- rep(NA,10);
  nevek <-
    c("year", "n.0","n.1","n.2","pred.0","pred.1","pred.2","starv.0",
      "starv.1","starv.2");
  names(r) <- nevek;
  if(file.exists(f.nev)) {
    command <- paste(.path.package("Rone"),"/scripts/census.awk ",
                     f.nev,sep="");
    zz1 <- pipe(command);
    r <- scan(zz1,what=list(year=numeric(0), n.0=numeric(0),n.1=numeric(0),
                    n.2=numeric(0), pred.0=numeric(0),pred.1=numeric(0),
                    pred.2=numeric(0), starv.0=numeric(0),starv.1=numeric(0),
                    starv.2=numeric(0)),
                    quiet=TRUE);
    close(zz1);
    r <- as.data.frame(r);
    if(!is.null(y)) {
      nevek <- colnames(r[,-1]);
      if (max(r$year) < y) {
        r <- rep(NA,length(nevek));
      } else {
        r <- as.matrix(r[r$year==y,-1]);
        r <- as.vector(r);
      }
      names(r) <- nevek;
    }
  } else {
    cat("ERROR: File '",f.nev,"' is not exists.\n");
  }
  r;
}


read.inis <- function(path=".", pattern="\\.ini$",y=2) {
  f.inik <- list.files(path=path,pattern=pattern,full.names=TRUE);
  start <- TRUE;
  for(f.nev in f.inik) {
    print(f.nev);
    zz1 <- pipe(paste("sed -e 's/		*/	/' ",f.nev,sep=""));
    ini1 <- read.table(file=zz1,header=FALSE,sep="\t");
    zz2 <- pipe(paste("grep -c ERR ",sub("\\.ini$",".sum",f.nev),sep=""));
    hiba <- scan(zz2,quiet=TRUE);
    close(zz2);
    zz3 <- pipe(paste("grep forw ",sub("\\.ini$",".txt",f.nev),
                      " |tail -1|tr \"=\" \" \"|tr -d \",\"|cut -d \" \" -sf 4",
                      sep=""));
    l.for <- scan(zz3,quiet=TRUE);
    close(zz3);
    zz4 <- pipe(paste("grep back ",sub("\\.ini$",".txt",f.nev),
                      " |tail -1|cut -d \" \" -sf 6",sep=""));
    l.bac <- scan(zz4,quiet=TRUE);
    close(zz4);
    
  print(hiba)
    if(start) {
      start <- FALSE;
      inik.df <- data.frame(matrix(as.vector(ini1[,2]), nrow=1),hiba,
                            l.bac,l.for);
      nevek <- c(gsub(" +","",gsub("_",".",as.character(ini1[,1]))),"hiba",
                 "la.back","la.forw");
      names(inik.df) <- make.names(nevek);
    } else {
      inik.df <- rbind(inik.df,c(ini1[,2],hiba,l.bac,l.for));
    }
  }
  f.inik <- gsub("\\.ini$","",f.inik);
  f.inik <- gsub("^.*/","",f.inik);
  rownames(inik.df) <- f.inik;
  invisible(inik.df);
}




moult.proc.dat <- function(path=".", pattern="\\.dat$",y=2,f.list=NULL) {
  if(is.null(f.list)) {
    f.inik <- list.files(path=path,pattern=pattern,full.names=TRUE);
  } else {
    f.inik <- f.list;
  }
  m.types <- list();
  start <- TRUE;
  for(f.nev in f.inik) {
    print(f.nev);
    m.desc <- moult.descr(y=y,nev=f.nev,ext=NULL);
    m.census <- moult.census(y=y,nev=f.nev,ext=NULL);
    m.types[[f.nev]] <- moult.types(y=y,nev=f.nev,ext=NULL);
    m.t.desc <- moult.typ.descr(m.types[[f.nev]]);
    m.class <- moult.classify(m.types[[f.nev]]);
    if(start) {
      inik.df <- data.frame(matrix(c(m.desc,m.t.desc,m.class,m.census),nrow=1));
      nevek <- c(names(m.desc),names(m.t.desc),names(m.class),names(m.census));
      names(inik.df) <- make.names(nevek);
      start <- FALSE;
    } else {
      inik.df <- rbind(inik.df,c(m.desc,m.t.desc,m.class,m.census));
    }
  }
  f.inik <- gsub("\\.dat$","",f.inik);
  f.inik <- gsub("^.*/","",f.inik);
  rownames(inik.df) <- f.inik;
  invisible(list(database=inik.df,type.list=m.types))
}



moult.policy <- function(x) {
  n.r <- c("C","I","BI","IB","CB","BC");
  if (length(x) != 0) {
    n <- sum(x,na.rm=TRUE);
    n.x <- names(x);
    n.x <- sub("BB*","B",n.x);
    n.x <- sub("^(sp|ps)$","C",n.x);
    n.x <- sub("^(p|s)$","I",n.x);
    n.x <- sub("^(Bp|Bs)$","BI",n.x);
    n.x <- sub("^(sB|pB)$","IB",n.x);
    n.x <- sub("^(spB|psB)$","CB",n.x);
    n.x <- sub("^(Bsp|Bps)$","BC",n.x);
    r <- rep(0,6);
    names(r) <- n.r;
    for(i in 1:length(x)) {
      if (n.x[i] %in% names(r)) {
        r[n.x[i]] <- r[n.x[i]] + x[i]/n;
      }
    }
  } else {
    r <- rep(NA,6);
    names(r) <- n.r;
  }
  r;
}



moult.policy.red <- function(x) {
  n.r <- c("m","Bm","mB");
  if (length(x) != 0) {
    n <- sum(x,na.rm=TRUE);
    n.x <- names(x);
    n.x <- gsub("[ad]","",n.x);
    n.x <- sub("BB*","B",n.x);
    n.x <- sub("^(sp|ps)$","m",n.x);
    n.x <- sub("^(p|s)$","m",n.x);
    n.x <- sub("^(Bp|Bs)$","Bm",n.x);
    n.x <- sub("^(sB|pB)$","mB",n.x);
    n.x <- sub("^(spB|psB)$","mB",n.x);
    n.x <- sub("^(Bsp|Bps)$","Bm",n.x);
    r <- rep(0,3);
    names(r) <- n.r;
    for(i in 1:length(x)) {
      if (n.x[i] %in% names(r)) {
        r[n.x[i]] <- r[n.x[i]] + x[i]/n;
      }
    }
  } else {
    r <- rep(NA,3);
    names(r) <- n.r;
  }
  r;
}


moult.order <- function(mt) {
  r <- numeric();
  if (length(mt)==0) {
    r <- c(sp=NA,ps=NA);
  } else {
    n.m <- names(mt);
    n.m <- gsub("ss*","s",n.m);
    n.m <- gsub("pp*","p",n.m);
    n.m <- gsub("[Bad]","",n.m);
    for(i in 1:length(mt)) {
      if(is.na(r[n.m[i]])) {r[n.m[i]] <- mt[i];}
      else {r[n.m[i]] <- r[n.m[i]] + mt[i];}
    }
    r <- r[c("sp","ps")]/sum(mt);
  }
  r;
}

