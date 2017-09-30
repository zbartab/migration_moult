
O.routine <- function(name="baseline",y=3,ext=".dat",type="adult") {
  ## read the routines individuals follow
  if(is.null(ext)) {
    f.name=name;
  } else {
    f.name <- paste(name,ext,sep="");
  }
  if(file.exists(f.name)) {
    if(type=="adult") {
      command <- paste(.path.package("RMigone"),
                       "/scripts/pre_process.awk -v y=",y," ",
                       f.name,sep="");
    } else {
      command <- paste(.path.package("RMigone"),
                       "/scripts/pre_all.awk -v y=",y," ",
                       f.name,sep="");
    }
    zz1 <- pipe(command);
    r <- scan(zz1,what=character(0),quiet=TRUE);
    close(zz1);
  } else {
    r <- NA;
  }
  if(type=="adult") {
    invisible(rou.split(r));
  } else {
    invisible(r);
  }
}

rou.split <- function(r) {
  r.s <- matrix(unlist(strsplit(sub("^\\|","",r),"\\|")),nrow=length(r),
                byrow=T);
  r.s;
}


rou.df.act <- function(r, act="B", loc=0) {
  a.w <- rou.w.act(r,act,loc);
  a.r <- rou.r.act(r,act,loc);
  a.f <- rou.f.act(r,act,loc);
  a.u <- rou.u.act(r,act,loc);
  rr <- cbind(a.w,a.r,a.f,a.u);
  colnames(rr) <- c("week","res","fea","u");
  rr;
}

rou.n.act <- function(r,act) {
  if(nchar(act) != 1) {
    stop("`act' should be one character");
  }
  r.b <- gsub(paste("[^",act,"]",sep=""),"",r);
  nchar(r.b);
}

rou.w.act <- function(r,act="B",loc=0) {
  if(nchar(act) != 1) {
    stop("`act' should be one character");
  }
  r.s <- strsplit(gsub("\\)","",sub("^\\(","",r)),"\\(");
  rr <- sapply(r.s, function(x) {
    y <- as.numeric(sub("^.,([0-9]+),.*$","\\1",
                        x[grep(paste("^",act,",[0-9]+,",loc,",.*$",sep=""),
                               x)[1]]));
    if(length(y) == 0) {NA}
    else {y}
  });
  rr;
}

rou.r.act <- function(r,act="B",loc=0) {
  if(nchar(act) != 1) {
    stop("`act' should be one character");
  }
  r.s <- strsplit(gsub("\\)","",sub("^\\(","",r)),"\\(");
  rr <- sapply(r.s, function(x) {
    y <- as.numeric(sub("^.,[0-9]+,[0-9]+,([0-9]+).*$","\\1",
                        x[grep(paste("^",act,",[0-9]+,",loc,",.*$",sep=""),
                               x)[1]]));
    if(length(y) == 0) {NA}
    else {y}
  });
  rr;
}

rou.f.act <- function(r,act="B",loc=0,MOL=8) {
  if(nchar(act) != 1) {
    stop("`act' should be one character");
  }
  r.s <- strsplit(gsub("\\)","",sub("^\\(","",r)),"\\(");
  rr <- sapply(r.s, function(x) {
    y <- as.numeric(sub("^.,[0-9]+,[0-9]+,[0-9]+,([0-9]+).*$","\\1",
                        x[grep(paste("^",act,",[0-9]+,",loc,",.*$",sep=""),
                               x)[1]]));
    if(length(y) == 0) {NA}
    else {y}
  });
  rr-MOL;
}

rou.u.act <- function(r,act="B",loc=0) {
  if(nchar(act) != 1) {
    stop("`act' should be one character");
  }
  r.s <- strsplit(gsub("\\)","",sub("^\\(","",r)),"\\(");
  rr <- sapply(r.s, function(x) {
    y <- as.numeric(sub("^.,[0-9]+,[0-9]+,[0-9]+,[0-9]+,([0-9.]+).*$","\\1",
                        x[grep(paste("^",act,",[0-9]+,",loc,",.*$",sep=""),
                               x)[1]]));
    if(length(y) == 0) {NA}
    else {y}
  });
  rr;
}

weekly.freq <- function(r) {
  h <- hist(r[,"week"],breaks=0:52,plot=FALSE);
  h$counts;
}

O.weekly <- function(name="baseline",y=3,ext=".dat") {
  ## read the routines individuals follow
  if(is.null(ext)) {
    f.name=name;
  } else {
    f.name <- paste(name,ext,sep="");
  }
  if(file.exists(f.name)) {
    command <- paste(.path.package("RMigone"),
                     "/scripts/weekly_ave.awk -v y=",y," ",
                     f.name,sep="");
    zz1 <- pipe(command);
    w <- read.table(file=zz1,header=FALSE,
                    col.names=c("week","exper","loc","res",
                      "fea","u","n","starv","pred"));
    w <- w[order(w$week,w$exper,w$loc),];
  } else {
    w <- NA;
  }
  invisible(w);
}

sig <- function(x,y)
{
 x1 <- floor(x*10^y+0.5)/10^y
 format(x1,nsmall=y,decimal.mark=".",trim=TRUE)
}




denz.2d <- function(x,y,rmode="ALL",dec.x=0,dec.y=0) {
  es <- !is.na(x) & !is.na(y)
  x <- sig(x,dec.x)
  y <- sig(y,dec.y)
  c.t <- table(x[es],y[es]);
  cc <- expand.grid(x=as.numeric(rownames(c.t)),y=as.numeric(colnames(c.t)));
  if(rmode=="ALL") {
    cc$count <- as.vector(c.t);
    cc$count <- cc$count/max(cc$count);
  } else {
    c.t <- c.t/rowSums(c.t)
    c.t <- t(apply(c.t,1,function(x) x/max(x)))
    cc$count <- as.vector(c.t);
  }
  cc <- cc[cc$count>0,];
  cc;
}

denz.plot <- function(x,y,rmode="ALL",dec.x=0,dec.y=0,size=3,lattice=TRUE,...){
  xx <- denz.2d(x,y,rmode=rmode,dec.x=dec.x,dec.y=dec.y);
  if(lattice) {
    l <- xyplot(y~x,xx,cex=sqrt(xx$count)*size,pch=16,...);
  } else {
    plot(xx$x,xx$y,cex=sqrt(xx$count)*size,pch=16,bty="l",...);
    l <- NA
  }
  invisible(l);
}
