O.txt <- function(name="baseline",ext=".txt",rajz=FALSE,offset=10) {
  if(is.null(ext)) {
    f.name <- name;
  } else {
    f.name <- paste(name,ext,sep="");
  }
  file.OK <- FALSE;
  if(file.exists(f.name)) {
    file.OK <- TRUE;
    zz <- file(f.name);
  } else {
    f.name <- paste(f.name,".gz",sep="");
    if(file.exists(f.name)) {
      file.OK <- TRUE;
      zz <- gzfile(f.name);
    } else {
      file.OK <- FALSE;
    }
  }
  if(file.OK) {
    txt <- scan(zz,what=character(0),sep="\n");
    close(zz);
    f.txt <- txt[grep("^for",txt)];
    f.txt <- gsub("[a-zA-Z,=]","",f.txt)
    f.txt <- gsub(" +"," ",f.txt)
    f.txt <- gsub("^ +","",f.txt)
    f.txt <- gsub(" +$","",f.txt)
    f.txt <- matrix(as.numeric(unlist(strsplit(f.txt," +"))),ncol=4,byrow=TRUE)
    f.txt <- f.txt[,1:2,drop=FALSE]
    b.txt <- txt[grep("^back",txt)];
    b.txt <- gsub("[a-zA-Z,=]","",b.txt)
    b.txt <- gsub(" +"," ",b.txt)
    b.txt <- gsub("^ +","",b.txt)
    b.txt <- gsub(" +$","",b.txt)
    b.txt <- matrix(as.numeric(unlist(strsplit(b.txt," +"))),ncol=4,byrow=TRUE)
    b.txt <- b.txt[b.txt[,3]<0.005,1:2,drop=FALSE]
    r <- list(backward=b.txt,forward=f.txt);
    a <- NA;
    if(rajz & nrow(b.txt) > 0) {
      if(f.txt[nrow(f.txt),1] > 10000) {
        f.txt[nrow(f.txt),1] <- f.txt[nrow(f.txt)-1,1]+10;
      }
      layout(matrix(1:4,ncol=2));
      opar <- par();
      par(mar=c(4,4,2,1));
      cim <- name;
      plot(b.txt,type="l",cex=0.5,xlab="year",ylab="lambda",main=cim);
      abline(h=1+c(-0.005,0.005))
      if(abs(b.txt[nrow(b.txt),2]-1) < 0.005) {cim <- "converged"}
      else {cim <- "not converged"}
      plot(f.txt,type="l",cex=0.5,xlab="year",ylab="population size",main=cim);
      a <- converg.bf(r,offset);
      if(is.matrix(a)) {
        plus <- sum(a[,2]>0);
        plot(a[,1],a[,2],type="o",cex=0.5,xlab="year",ylab="coeff",ylim=c(-1,1),
             main=paste("No. of pluses:",plus));
##        abline(v=a[a[,2]>0,1],lty=3,col=2);
        points(a[a[,2]>0,],pch=16,col="red");
        abline(h=0);
        aa <- list();
        if(nrow(a) > offset) {
          for(i in (offset+1):nrow(a)) {
            e <- (i-offset):i;
            aa[[i]] <- c(a[i,1],sum(a[e,2])/length(e));
          }
          if(length(aa) > 0) {
            aa <- matrix(unlist(aa),ncol=2,byrow=TRUE);
            plot(aa,xlab="year",ylab="dev",type="o",cex=0.5,ylim=c(-1,1))
            abline(h=0);
            points(aa[aa[,2]>0,],pch=16,col="red")
          }
        }
      }
      layout(1);
      par(opar);
    }
    invisible(list(backward=b.txt,forward=f.txt,a=a));
  }
}


converg.bf <- function(bf,offset=10) {
  b.txt <- bf$backward[-1,,drop=FALSE];
  r <- NA;
  if(nrow(b.txt) > offset) {
    b.txt <- cbind(b.txt[,1],abs(1-b.txt[,2]))
    r <- list();
    for(i in offset:nrow(b.txt)) {
      e <- (i-offset+1):i
      a <- cor(b.txt[e,2],1:offset);##lm(b.txt[e,2]~b.txt[e,1]);
      r[[i]] <- c(b.txt[i,1],a);##coefficients(a)[2]);
##      a <- lm(b.txt[e,2]~b.txt[e,1]);
##      r[[i]] <- c(b.txt[i,1],coefficients(a)[2]);
    }
    r <- matrix(unlist(r),ncol=2,byrow=TRUE);
  } else {
    cat("Too few observation\n");
  }
  invisible(r);
}


bf.ser <- function(y=3,m.path=".",patt="\\.txt.*$",f.list=NULL,ML=4, ask=TRUE,
                   offset=10) {
  if(is.null(f.list)) {
    f.list <- list.files(path=m.path,pattern=patt,full.names=TRUE);
    f.list <- gsub("\\.txt.*$","",f.list);
  }
  r <- list();
  for(i in f.list){
    print(i);
    r[[i]] <- O.txt(name=i,rajz=TRUE,offset=offset);
    if (ask) {readline("...");}
  }
  invisible(r);
}


bf.converged <- function(name) {
  r <- O.txt(name);
  rr <- r$backward[nrow(r$backward),];
  if(is.null(rr)) {
    rr <- rep(NA,2)
  }
  rr <- data.frame(sz=rr[1],b.lambda=rr[2]);
  rr;
}

