beh.timing.ser <- function(y=2,m.path=".",patt="\\.dat$",f.list=NULL){
  
  ## data processing for dat files
  
  if(is.null(f.list)) {
    f.list <- list.files(path=m.path,pattern=patt,full.names=TRUE);
    f.list <- gsub("\\.dat$","",f.list);
  }
  START <- TRUE;
  for(i in f.list){
    cat("-----------------------------------------\n");
    cim <- paste("File:", i);
    cat(cim,"\n");
    system(paste("zcat ", i,".sum.gz|", "grep '^ERR' ",sep=""));
    if(START) {
      r <- beh.timing(name=i,y=y);
      START <- FALSE;
    } else {
      r <- rbind(r,beh.timing(name=i,y=y));
    }
  }
  row.names(r) <- 1:nrow(r);
  invisible(r);
}


beh.timing.read <- function(name="baseline",y=3) {
  f.name <- paste(name,".dat",sep="");
  zz <- pipe(paste(.path.package("RMigone"),"/scripts/behav_time.awk ",
                   f.name,sep=""));
  r <- read.table(zz,col.names=c("behav","week","exper","loc","count"));
  r;
}

beh.timing <- function(name="baseline",y = 3) {
  r <- beh.timing.read(name,y);
  exper <- min(r$exper):max(r$exper);
  loc <- min(r$loc):max(r$loc);
  behav <- unique(r$behav);
  s <- expand.grid(behav=behav,exper=exper,loc=loc);
  p <- apply(s,1,function(x) {
    beh.timing.mean(r,r$behav==x[1] & r$exper==x[2] & r$loc==x[3]);
  })
  r <- cbind(s,t(p));
  r$file <- rep(name,nrow(r));
  r;
}

beh.timing.mean <- function(df,eset=TRUE) {
  x <- df$week[eset];
  x <- x*(pi/26);
  w <- df$count[eset];
  sinr <- sum(w*sin(x));
  cosr <- sum(w*cos(x));
  circmean <- atan(sinr,cosr);
  r <- (circmean/pi)*26;
  if(r < 0) {r <- 52+r};
  r <- c(r,sum(w));
  names(r) <- c("circ.mean","n");
  r;
}
