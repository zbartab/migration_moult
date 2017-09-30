rand.ini <- function(nev="moult",n=100,tart=0.2){
  ini.start <- read.table(paste(nev,".ini",sep=""),col.names=c("var","value"));
  for(j in 1:n) {
    uj.nev <- paste("rand",formatC(j,width=4,flag="0"),".ini",sep="");
    cat(file=uj.nev,append=FALSE);
    for(i in 1:nrow(ini.start)) {
      if (i<36){
        cat(as.character(ini.start[i,"var"]),"\t\t",
            runif(1,min=ini.start[i,"value"]*(1.0-tart),
                  max=ini.start[i,"value"]*(1.0+tart)),
            "\n",file=uj.nev,append=TRUE,sep="");
      }
      else {
        cat(as.character(ini.start[i,"var"]),"\t\t",ini.start[i,"value"],
            "\n",file=uj.nev,append=TRUE,sep="");
      }
    }
  }
}

rand.range.ini <- function(nev="moult",n=100){
  ini.start <- read.table(paste(nev,".rini",sep=""),
                          col.names=c("var","low","high"));
##  ini.start$low <- as.numeric(as.character(ini.start$low));
##  ini.start$high <- as.numeric(as.character(ini.start$high));
  for(j in 1:n) {
    uj.nev <- paste("rr",formatC(j,width=4,flag="0"),".ini",sep="");
    cat(file=uj.nev,append=FALSE);
    for(i in 1:nrow(ini.start)) {
      if (i<36){
        cat(as.character(ini.start[i,"var"]),"\t",
            runif(1,min=ini.start[i,"low"],
                  max=ini.start[i,"high"]),
            "\n",file=uj.nev,append=TRUE,sep="");
      }
      else {
        cat(as.character(ini.start[i,"var"]),"\t\t",ini.start[i,"low"],
            "\n",file=uj.nev,append=TRUE,sep="");
      }
    }
  }
}


ini.gen <- function(a.nev="baseline") {

### read a base config file (.bini) and generate a series of ini files
### for the moult project. The bini file is like an ini file for the
### moult program except of two things. First, its first line contains
### two parameters for this function, namely basename of the
### generating ini files and whether the feathers should be symmetric
### or not (symmetric feathers mean that both type of feather has the
### same parameter value. The second exception is that several
### variable can have multiply values after the variable
### names. 'ini.gen' will generate ini files for all combination for
### these multiply values. 

  a.f.nev <- paste(a.nev,".bini",sep="");
  gen.param <-  scan(a.f.nev,what=character(0),sep="\n",quiet=TRUE,n=1);
  gen.param <- unlist(strsplit(gen.param,"[\t \n]+"));
  if(length(gen.param)!=2) {stop("Not proper parametrization of ini.gen");}
  names(gen.param) <- c("basename","fixed");
  base.config <- scan(a.f.nev,what=character(0),sep="\n",
                      quiet=TRUE,skip=1);
  base.config <- strsplit(base.config,"[\t \n]+");

  V <- numeric(38);
  
  for(V1 in base.config[[1]][-1]){
    V[1] <- V1;
  for(V2 in base.config[[2]][-1]){
    V[2] <- V2;
  for(V3 in base.config[[3]][-1]){
    V[3] <- V3;
  for(V4 in base.config[[4]][-1]){
    V[4] <- V4;
  for(V5 in base.config[[5]][-1]){
    V[5] <- V5;
  for(V6 in base.config[[6]][-1]){
    V[6] <- V6;
  for(V7 in base.config[[7]][-1]){
    V[7] <- V7;
  for(V8 in base.config[[8]][-1]){
    V[8] <- V8;
  for(V9 in base.config[[9]][-1]){
    V[9] <- V9;
  for(V10 in base.config[[10]][-1]){
    V[10] <- V10;
  for(V11 in base.config[[11]][-1]){
    V[11] <- V11;
  for(V12 in base.config[[12]][-1]){
    V[12] <- V12;
  for(V13 in base.config[[13]][-1]){
    V[13] <- V13;
  for(V14 in base.config[[14]][-1]){
    V[14] <- V14;
  for(V15 in base.config[[15]][-1]){
    V[15] <- V15;
  for(V16 in base.config[[16]][-1]){
    V[16] <- V16;
  for(V17 in base.config[[17]][-1]){
    V[17] <- V17;
  for(V18 in base.config[[18]][-1]){
    V[18] <- V18;
  for(V19 in base.config[[19]][-1]){
    V[19] <- V19;
  for(V20 in base.config[[20]][-1]){
    V[20] <- V20;
  for(V21 in base.config[[21]][-1]){
    V[21] <- V21;
  for(V22 in base.config[[22]][-1]){
    V[22] <- V22;
  for(V23 in base.config[[23]][-1]){
    V[23] <- V23;
  for(V24 in base.config[[24]][-1]){
    V[24] <- V24;
  for(V25 in base.config[[25]][-1]){
    V[25] <- V25;
  for(V26 in base.config[[26]][-1]){
    V[26] <- V26;
  for(V27 in base.config[[27]][-1]){
    V[27] <- V27;
  for(V28 in base.config[[28]][-1]){
    V[28] <- V28;
  for(V29 in base.config[[29]][-1]){
    V[29] <- V29;
  for(V30 in base.config[[30]][-1]){
    V[30] <- V30;
  for(V31 in base.config[[31]][-1]){
    V[31] <- V31;
  for(V32 in base.config[[32]][-1]){
    V[32] <- V32;
  for(V33 in base.config[[33]][-1]){
    V[33] <- V33;
  for(V34 in base.config[[34]][-1]){
    V[34] <- V34;
  for(V35 in base.config[[35]][-1]){
    V[35] <- V35;
  for(V36 in base.config[[36]][-1]){
    V[36] <- V36;
  for(V37 in base.config[[37]][-1]){
    V[37] <- V37;
  for(V38 in base.config[[38]][-1]){
    V[38] <- V38;
    out.nev <- paste(gen.param["basename"],
                     formatC(as.integer(runif(1,min=0,max=100000)),
                             width=6,flag="0"),".ini",sep="");
    out.f <- file(out.nev,"w");
    cat("",file=out.f,append=FALSE);
    for(i in 1:38) {
      cat(base.config[[i]][1],"\t",V[i],"\n",file=out.f,append=TRUE);
    }
    close(out.f);
  }}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}
}


ini.ser <- function(path=".", pattern="\\.bini$") {
  f.list <- list.files(path=path,pattern=pattern,full.names=TRUE);
  f.list <- gsub("\\.bini$","",f.list);
  for(i in f.list){
    cat(i,"...")
    ini.gen(a.nev=i);
    cat("done.\n");
  }

}
