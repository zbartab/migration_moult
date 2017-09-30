
loc.ab <- function(nev="baseline", y=3) {

  ## Return the number of birds at each location at mid summer and at
  ## midwinter of year `y'. Only those birds in summer counted who are
  ## still alive at the winter count.


  f.name <- paste(nev,".dat",sep=""); 
  command <- paste(.path.package("RMigone"),"/scripts/count.awk -v y=",y," ",
                   f.name,sep="");
  zz1 <- pipe(command);
  df <- read.table(file=zz1,col.names=c("loc","summ","wint"));
  df;
}

loc.ab.ser <- function(path=".",patt="\\.dat$",f.list=NULL,y=3) {
  if(is.null(f.list)) {
    f.list <- list.files(path=path,pattern=patt,full.names=TRUE);
    f.list <- gsub("\\.dat$","",f.list);
  }
  r <- list();
  for(i in f.list){
    r[[i]] <- loc.ab(nev=i, y=y);
  }
  r;
}

loc.ab.plot <- function(r) {
  
}
