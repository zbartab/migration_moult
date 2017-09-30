
s.runner <- function(t1, res = 0, f1 = 0, f2 = 0, a.o = 0, a.t = 0,food = 1,
                     t2 = t1, nev = "baseline", ext = NULL, y1 = 2, y2 = y1){
  pars <- "";
  if(is.matrix(t1)) {
    for(i in 1:nrow(t1)) {
      pars <- paste(pars,"-t",paste(t1[i,],collapse=" "));
    }
  } else {
    pars <- paste("-t",y1*52+t1,y2*52+t2,res,f1,f2,a.o,a.t,food);
  }
  if (is.null(ext)) {
    ext <- formatC(as.integer(runif(1,min=0,max=100000)),width=6,flag="0");
  }
  comm <- paste("./simrunner",pars,nev,ext);
  print(comm);
  system(comm);
  ext;
}
