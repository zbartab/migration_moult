

create.P <- function(P=NULL,...) {
  n.valt <- names(valt <- list(...));
  if(is.null(P)) {
    P <- list(F.fixed=0,
              F.fora=0.03,
              F.flight=0.6,
              F.worn=0.2,
              F.NU=0.15,
              M.feather=0.05,
              M.energy=0.3,
              M.preda=0.001,
              C.fora=0.2,
              C.propF=0.5,
              C.def=0.75,
              C.birtok=1.0,
              C.start=1.1,
              C.keep=1.1,
              C.moult=0.85,
              C.mass=0.04,
              C.baseM=0.3,
              C.baseS=0.0,
              P.base=0.0005,
              P.mass=0.002,
              L.THETA=0.7,
              L.MAXEPS=0.95,
              L.FOOD=1.6,
              L.p=0.5,
              L.omega=0.01,
              L.n=3,
              L.N0=1000,
              L.DD=0.1,
              L.pe=0.025,
              L.nbirtok=500,
              L.K=0.9,
              L.al=1.0,
              Mtol=0.005,
              Malpha=0.25,
              MK=0.0005,
              Msimit=0.75,
              MR=12,
              MT=52,
              MF=18,
              MOL=8,
              MA=8,
              ME=2,
              ML=4,
              MNP=7);
  }
  for(i in n.valt) {
    if(i %in% names(P)) {
      P[[i]] <- valt[[i]];
    } else {
      stop(paste("Setting unused variable: ",i));
    }
  }
  for(i in names(P)) {
    if(regexpr("^M[a-zA-Z]",names(P[i]))==-1) {
      P[[i]] <- rep(P[[i]],length.out=P[["ML"]]);
    } 
  }
  invisible(P);
}


write.P <- function(P,nev="baseline") {
  f.nev <- paste(nev,".ini",sep="");
  cat("",file=f.nev);
  for(i in 1:length(P)) {
    cat(names(P[i]),"\t\t",paste(P[[i]],collapse="\t"),"\n",
        sep="",file=f.nev,append=TRUE);
  }
}
