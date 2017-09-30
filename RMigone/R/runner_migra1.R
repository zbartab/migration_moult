

run.migra1 <- function(P=NULL,do.plot=TRUE,f.name="run_migra1",...) {

  ## run migra1 one with arbitary parameters
  
  if(!require(RMigone,quietly=TRUE)) {
    stop("No RMigone package can be loaded");
  }

  if(is.null(P)) {
    P <- create.P(F.fixed=0,
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
                  L.MAXEPS=c(1.25,0.7,0,-0.3),
                  L.FOOD=1.5,
                  L.p=0.5,
                  L.omega=0.01,
                  L.n=3,
                  L.N0=1000,
                  L.pe=0.025,
                  L.nbirtok=500,
                  L.K=0.9,
                  L.al=2.0,
                  Mtol=0.005,
                  Malpha=0.25,
                  MK=0.0005,
                  Msimit=0.75,
                  MR=12,
                  MT=52,
                  MF=23,
                  MOL=13,
                  MA=8,
                  ME=2,
                  ML=4,
                  MNP=7);
  }
  P <- create.P(P,...);
  write.P(P,f.name);
  command <- paste("./migra1 ",f.name," > ",f.name,".txt",sep="");
  system(command);
  if(do.plot) {plot.sim1(nev=f.name);}
}
