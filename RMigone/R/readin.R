O.bin <- function(nev="baseline",typ="uv"){

  ## this function read a binary dump of the migra1 program
  ## `nev' is the name of the file WITHOUT extension
  ## `typ' is the type of the binary file. Basicly it is the extension without
  ## the terminating `b' character

  f.nev <- paste(nev,".",typ,"b",sep="");
  zz <- file(f.nev,"rb");
  n <- readBin(zz,what="int",n=1);
  my.what <- "double";
  uv <- readBin(zz,what=my.what,n=n);
  if (typ %in% c("uv","ac")) {
    uv <- matrix(uv,ncol=4,byrow=TRUE);
  }
  close(zz);
  invisible(uv);
}

O.food <- function(nev="baseline") {

  ## this function read the amount of available food. Note this is scaled up by
  ## MR
  ## `nev' is the name of the file WITHOUT extension
  
#  i.nev <- sub("-.*$","",nev);
  i.nev <- nev;
  f.nev <- paste(i.nev,".ini",sep="");
  zz <- pipe(paste("sed -n '/^M[^.]/p' ",f.nev));
  f.ini <- read.table(zz,col.names=c("Var","Val"));
  ss <- list();
  ss[["L"]] <- 0:(f.ini[f.ini$Var=="ML","Val"]-1);
  ss[["E"]] <- 0:f.ini[f.ini$Var=="ME","Val"];
  ss[["T"]] <- 0:f.ini[f.ini$Var=="MT","Val"];
  states <- gen.states(ss);
  G <- scan(paste(i.nev,".food",sep=""));
  food <- data.frame(states,G=G);
  invisible(food);
}


G.sta <- function(nev="baseline",typ="st"){

  ## function used by others in this file
  
  i.nev <- sub("-.*$","",nev);
  f.nev <- paste(i.nev,".ini",sep="");
  zz <- pipe(paste("sed -n '/^M[^.]/p' ",f.nev));
  f.ini <- read.table(zz,col.names=c("Var","Val"));
  ss <- list();
  ss[["L"]] <- 0:(f.ini[f.ini$Var=="ML","Val"]-1);
  ss[["A"]] <- 0:f.ini[f.ini$Var=="MA","Val"];
  ss[["E"]] <- 0:f.ini[f.ini$Var=="ME","Val"];
  ss[["F"]] <- -(mol <- f.ini[f.ini$Var=="MOL","Val"]):(
                                f.ini[f.ini$Var=="MF","Val"]-mol);
  ss[["R"]] <- 0:f.ini[f.ini$Var=="MR","Val"];
  if(typ %in% c("st","rv")) {
    ss[["T"]] <- 0:f.ini[f.ini$Var=="MT","Val"];
  } else {
    ss[["T"]] <- 0:(f.ini[f.ini$Var=="MT","Val"]-1);
  }
  states <- gen.states(ss);
  invisible(states);
}

gen.states <- function (...) {
  
  ## function used by others in this file
  
  nargs <- length(args <- list(...))
  if (nargs == 1 && is.list(a1 <- args[[1]])) 
    nargs <- length(args <- a1)
  cargs <- args
  nmc <- paste("Var", 1:nargs, sep = "")
  nm <- names(args)
  if (is.null(nm)) 
    nm <- nmc
  if (any(ng0 <- nchar(nm) > 0)) 
    nmc[ng0] <- nm[ng0]
  names(cargs) <- nmc
  rep.fac <- 1
  orep <- final.len <- prod(sapply(args, length))
  for (i in 1:nargs) {
    x <- args[[i]]
    nx <- length(x)
    orep <- orep/nx
    x <- rep(rep(x, rep(rep.fac, nx)), orep)
    if (!is.factor(x) && is.character(x)) 
      x <- factor(x, levels = unique(x))
    cargs[[i]] <- x
    rep.fac <- rep.fac * nx
  }
  invisible(cargs);
}

O.all <- function(nev="baseline") {

  ## read all binary dump file created by migra1
  ## it returns a list of those supplemented by the state variables
  ## `nev' is the name of the file WITHOUT extension
  
  sr <- G.sta(nev=nev,typ="rv");
  sa <- G.sta(nev=nev,typ="ac");
  st <- O.bin(nev=nev,typ="st");
  rv <- O.bin(nev=nev,typ="rv");
  st <- st[sr$T!=52];
  rv <- rv[sr$T!=52];
  rm(sr);
  ac <- O.bin(nev=nev,typ="ac");
  uv <- O.bin(nev=nev,typ="uv");
  mo <- O.bin(nev=nev,typ="mo");
  sa$uv <- uv;
  sa$ac <- ac;
  sa$mo <- mo;
  sa$rv <- rv;
  sa$st <- st;
  invisible(sa);
}

