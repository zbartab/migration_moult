moult.brood.man <- function(path=".", pattern="\\.dat$",y=2,f.list=NULL,MA=7) {
  if(is.null(f.list)) {
    f.inik <- list.files(path=path,pattern=pattern,full.names=TRUE);
  } else {
    f.inik <- f.list;
  }
  r <- list(0);
  start <- TRUE;
  for(f.nev in f.inik) {
    print(f.nev);
    ds <- O.sim(nev=sub("\\.dat$","",f.nev),y=2);
    y.man <- as.numeric(gsub("^.*-\\(.*\\)\\.dat$","\\1",f.nev));
    b.bird <- ds$bird %in% ds$bird[ds$action=="START"];
    m.bird <- ds$bird %in% ds$bird[ds$action=="START" & ds$week<y.man
                                   & ds$week >= (y.man-MA)];
    r.bird <- ds$bird %in% ds$bird[m.bird & ds$action=="START"
                                   & ds$week>=y.man];
    pr.manip <- length(unique(ds$bird[m.bird]))/length(unique(ds$bird[b.bird]));
    pr.renest <- length(unique(ds$bird[r.bird]))/
      length(unique(ds$bird[m.bird]));
    if(sum(m.bird) > 0) {
      ds <- data.frame(ds);
      m.typ <- sapply(unique(ds$bird[m.bird]),
                      function(i) moult.types2(ds[ds$bird==i,]));
      m1.time <- mean(ds$week[m.bird&ds$moult1=="START"])
      m2.time <- mean(ds$week[m.bird&ds$moult2=="START"])
      br.time <- sapply(unique(ds$bird[m.bird]),
                             function(i) {
                               min(ds$week[ds$bird==i & ds$week>=y.man &
                                           ds$action!="START" &
                                           ds$action!="KEEP"])
                             });
      br.time <- mean(br.time[is.finite(br.time)]);
    } else {
      m.typ <- list();
      m1.time <- NA;
      m2.time <- NA;
      br.time <- NA;
    }
    if (length(m.typ)==0) {
      pr.mo <- c(sp=NA,ps=NA);
    } else {
      pr.mo <- moult.order(table(m.typ));
    }
    if(start) {
      start <- FALSE;
      r <- data.frame(pr.manip=pr.manip,pr.renest=pr.renest,ps=pr.mo["ps"],
                      sp=pr.mo["sp"],m1.time=m1.time,m2.time=m2.time,
                      br.time=br.time);
    } else {
      r <- rbind(r,c(pr.manip,pr.renest,pr.mo["ps"],pr.mo["sp"],m1.time,
                     m2.time,br.time));
    }
  }
  rownames(r) <- f.inik;
  r;
}


