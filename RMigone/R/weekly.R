
########
## routines to process the results of the weekly_ave.awk script
########


O.weekly <- function(name="baseline",ext=".week") {
  ## read the weekly average of state variables and counts
  if(is.null(ext)) {
    f.name=name;
  } else {
    f.name <- paste(name,ext,sep="");
  }
  if(file.exists(f.name)) {
    r <- read.table(file=f.name,col.names=c("year","week","exper","loc","res",
                               "fea","u","n","starv","pred","m.starv",
                               "m.pred"));
    r <- r[order(r$year,r$week,r$loc,r$exper),];
  } else {
    r <- NA;
  }
  invisible(r);
}

  
