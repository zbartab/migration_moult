
table.item <- function(papd,bl.papd,var,case,bl.val,numbers=FALSE,
                       limits=c(0.1,0.5,Inf),p.marks=c("","+","++"),
                       n.marks=c("","$-$","$--$"))
{
  val.str <- sub("\\.","_",as.character(bl.val))
  bl.papd$file <- sub("baseline-sim",
                      paste(case,"_baseline-",val.str,"-sim",sep=""),
                      bl.papd$file)
  p <- rbind(papd,bl.papd)
  r.bl <- ex.papd(var,bl.papd,"baseline")
  r <- ex.papd(var,p,case)
  if(sum(!is.na(r$r)) > 0) {
    lm.lo <- try(lm(r$r~r$v,subset=r$v <= bl.val),TRUE)
    if(!inherits(lm.lo,"try-error")) {
      pr.lo <- predict(lm.lo)
      pr.lo <- pr.lo[length(pr.lo)] - pr.lo[1] 
      for(i in 1:length(limits)) {
        if (abs(pr.lo) < limits[i]) {
          if(pr.lo<0) m.lo <- n.marks[i] else m.lo <- p.marks[i]
          break
        }
      }
    } else {
      m.lo <- "NA"
    }
    lm.hi <- try(lm(r$r~r$v,subset=r$v >= bl.val),TRUE)
    if(!inherits(lm.hi,"try-error")) {
      pr.hi <- predict(lm.hi)
      pr.hi <- pr.hi[length(pr.hi)] - pr.hi[1]
      for(i in 1:length(limits)) {
        if (abs(pr.hi) < limits[i]) {
          if(pr.hi<0) m.hi <- n.marks[i] else m.hi <- p.marks[i]
          break
        }
      }
    } else {
      m.hi <- "NA"
    }
    if(numbers & regexpr("(\\.br|\\.mo)",var) > -1) {
      t.item <- paste(m.lo,round(r.bl$r,2),m.hi,sep=" & ")
    } else {
      t.item <- paste(m.lo,round(100*r.bl$r,0),m.hi,sep=" & ")
    }
  } else {
    t.item <- " & NA & "
  }
#  list(lo=pr.lo,hi=pr.hi,data=r,t.item=t.item)
  t.item
}


f.table.item <- function(p,var,case,index,numbers=FALSE,
                         limits=c(0.1,0.5,Inf),p.marks=c("","+","++"),
                         n.marks=c("","$-$","$--$"))
{
  m.hi <- NA
  r.bl <- ex.papd(var,p,paste(case,"-0-sim",sep=""))
  r <- ex.papd(var,p,case)
  r$v <- r$v[index]
  r$r <- r$r[index]
  lm.hi <- lm(r$r~r$v)
  pr.hi <- predict(lm.hi)
  pr.hi <- pr.hi[length(pr.hi)] - pr.hi[1]
  for(i in 1:length(limits)) {
    if (abs(pr.hi) < limits[i]) {
      if(pr.hi<0) m.hi <- n.marks[i] else m.hi <- p.marks[i]
      break
    }
  }
  t.item <- paste(round(100*r.bl$r,0),m.hi,sep=" & ")
  t.item
#  list(hi=pr.hi,data=r,t.item=t.item)
}

n.table.item <- function(p,var,case,index,cut.n=150,
                         limits=c(0.1,0.5,Inf),p.marks=c("","+","++"),
                         n.marks=c("","$-$","$--$"))
{
  if(regexpr("nwb",case) != -1) {m.hi <- NA} else {m.hi <- ""}
  r <- ex.papd(var,p,case)
  n <- ex.papd("^n$",p,case)
  r$v <- r$v[index]
  r$r <- r$r[index]
  n$v <- n$v[index]
  n$r <- n$r[index]
  r$v <- r$v[n$r >= cut.n]
  r$r <- r$r[n$r >= cut.n]
  if(sum(!is.na(r$r)) > 0) {
    lm.hi <- lm(r$r~r$v)
    pr.hi <- predict(lm.hi)
    pr.hi <- pr.hi[length(pr.hi)] - pr.hi[1]
    for(i in 1:length(limits)) {
      if (abs(pr.hi) < limits[i]) {
        if(pr.hi<0) m.hi <- n.marks[i] else m.hi <- p.marks[i]
        break
      }
    }
  }
  m.hi
}

short.item <- function(p,var,case, symb,sign = "+",
                       index=TRUE, limits=c(0.1,0.5,Inf),
                       b.mark=c("","$","$\\\\mathbf{"),
                       e.mark=c("","$","}$")
                       )
{
  m.hi = ""
  symb = c("",rep(symb,length(limits)-1))
  r <- ex.papd(var,p,case)
  r$v <- r$v[index]
  r$r <- r$r[index]
  if(sum(!is.na(r$r)) > 0) {
    lm.hi <- lm(r$r~r$v)
    pr.hi <- predict(lm.hi)
    pr.hi <- pr.hi[length(pr.hi)] - pr.hi[1]
    for(i in 1:length(limits)) {
      if (abs(pr.hi) < limits[i]) {
        if(sign == "+" & pr.hi > 0) {
          m.hi = paste(b.mark[i],symb[i],e.mark[i],sep="")
        } else if (sign == "-" & pr.hi < 0) {
          m.hi = paste(b.mark[i],symb[i],e.mark[i],sep="")
        }
        break
      }
    }
  }
  m.hi
}

sh.cell <- function(papd,v,prefix="",sign="+",index=TRUE) {
  cases <- c("C_mass","C_moult","F_flight","F_fora","F_worn",
             "M_energy","M_feather","M_preda","P_base","P_mass")
  cases <- paste(prefix,cases,sep="")
  symb <- c("c_r","\\\\kappa","\\\\alpha","f_f","m_A","C_{mig}",
            "f_m","M_m","M_b","M_f")
  r <- ""
  for(i in 1:length(cases)) {
    s <- short.item(papd,v,cases[i],symb[i],sign=sign,index=index)
    if(s != "") {
      if(r == "") {
        r <- s
      } else {
        r <- paste(r,s,sep=", ")
      }
    }
  }
  paste("\\\\parbox[c]{0.18\\\\textwidth}{\\\\rule{0pt}{3ex}\\\\centering ",r,"}",sep="")
}

f.sh.cell <- function(papd,v,prefix="",sign="+",index=TRUE) {
  
  cases <- paste("food_",1:4,sep="")
  symb <- paste("S",1:4,sep="")
  if(regexpr("w_",prefix) == -1) {
    cases <- c(cases,"food_all")
    symb <- c(symb,"A")
  }
  cases <- paste(prefix,cases,sep="")
  b.mark <- c("","","\\\\textbf{")
  e.mark <- c("","","}")
  r <- ""
  for(i in 1:length(cases)) {
    if(regexpr("all",cases[i]) > -1) {
      index <- 1:21
    } else {
      index <- 1:11
    }
    s <- short.item(papd,v,cases[i],symb[i],sign=sign,index=index,b.mark=b.mark,
                    e.mark=e.mark)
    if(s != "") {
      if(r == "") {
        r <- s
      } else {
        r <- paste(r,s,sep=", ")
      }
    }
  }
  paste("\\\\parbox[c]{0.18\\\\textwidth}{\\\\rule{0pt}{3ex}\\\\centering ",r,"}",sep="")
}

table.row <- function(papd,case,numbers=TRUE,index=TRUE,cut.n=150) {
  vars <- c("p.mo0","p.mo3","mo.M","mo.B","p.br0","p.br3","p.s0","p.w3")
  if(numbers) {
    vars <- sub("^p\\.([b])","na.\\1",vars) # change [b] -> [bm] if
                                        # you want average numbers for moult too
  }
  r <- ""
  for(v in vars) {
    r <- paste(r,n.table.item(papd,v,case,index,cut.n),sep=" & ")
  }
  r <- paste(r,"\\\\\\\\")
  r
}

f.table.row <- function(papd,case,numbers=FALSE) {
  if(regexpr("all",case) > -1) i <- 1:21 else i <- 1:11
  vars <- c("p.s0","p.w3","p.br0","p.br3","p.mo0","p.mo3","mo.B","mo.M")
  if(numbers) {
    vars <- sub("^p\\.","na.",vars)
  }
  r <- ""
  for(v in vars) {
    r <- paste(r,f.table.item(papd,v,case,i,numbers),sep=" & ")
  }
  r <- paste(r,"\\\\\\\\")
  r
}
