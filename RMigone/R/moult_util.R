f.nev <- function(val, var, df) {
  unique(rownames(df[df[,var]==val,]));
}


flight.eff <- function(f1,f2,P) {
  dfe1 <- P$f1$worn + f1^P$f1$flight*(1-P$f1$worn);
  dfe2 <- P$f2$worn + f2^P$f2$flight*(1-P$f2$worn);
  dfe <- P$f1$def.coeff*dfe1 + P$f2$def.coeff*dfe2 + P$dfeint*dfe1*dfe2;
  dfe <- dfe/(P$f1$def.coeff+P$f2$def.coeff+P$dfeint);
  dfe;
}

m.cost <- function(r,u,f1,f2,P) {
  fora.dfe <- (1-P$propF)*P$fora + P$propF*(P$def/flight.eff(f1,f2,P));
  mc <- P$base + P$mass * r^2 + u^2*(1+0.1*r^2)*fora.dfe;
  mc;
}

r.P <- function(P.worn=0.18,S.worn=0.18,P.flight=0.6,S.flight=0.6,
                P.def.coeff=1.5,S.def.coeff=1.5,dfeint=1,propF=0.5,fora=0.2,
                def=0.75,base=0.3,mass=0.03,FOOD=1.55,MAXEPS=0.85) {
  P <- list(f1=list(worn=P.worn,flight=P.flight,def.coeff=P.def.coeff),
            f2=list(worn=S.worn,flight=S.flight,def.coeff=S.def.coeff),
            dfeint=dfeint,propF=propF,fora=fora,def=def,base=base,mass=mass,
            FOOD=FOOD,MAXEPS=MAXEPS);
  P;
}

m.gain <- function(u, tt=13,P) {
  g <- P$FOOD+(P$MAXEPS*sin(2*pi*tt/52-pi/2));
  g*u;
}

p.metabol <- function(P){
  uu <- seq(0,1,0.1);
  plot(1,1,xlim=c(0,1.2),ylim=c(0,4),xlab="u",
       ylab="cost or gain",type="n");
  lines(uu,m.gain(uu,0,P),lty=2,col="red",lwd=0.5)
  lines(uu,m.gain(uu,13,P),lty=2,col="red",lwd=1)
  lines(uu,m.gain(uu,26,P),lty=2,col="red",lwd=2)
  lines(uu,m.cost(0,uu,0,0,P),lty=1,col="black",lwd=1)
  lines(uu,m.cost(1,uu,0,0,P),lty=1,col="black",lwd=1)
  lines(uu,m.cost(0,uu,0.1,0.1,P),col="green",lwd=2,lty=2)
  lines(uu,m.cost(1,uu,0.1,0.1,P),col="green",lwd=2,lty=2)
  lines(uu,m.cost(0,uu,1,1,P),col="blue",lwd=2,lty=3)
  lines(uu,m.cost(1,uu,1,1,P),col="blue",lwd=2,lty=3)
  lines(c(1.05,1.1),c(0.75,0.75),col="red",lty=2)
  lines(c(1.05,1.1),c(1.5,1.5),col="black",lty=1)
  lines(c(1.05,1.1),c(2.25,2.25),col="green",lwd=2,lty=2)
  lines(c(1.05,1.1),c(3,3),col="blue",lwd=2,lty=3)
  text(1.1,0.75,labels="gain",pos=4);
  text(1.1,1.5,labels="bad fea",pos=4);
  text(1.1,2.25,labels="poor fea",pos=4);
  text(1.1,3,labels="good fea",pos=4);
}
