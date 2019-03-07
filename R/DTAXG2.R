DTAXG2=function(
  group1,
  group2,
  prior.se.group1,
  prior.sp.group1,
  prior.se.group2,
  prior.sp.group2,
  prior.pi,
  n.sample,
  n.burnin,
  SUM=TRUE){

  group=data.frame(group1,group2)
  U1=nrow(group[which(group1==1 & group2==1),])
  U2=nrow(group[which(group1==1 & group2==0),])
  U3=nrow(group[which(group1==0 & group2==1),])
  U4=nrow(group[which(group1==0 & group2==0),])

  if(U1==0){U1=1}
  if(U2==0){U2=1}
  if(U3==0){U3=1}
  if(U4==0){U4=1}

  muPI<-sum(prior.pi)/2
  varPI<-(diff(prior.pi)/4)^2
  alphaPI<-muPI*((1-muPI)*muPI/varPI-1)
  betaPI<-(1-muPI)*((1-muPI)*muPI/varPI-1)

  muS1<-sum(prior.se.group1)/2
  varS1<-(diff(prior.se.group1)/4)^2
  alphaS1<-muS1*((1-muS1)*muS1/varS1-1)
  betaS1<-(1-muS1)*((1-muS1)*muS1/varS1-1)

  muS2<-sum(prior.se.group2)/2
  varS2<-(diff(prior.se.group2)/4)^2
  alphaS2<-muS2*((1-muS2)*muS2/varS2-1)
  betaS2<-(1-muS2)*((1-muS2)*muS2/varS2-1)

  muC1<-sum(prior.sp.group1)/2
  varC1<-(diff(prior.sp.group1)/4)^2
  alphaC1<-muC1*((1-muC1)*muC1/varC1-1)
  betaC1<-(1-muC1)*((1-muC1)*muC1/varC1-1)

  muC2<-sum(prior.sp.group2)/2
  varC2<-(diff(prior.sp.group2)/4)^2
  alphaC2<-muC2*((1-muC2)*muC2/varC2-1)
  betaC2<-(1-muC2)*((1-muC2)*muC2/varC2-1)

  inits=c(X1=U1,X2=U2,X3=U3,X4=U4)
  parameters=c("PI","S1","S2","C1","C2")

  mchain=data.frame(matrix(rep(NA,(n.sample-n.burnin)*(length(inits)+length(parameters))),nrow=(n.sample-n.burnin)))
  names(mchain)=c(names(inits),parameters)
  mchain[1,names(inits)]=inits

  for(i in 1:(n.sample-n.burnin)){
    X1=mchain[i,"X1"]
    X2=mchain[i,"X2"]
    X3=mchain[i,"X3"]
    X4=mchain[i,"X4"]

    dbeta1PI<-X1+X2+X3+X4+alphaPI
    dbeta2PI<-(U1+U2+U3+U4)-(X1+X2+X3+X4)+betaPI
    dbeta1S1<-X1+X2+alphaS1
    dbeta2S1<-X3+X4+betaS1
    dbeta1C1<-U1+U2-(X3+X4)+alphaC1
    dbeta2C1<-U3+U4-(X1+X2)+betaC1
    dbeta1S2<-X1+X3+alphaS2
    dbeta2S2<-X2+X4+betaS2
    dbeta1C2<-U1+U3-(X2+X4)+alphaC2
    dbeta2C2<-U2+U4-(X1+X3)+betaC2

    mchain[i,"PI"]<-suppressWarnings(rbeta(n=1,abs(dbeta1PI),abs(dbeta2PI)))
    mchain[i,"S1"]<-suppressWarnings(rbeta(n=1,abs(dbeta1S1),abs(dbeta2S1)))
    mchain[i,"C1"]<-suppressWarnings(rbeta(n=1,abs(dbeta1C1),abs(dbeta2C1)))
    mchain[i,"S2"]<-suppressWarnings(rbeta(n=1,abs(dbeta1S2),abs(dbeta2S2)))
    mchain[i,"C2"]<-suppressWarnings(rbeta(n=1,abs(dbeta1C2),abs(dbeta2C2)))

    PI=mchain[i,"PI"]
    S1=mchain[i,"S1"]
    C1=mchain[i,"C1"]
    S2=mchain[i,"S2"]
    C2=mchain[i,"C2"]

    dbinX1<-PI*S1*S2/(PI*S1*S2+(1-PI)*(1-C1)*(1-C2))
    dbinX2<-PI*S1*(1-S2)/(PI*S1*(1-S2)+(1-PI)*(1-C1)*C2)
    dbinX3<-PI*(1-S1)*S2/(PI*(1-S1)*S2+(1-PI)*C1*(1-C2))
    dbinX4<-PI*(1-S1)*(1-S2)/(PI*(1-S1)*(1-S2)+(1-PI)*C1*C2)

    if(i+1<nrow(mchain)){
      mchain[i+1,"X1"]<-rbinom(n=1,U1,dbinX1)
      mchain[i+1,"X2"]<-rbinom(n=1,U2,dbinX2)
      mchain[i+1,"X3"]<-rbinom(n=1,U3,dbinX3)
      mchain[i+1,"X4"]<-rbinom(n=1,U4,dbinX4)
    }
  }
  quantile2=function(x){quantile(x,probs=c(0.5,0.025,0.975),na.rm=T)}
  quantiles=apply(mchain[,parameters],2,"quantile2")
  if(SUM){return(quantiles)}else(return(mchain))
}

