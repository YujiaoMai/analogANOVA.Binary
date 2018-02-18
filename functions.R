############################################################################
tlog <- function(V){
  logs <- rep(NA,length(V))
  idx0 <- which(V==0)
  if(length(idx0)>=1) {
    logs[idx0]=0 
    idx1 <- which(is.na(logs))
    if(length(idx1)>=1) logs[idx1]<-log(V[idx1])
  } else logs <-log(V)
  logs
}
## begin power function
chisq.power<-function(chi, df=1, alpha=.05){
  if(missing(chi))stop()
  crit.val <-  qchisq(p=1-alpha,df=df,ncp=0)
  power <- 1-pchisq(crit.val,df=df,ncp=chi)
  power
}
effect.size.phi.calculator <- function(mean, weight){
  Pi <- mean
  W <-weight
  if(length(Pi)!=length(W))stop()
  p <- sum(diag(W)%*%Pi)
  temp <- (Pi*tlog(p)+(1-Pi)*tlog(1-p)) - (Pi*tlog(Pi)+(1-Pi)*tlog(1-Pi))
  phi.square <-  -2*sum(diag(W)%*%temp)
  if(phi.square<0) stop(paste("Phi.square negative : ", phi.square, " p1=",Pi[1]," p2=",Pi[2], " p3=",Pi[3], sep='' ) )
  sqrt(phi.square)
}
onewayANOVA.binary.power.phi <- function(effect.size, sample.size, ngroup, alpha=.05){
  phi <- effect.size
  n <- sample.size 
  SS1 <- chi <- (phi)^2*n
  df <- ngroup-1
  chisq.power(chi=SS1,df=df,alpha=alpha)
}
## end of power funtion
############################################################################
## begin ANOVA function
onewayANOVA.binary <- function(y,A,alpha=.05){
  if(length(y)!=length(A))stop()
  n <- length(y)
  p1 <- mean(y)
  N <- table(A)
  ngroup <- length(N)
  ps <- tapply(y, INDEX=A, FUN = mean); if(length(ps)!=ngroup)stop()
  P <- ps
  SS0 <- -2*n*(p1*tlog(p1)+(1-p1)*tlog(1-p1))
  SS1 <- -2*sum(diag(N)%*%((P*tlog(p1)+(1-P)*tlog(1-p1)) - (P*tlog(P)+(1-P)*tlog(1-P))) )
  df1 <- ngroup-1
  crit.val1 <- qchisq(p=1-alpha,df=df1)
  sig <- SS1 >crit.val1
  power <- chisq.power(chi=SS1,df=df1,alpha=alpha)
  phi <- sqrt(SS1/n)
  R.square <- SS1/SS0
  SS2 <- SS0-SS1
  ANOVA.table <- cbind(Model=c('M0-M1','Mf-M1','M0'),
                       Source=c('Between-group','Within-group','Intercept'),
                       SSS=c(SS1,SS2,SS0),
                       df=c(ngroup-1,n-ngroup,n-1))
  list(ANOVA.table=ANOVA.table, sig=sig, power=power,effect.size=phi, R.square =R.square)
}
####################################
## begin essential functions           
####################################
## variation between x and pi , x =1 or 0, pi from 0 to 1
## a variation function S(x,pi)
## S1()  loglikelihood variation
S1 <- function(x,pi){
  if(!(pi>=0 & pi<=1))stop("pi is out of its domain")
  if(pi==0) return(0) 
  if(x==1) -2*log(pi) else if(x==0) -2*log(1-pi)
}
S <- S1
## total variation of X conditional on Pi: 
SS <- function(X,Pi){
  sum(apply(cbind(x=X,pi=Pi),1,function(r){S(r['x'],r['pi'])}))
}
####################################
## variation between different p/pi,  p from 0 to 1, pi from 0 to 1
s <- function(pi){
  S(1,pi)
}
fei <- function(pi){
  if(!(pi>=0 & pi<=1))stop("pi is out of its domain")
  1-pi
}
sigma <- function(p){
  p*s(p)+fei(p)*s(fei(p))
}
sigma.pi <- function(p){
  p*s(pi)+fei(p)*s(fei(pi))
  # equivalently : S(p,pi)  
}
S.p.pi <- function(p,pi){
  p*s(pi)+fei(p)*s(fei(pi))-sigma(p)
}
SS.P.Pi <- function(P,Pi){
  sum(unlist( apply(cbind(p=P,pi=Pi),1,function(r){S.p.pi(r['p'],r['pi'])})) )
}
####################################
## end of essential functions           
####################################
## end of ANOVA function
############################################################################
#  logistic regression for binary response and categorical predictor
logistic.ANOVA<-function(y,A){
  A.factor <- as.factor(A)
  m1 <- glm(y ~ A.factor, binomial(link = "logit"))
  m0 <- glm(y ~ 1, binomial(link = "logit"))
  anova(m0,m1, test="Chisq")
}

#####################################
## simulation functions
MCpower.single<- function(ngroup,N,Pi,factor.A,R=1000,alpha=.05){
  if(missing(ngroup)) stop("missing argument ngroup")
  if(missing(N)) stop("missing arg sample size N")
  if(missing(Pi)) stop("missing arg Pi")
  if(length(N)!=ngroup) stop("the number of group sample sizes N not equal to the arg ngroup")
  if(length(Pi)!=ngroup) stop("the number of group means Pi not equal to the arg ngroup")
  if(missing(factor.A)) factor.A <- seq(1,ngroup,by=1)
  n <- sum(N)
  # MC power
  sig <- rep(NA,R)
  sig.logistic <- rep(NA,R)
  Dev.mat <- matrix(NA,nrow=R,ncol=4)
  colnames(Dev.mat) <- c("ANOVA.Dev","ANOVA.Df","Logistic.Dev","Logistic.Df")
  for(i in 1:R){
    ####First Step: Generate data
    ## Dset when Alternative Model is true
    dset <- NULL
    temp <- apply(cbind(Pi,N,factor.A),1,function(r){
      p.j <- r[1]; n.j <- r[2]; a.j <- r[3]
      list(y=rbinom(n=n.j,size=1,prob=p.j), A =rep(a.j,n.j))
    })
    temp <- do.call(rbind, temp)
    dset <- cbind(y = unlist( temp[,'y']), A = unlist( temp[,'A']))
    ####Second Step: Fit model
    res <- onewayANOVA.binary(dset[,'y'],dset[,'A'],alpha=.05)
    sig[i] <- res$sig
    Dev.mat[i,'ANOVA.Dev'] <- res$ANOVA.table[1,'SSS']
    Dev.mat[i,'ANOVA.Df']  <- res$ANOVA.table[1,'df']
    res.logistic <- logistic.ANOVA(dset[,'y'],dset[,'A']) 
    sig.logistic[i] <- res.logistic[2,'Deviance']>qchisq(1-alpha,df=res.logistic[2,'Df'])
    Dev.mat[i,'Logistic.Dev']  <- res.logistic[2,'Deviance']
    Dev.mat[i,'Logistic.Df'] <-res.logistic[2,'Df'] 
  }
  
  MC.power.ANOVA<- sum(sig)/R
  MC.power.logistic <- sum(sig.logistic)/R  
  list(MC.power.ANOVA=MC.power.ANOVA,MC.power.logistic=MC.power.logistic,Dev.mat=Dev.mat)
}
##  end of simulation functions
#####################################

