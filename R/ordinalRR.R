## Control Function
ordinalRR.control<-function (mu.mu.alpha = 0.8, tau.mu.alpha = 0.4, mu.tau.alpha = 4,
tau.tau.alpha = 0.4, mu.lambda = 2, tau.lambda = 0.2, rjags.B = 10000L,
rjags.Burn = 1000L, rjags.n.chains = 1L, rjags.n.adapt = 5000L,r.seed=10L,rjags.seed=10L)
{
  if (!(is.numeric(mu.mu.alpha)    && mu.mu.alpha > 0))     {stop("mu.mu.alpha>0")}
  if (!(is.numeric(tau.mu.alpha)   && tau.mu.alpha > 0))    {stop("tau.mu.alpha>0")}
  if (!(is.numeric(mu.tau.alpha)   && mu.tau.alpha > 0))    {stop("mu.tau.alpha>0")}
  if (!(is.numeric(tau.tau.alpha)  && tau.tau.alpha > 0))   {stop("tau.tau.alpha>0")}
  if (!(is.numeric(mu.lambda)      && mu.lambda > 0))       {stop("mu.lambda>0")}
  if (!(is.numeric(tau.lambda)     && tau.lambda > 0))      {stop("tau.lambda>0")}
  if (!(is.numeric(rjags.B)        && rjags.B >= 1))        {stop("rjags.B>=1")}
  if (!(is.numeric(rjags.Burn)     && rjags.Burn > 0))      {stop("rjags.Burn>0")}
  if (!(is.numeric(rjags.n.chains) && rjags.n.chains == 1)) {stop("rjags.n.chains==1")}
  if (!(is.numeric(rjags.n.adapt)  && rjags.n.adapt >= 1))  {stop("rjags.n.adapt>=1")}
  if (!(is.numeric(r.seed)  && r.seed > 0))  {stop("r.seed>0")}
  if (!(is.numeric(rjags.seed)  && rjags.seed > 0))  {stop("rjags.seed>0")}

  return(list(mu.mu.alpha = mu.mu.alpha, tau.mu.alpha = tau.mu.alpha,
  mu.tau.alpha = mu.tau.alpha, tau.tau.alpha = tau.tau.alpha,
  mu.lambda = mu.lambda, tau.lambda = tau.lambda, rjags.B = rjags.B,
  rjags.Burn = rjags.Burn, rjags.n.chains = rjags.n.chains,
  rjags.n.adapt = rjags.n.adapt,r.seed=r.seed,rjags.seed=rjags.seed))
}

## Main Function
ordinalRR<-function (x, random = TRUE, control = ordinalRR.control())
{
  cl <- match.call(expand.dots = TRUE)
  cl[[1]] <- as.name("ordinalRR")
  if (!is.list(x)) {stop("x must be a list")}
  if (is.null(x$preprocess)) {stop("x must be from call .....")}
  if (!x$preprocess) {stop("x must be from call .....")}
  dat <- list()
  dat[[1]]  = x$I
  dat[[2]]  = x$J
  dat[[3]]  = x$K
  dat[[4]]  = x$H
  dat[[5]]  = x$R
  if(random){
   dat[[6]]  = control$mu.mu.alpha
   dat[[7]]  = control$tau.mu.alpha
   dat[[8]]  = control$mu.tau.alpha
   dat[[9]]  = control$tau.tau.alpha
   dat[[10]] = control$mu.lambda
   dat[[11]] = control$tau.lambda
   names(dat) <- c("I", "J", "K", "H", "R", "mu.mu.alpha", "tau.mu.alpha",
                   "mu.tau.alpha", "tau.tau.alpha", "mu.lambda", "tau.lambda")
  }else{names(dat) <- c("I", "J", "K", "H", "R")}
  
  #setwd(tempdir()) #fix warning March 2018
  if (random) {#random
    modelString = "
    model{#random
        for(j in 1:J){
            alpha[j]~dlnorm(mu.alpha,tau.alpha) #RANDOM
            #alphainv[j]~dgamma(.001,.5) #FIXED
            #alpha[j]<-1/alphainv[j] #FIXED
            pi[j,1:H]~ddirch(lambda)
            for(h in 1:(H-1)){delta[j,h]<-qnorm(sum(pi[j,1:h]),0,1)}
        }
        for(i in 1:I){
            X[i]~dnorm(0,1)
            for(j in 1:J){
                p[i,j,1]<-1
                for(h in 2:H){p[i,j,h]<-exp(sum(alpha[j]*(X[i]-delta[j,1:(h-1)])))}
                R[i,j,1:H]~dmulti(p[i,j,1:H]/sum(p[i,j,1:H]),K)
            }
        }
        mu.alpha ~dnorm( mu.mu.alpha, tau.mu.alpha)  #RANDOM
        tau.alpha~dlnorm(mu.tau.alpha,tau.tau.alpha) #RANDOM
        for(h in 1:H){lambda[h]~dlnorm(mu.lambda,tau.lambda)} #RANDOM
        #for(h in 1:H){lambda[h]<-1/2} #FIXED
    }"
  }
    else {#fixed
        modelString = "
        model{#fixed
            for(j in 1:J){
                #alpha[j]~dlnorm(mu.alpha,tau.alpha) #RANDOM
                alphainv[j]~dgamma(.001,.5) #FIXED
                alpha[j]<-1/alphainv[j] #FIXED
                pi[j,1:H]~ddirch(lambda)
                for(h in 1:(H-1)){delta[j,h]<-qnorm(sum(pi[j,1:h]),0,1)}
            }
            for(i in 1:I){
                X[i]~dnorm(0,1)
                for(j in 1:J){
                    p[i,j,1]<-1
                    for(h in 2:H){p[i,j,h]<-exp(sum(alpha[j]*(X[i]-delta[j,1:(h-1)])))}
                    R[i,j,1:H]~dmulti(p[i,j,1:H]/sum(p[i,j,1:H]),K)
                }
            }
            #mu.alpha ~dnorm( mu.mu.alpha, tau.mu.alpha)  #RANDOM
            #tau.alpha~dlnorm(mu.tau.alpha,tau.tau.alpha) #RANDOM
            #for(h in 1:H){lambda[h]~dlnorm(mu.lambda,tau.lambda)} #RANDOM
            for(h in 1:H){lambda[h]<-1/2} #FIXED
        }"
    }
    #writeLines(modelString, con = "jags.bug")#fix warning March 2018
  temp=textConnection(modelString)#fix warning March 2018
  jfit = jags.model(temp, data = dat, n.chains = control$rjags.n.chains,
  n.adapt = control$rjags.n.adapt,
  inits=list(.RNG.name = "base::Mersenne-Twister",.RNG.seed=control$rjags.seed))
  close(temp)#fix warning March 2018
  update(jfit, control$rjags.Burn)
  obj <- NULL
  obj$dat=x
  obj$call <- cl
  obj$control <- control
  obj$random <- random
  if (random) {
    obj$post <- coda.samples(jfit, c("alpha", "delta", "mu.alpha",
    "tau.alpha", "lambda", "X"), n.iter = control$rjags.B)
    
    obj$x=obj$post[[1]][,1:x$I]
    obj$a=obj$post[[1]][,x$I+1:x$J]
    obj$d=obj$post[[1]][,x$I+x$J+1:((x$H-1)*x$J)]
    temp=seq(from=0,by=x$J,length=x$H-1)
    permute=NULL
    for(i in 1:x$J) permute=c(permute,temp+i)
    obj$d=obj$d[,permute]
    obj$lambda=obj$post[[1]][,x$I+x$J*x$H+1:x$H]
    hyper.a=obj$post[[1]][,x$I+(x$J+1)*x$H+1:2]
    obj$mu.a=hyper.a[,1]
    obj$sigma.a=1/sqrt(hyper.a[,2])

set.seed(control$r.seed)
obj$dnew=cbind(rdelta(obj$lambda),rdelta(obj$lambda))
obj$anew=cbind(rlnorm(control$rjags.B,obj$mu.a,obj$sigma.a),
               rlnorm(control$rjags.B,obj$mu.a,obj$sigma.a))

  }
  else {
      obj$post <- coda.samples(jfit, c("alpha","delta","X"),n.iter=control$rjags.B)
      obj$x=obj$post[[1]][,1:x$I]
      obj$a=obj$post[[1]][,x$I+1:x$J]
      obj$d=obj$post[[1]][,x$I+x$J+1:((x$H-1)*x$J)]
      temp=seq(from=0,by=x$J,length=x$H-1)
      permute=NULL
      for(i in 1:x$J) permute=c(permute,temp+i)
      obj$d=obj$d[,permute]
  }

  structure(obj, class = "ordinalRR")
}
## preProcessData
preprocess<-function (x, J = 3, K = 2, H = 4)
{
  if (!is.data.frame(x)) {stop("x must be a data frame.")}
  if (!(is.numeric(J) && J>=1)) {stop("Number of operators must be J>=1.")}
  if (!(is.numeric(K) && K>=1)) {stop("Number of repetitions must be K>=1.")}
  if (!(is.numeric(H) && H>=2)) {stop("Number of ordinal categories must be H>=2.")}
  if (J * K != ncol(x)){stop("The number of columns in x must be J*K.")}
  I = nrow(x)
  R = array(0, c(I, J, H))
  for (i in 1:I) {
    for (j in 1:J) {
      for (k in 1:K) {
        if(!sum(x[i,(j - 1) * K + k]==1:H)){stop("Entries of x must be from {1,...,H}.")}
        R[i, j, x[i, (j - 1) * K + k]] = R[i, j, x[i,(j - 1) * K + k]] + 1
      }
    }
  }
  list(I = I, J = J, K = K, H = H, x=x, R = R, preprocess = TRUE)
}

ordinalRR.sim=function(H=4L,I=30L,J=3L,K=2L,mu.a=2.6,sigma.a=.2,lambda=c(11,44,29,40),seed=10L)
{
  set.seed(seed)
  dataset=matrix(0,nrow=I,ncol=J*K)
  a=rlnorm(J,mu.a,sigma.a) #alpha model parameters
  d=matrix(0,nrow=J,ncol=H-1)
  for(j in 1:J)
  {
      temp=rgamma(H,lambda)
      temp=temp/sum(temp)
      temp=cumsum(temp)[-H]
      d[j,]=qnorm(temp)
  }
  temp=paste("Simulated parameters for operator 1 are: alpha1=", round(a[1],2), " and delta1=(", sep="")
  for(j in 1:(H-2)) temp=paste(temp,round(d[1,j],2),",",sep="")
  temp=paste(temp,round(d[1,H-1],2),").",sep="")
  print(temp)
  
  x=rnorm(I)
  for(i in 1:I)
  for(j in 1:J)
  {
      p=1
      for(h in 2:H)
      p=c(p,exp(sum(a[j]*(x[i]-d[j,1:(h-1)]))))
      p=p/sum(p)
      dataset[i,(j-1)*K+1:K]=sample(1:H,K,replace=T,p)
  }
  preprocess(as.data.frame(dataset),J,K,H)
}

make.rater=function(alpha,cutpoints)
{
obj=list(alpha,cutpoints)
names(obj)=c("alpha","cutpoints")
structure(obj, class = "rater")
}

compute.q=function(rater,x)
{
    if (class(rater) != "rater") stop("Object must be of class `rater'")
    a=rater$alpha
    d=rater$cutpoints
    H=length(d)+1
    p=sapply(1:(H-1),function(h)a*(x-d[h]))
    if(is.vector(p)) p=matrix(p,nrow=1)
    p=t(apply(p,1,cumsum))
    p=cbind(0,p)
    p=exp(p-apply(p,1,max))
    p=p/as.vector(apply(p,1,sum))
    
    p
}



## S3-Generics

#plot.ordinalRR=function()

plot.rater=function(x,y,plt.type=c("rater","measure"),m=0,lwd=1.2,...)
{
    glen=10^3
    if(plt.type=="rater")
    {
        if (class(x) != "rater") stop("Object x must be of class `rater'")
        alpha=rep(x$alpha,glen)
        cuts=x$cutpoints
        delta=matrix(rep(cuts,glen),nrow=glen,byrow=T)
        xgrid=seq(-3,3,length=glen)
        p=computep(alpha,xgrid,delta)
        plot(0,.5,ylim=c(0,1),xlim=c(-3,3),xaxt="n",yaxt="n",type="n",...)
        box(lwd=lwd)
        axis(1,-3:3,lwd=lwd)
        axis(2,(0:4)/4,c("0","","0.5","","1"),lwd=lwd)
        for(i in 1:4) lines(xgrid,p[,i],lwd=lwd)
        abline(v=cuts,lty=2)
    }
    if(plt.type=="measure")
    {
        if (class(x) != "rater") stop("Object x must be of class `rater'")
        if (class(y) != "rater") stop("Object y must be of class `rater'")
        if (length(x$cutpoints)!=length(y$cutpoints)) stop("Rater must have same number of cutpoints.")
        H=length(x$cutpoints)+1
        m=m+1
        Bm=toeplitz(c(rep(1,m),rep(0,H-m)))
        
        xgrid=matrix(seq(-3,3,length=glen),ncol=1)
        
        alpha=rep(x$alpha,glen)
        delta=matrix(rep(x$cutpoints,glen),nrow=glen,byrow=T)
        p1=computep(alpha,xgrid,delta)

        alpha=rep(y$alpha,glen)
        delta=matrix(rep(y$cutpoints,glen),nrow=glen,byrow=T)
        p2=computep(alpha,xgrid,delta)

        repeat1=rowSums(p1%*%Bm*p1)
        repeat2=rowSums(p2%*%Bm*p2)
        rr=rowSums(p1%*%Bm*p2)
        prop=rr^2/(repeat1*repeat2)
        
        plot(0,.5,ylim=c(0,1),xlim=c(-3,3),xaxt="n",yaxt="n",type="n",...)
        box(lwd=lwd)
        axis(1,-3:3,lwd=lwd)
        axis(2,(0:4)/4,c("0","","0.5","","1"),lwd=lwd)
        lines(xgrid,rr,lwd=lwd)
        lines(xgrid,prop,lty=2,lwd=lwd)
    }
    
}

hist.ordinalRR=function(x,x.low=-4,x.high=4,col="grey",...)
{
    if (class(x) != "ordinalRR") stop("Object must be of class `ordinalRR'")
    den=list()
    denmax=0
    n<-x$dat$I+1
    
    for(i in 1:n)
    {
        if(i<=x$dat$I) den[[i]]=density(x$x[,i],from=x.low,to=x.high)
        if(i>x$dat$I){den[[i]]=den[[i-1]]; den[[i]]$x=seq(x.low,x.high,length=10^3); den[[i]]$y=dnorm(den[[i]]$x)}
        denmax=max(c(denmax,den[[i]]$y))
    }
    plot(den[[i]]$x,den[[i]]$y,  type="n",xlim=c(x.low,x.high),ylim=c(0,denmax),col=col,...)
    if(length(col)==1)col<-rep(col,n)
    if(length(col)<n){
        col<-c(col,rep("grey",n-length(col)+1))
        warning("Color length mismatch, grey's added")
    }
    for(i in 1:(x$dat$I+1))lines(den[[i]]$x,den[[i]]$y,col=col[i])
}


density.ordinalRR<-function(x,plt.type=c("repeat","rr","prop","all"),m=0,...){
    if(class(x)!="ordinalRR")stop("error")
    I=ncol(x$x)
    J=ncol(x$a)
    B=nrow(x$a)
    if(missing(plt.type)){
        plt.type="repeat"
    }
    final<-list()
    k1<-1
    for(j1 in 1:(J-1)){
        for(j2 in (j1+1):J){
            RRavg=matrix(0,nrow=B,ncol=4)
            for(i in 1:I){
                RRavg=RRavg+RR(x,j1,j2,i=i,m=m)
            }
            final[[k1]]<-RRavg/I
            k1<-k1+1
        }
    }
    if(x$random)
    {
     RRnew=matrix(0,nrow=B,ncol=4)
     for(i in 1:I) RRnew=RRnew+RR(x,J+1,J+2,i=i,m=m)
     RRnew=RRnew/I
    }
    
    if(plt.type=="all")par(mfrow=c(3,1))
    if(plt.type=="repeat"||plt.type=="all"){
        temp<-list()
        temp[[1]]<-density(final[[1]][,1],from=0,to=1)
        ymax<-max(temp[[1]]$y)
        for(i in 2:J){
            temp[[i]]<-density(final[[i-1]][,2],from=0,to=1)
            ymax=max(c(ymax,temp[[i]]$y))
        }
        if(plt.type=="all") plot(temp[[1]]$x,temp[[1]]$y,ylim=c(0,ymax),xlab="Repeatability",ylab="",...)
        else plot(temp[[1]]$x,temp[[1]]$y,ylim=c(0,ymax),...)
        lapply(temp,function(i){lines(i,col="grey")})
        if(x$random)lines(density(RRnew[,1],from=0,to=1),lwd=2)
    }
    if(plt.type=="rr"||plt.type=="all"){
        temp<-list()
        temp[[1]]<-density(final[[1]][,3],from=0,to=1)
        ymax<-max(temp[[1]]$y)
        for(i in 2:(choose(J,2))){
            temp[[i]]<-density(final[[i]][,3],from=0,to=1)
            ymax=max(c(ymax,temp[[i]]$y))
        }
        if(plt.type=="all") plot(temp[[1]]$x,temp[[1]]$y,ylim=c(0,ymax),xlab="R&R",ylab="Density",...)
        else plot(temp[[1]]$x,temp[[1]]$y,ylim=c(0,ymax),...)
        lapply(temp,function(i){lines(i,col="grey")})
        if(x$random)lines(density(RRnew[,3],from=0,to=1),lwd=2)
    }
    if(plt.type=="prop"||plt.type=="all"){
        temp<-list()
        temp[[1]]<-density(final[[1]][,4],from=0,to=1)
        ymax<-max(temp[[1]]$y)
        for(i in 2:(choose(J,2))){
            temp[[i]]<-density(final[[i]][,4],from=0,to=1)
            ymax=max(c(ymax,temp[[i]]$y))
        }
        if(plt.type=="all") plot(temp[[1]]$x,temp[[1]]$y,ylim=c(0,ymax),xlab="Proportion",ylab="",...)
        else plot(temp[[1]]$x,temp[[1]]$y,ylim=c(0,ymax),...)
        lapply(temp,function(i){lines(i,col="grey")})
        if(x$random)lines(density(RRnew[,4],from=0,to=1),lwd=2)
    }
    
    invisible(x)
}
## print
print.ordinalRR<-function (x, ...){
  if (class(x) != "ordinalRR")
  stop("Object must be of class `ordinalRR'.")
  if (!is.null(cl <- x$call)) {
    names(cl)[2] <- ""
    cat("Call:\n")
    dput(cl)
  }
  cat("\nData:", x$dat$I, "parts,", x$dat$J, "operators,", x$dat$K, "repetitions with", x$dat$H, "ordinal categories.\n")
txt <- "fixed-effects"; if (x$random) {txt <- "random-effects"}
cat("\nA single MCMC chain of the", txt, "ordinal model was fit:", x$control$rjags.Burn, "burn-in and", x$control$rjags.B, "retained.\n")

invisible(x)
}

##summary
summary.ordinalRR<-function (object, decimals=1,...){
    if (class(object) != "ordinalRR")
    stop("Object must be of class `ordinalRR'.")
    if (!is.null(cl <- object$call)) {
        names(cl)[2] <- ""
        cat("Call:\n")
        dput(cl)
    }
    cat("\nData:", object$dat$I, "parts,", object$dat$J, "operators,", object$dat$K, "repetitions with", object$dat$H, "ordinal categories.\n")
    txt <- "Fixed-effects"; if (object$random) {txt <- "Random-effects"}
    cat(txt, "model MCMC chain:", object$control$rjags.Burn, "burn-in and", object$control$rjags.B, "retained.\n")
    
    I=object$dat$I
    J=object$dat$J
    K=object$dat$K
    H=object$dat$H
    x=object$dat$x
    a=object$a
    d=object$d

    matches=rep(0,J) #Repeatability
    for(j in 1:J)
    for(k1 in 1:(K-1))
    for(k2 in (k1+1):K)
    matches[j]=matches[j]+sum(x[,K*(j-1)+k1]==x[,K*(j-1)+k2])
    matches=round(matches/(I*choose(K,2)),decimals+2)

    dat=cbind(apply(a,2,median),t(matrix(apply(d,2,median),3)))
    dat=cbind(1:dim(dat)[1],matches,round(dat,decimals))
    dimnames(dat)[[2]]=c("Rater j","Repeatability","a_j",paste0("d_{j,",1:(H-1),"}"))
    cat("\nSimple repeatability and model parameter estimates by rater:\n")
    print(as.data.frame(dat),row.names=FALSE)
    
    matches=rep(0,choose(J,2)) #reproducibility
    r1=rep(0,choose(J,2))
    r2=rep(0,choose(J,2))
    b=0
    for(j1 in 1:(J-1))
    for(j2 in (j1+1):J)
    {
        b=b+1
        r1[b]=j1
        r2[b]=j2
        for(k1 in 1:K)
        for(k2 in 1:K)
        matches[b]=matches[b]+sum(x[,K*(j1-1)+k1]==x[,K*(j2-1)+k2])
    }
    matches=round(matches/(K^2*I),decimals+2)
    
    matches=cbind(r1,r2,matches)
    dimnames(matches)[[1]]=NULL
    dimnames(matches)[[2]]=c("Rater j","Rater j\'","(R&R)_{j,j\'}")
    cat("\nSimple repeatability and reproducibility (R&R) point estimates for pairs of raters:\n")
    print(as.data.frame(matches),row.names=FALSE)
    
    invisible(object)
}

### Internal

rdelta=function(lambda)
{
    B=nrow(lambda)
    H=ncol(lambda)
    p=matrix(rgamma(B*H,lambda),nrow=B,ncol=H)
    p=p/apply(p,1,sum)
    
    delta=t(apply(p,1,cumsum))[,-H]
    delta[delta>1]=1
    delta=qnorm(delta)
    
    delta
}


computep=function(a,x,d)
{
    H=ncol(d)+1
    p=sapply(1:(H-1),function(h)a*(x-d[,h]))
    if(is.vector(p)) p=matrix(p,nrow=1)
    p=t(apply(p,1,cumsum))
    p=cbind(0,p)
    p=exp(p-apply(p,1,max))
    p=p/as.vector(apply(p,1,sum))
    
    p
}



RR=function(g,j1,j2,i=1,m=0)
{
    if(class(g)!="ordinalRR"){stop("some error")}
   
    ## get stuff from g$post
    I=ncol(g$x)
    B=nrow(g$a)
    J=ncol(g$a)
    a=g$a
    d=g$d
    x=g$x
    H=ncol(d)/J+1
    
    m=m+1
    Bm=toeplitz(c(rep(1,m),rep(0,H-m)))
    
    if(j1>j2){temp=j1; j1=j2; j2=temp}
    same=0; if(j1==j2) same=1
    if(j1>J){J=J+1; j1=J; d=cbind(d,g$dnew); a=cbind(a,g$anew)}
    if(same==1) j2=j1
    if(j2>J){J=J+1; j2=J;}
    
    p1=computep(a[,j1],x[,i],d[,(j1-1)*(H-1)+1:(H-1)])
    if(same==1){p2=p1} else {p2=computep(a[,j2],x[,i],d[,(j2-1)*(H-1)+1:(H-1)])}
    
    repeat1=rowSums(p1%*%Bm*p1)
    repeat2=rowSums(p2%*%Bm*p2)
    rr=rowSums(p1%*%Bm*p2)
    prop=rr^2/(repeat1*repeat2)
    
    cbind(repeat1,repeat2,rr,prop)
}
