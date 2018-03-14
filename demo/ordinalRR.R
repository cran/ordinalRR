library(ordinalRR)

#Part I of Demo: Figure 1 from Culp, Ryan, Chen, and Hamada (2018)
r1=make.rater(1,c(-1.7,-0.5,1.6))
r2=make.rater(1,c(-1.4,0.0,1.2))
r3=make.rater(3,c(-2.6,-1.6,0.3))
r4=make.rater(3,c(-1.7,1.2,2.6))

par(mfrow=c(3,2))
plot(r1,plt.type="rater",xlab="",ylab="q1",main="alpha=1")

#Display multinomial probabilities in upper left panel when x=0.1 and plot them with circles
(probs=compute.q(r1,.1))
sum(probs)
points(rep(.1,4),probs)

plot(r3,plt.type="rater",xlab="",ylab="",main="alpha=3")
plot(r2,plt.type="rater",xlab="",ylab="q2")
plot(r4,plt.type="rater",xlab="",ylab="")
plot(r1,r2,plt.type="measure",xlab="x",ylab="Measure")
plot(r3,r4,plt.type="measure",xlab="x",ylab="")
title("Figure 1: Culp, Ryan, Chen, and Hamada (2018)",outer=T,line=-1)

invisible(readline(prompt="Press [enter] to continue."))

#Part II of Demo: Analyze followup data
data(followup)
x=preprocess(followup)
g.random<-ordinalRR(x)
g.fixed<-ordinalRR(x,random=F)

#Table 2 from Culp, Ryan, Chen, and Hamada (2018)
random=round(apply(g.random$a,2,median),1)
fixed=round(apply(g.fixed$a,2,median),1)
rbind(fixed,random)

#Compare 95% Bayes intervals for alpha1 (random interval is shorter)
fixed=round(quantile(g.fixed$a[,1],c(.025,.975)),1)
random=round(quantile(g.random$a[,1],c(.025,.975)),1)
rbind(fixed,random)

#Figure 4 from Culp, Ryan, Chen, and Hamada (2018)
par(mfrow=c(1,1))
hist(g.random,xlab="x",ylab="Density",main="Figure 4: Culp, Ryan, Chen, and Hamada (2018)")

#Figure 5 from Culp, Ryan, Chen, and Hamada (2018)
density(g.random,plt.type="all",xlim=c(.6,1),type="l",col="grey")
title("Figure 5: Culp, Ryan, Chen, and Hamada (2018)",outer=T,line=-2)

invisible(readline(prompt="Press [enter] to continue."))

#Part III of Demo: Simulation
#Simulate data sets (note parameters for operator 1 are printed to screen)
data.small=ordinalRR.sim(I=10,J=1)
data.bigger=ordinalRR.sim(I=50,J=5)

#Fit random-effects model
g.small =ordinalRR(data.small)
g.bigger=ordinalRR(data.bigger)

#Compare 95% Bayes intervals for alpha1 to each other and to the true value
#Expect big data to have the shorter interval
small=round(quantile(g.small$a,c(.025,.975)),1)
big  =round(quantile(g.bigger$a[,1],c(.025,.975)),1)
rbind(small,big)
