# sample size calculation and power analysis in R

```
library(pwr)
```

```
effect_size<-1.1
power.t.test(delta=effect_size,sd=1,sig.level=0.05,power=.80,alternative="two.sided")
## equivalent to
uniroot(function(n)pt(qt(0.025,df=(n-1)*2,lower.tail=F),df=(n-1)*2,ncp=sqrt(n/2)*effect_size,lower.tail=F)-0.8,c(3,20))$root
```

# linear regression

https://advstats.psychstat.org/book/power/index.php

```
# Eg: GPA_college ~ SAT + GPA_high_school
# SAT + GPA_high_school explains 50% of variance
# Test: SAT and GPA_high_school
# R2full=0.5
# R2reduced=0.0

R2full<-0.5
R2reduced<-0.0
f2<-(R2full-R2reduced)/(1-R2full)

wp.regression(p1=2,p2=0,f2=f2,alpha=0.05,power=0.8)

# Eg: GPA_college ~ SAT + GPA_high_school + letter
# SAT + GAP_high_school explains 50% of variance
# letter additionally explains 5%
# Test: letter
# R2full=0.55
# R2reduced=0.5

R2full<-0.55
R2reduced<-0.5
f2<-(R2full-R2reduced)/(1-R2full)

wp.regression(p1=3,p2=2,f2=f2,alpha=0.05,power=0.8)
```

## Survival Continuous

```

library(powerSurvEpi)
logHR<-2
ssizeEpiCont.default(
  power = 0.8,
  theta = 1.55,
  sigma2 = 1,
  psi = 0.5,
  rho2 = 0.0,
  alpha = 0.05)

ssizeEpiCont.default(
  power = 0.8,
  theta = 1.85,
  sigma2 = 1,
  psi = 0.5,
  rho2 = 0.0,
  alpha = 0.05)
```

```
library(survival)

nsim<-10000
reject<-0
for(sim in 1:nsim){
  if(sim%%100==0)cat(sim,"|")
  nind<-82
  xx<-rnorm(nind)
  HR<-1.56
  beta<-logHR<-log(HR)
  ft<-rexp(nind)*exp(-beta*xx)
  ct<-quantile(ft,0.5)
  ot<-pmin(ft,ct)
  delta<-ft<=ct
  CI<-confint(coxph(Surv(ot,delta)~xx))["xx",]
  if(all(CI>=0))reject<-reject+1
}
reject/nsim
```

## correlation and R2

```
library(MASS)

allR<-c()
for(sim in 1:10000){
  Sigma<-matrix(c(1,0.3,0.3,1),2,2)
  xx<-mvrnorm(100,c(0,0),Sigma)
  x1<-xx[,1]
  x2<-xx[,2]
  allR[sim]<-summary(lm(x1~x2))$r.squared
}

summary(lm(x1~x2))$r.squared
cor(x1,x2)^2
```

## Two-sample variance

```
ratio<-1.48
nsim<-100000
reject<-0
for(sim in 1:nsim){
  if(sim%%1000==0)print(sim)
  x1<-rnorm(82,0,sd=1*ratio)
  x2<-rnorm(42,0,sd=1)
  if(var.test(x1,x2,alternative="greater")$p.value<0.025)reject<-reject+1
}
reject/nsim

```
