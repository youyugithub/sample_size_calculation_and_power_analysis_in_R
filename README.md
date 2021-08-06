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

```
library(seqtest)
size.cor(rho=0.6, 0.2, "two.sided", alpha = 0.05, beta = 0.1, output = TRUE)
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

## Survival two groups

Number of Events needed $d=(z_\beta+z_\alpha)^2/[P_A P_B (\log HR)^2]$

#### Freedman's

TABLES OF THE NUMBER OF PATIENTS REQUIRED IN CLINICAL TRIALS USING THE LOGRANK TEST 

https://eclass.uoa.gr/modules/document/file.php/MATH301/PracticalSession2/Freedman2006.pdf

Allocation ratio: $\phi:1$
Hazard ratio: $\theta:1$

Number of Events needed: $d=(z_\beta+z_\alpha)^2(1+\theta\phi)^2/(1-\theta)/\phi$
Total numbers needed: $d*(1+\phi)/(\phi*(1-P_1)+(1-P_2))$

```
# sample size

# Control: 2-year survival 0.5
# Treatment: HR 0.5 less event

# Freedman's
HR<-0.5
S0_2y<-0.6
S1_2y<-exp(HR*log(S0_2y))
hazard0<--log(S0_2y)/2
hazard1<--log(S1_2y)/2
# because HR==log(S1_2y)/log(S0_2y)
ratio<-2
d1<-(1+HR*ratio)^2/(1-HR)^2/ratio*(qnorm(0.975)+qnorm(0.8))^2
d1*(1+ratio)/(ratio*(1-S0_2y)+(1-S1_2y))

# Schoenfeld's
HR<-0.5
d2<-(qnorm(0.975)+qnorm(0.8))^2/(1/3*2/3*log(HR)^2)
d2*(1+ratio)/(ratio*(1-S0_2y)+(1-S1_2y))

library(powerSurvEpi)
powerCT.default(
  nE=122,nC=61,pE=1-S0_2y,pC=1-S1_2y,
  RR=0.5,alpha=0.05)
```

## For two sample survival, information can be approximated by ndeath*P0*P1

```
library(survival)

n_pl<-60
n_tr<-60

lambda_pl<--log(0.5)/2
lambda_tr<-0.43*lambda_pl

allvar<-rep(NA,10000)
for(sim in 1:10000){
  if(sim%%100==0)print(sim)
  Fi_pl<-rexp(n_pl,lambda_pl)
  Fi_tr<-rexp(n_tr,lambda_tr)
  
  Fi<-c(Fi_pl,Fi_tr)
  Ci<-sort(Fi)[50]
  Ti<-pmin(Fi,Ci)
  deltai<-Fi<=Ci
  treat<-rep(c(0,1),c(n_pl,n_tr))
  
  fit<-survreg(Surv(Ti,deltai)~treat)
  allvar[sim]<-fit$var["treat","treat"]
}

mean(1/allvar)
50*(1/3)*(2/3)
```

# F test related
```
set.seed(0)

y<-rnorm(100)
x<-rep(LETTERS[1:4],each=25)
fit<-lm(y~x)
pred<-predict(fit)

anova(fit)$`F value`[1]
anova(fit)$`Pr(>F)`[1]

numerator<-sum((pred-mean(y))^2)
df_numerator<-4-1
denominator<-sum((y-pred)^2)
df_denominator<-100-4
my_f<-(numerator/df_numerator)/(denominator/df_denominator)
numerator+denominator==sum((y-mean(y))^2)

library(pwr)
pwr.anova.test(k=4,n=22,power=0.80,sig.level=0.05)$f
sqrt(pwr.f2.test(u=4-1,v=88-4,power=0.80,sig.level=0.05)$f2)
library(WebPower)
wp.anova(k=4,n=88,power=0.8,alpha=0.05)
```

# t-test

```
mu0_log<-0.36
mu0_exp<-exp(mu0_log)-1
mu1_exp<-1.5*mu0_exp
mu1_log<-log(mu1_exp+1)
sigma<-0.167

reject_B<-rep(0,1000000)
reject_C<-rep(0,1000000)
signif_B<-rep(0,1000000)
signif_C<-rep(0,1000000)
se<-sigma/sqrt(23/2)
for(iter in 1:1000000){
  A<-rnorm(1,mu0_log,se)
  B<-rnorm(1,mu1_log,se)
  C<-rnorm(1,mu1_log,se)
  if(A>B)reject_B[iter]<-1
  if(A>C)reject_C[iter]<-1
  if((B-A)>qt(0.95,2*23/2-2)*se)signif_B[iter]<-1
  if((C-A)>qt(0.95,2*23/2-2)*se)signif_C[iter]<-1
}
table(reject_B)/1000000
table(reject_B+reject_C>0)/1000000
table(signif_B+signif_C>0)/1000000

table(signif_B)/1000000

reject_B<-rep(0,100000)
reject_C<-rep(0,100000)
for(iter in 1:100000){
  A<-rnorm(1)
  B<-rnorm(1)
  C<-rnorm(1)
  if(A>B)reject_B[iter]<-1
  if(A>C)reject_C[iter]<-1
}
table(reject_B)
table(reject_B+reject_C)
```
