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

## Survival

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
