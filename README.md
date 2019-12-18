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

