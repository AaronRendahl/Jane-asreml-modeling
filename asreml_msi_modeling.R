library(asreml)

OST_R_asreml_NS_hr <- read.csv("OST_R_asreml_NS_hr.csv", check.names=FALSE)

m1 <- asreml(glucose ~ time + time2, data=OST_R_asreml_NS_hr, random=~str(~Number_ID + time:Number_ID + time2:Number_ID, ~us(3):id(Number_ID)) + ~Owner, rcov=~units, maxit=400,  Cfixed=TRUE)
m1 <- asreml(glucose ~ time + time2, data=OST_R_asreml_NS_hr, random=~str(~Number_ID + time:Number_ID + time2:Number_ID, ~us(3):id(Number_ID)))

## previous output
## glucose works with farm random intercept and breed fixed effect
#
# > m1 <- asreml(glucose ~ hour + hour2 + Breed + Breed:hour + Breed:hour2, data=df2, random=~str(~Number_ID + hour:Number_ID + hour2:Number_ID, ~us(3):id(Number_ID)) + ~Owner, rcov=~units, maxit=400,  Cfixed=TRUE)
# ASReml: Mon Oct 31 17:37:13 2016
# 
#  US matrix updates modified 1 times to remain positive definite.
#      LogLik         S2      DF      wall     cpu
#   -2206.1541     95.3231   723  17:37:13     0.0 (6 restrained)
#   -2151.5520     94.8190   723  17:37:13     0.0 (3 restrained)
#  US matrix updates modified 1 times to remain positive definite.
#   -2142.4206     90.3250   723  17:37:13     0.0 (4 restrained)
#  US matrix updates modified 1 times to remain positive definite.
#   -2127.6595     87.2655   723  17:37:13     0.0 (6 restrained)
#  US matrix updates modified 1 times to remain positive definite.
#   -2116.6919     83.3261   723  17:37:13     0.0 (6 restrained)
#  US matrix updates modified 1 times to remain positive definite.
#   -2107.3980     79.8055   723  17:37:13     0.0 (6 restrained)
#  US matrix updates modified 1 times to remain positive definite.
#   -2099.4692     76.5959   723  17:37:13     0.0 (6 restrained)
#  US matrix updates modified 1 times to remain positive definite.
#   -2092.5992     73.8339   723  17:37:13     0.0 (6 restrained)
#   -2086.6355     71.4900   723  17:37:13     0.0 (1 restrained)
#   -2063.4866     63.8633   723  17:37:13     0.0
#   -2048.7861     61.2144   723  17:37:13     0.0
#   -2044.8604     57.7331   723  17:37:13     0.0
#   -2044.5356     56.9646   723  17:37:13     0.0
#   -2044.5296     56.8785   723  17:37:13     0.0
#   -2044.5294     56.8769   723  17:37:13     0.0
#   -2044.5293     56.8768   723  17:37:13     0.0
# US variance structures were modified in 7 instances to make them positive definite
# 
# Finished on: Mon Oct 31 17:37:13 2016
#  
# LogLikelihood Converged 



m1 <- asreml(logINS ~ hour + hour2 + Breed + Breed:hour + Breed:hour2, data=df2, random=~str(~Number_ID + hour:Number_ID + hour2:Number_ID, ~us(3):id(Number_ID)) , rcov=~units, maxit=400,  Cfixed=TRUE)
> 


##insulin model works with log transformed insulin



> m1 <- asreml(logINS ~ hour + hour2, data=df2, random=~str(~Number_ID + hour:Number_ID + hour2:Number_ID, ~us(3):id(Number_ID)) + ~Owner, rcov=~units, maxit=400,  Cfixed=TRUE)
ASReml: Mon Oct 31 17:38:22 2016

     LogLik         S2      DF      wall     cpu
    446.1769      0.0665   735  17:38:22     0.0 (3 restrained)
    536.4241      0.0586   735  17:38:22     0.0 (1 restrained)
 US matrix updates modified 1 times to remain positive definite.
    611.1088      0.0387   735  17:38:22     0.0 (6 restrained)
    662.2820      0.0333   735  17:38:22     0.0
    682.7249      0.0277   735  17:38:22     0.0
    685.6366      0.0263   735  17:38:22     0.0
    685.8624      0.0260   735  17:38:22     0.0
    685.8717      0.0260   735  17:38:22     0.0
    685.8717      0.0260   735  17:38:22     0.0
    685.8717      0.0260   735  17:38:22     0.0
US variance structures were modified in 1 instances to make them positive definite

Finished on: Mon Oct 31 17:38:22 2016
 
LogLikelihood Converged 
> 
> 

## models for insulin and glucose trajectories


m1 <- asreml(logINS ~ hour + hour2 + Breed + Breed:hour + Breed:hour2, data=df2, random=~str(~Number_ID + hour:Number_ID + hour2:Number_ID, ~us(3):id(Number_ID)), rcov=~units, maxit=400,  Cfixed=TRUE)


m1 <- asreml(glucose ~ hour + hour2 + Breed + Breed:hour + Breed:hour2, data=df2, random=~str(~Number_ID + hour:Number_ID + hour2:Number_ID, ~us(3):id(Number_ID)), rcov=~units, maxit=400,  Cfixed=TRUE)


##Conditional model (estimate effects of age)
m1_plusAge<-update(m1, ~.+ AGE + time:AGE + time2:AGE)
x<-as.data.frame(summary(m1_plusAge, all=TRUE)$coef.fi)
pval<-as.vector(unname(2*pnorm(-abs(x[,'z ratio']))))
cbind(x, pval)


##Conditional model (estimate effects of breed)
m1_plusBreed<-update(m1, ~.+ Breed + time:Breed + time2:Breed)
x<-as.data.frame(summary(m1_plusBreed, all=TRUE)$coef.fi)
pval<-as.vector(unname(2*pnorm(-abs(x[,'z ratio']))))
cbind(x, pval) 
(summary(m1_plusBreed, all=TRUE)$coef.fi)
> pval<-as.vector(unname(2*pnorm(-abs(x[,'z ratio']))))
> cbind(x, pval)
                   solution std error      z ratio         pval
Breed_A:time2     0.0000000        NA           NA           NA
Breed_M:time2     0.0000000        NA           NA           NA
Breed_QH:time2    0.0000000        NA           NA           NA
Breed_T:time2     0.0000000        NA           NA           NA
Breed_WP :time2   0.0000000        NA           NA           NA
Breed_A:time      0.0000000        NA           NA           NA
Breed_M:time      0.0000000        NA           NA           NA
Breed_QH:time     0.0000000        NA           NA           NA
Breed_T:time      0.0000000        NA           NA           NA
Breed_WP :time    0.0000000        NA           NA           NA
hour2:Breed_A     0.0000000        NA           NA           NA
hour2:Breed_M    -4.7663142  2.126516  -2.24137219 2.500198e-02
hour2:Breed_QH    7.1911018  2.102212   3.42073048 6.245319e-04
hour2:Breed_T     3.2156348  3.189774   1.00810733 3.134029e-01
hour2:Breed_WP   -5.2005287  2.493561  -2.08558291 3.701642e-02
hour:Breed_A      0.0000000        NA           NA           NA
hour:Breed_M     18.1633158  6.351715   2.85959257 4.241856e-03
hour:Breed_QH   -16.1663662  6.279121  -2.57462242 1.003496e-02
hour:Breed_T      5.8858334  9.527572   0.61776846 5.367280e-01
hour:Breed_WP    20.7393516  7.448046   2.78453609 5.360436e-03
Breed_A           0.0000000        NA           NA           NA
Breed_M          -0.8801508  2.218268  -0.39677381 6.915343e-01
Breed_QH         -2.4838650  2.192916  -1.13267679 2.573500e-01
Breed_T           5.6018091  3.327403   1.68353807 9.227098e-02
Breed_WP          0.2399683  2.601150   0.09225469 9.264957e-01
hour2           -15.4763315  1.503674 -10.29234492 7.629533e-25
hour             45.3175308  4.491340  10.08997898 6.118262e-24
(Intercept)      85.4192032  1.568553  54.45733806 0.000000e+00
> ####overall Wald chi square tests of breed effect on entire trajectory
> covariance_matrix<-diag(length(m1_plusBreed$coefficients$fixed))
> covariance_matrix[upper.tri(covariance_matrix, diag=TRUE)]<-m1_plusBreed$Cfixed
> covariance_matrix[lower.tri(covariance_matrix)]<-t(covariance_matrix)[lower.tri(t(covariance_matrix))]
> colnames(covariance_matrix)<-names(m1_plusBreed$coefficients$fixed)
> rownames(covariance_matrix)<-names(m1_plusBreed$coefficients$fixed)
> terms<-grep("Breed", rownames(covariance_matrix))[rowSums(covariance_matrix[grep("Breed", rownames(covariance_matrix)),])!=0]
> wald.test(b = m1_plusBreed$coefficients$fixed[terms], Sigma =  covariance_matrix[terms, terms], Terms=1:length(terms))$result$chi2
    chi2       df        P 
130.4056  12.0000   0.0000 
>    
##Conditional model (estimate effects of breed)
m1_plusBreed<-update(m1, ~.+ Breed + time:Breed + time2:Breed)
x<-as.data.frame(summary(m1_plusBreed, all=TRUE)$coef.fi)
pval<-as.vector(unname(2*pnorm(-abs(x[,'z ratio']))))
cbind(x, pval) 


####overall Wald chi square tests of breed effect on entire trajectory
covariance_matrix<-diag(length(m1_plusBreed$coefficients$fixed))
covariance_matrix[upper.tri(covariance_matrix, diag=TRUE)]<-m1_plusBreed$Cfixed
covariance_matrix[lower.tri(covariance_matrix)]<-t(covariance_matrix)[lower.tri(t(covariance_matrix))]
colnames(covariance_matrix)<-names(m1_plusBreed$coefficients$fixed)
rownames(covariance_matrix)<-names(m1_plusBreed$coefficients$fixed)
terms<-grep("Breed", rownames(covariance_matrix))[rowSums(covariance_matrix[grep("Breed", rownames(covariance_matrix)),])!=0]
wald.test(b = m1_plusBreed$coefficients$fixed[terms], Sigma =  covariance_matrix[terms, terms], Terms=1:length(terms))$result$chi2



m1 <- asreml(logINS ~ hour + hour2 + Breed + Breed:hour + Breed:hour2, data=df2, random=~str(~Number_ID + hour:Number_ID + hour2:Number_ID, ~us(3):id(Number_ID)), rcov=~units, maxit=400,  Cfixed=TRUE)
##Conditional model (estimate effects of breed)
m1_plusBreed<-update(m1, ~.+ Breed + time:Breed + time2:Breed)
x<-as.data.frame(summary(m1_plusBreed, all=TRUE)$coef.fi)
pval<-as.vector(unname(2*pnorm(-abs(x[,'z ratio']))))
cbind(x, pval) 

##Conditional model (estimate effects of age)
m1_plusAge<-update(m1, ~.+ AGE + time:AGE + time2:AGE)
x<-as.data.frame(summary(m1_plusAge, all=TRUE)$coef.fi)
pval<-as.vector(unname(2*pnorm(-abs(x[,'z ratio']))))
cbind(x, pval) #estimates of age on the intercept (AGE), initial rate of change (time_min:AGE), and deceleration rate (time_min2:AGE2)
  
##conditional model adiponectin effects
m1_plusAdiponectin<-update(m1, ~.+ Adiponectin + time:Adiponectin + time2:Adiponectin)
x<-as.data.frame(summary(m1_plusAdiponectin, all=TRUE)$coef.fi)
pval<-as.vector(unname(2*pnorm(-abs(x[,'z ratio']))))
cbind(x, pval) 

####overall Wald chi square tests of Adiponectin effect on entire trajectory
covariance_matrix<-diag(length(m1_plusAdiponectin$coefficients$fixed))
covariance_matrix[upper.tri(covariance_matrix, diag=TRUE)]<-m1_plusAdiponectin$Cfixed
covariance_matrix[lower.tri(covariance_matrix)]<-t(covariance_matrix)[lower.tri(t(covariance_matrix))]
colnames(covariance_matrix)<-names(m1_plusAdiponectin$coefficients$fixed)
rownames(covariance_matrix)<-names(m1_plusAdiponectin$coefficients$fixed)
terms<-grep("Adiponectin", rownames(covariance_matrix))[rowSums(covariance_matrix[grep("Adiponectin", rownames(covariance_matrix)),])!=0]
wald.test(b = m1_plusAdiponectin$coefficients$fixed[terms], Sigma =  covariance_matrix[terms, terms], Terms=1:length(terms))$result$chi2

##conditional model TG effects
m1_plusTG<-update(m1, ~.+ TG + time:TG + time2:TG)
x<-as.data.frame(summary(m1_plusTG, all=TRUE)$coef.fi)
pval<-as.vector(unname(2*pnorm(-abs(x[,'z ratio']))))
cbind(x, pval) 

####overall Wald chi square tests of TG effect on entire trajectory
covariance_matrix<-diag(length(m1_plusTG$coefficients$fixed))
covariance_matrix[upper.tri(covariance_matrix, diag=TRUE)]<-m1_plusTG$Cfixed
covariance_matrix[lower.tri(covariance_matrix)]<-t(covariance_matrix)[lower.tri(t(covariance_matrix))]
colnames(covariance_matrix)<-names(m1_plusTG$coefficients$fixed)
rownames(covariance_matrix)<-names(m1_plusTG$coefficients$fixed)
terms<-grep("TG", rownames(covariance_matrix))[rowSums(covariance_matrix[grep("TG", rownames(covariance_matrix)),])!=0]
wald.test(b = m1_plusTG$coefficients$fixed[terms], Sigma =  covariance_matrix[terms, terms], Terms=1:length(terms))$result$chi2

##conditional model NEFA effects
m1_plusNEFA<-update(m1, ~.+ NEFA + time:NEFA + time2:NEFA)
x<-as.data.frame(summary(m1_plusNEFA, all=TRUE)$coef.fi)
pval<-as.vector(unname(2*pnorm(-abs(x[,'z ratio']))))
cbind(x, pval) 

####overall Wald chi square tests of NEFA effect on entire trajectory
covariance_matrix<-diag(length(m1_plusNEFA$coefficients$fixed))
covariance_matrix[upper.tri(covariance_matrix, diag=TRUE)]<-m1_plusNEFA$Cfixed
covariance_matrix[lower.tri(covariance_matrix)]<-t(covariance_matrix)[lower.tri(t(covariance_matrix))]
colnames(covariance_matrix)<-names(m1_plusNEFA$coefficients$fixed)
rownames(covariance_matrix)<-names(m1_plusNEFA$coefficients$fixed)
terms<-grep("NEFA", rownames(covariance_matrix))[rowSums(covariance_matrix[grep("NEFA", rownames(covariance_matrix)),])!=0]
wald.test(b = m1_plusNEFA$coefficients$fixed[terms], Sigma =  covariance_matrix[terms, terms], Terms=1:length(terms))$result$chi2

##conditional model Leptin effects
m1_plusLeptin<-update(m1, ~.+ Leptin + time:Leptin + time2:Leptin)
x<-as.data.frame(summary(m1_plusLeptin, all=TRUE)$coef.fi)
pval<-as.vector(unname(2*pnorm(-abs(x[,'z ratio']))))
cbind(x, pval) 

####overall Wald chi square tests of Leptin effect on entire trajectory
covariance_matrix<-diag(length(m1_plusLeptin$coefficients$fixed))
covariance_matrix[upper.tri(covariance_matrix, diag=TRUE)]<-m1_plusLeptin$Cfixed
covariance_matrix[lower.tri(covariance_matrix)]<-t(covariance_matrix)[lower.tri(t(covariance_matrix))]
colnames(covariance_matrix)<-names(m1_plusLeptin$coefficients$fixed)
rownames(covariance_matrix)<-names(m1_plusLeptin$coefficients$fixed)
terms<-grep("Leptin", rownames(covariance_matrix))[rowSums(covariance_matrix[grep("Leptin", rownames(covariance_matrix)),])!=0]
wald.test(b = m1_plusLeptin$coefficients$fixed[terms], Sigma =  covariance_matrix[terms, terms], Terms=1:length(terms))$result$chi2

##conditional model SEX effects
m1_plusSEX<-update(m1, ~.+ SEX + time:SEX + time2:SEX)
x<-as.data.frame(summary(m1_plusSEX, all=TRUE)$coef.fi)
pval<-as.vector(unname(2*pnorm(-abs(x[,'z ratio']))))
cbind(x, pval) 

####overall Wald chi square tests of SEX effect on entire trajectory
covariance_matrix<-diag(length(m1_plusSEX$coefficients$fixed))
covariance_matrix[upper.tri(covariance_matrix, diag=TRUE)]<-m1_plusSEX$Cfixed
covariance_matrix[lower.tri(covariance_matrix)]<-t(covariance_matrix)[lower.tri(t(covariance_matrix))]
colnames(covariance_matrix)<-names(m1_plusSEX$coefficients$fixed)
rownames(covariance_matrix)<-names(m1_plusSEX$coefficients$fixed)
terms<-grep("SEX", rownames(covariance_matrix))[rowSums(covariance_matrix[grep("SEX", rownames(covariance_matrix)),])!=0]
wald.test(b = m1_plusSEX$coefficients$fixed[terms], Sigma =  covariance_matrix[terms, terms], Terms=1:length(terms))$result$chi2

##conditional model TBFM effects
m1_plusTBFM<-update(m1, ~.+ TBFM + time:TBFM+ time2:TBFM)
x<-as.data.frame(summary(m1_plusTBFM, all=TRUE)$coef.fi)
pval<-as.vector(unname(2*pnorm(-abs(x[,'z ratio']))))
cbind(x, pval) 

####overall Wald chi square tests of TBFM effect on entire trajectory
covariance_matrix<-diag(length(m1_plusTBFM$coefficients$fixed))
covariance_matrix[upper.tri(covariance_matrix, diag=TRUE)]<-m1_plusTBFM$Cfixed
covariance_matrix[lower.tri(covariance_matrix)]<-t(covariance_matrix)[lower.tri(t(covariance_matrix))]
colnames(covariance_matrix)<-names(m1_plusTBFM$coefficients$fixed)
rownames(covariance_matrix)<-names(m1_plusTBFM$coefficients$fixed)
terms<-grep("TBFM", rownames(covariance_matrix))[rowSums(covariance_matrix[grep("TBFM", rownames(covariance_matrix)),])!=0]
wald.test(b = m1_plusTBFM$coefficients$fixed[terms], Sigma =  covariance_matrix[terms, terms], Terms=1:length(terms))$result$chi2


##moving to glucose trajectory analysis again....



