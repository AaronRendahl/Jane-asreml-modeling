library(asreml)
library(aod)

OST_R_asreml_NS_hr <- read.csv("OST_R_asreml_NS_hr.csv", check.names=FALSE)
OST_R_asreml_NS_hr$Number_ID <- factor(OST_R_asreml_NS_hr$Number_ID)
OST_R_asreml_NS_hr$logINS <- log(OST_R_asreml_NS_hr$insulin)

## glucose, without and with Owner
glu <- asreml(glucose ~ time + time2, data=OST_R_asreml_NS_hr, random=~str(~Number_ID + time:Number_ID + time2:Number_ID, ~us(3):id(Number_ID)), rcov=~units, maxit=400,  Cfixed=TRUE)
glu_own <- asreml(glucose ~ time + time2, data=OST_R_asreml_NS_hr, random=~str(~Number_ID + time:Number_ID + time2:Number_ID, ~us(3):id(Number_ID)) + ~Owner, rcov=~units, maxit=400,  Cfixed=TRUE)

## logIns, without and with Owner
ins <- asreml(logINS ~ time + time2, data=OST_R_asreml_NS_hr, random=~str(~Number_ID + time:Number_ID + time2:Number_ID, ~us(3):id(Number_ID)), rcov=~units, maxit=400,  Cfixed=TRUE)
ins_own <- asreml(logINS ~ time + time2, data=OST_R_asreml_NS_hr, random=~str(~Number_ID + time:Number_ID + time2:Number_ID, ~us(3):id(Number_ID)) + ~Owner, rcov=~units, maxit=400,  Cfixed=TRUE)


## set m1 to desired model before continuing


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

####overall Wald chi square tests of breed effect on entire trajectory
covariance_matrix<-diag(length(m1_plusBreed$coefficients$fixed))
covariance_matrix[upper.tri(covariance_matrix, diag=TRUE)]<-m1_plusBreed$Cfixed
covariance_matrix[lower.tri(covariance_matrix)]<-t(covariance_matrix)[lower.tri(t(covariance_matrix))]
colnames(covariance_matrix)<-names(m1_plusBreed$coefficients$fixed)
rownames(covariance_matrix)<-names(m1_plusBreed$coefficients$fixed)
terms<-grep("Breed", rownames(covariance_matrix))[rowSums(covariance_matrix[grep("Breed", rownames(covariance_matrix)),])!=0]
wald.test(b = m1_plusBreed$coefficients$fixed[terms], Sigma =  covariance_matrix[terms, terms], Terms=1:length(terms))$result$chi2

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
