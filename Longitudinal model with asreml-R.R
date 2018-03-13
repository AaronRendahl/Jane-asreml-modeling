library(asreml)
library(abind)
library(MIfuns)
library(gdata)
library(haplo.stats)
library(HapEstXXR)
library(synbreed)
library(aod)

pin <-function (object, transform)
{
pframe <- as.list(object$gammas)
names(pframe) <- paste("V", seq(1, length(pframe)), sep = "")
tvalue <- eval(deriv(transform[[length(transform)]], names(pframe)),
pframe)
X <- as.vector(attr(tvalue, "gradient"))
X[object$gammas.type == 1] <- 0
tname <- transform[[2]]
n <- length(pframe)
i <- rep(1:n, 1:n)
j <- sequence(1:n)
k <- 1 + (i > j)
Vmat <- object$ai
se <- sqrt(sum(Vmat * X[i] * X[j] * k))
Z<-tvalue/se
pval<-as.vector(unname(2*pnorm(-abs(Z))))
a<-data.frame(row.names = tname, Estimate = tvalue, SE = se, Z=Z, pval)
print(c(a$Estimate, a$SE, a$Z, a$pval))
}


#########################
####Incretin analysis####
#########################

#read in data
d = read.table("/Users/Nichol/Google Drive/backup April 1 2013/phenotype manuscript/EMS_final_data_2014.txt", header=T)
names(d)

rownames(d)<-d$Number_ID
d$INS_75[d$ID=='423']<-d$INS_OGT[d$ID=='423']
d<-subset (d, d$Incretin_sampling!='n')
d<-subset (d, d$Breed_new!='LR')
d<-droplevels(d)
d$Breed<-d$Breed_new
levels(d$Breed) <-c('Morgan', 'Welsh Pony')

collection2varA<-names(d[,which(names(d)=='Collection2_Mcal'):which(names(d)=='Collection2_ACTH')])
collection2varB<-list()
collection2varB <-abind(lapply(strsplit(collection2varA, '2_'), `[[`, 2))
for (i in 1:length(collection2varA)){
d[,collection2varB[i]][!is.na(d[,collection2varA[i]])]<-d[,collection2varA[i]][!is.na(d[,collection2varA[i]])]
}
d$Month<-relevel(d$Month, ref='DEC')
d$Group<-NULL
d$Group[d$BCS<7 & d$LAM=='n']<-'(-)O(-)PL'
d$Group[d$BCS>=7 & d$LAM=='n']<-'(+)O(-)PL'
d$Group[d$BCS<7 & d$LAM=='y']<-'(-)O(+)PL'
d$Group[d$BCS>=7 & d$LAM=='y']<-'(+)O(+)PL'
d$Group <- relevel(as.factor(d$Group), ref="(-)O(-)PL")
d$GUE_5<-relevel(d$GUE_5, ref='G')
d$INS_Baseline[d$INS_Baseline<1.5]<-1.5
d$INS_30[d$INS_30<1.5]<-1.5
d$GLP1a_30[d$GLP1a_30<1.5]<-1.5
d$GLP1tot_75[d$GLP1tot_75<0.36]<-0.36
d$GLP1tot_120[d$GLP1tot_120<0.29]<-0.29
d$NH_s<-scale(d$NH)[1:nrow(d)]
d$GH_s<-scale(d$GH)[1:nrow(d)]
d$TG_s<-scale(d$TG)[1:nrow(d)]
d$NEFA_s<-scale(d$NEFA)[1:nrow(d)]
d$ACTH_s<-scale(d$ACTH)[1:nrow(d)]
d$Leptin_s<-scale(d$Leptin)[1:nrow(d)]
d$Adiponectin_s<-scale(d$Adiponectin)[1:nrow(d)]
d$AGE_s<-scale(d$AGE)[1:nrow(d)]
d$BCS_s<-scale(d$BCS)[1:nrow(d)]
d$DPP4_activity<-d$DPP4.1
d$DPP4_activity_s<-scale(d$DPP4_activity)
d$Mcal_s<-scale(d$Mcal)
d$CP_s<-scale(d$CP)
d$NDF_s<-scale(d$NDF)
d$Starch_s<-scale(d$Starch)
d$WSC_s<-scale(d$WSC)
d$Owner<-reorder(d$Owner, new.order=c('Miner', 'UCONN', 'RisingSunMorgans_Cindy_Lundgren', 'R2_Ranch_Sharlene_Anderson', 'Thomson_Gail','Wilburn_Ruth', 'Farnley'))
d$Owner<-factor(d$Owner, labels=c('MorganFarm1','MorganFarm2','MorganFarm3','MorganFarm4','WelshFarm1','WelshFarm2','WelshFarm3'))
d$obese<-0
d$obese[d$BCS>=7]<-1
d$obese_factor<-'n'
d$obese_factor[d$BCS>=7]<-'y'
d$LAM_01<-0
d$LAM_01[d$LAM=='y']<-1
d$DPP4_012<-0
d$DPP4_012[d$DPP4=='AG']<-1
d$DPP4_012[d$DPP4=='GG']<-2
d$GCG_012<-0
d$GCG_012[d$GCG=='CG']<-1
d$GCG_012[d$GCG=='GG']<-2
d$GUE1_012<-as.numeric(reorder(d$GUE_1, new.order=c('C','S','G')))-1
d$GUE2_012<-as.numeric(reorder(d$GUE_2, new.order=c('C','M','A')))-1
d$GUE3_012<-as.numeric(reorder(d$GUE_3, new.order=c('A','R')))-1
d$GUE4_012<-as.numeric(reorder(d$GUE_4, new.order=c('G','K','T')))-1
d$GUE5_012<-as.numeric(reorder(d$GUE_5, new.order=c('G','S','C')))-1
d$stallion<-0
d$stallion[d$SEX=='s']<-1
d$mare<-0
d$mare[d$SEX=='m']<-1
d$welsh<-0
d$welsh[d$Breed=='Welsh Pony']<-1
d$Number_ID<-as.factor(d$Number_ID)

#recode proglucagon upstream enhancer genotypes
h_mat<-as.matrix(d[,c('GUE_1','GUE_2','GUE_3','GUE_4','GUE_5', 'GCG')])
h_mat[,'GUE_1']<-replace(h_mat[,'GUE_1'], h_mat[,'GUE_1']=='S', 'C/G')
h_mat<-replace(h_mat, h_mat=='M', 'C/A')
h_mat<-replace(h_mat, h_mat=='R', 'A/G')
h_mat<-replace(h_mat, h_mat=='K', 'G/T')
h_mat[,'GUE_5']<-replace(h_mat[,'GUE_5'], h_mat[,'GUE_5']=='S', 'G/C')
h_mat<-replace(h_mat, h_mat=='AG', 'A/G')
h_mat<-replace(h_mat, h_mat=='AA', 'A/A')
h_mat<-replace(h_mat, h_mat=='CG', 'C/G')
h_mat<-replace(h_mat, h_mat=='CC', 'C/C')
h_mat<-replace(h_mat, h_mat=='GG', 'G/G')
h_mat<-replace(h_mat, h_mat=='A', 'A/A')
h_mat<-replace(h_mat, h_mat=='T', 'T/T')
h_mat<-replace(h_mat, h_mat=='G', 'G/G')
h_mat<-replace(h_mat, h_mat=='C', 'C/C')


#generate haplotypes
gp <- create.gpData(geno=h_mat)
h_mat3<-codeGeno(gp, label.heter="alleleCoding")
geno<-allele1to2(h_mat3$geno)
label<- c('GUE_1','GUE_2','GUE_3','GUE_4','GUE_5', 'GCG')

save.em<-haplo.em(geno=geno, locus.label=label, miss.val=NA)
summary(save.em)
df<-summary(save.em)[1:4]
df<-df[order(-df[,4]),]
df<-df[unique(df[,1]),]
df<-df[order(df[,1]),]
rownames(df)<-rownames(h_mat)
df<-df[,-1]
df_2<-cbind(save.em$haplotype, save.em$hap.prob)
df_2<-df_2[df_2[,'save.em$hap.prob']>=0.03,]
d$hap1<-0
d$hap1[which(df[,1]==1, arr.ind=TRUE)]<-1
d$hap1[which(df[,2]==1, arr.ind=TRUE)]<-1
d$hap1[which(df[,1]==1, arr.ind=TRUE)[which(df[,1]==1, arr.ind=TRUE)%in%which(df[,2]==1, arr.ind=TRUE)]]<-2
d$hap2<-0
d$hap2[which(df[,1]==3, arr.ind=TRUE)]<-1
d$hap2[which(df[,2]==3, arr.ind=TRUE)]<-1
d$hap2[which(df[,1]==3, arr.ind=TRUE)[which(df[,1]==3, arr.ind=TRUE)%in%which(df[,2]==3, arr.ind=TRUE)]]<-2
d$hap3<-0
d$hap3[which(df[,1]==8, arr.ind=TRUE)]<-1
d$hap3[which(df[,2]==8, arr.ind=TRUE)]<-1
d$hap3[which(df[,1]==8, arr.ind=TRUE)[which(df[,1]==8, arr.ind=TRUE)%in%which(df[,2]==8, arr.ind=TRUE)]]<-2
d$hap4<-0
d$hap4[which(df[,1]==11, arr.ind=TRUE)]<-1
d$hap4[which(df[,2]==11, arr.ind=TRUE)]<-1
d$hap4[which(df[,1]==11, arr.ind=TRUE)[which(df[,1]==11, arr.ind=TRUE)%in%which(df[,2]==11, arr.ind=TRUE)]]<-2
d$hap5<-0
d$hap5[which(df[,1]==6, arr.ind=TRUE)]<-1
d$hap5[which(df[,2]==6, arr.ind=TRUE)]<-1
d$hap5[which(df[,1]==6, arr.ind=TRUE)[which(df[,1]==6, arr.ind=TRUE)%in%which(df[,2]==6, arr.ind=TRUE)]]<-2

hap.dat<-haplo.design(save.em, haplo.effect="additive", hapcodes=NA, min.count=round(nrow(d)*2*0.00), haplo.base=NA)
hap.dat<-as.data.frame(hap.dat)
head(hap.dat)
colSums(hap.dat)
rowSums(hap.dat)
save.em$hap.prob
unname(which((rowSums(hap.dat))>1.8))
hap.dat2<-hap.dat[,which(save.em$hap.prob>=0.05)]
colnames(hap.dat2)<-c('hapA', 'hapB', 'hapC','hapD')
hap.dat2$hapE<-rowSums(hap.dat[,which(save.em$hap.prob<0.05)])  ##probability of rare haplotype
d<-cbind(d, hap.dat2)
table(d$hapA, d$Owner)




#generate AUC
xx<-as.data.frame(cbind(c(d$GLU_Baseline, d$GLU_15, d$GLU_30, d$GLU_60, d$GLU_75, d$GLU_90, d$GLU_120),
rep(c(0,15,30,60,75,90,120), each=nrow(d)),
rep(paste('s',d$Number_ID,sep="_"), 7)))
colnames(xx)<-c('trait', 'time', 'ID')
xx[,1]<-as.numeric(as.character(xx[,1]))
xx[,2]<-as.numeric(as.character(xx[,2]))
xx[,3]<-as.factor(xx[,3])
d$GLU_AUC<-c(AUC(xx, time = 'time', id = 'ID', dv = 'trait')[2])$AUC
xx<-as.data.frame(cbind(c(d$INS_Baseline, d$INS_15, d$INS_30, d$INS_60, d$INS_75, d$INS_90, d$INS_120),
rep(c(0,15,30,60,75,90,120), each=nrow(d)),
rep(paste('s',d$Number_ID,sep="_"), 7)))
colnames(xx)<-c('trait', 'time', 'ID')
xx[,1]<-as.numeric(as.character(xx[,1]))
xx[,2]<-as.numeric(as.character(xx[,2]))
xx[,3]<-as.factor(xx[,3])
d$INS_AUC<-c(AUC(xx, time = 'time', id = 'ID', dv = 'trait')[2])$AUC
xx<-as.data.frame(cbind(c(d$GLP1a_Baseline, d$GLP1a_15, d$GLP1a_30, d$GLP1a_60, d$GLP1a_75, d$GLP1a_90, d$GLP1a_120),
rep(c(0,15,30,60,75,90,120), each=nrow(d)),
rep(paste('s',d$Number_ID,sep="_"), 7)))
colnames(xx)<-c('trait', 'time', 'ID')
xx[,1]<-as.numeric(as.character(xx[,1]))
xx[,2]<-as.numeric(as.character(xx[,2]))
xx[,3]<-as.factor(xx[,3])
d$GLP1a_AUC<-c(AUC(xx, time = 'time', id = 'ID', dv = 'trait')[2])$AUC
xx<-as.data.frame(cbind(c(d$GLP1tot_Baseline, d$GLP1tot_15, d$GLP1tot_30, d$GLP1tot_60, d$GLP1tot_75, d$GLP1tot_90, d$GLP1tot_120),
rep(c(0,15,30,60,75,90,120), each=nrow(d)),
rep(paste('s',d$Number_ID,sep="_"), 7)))
colnames(xx)<-c('trait', 'time', 'ID')
xx[,1]<-as.numeric(as.character(xx[,1]))
xx[,2]<-as.numeric(as.character(xx[,2]))
xx[,3]<-as.factor(xx[,3])
d$GLP1tot_AUC<-c(AUC(xx, time = 'time', id = 'ID', dv = 'trait')[2])$AUC


#reshape data from wide to long format
phenotype_levels<-c("GLU_Baseline", "GLU_15", "GLU_30", "GLU_60", "GLU_75", "GLU_90", "GLU_120", "INS_Baseline", "INS_15", "INS_30", "INS_60", "INS_75", "INS_90", "INS_120", "GLP1a_Baseline", "GLP1a_15", "GLP1a_30", "GLP1a_60", "GLP1a_75", "GLP1a_90", "GLP1a_120", "GLP1tot_Baseline", "GLP1tot_15", "GLP1tot_30", "GLP1tot_60", "GLP1tot_75", "GLP1tot_90", "GLP1tot_120")
df<-d[,c('Number_ID', 'AGE', 'Breed', 'SEX', 'LAM', 'obese_factor', 'Owner', 'DPP4_activity', 'Mcal', 'CP', 'NDF', 'Starch', 'WSC', 'DPP4_012', 'GCG_012', 'DPP4', 'GCG', 'GUE1_012','GUE2_012','GUE3_012','GUE4_012','GUE5_012', 'hapA', 'hapB', 'hapC', 'hapD','hapE', "NH", "GH", "TG", "NEFA" , "ACTH" , "Leptin", "Adiponectin", "Group", phenotype_levels)]
df2<-reshape(df, idvar='Number_ID', timevar='time',direction='long', varying=colnames(df)[(which(colnames(df)=='GLU_Baseline')):(which(colnames(df)=='GLP1tot_120'))], sep='_')
rownames(df2)<-c(1:nrow(df2))
df2$Number_ID<-as.factor(df2$Number_ID)
df2$time <-as.factor(df2$time)
df2$time <-reorder(df2$time, new.order=c('Baseline','15','30','60','75','90','120'))
df2$time_min<-as.numeric(gsub('Baseline', '0', df2$time))
df2$time_min2<-df2$time_min^2
df2$time2<-as.numeric(gsub('Baseline', '0', df2$time))
df2$Breed<-as.factor(df2$Breed)
df2$obese_factor<-as.factor(df2$obese_factor)
df2<-df2[order(df2$Number_ID),]
df2<-droplevels(df2)
df2$minutes<-0
df2$minutes[df2$time!='Baseline']<-as.numeric(as.character(df2$time[df2$time!='Baseline']))
df2$minutes<-as.numeric(df2$minutes/15)
df2$minutes2<-df2$minutes^2
df2$GLU_s<-scale(df2$GLU)
df2$INS_s<-scale(df2$INS)
df2$logINS_s<-scale(log10(df2$INS))
df2$GLP1a_s<-scale(df2$GLP1a)
df2$GLP1tot_s<-scale(df2$GLP1tot)
ID_list<-unique(df2$Number_ID)
time_list<-unique(df2$time)
time_list<-gsub('Baseline', '0', df2$time)
df2$time2_s<-scale(df2$time2)
df2$time_squared<-(df2$time2^2)
df2$time_squared_s<-scale(df2$time_squared)
df2$time_cubed<-(df2$time2^3)
df2$Mcal_s<-scale(df2$Mcal)
df2$CP_s<-scale(df2$CP)
df2$NDF_s<-scale(df2$NDF)
df2$Starch_s<-scale(df2$Starch)
df2$WSC_s<-scale(df2$WSC)
df2$AGE_s<-scale(df2$AGE)
df2$TG_s<-scale(df2$TG)
df2$NEFA_s<-scale(df2$NEFA)
df2$ACTH_s<-scale(df2$ACTH)
df2$Leptin_s<-scale(df2$Leptin)
df2$Adiponectin_s<-scale(df2$Adiponectin)
df2$DPP4_activity_s<-scale(df2$DPP4_activity)
df2$GLP1_diff<-df2$GLP1tot-df2$GLP1a

##Unconditional model (quadratic fixed and random time effects, random farm intercept)
m1 <- asreml(GLU ~ time_min + time_min2, data=df2, random=~str(~Number_ID + time_min:Number_ID + time_min2:Number_ID, ~us(3):id(Number_ID)) + ~Owner, rcov=~units, maxit=400,  Cfixed=TRUE)

##example of first horse timepoint 0, = horse's value = grand mean + individual random intercept + owner intercept + residual 
i<-1
m1$coefficients$fixed['(Intercept)'] + m1$coefficients$random['Number_ID_190'] + m1$coefficients$random['Owner_MorganFarm2'] + m1$residual[i]
df2$GLU[i]
##example of first horse timepoint 15
i<-2
m1$coefficients$fixed['(Intercept)'] + m1$coefficients$random['Number_ID_190']+ ((unique(df2$time_min)[i])*(m1$coefficients$fixed['time_min']+ m1$coefficients$random['time_min:Number_ID_190'])) + ((unique(df2$time_min2)[i])*(m1$coefficients$fixed['time_min2']+ m1$coefficients$random['time_min2:Number_ID_190'] ))+ m1$coefficients$random['Owner_MorganFarm2'] + m1$residual[i]
df2$GLU[i]
m1$fitted[i]+ m1$residual[i]
##example of first horse timepoint 30
i<-3
m1$coefficients$fixed['(Intercept)'] + m1$coefficients$random['Number_ID_190']+ ((unique(df2$time_min)[i])*(m1$coefficients$fixed['time_min']+ m1$coefficients$random['time_min:Number_ID_190'])) + ((unique(df2$time_min2)[i])*(m1$coefficients$fixed['time_min2']+ m1$coefficients$random['time_min2:Number_ID_190'] ))+ m1$coefficients$random['Owner_MorganFarm2'] + m1$residual[i]
df2$GLU[i]
m1$fitted[i]+ m1$residual[i]

##Conditional model (estimate effects of age)
m1_plusAge<-update(m1, ~.+ AGE + time_min:AGE + time_min2:AGE)
x<-as.data.frame(summary(m1_plusAge, all=TRUE)$coef.fi)
pval<-as.vector(unname(2*pnorm(-abs(x[,'z ratio']))))
cbind(x, pval) #estimates of age on the intercept (AGE), initial rate of change (time_min:AGE), and deceleration rate (time_min2:AGE2)

####overall Wald chi square tests of Age effect on entire trajectory
covariance_matrix<-diag(length(m1_plusAge$coefficients$fixed))
covariance_matrix[upper.tri(covariance_matrix, diag=TRUE)]<-m1_plusAge$Cfixed
covariance_matrix[lower.tri(covariance_matrix)]<-t(covariance_matrix)[lower.tri(t(covariance_matrix))]
colnames(covariance_matrix)<-names(m1_plusAge$coefficients$fixed)
rownames(covariance_matrix)<-names(m1_plusAge$coefficients$fixed)
terms<-grep("AGE", rownames(covariance_matrix))[rowSums(covariance_matrix[grep("AGE", rownames(covariance_matrix)),])!=0]
wald.test(b = m1_plusAge$coefficients$fixed[terms], Sigma =  covariance_matrix[terms, terms], Terms=1:length(terms))$result$chi2



##Conditional model (estimate effects of breed)
m1_plusBreed<-update(m1, ~.+ Breed + time_min:Breed + time_min2:Breed)
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



##Obtain and plot predicted trajectories 
time1<-predict(m1_plusBreed, classify="time_min2:time_min:Breed", levels = list(time_min=df2$time_min[1], time_min2=df2$time_min2[1]))$predictions$pvals
time2<-predict(m1_plusBreed, classify="time_min2:time_min:Breed", levels = list(time_min=df2$time_min[2], time_min2=df2$time_min2[2]))$predictions$pvals
time3<-predict(m1_plusBreed, classify="time_min2:time_min:Breed", levels = list(time_min=df2$time_min[3], time_min2=df2$time_min2[3]))$predictions$pvals
time4<-predict(m1_plusBreed, classify="time_min2:time_min:Breed", levels = list(time_min=df2$time_min[4], time_min2=df2$time_min2[4]))$predictions$pvals
time5<-predict(m1_plusBreed, classify="time_min2:time_min:Breed", levels = list(time_min=df2$time_min[5], time_min2=df2$time_min2[5]))$predictions$pvals
time6<-predict(m1_plusBreed, classify="time_min2:time_min:Breed", levels = list(time_min=df2$time_min[6], time_min2=df2$time_min2[6]))$predictions$pvals
time7<-predict(m1_plusBreed, classify="time_min2:time_min:Breed", levels = list(time_min=df2$time_min[7], time_min2=df2$time_min2[7]))$predictions$pvals
res<-as.data.frame(rbind(time1, time2, time3,time4, time5, time6,time7))
res$lower<-res$predicted.value - (res$standard.error)
res$upper<-res$predicted.value + (res$standard.error)


res2<-res[res$Breed=='Morgan',]
ylim_min<-min(res$lower)-3
ylim_max<-max(res$upper)+3
x <-  c( 0,  15,  30,  60,  75,  90, 120)
xl <- seq(min(x),max(x), (max(x) - min(x))/1000)
plot(x=x, y=res2$predicted.value, col='white', ylim=c(ylim_min,ylim_max), xlab='minutes post OST', ylab="glucose mg/dl", bty='n')
matlines(x=xl, y=cbind(predict(loess(res2$predicted.value~x),xl), predict(loess(res2$lower~x),xl),predict(loess(res2$upper~x),xl)), col='bisque4', lty=c(1,2,2), lwd=c(3,1,1))
points(x=c(0,15,30,45,60,75,90,120), y=predict(loess(res2$predicted.value~x),xl)[match(c(0,15,30,45,60,75,90,120), xl)], pch=0, col='bisque4', cex=1.5)
res2<-res[res$Breed=='Welsh Pony',]
matlines(x=xl, y=cbind(predict(loess(res2$predicted.value~x),xl), predict(loess(res2$lower~x),xl),predict(loess(res2$upper~x),xl)), col='deepskyblue1', lty=c(1,2,2), lwd=c(3,1,1))
points(x=c(0,15,30,45,60,75,90,120), y=predict(loess(res2$predicted.value~x),xl)[match(c(0,15,30,45,60,75,90,120), xl)], pch=15, col='deepskyblue1', cex=1.5)
legend("topleft", levels(res$Breed), col=c('bisque4','deepskyblue1'), lty=1,lwd=2,pch=c(0,15), xpd = TRUE, inset=c(0,-0), horiz = FALSE, bty = "n",  cex =1)


##########################################################################################################################################
###Unconditional model multivariate response correlation#################
##########################################################################################################################################
m1 <- asreml(cbind(GLU, INS) ~ trait + trait:time_min + trait:time_min2, data=df2, random=~str(~trait:Number_ID + trait:time_min:Number_ID + trait:time_min2:Number_ID, ~us(6):id(Number_ID))+ us(trait):Owner, rcov=~Number_ID:time:us(trait), maxit=400) 
m1<-update(m1)

##converting time to number of 15 minute intervals and a scaled response improves convergence speed
m1 <- asreml(cbind(GLU_s, INS_s) ~ trait + trait:minutes + trait:minutes2, data=df2, random=~str(~trait:Number_ID + trait:minutes:Number_ID + trait:minutes2:Number_ID, ~us(6):id(Number_ID))+ us(trait):Owner, rcov=~Number_ID:time:us(trait), maxit=400) 


##covarTOcorr_matrix function to obtain correlation estimates from covariance estimates
# x= #number of parameter estimates
# y<- indicator for farm, individual, and residual 
# z<- asreml model object
# trait_names<-c('GLU', 'INS')

covarTOcorr_matrix<-function (x, y, z, trait_names){	
	diag_size<-x
	T<-diag(diag_size)
	length_param<-length(grep(y, names(z$gammas)))
	T[upper.tri(T,diag=T)] <-seq(1:length_param)
	T<-t(T)
	list1<-T[lower.tri(T)]
	list2<-list()
	for (i in 1:length(list1)){list2[i]<-unname(which(T==list1[i], arr.ind=TRUE)[1,2])}
	list2<-abind(list2)
	list4<-diag(T)[list2]
	list3<-list()
	for (i in 1:length(list1)){list3[i]<-unname(which(T==list1[i], arr.ind=TRUE)[1,1])}
	list3<-abind(list3)
	list7<-diag(T)[list3]
	res<-mat.or.vec(length(list1),4)
	pframe <- as.list(z$gammas[grep(y, names(z$gammas))])
	pf_names<-names(pframe)[list1]
	if (length(grep(x,pf_names,fixed=TRUE))>1){
		AA<-as.numeric(sapply(strsplit(sapply(strsplit(pf_names, ').',fixed=TRUE),"[[",2),":"),"[[",1))
		BB<-as.numeric(sapply(strsplit(sapply(strsplit(pf_names, ').',fixed=TRUE),"[[",2),":"),"[[",2))
		pf_names<-paste(paste(rep(trait_names, 3),rep(c('Intercept','Linear','Quad'), each=length(trait_names)),sep='_')[AA], paste(rep(trait_names, 	3),rep(c('Intercept','Linear','Quad'), each=length(trait_names)),sep='_')[BB],sep=":")
	}
	names(pframe) <- paste("V", grep(y, names(z$gammas)), sep = "")	

	for (ii in 1:length(list1)){
		i<-list1[ii]
		j<-list4[ii]
		k<-list7[ii]
		transform<-as.formula(paste("gen.cor ~ ",as.name(names(pframe)[i]), "/sqrt(", as.name(names(pframe)[j]), "*", as.name(names(pframe)[k]), ")" ))
		res[ii,]<-pin(z, transform)
		}

	res<-as.data.frame(abind(res))
	rownames(res)<-pf_names
	res<-as.data.frame(res)
	res$star<-''
	res$star[res[,4]<0.05]<-'*'
	res$star[res[,4]<0.01]<-'**'
	res$star[res[,4]<0.001]<-'***' 
	res$result<-gsub("\\( ", "(", paste(format(round(res[,1],2),nsmall=2), '(', format(round(res[,2],2),nsmall=2), ')', res$star, sep=""))
	colnames(res)<-c('correlation_estimate', 'SE', 'Z', 'pval', 'star', 'result')
	return(res)
}

Owner_covar<-covarTOcorr_matrix(2,'trait:Owner!trait.', m1, trait_names=c('GLU', 'INS'))
Owner_covar

Individual_between_covar<-covarTOcorr_matrix(6,'minutes', m1, trait_names=c('GLU', 'INS'))
Individual_between_covar

Individual_within_covar<-covarTOcorr_matrix(2,'R!trait.', m1, trait_names=c('GLU', 'INS'))
Individual_within_covar
