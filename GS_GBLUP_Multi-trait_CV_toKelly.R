###PART1 The main job is to prepare data format as required by ASreml 

# This is loading my data
#Y: phenotypes, nxt 
#A: marker relationship matrix estimated from A.mat(), nxn
#create a dummy pedigree

install.packages("matrixcalc")
load("Africa_X,A,y_Data.rdata")
X.t<-t(X.africa)  ## marker data
Pop<-"Afri"    ## Study name
Y<-Y            ## phenotypic data

colnames(Y)[1]<-"Taxa"
Pop.asr<-merge(Y,X.t,by.x="Taxa",by.y="row.names",all=TRUE)
Pop.asr.sort<-Pop.asr[order(match(Pop.asr$Taxa,rownames(Amat))),]
colnames(Pop.asr.sort)[1]<-"Genot"
N<-nrow(Pop.asr.sort)

library(matrixcalc)
is.positive.definite(Amat)

###Function to make matrix into sparse format
mat2sparse <- function (X, rowNames = dimnames(X)[[1]])
{
  which <- (X != 0 & lower.tri(X, diag = TRUE))
  df <- data.frame(row = t(row(X))[t(which)], col = t(col(X))[t(which)],
                   val = t(X)[t(which)])
  if (is.null(rowNames))
    rowNames <- as.character(1:nrow(X))
  attr(df, "rowNames") <- rowNames
  df
}


#### Preparing AHAT matrix, making positive definite
diag(Amat) <- diag(Amat) + 1e-6
is.positive.definite(Amat)

AHAT.inv<-solve(Amat)
AHAT.inv[1:5,1:5]
det(AHAT.inv)
AHAT.inv.sparse<-mat2sparse(AHAT.inv)  #lower diag sparse matrix
colnames(AHAT.inv.sparse)<-c('Row','Column','Ainverse')
head(AHAT.inv.sparse)
write.table(AHAT.inv.sparse,file=paste(Pop,"_Amat-Inv-sparse.txt",sep=""))


##### Creating a dummy pedigree file and save it
peddummy<-matrix(ncol=3,nrow=N)
colnames(peddummy)<-c("Individual","Female","Male")
peddummy[,1]<-rownames(Amat)
peddummy[,2]<-rep(0,length=N)
peddummy[,3]<-rep(0,length=N)
peddummy<-as.data.frame(peddummy)
rownames(peddummy)<-rownames(Amat)
write.table(peddummy,file=paste(Pop,"_peddummy.txt",sep=""))


gmatrix<-data.frame(AHAT.inv.sparse)
ainvped<-asreml.Ainverse(peddummy)$ginv
attr(gmatrix,"rowNames")<-attr(ainvped,"rowNames")
colnames(Y)[1]<-"Genot"
Y$Genot<-as.factor(Y$Genot)


###PART2 decide your own model. Onetime predictions:
### Multi traits model
modelGBLUP<-asreml(fixed=cbind(TZ_01,TZ_12,UG_05,UG_11,KE_37,BF_07)~trait,rcov=~units:us(trait),random=~corh(trait):giv(Genot),ginverse=list(Genot=gmatrix),workspace=128e06,na.method.Y="include",data=Y)

GenVar<-summary(modelGBLUP)$varcomp
write.csv(GenVar,paste(Pop,"_GBLUP_MultiTrait_Genetic_Variance_Components_Onetime.csv",sep=""))

predGBLUP<-predict(modelGBLUP,classify="Genot")$predictions$pvals
predGBLUP<-predGBLUP[order(match(predGBLUP$Genot,Y$Genot)),]
cor.onetime<-matrix(nrow=traits,ncol=1)
rownames(cor.onetime)<-colnames(Y[,-1])

for (i in 1:traits){
    cor.onetime[i,]<-cor(predGBLUP[,2],Y[,-1][,i],use="complete")
}
    round(cor.onetime,digit=2)

write.csv(predGBLUP,paste0("GBLUP_MultiTrait_PredGEBV_Onetime.csv"))
write.csv(cor.onetime,paste0("GBLUP_MultiTrait_cor_Onetime.csv"))


#### Additional Part
#### to do Cross-validation, read in Sampling order file 
sample<-read.csv(paste("sample_5000cycles_sets_order_",Pop,".csv",sep=""),sep=",",header=T,row.names=1)

cycles=50 ## U decide !!!
traits= ncol(Y)-1 ## U decide !!!
sample2<-sample[,1:cycles]
cor<-matrix(nrow=cycles,ncol=traits)
y_hat<-matrix(nrow=nrow(Y),ncol=cycles)
rownames(y_hat)<-Y[,1]
colnames(cor)<-colnames(Y[,-1])

folds=10 ## U decide !!!
dim(sample2)

library(asreml)
for (i in 1:cycles){
  tmp<-NULL
  for(fold in 1:folds){
    tst<-which(sample2[,i]==fold)
    Y2<-Y
    Y2[tst,c(2:ncol(Y2))]=NA  
    modelGBLUP<-asreml(fixed=cbind(TZ_01,TZ_12,UG_05,UG_11,KE_37,BF_07)~trait,rcov=~units:us(trait),random=~corh(trait):giv(Genot),ginverse=list(Genot=gmatrix),workspace=128e06,na.method.Y="include",data=Y2)
    
    predGBLUP<-predict(modelGBLUP,classify="Genot")$predictions$pvals ### Get the predGEBVs
    preds<-as.vector(predGBLUP[,2][tst])	
    names(preds)<-predGBLUP[,1][tst]
    tmp<-c(tmp,preds)
  } 
  for (j in 1:traits){
    y<-Y[,-1][,j]     
    names(y)<-Y[,1]   
    tmp2<-tmp[names(y)]
    y_hat[,i]<-tmp2          ### Each ith cycle/col is a new cycle
    cor[i,j]<-cor(tmp2,y,use="complete")    ### Each ith row is a cycle, jth trait
  }
  print (i)
}

cor.Model<-cor
colnames(cor.Model)<-colnames(Y[,-1])
write.csv(cor.Model,file=paste0("cor.Model_use_same_order_",cycles,"cyc.csv"))

cor.mean<-colMeans(cor.Model)
write.csv(cor.mean,paste0("cor.Model_",cycles,"cyc_mean.csv"))

r<-as.data.frame(rowMeans(y_hat)) 
write.csv(r,file=paste0("y_hat_",cycles,"cycle_mean_all_traits_ASReml.csv"))
######## CV finish


###Try different error structure
#modelGBLUP<-asreml(fixed=cbind(TZ_01,TZ_12,UG_05,UG_11,KE_37,BF_07)~trait,rcov=~units:us(trait),random=~corh(trait):giv(Genot,var=T),ginverse=list(Genot=gmatrix),workspace=128e06,na.method.Y="include",data=Y2)
#modelGBLUP<-asreml(fixed=cbind(TZ_01,TZ_12,UG_05,UG_11,KE_37,BF_07)~trait,rcov=~units:diag(trait),random=~giv(Genot):corgh(trait),ginverse=list(Genot=gmatrix),workspace=128e06,na.method.Y="include",data=Y)

###Single trait model, the onetime prediction is the same as the GBLUP model specified in rrBLUP()
modelGBLUP<-asreml(fixed=Pheno.PC1~1,random=~giv(Genot),ginverse=list(Genot=gmatrix),workspace=128e06,na.method.Y="include",data=Y)
GBLUP<-summary(modelGBLUP,all=TRUE)$coef.random


#### Estimate repeatability of the GEBVs
summary(modelGBLUP)$varcomp
(h2_GBLUP<-nadiv:::pin(modelGBLUP,h2_ind~V1/(V1+V2)))