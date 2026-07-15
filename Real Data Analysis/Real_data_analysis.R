setwd('C:/Dropbox Backup/Sequential mediation project/2024/Real_Data')
set.seed(1234)

jaccard_index <- function(set1, set2) {
  if(length(union(set1, set2))!=0){
    length(intersect(set1, set2)) / length(union(set1, set2))
  }else{
    0
  }
}

otsuka_ochiai <- function(set1, set2) {
  intersection_size <- length(intersect(set1, set2))
  sqrt_prod <- sqrt(length(set1) * length(set2))
  
  if(sqrt_prod!=0){
    intersection_size / sqrt_prod
  }else{
    0
  }
  
}

# Function to calculate Sørensen–Dice coefficient
sorensen_dice <- function(set1, set2) {
  intersection_size <- length(intersect(set1, set2))
  if((length(set1) + length(set2))!=0){
    2 * intersection_size / (length(set1) + length(set2))
  }else{
    0
  }
}

Pairwise_index <- function(B) {
  n = dim(B)[1]
  JI = c()
  SD = c()
  OO = c()
  
  for(i in 1:n){
  
    for(j in setdiff(i:n,i)){
      JI = c(JI,jaccard_index(which(B[i,]!=0),which(B[j,]!=0)))
      SD = c(SD,sorensen_dice(which(B[i,]!=0),which(B[j,]!=0)))
      OO = c(OO,otsuka_ochiai(which(B[i,]!=0),which(B[j,]!=0)))
    }
  
  }
  return(data.frame(JI,SD,OO))
}
#### LOAD DATA ######################

## load ID
ID=read.table(file='sampleID.csv',header = TRUE,sep=',')
dim(ID) # 75 subjects

## load microbiome
OTU=read.table(file='otu.good.csv',header = TRUE,sep=',')
OTU[1:5,1:5] # read counts
dim(OTU) # 135 samples (some subj have multiple) by 140 OTUs

taxonomy=read.table(file='tax.good.csv',header = TRUE,sep=',')

OTU1 = as.data.frame(t(OTU[,-1]))
OTU1$X = rownames(OTU1)
OTU2 = merge(taxonomy,OTU1,by='X')

OTU_family = aggregate(. ~ Family, data = OTU2[, -c(1:5,7)], FUN = sum)
#OTU_family = aggregate(. ~ Genus, data = OTU2[, -c(1:6)], FUN = sum)

OTU3 = as.data.frame(t(OTU_family[,-1])); colnames(OTU3) = OTU_family$Family
OTU = cbind(X = OTU$X, OTU3)

## load metabolome
metabo.neg=read.table(file='Huang-data_IMPUTED-LOG-UNSCALED-ADJUSTED_EX01199_RP-NEG_4BINNER_20220313_173830-sal-rem_Report_PP-Named.csv',header = TRUE,sep=',')
metabo.pos=read.table(file='Huang-data_IMPUTED-LOG-UNSCALED-ADJUSTED_EX01199_RP-POS_4BINNER_20220313_173119-sal-rem_Report_PP-Named.csv',header = TRUE,sep=',')
head(metabo.neg)
dim(metabo.neg) # 75*404
dim(metabo.pos) # 75*589

## load meta data
meta=read.table(file='Final_Samples_PlasmaMetabUntarg_clinicalmd_YH.csv',header = TRUE,sep=',')
dim(meta) # 75*17
colnames(meta)

##### Process data 
# match samples
commonID=ID$Gut_16S_data_identifer[ID$Gut_16S_data_identifer!=""]
length(commonID) # n=55 
mch1=match(commonID,OTU$X)
myOTU=OTU[mch1,-1]
rownames(myOTU)=commonID

caarsID=as.numeric(gsub("\\D", "", metabo.neg$CAARS_ID)) # CAARS ID
mch2=match(ID$CAARS_ID[match(commonID,ID$Gut_16S_data_identifer)],caarsID)
mymetabo.neg=metabo.neg[mch2,4:ncol(metabo.neg)]
rownames(mymetabo.neg)=commonID

caarsID=as.numeric(gsub("\\D", "", metabo.pos$CAARS_ID)) # CAARS ID
mch3=match(ID$CAARS_ID[match(commonID,ID$Gut_16S_data_identifer)],caarsID)
mymetabo.pos=metabo.pos[mch3,4:ncol(metabo.pos)]
rownames(mymetabo.pos)=commonID

mch4=match(ID$CAARS_ID[match(commonID,ID$Gut_16S_data_identifer)],meta$CAARS_ID)
mymeta=meta[mch4,]

mymetabo = cbind(mymetabo.pos,mymetabo.neg)

clrX=as.matrix(myOTU)
for(i in 1:nrow(clrX)){
  temp=as.numeric(myOTU[i,])
  temp[temp==0]=1
  clrX[i,]=scale(log(temp),center=TRUE,scale=FALSE)
}
clr.OTU=as.data.frame(clrX)

# center/scale/combine metabolome data
variances <- apply(mymetabo.neg, 2, var, na.rm = TRUE)
top50_indices <- order(variances, decreasing = TRUE)[1:100]
mymetabo.neg_top50 <- mymetabo.neg[, top50_indices]

variances <- apply(mymetabo.pos, 2, var, na.rm = TRUE)
top50_indices <- order(variances, decreasing = TRUE)[1:100]
mymetabo.pos_top50 <- mymetabo.pos[, top50_indices]


#mymetabo=cbind(scale(mymetabo.neg_top50,center=TRUE,scale=TRUE),
#               scale(mymetabo.pos_top50,center=TRUE,scale=TRUE))

mymetabo=cbind(scale(mymetabo.neg_top50,center=TRUE,scale=TRUE),
               scale(mymetabo.pos_top50,center=TRUE,scale=TRUE))

mymetabo = as.matrix(mymetabo)
clr.OTU = as.matrix(clr.OTU)

Leave_One_Out_Analysis<-function(M1,M2=NULL,Y,Integfactor=F,Factor=F,MV=F){
  X = if(is.null(M2)) M1 else cbind(M1,M2)
  n = length(Y); p = dim(X)[2]
  B = matrix(0,1,p)
  Y_pre = c()
  Var = c()
  coefM = if(MV) matrix(0,dim(M1)[2],dim(M2)[2]) else 0
  
  nusiance = c()
  for(i in 1:n){
    Y1 = Y[-i]; X1 = X[-i,]
    Xtest = X[i,]
    
    if(Factor){
      fc2 = Factor_decompose(list(cbind(M1[-i,],M2[-i,])))
      X1 = cbind(as.matrix(fc2$U[[1]]),as.matrix(fc2$Fm[[1]]),1)
      nusiance = c(nusiance,dim(fc2$Fm[[1]])[2])
      
    }else if(Integfactor){
      fc1 = Factor_decompose(list(M1[-i,],M2[-i,]))
      F1 = do.call(cbind,fc1$Fm)
      U1 = do.call(cbind,fc1$U)
      X1 = cbind(U1,F1,1)
      
      nusiance = c(nusiance,dim(F1)[2])
    }else if(MV){
      rfit <- with(list(otu = M1[-i,],meta = M2[-i,]), rrr(meta, otu, maxrank = 10))
      coef1 = coef(rfit)
      rank1 = rfit$rank
      
      coefM = coefM + coef1
      
      Var = c(Var,Frobenius_Norm(M1[-i,] %*% coef1)/Frobenius_Norm(M2[-i,]))
      
      s = svd(as.matrix(M1[-i,]) %*% coef1)
      Res = as.matrix(M2[-i,]) - as.matrix(M1[-i,]) %*% coef1
      Pro = (diag(dim(M1[-i,])[1]) - s$u[,1:rank1] %*% t(s$u[,1:rank1])) %*% as.matrix(M1[-i,])
      X1 = cbind(Pro,Res,s$u[,1:rank1],1)
      
      nusiance = c(nusiance,rank1)
    }else{
      X1 = cbind(X,1)[-i,]
    }
    r = dim(X1)[2] - p
    
    lasso = glmnet::cv.glmnet(x = X1,y = Y1,penalty.factor = c(rep(1,p),rep(0,r)))
    beta1 = glmnet::glmnet(x = X1,y = Y1,penalty.factor = c(rep(1,p),rep(0,r)),
                           lambda = lasso$lambda.min)
    B = rbind(B,as.numeric(beta1$beta)[1:p])
    
    pre = Xtest %*% beta1$beta[1:p] + as.numeric(beta1$beta)[dim(X1)[2]]
    
    Y_pre = c(Y_pre,as.numeric(pre))
  }
  result = list()
  result$pre = Y_pre
  result$B = B[-1,]
  result$nusiance = nusiance
  result$Var = Var
  result$CoefM = coefM/n
  return(result)
}

##################################################################################
Y = sqrt(mymeta$absol_eos)#mymeta$absol_eos
Xdata = cbind(clr.OTU,mymetabo)
pdim = dim(Xdata)[2]

MVFit = Leave_One_Out_Analysis(M1 = clr.OTU,M2 = mymetabo,Y = Y,MV = T)
FCFit = Leave_One_Out_Analysis(M1 = clr.OTU,M2 = mymetabo,Y = Y,Integfactor = T)
ORFit = Leave_One_Out_Analysis(M1 = clr.OTU,M2 = mymetabo,Y = Y)
FCFit1 = Leave_One_Out_Analysis(M1 = clr.OTU,M2 = mymetabo,Y = Y,Factor = T)

#######################################
AverageCoefM = MVFit$CoefM
colnames(AverageCoefM) = 1:200
pheatmap::pheatmap(AverageCoefM,cluster_rows = F,cluster_cols = F,color = colorRampPalette(c("blue", "white", "red"))(100))


BI1 = MVFit$B!=0
BI2 = FCFit$B!=0
BI3 = ORFit$B!=0
BI4 = FCFit1$B!=0

MV = colSums(BI1)[1:pdim]
FC = colSums(BI2)[1:pdim]
OR = colSums(BI3)[1:pdim]
FC1 = colSums(BI4)[1:pdim]

sum((Y-MVFit$pre)**2)
sum((Y-FCFit$pre)**2)
sum((Y-ORFit$pre)**2)
sum((Y-FCFit1$pre)**2)

MV = colSums(BI1)[1:pdim]
FC = colSums(BI2)[1:pdim]
OR = colSums(BI3)[1:pdim]
FC1 = colSums(BI4)[1:pdim]

par(mfrow = c(2,2))
plot(MV,main = 'Multi-view Regression',ylim=c(0,55),
     ylab='Selection Frequency',xlab='Feature ID')
plot(FC,main = 'Integrative Factor Regression',ylim=c(0,55),
     ylab='Selection Frequency',xlab='Feature ID')
plot(OR,main = 'Lasso Regression',ylim=c(0,55),
     ylab='Selection Frequency',xlab='Feature ID')
plot(FC1,main = 'Factor-Adjusted Regularized Regression',ylim=c(0,55),
     ylab='Selection Frequency',xlab='Feature ID')

sum(MV>=49)
sum(FC>=49)
sum(OR>=49)
sum(FC1>=49)

Meta1 = which(MV>=40)-31

MetaD = AverageCoefM[,Meta1[-1]]
colnames(MetaD) = colnames(mymetabo[,Meta1[-1]])
pheatmap::pheatmap(t(MetaD),cluster_rows = F,cluster_cols = F)

sum((Y-MVFit$pre)**2)/length(Y)
sum((Y-FCFit$pre)**2)/length(Y)
sum((Y-ORFit$pre)**2)/length(Y)
sum((Y-FCFit1$pre)**2)/length(Y)

B1 = reshape2::melt(MVFit$B[,1:pdim])
B2 = reshape2::melt(FCFit$B[,1:pdim])
B3 = reshape2::melt(ORFit$B[,1:pdim])
B4 = reshape2::melt(FCFit1$B[,1:pdim])

p1 = ggplot(B1,aes(x=Var2,y=value,group=Var2))+
  geom_boxplot()+ggtitle('Multi-view Regression')
p2 = ggplot(B2,aes(x=Var2,y=value,group=Var2))+
  geom_boxplot()+ggtitle('Integrative Factor regression')
p3 = ggplot(B3,aes(x=Var2,y=value,group=Var2))+
  geom_boxplot()+ggtitle('Lasso Regression')
p4 = ggplot(B4,aes(x=Var2,y=value,group=Var2))+
  geom_boxplot()+ggtitle('Factor-Adjusted Regularized Regression')

ggpubr::ggarrange(p1,p2,p3,p4,ncol = 2,nrow = 2)

##########################################################
B1 = reshape2::melt(MVFit$B[,which(abs(colSums(MVFit$B[,1:pdim]))>0)])
B2 = reshape2::melt(FCFit$B[,which(abs(colSums(FCFit$B[,1:pdim]))>0)])
B3 = reshape2::melt(ORFit$B[,which(abs(colSums(ORFit$B[,1:pdim]))>0)])
B4 = reshape2::melt(FCFit1$B[,which(abs(colSums(FCFit1$B[,1:pdim]))>0)])


B1 = reshape2::melt(MVFit$B[,1:pdim])
B2 = reshape2::melt(FCFit$B[,1:pdim])
B3 = reshape2::melt(ORFit$B[,1:pdim])
B4 = reshape2::melt(FCFit1$B[,1:pdim])

colnames(B1) = c('ID','Loc','Coef')
colnames(B2) = c('ID','Loc','Coef')
colnames(B3) = c('ID','Loc','Coef')
colnames(B4) = c('ID','Loc','Coef')

data1 <- data_summary(B1, varname="Coef", 
                      groupnames=c("Loc"))
data2 <- data_summary(B2, varname="Coef", 
                      groupnames=c("Loc"))
data3 <- data_summary(B3, varname="Coef", 
                      groupnames=c("Loc"))
data4 <- data_summary(B4, varname="Coef", 
                      groupnames=c("Loc"))

p1 = ggplot(data1, aes(x=Loc, y=as.numeric(Coef))) + 
  geom_errorbar(aes(ymin=Coef-sd, ymax=Coef+sd), width=5,
                position=position_dodge(0.5))+ggtitle('Multi-view Regression')+ 
  theme(legend.position="none")+ylab('Coefficient estimations')+xlab('Feature ID')
p2 = ggplot(data2, aes(x=Loc, y=as.numeric(Coef))) + 
  geom_errorbar(aes(ymin=Coef-sd, ymax=Coef+sd), width=5,
                position=position_dodge(0.5))+ggtitle('Integrative Factor Regression')+ 
  theme(legend.position="none")+ylab('Coefficient estimations')+xlab('Feature ID')
p3 = ggplot(data3, aes(x=Loc, y=as.numeric(Coef))) + 
  geom_errorbar(aes(ymin=Coef-sd, ymax=Coef+sd), width=5,
                position=position_dodge(0.5))+ggtitle('Lasso Regression')+ 
  theme(legend.position="none")+ylab('Coefficient estimations')+xlab('Feature ID')
p4 = ggplot(data4, aes(x=Loc, y=as.numeric(Coef))) + 
  geom_errorbar(aes(ymin=Coef-sd, ymax=Coef+sd), width=5,
                position=position_dodge(0.5))+ggtitle('Factor-Adjusted Regularized Regression')+ 
  theme(legend.position="none")+ylab('Coefficient estimations')+xlab('Feature ID')
ggpubr::ggarrange(p1,p2,p3,p4,ncol = 2,nrow = 2)

p1 = ggplot(B1,aes(x=Var2,y=value,group=Var2))+
  geom_boxplot()+ggtitle('Multi-view regression')
p2 = ggplot(B2,aes(x=Var2,y=value,group=Var2))+
  geom_boxplot()+ggtitle('Integrative factor regression')
p3 = ggplot(B3,aes(x=Var2,y=value,group=Var2))+
  geom_boxplot()+ggtitle('Adaptive lasso regression')
p4 = ggplot(B4,aes(x=Var2,y=value,group=Var2))+
  geom_boxplot()+ggtitle('Latent factor regression')

ggpubr::ggarrange(p1,p2,p3,p4,ncol = 2,nrow = 2)

##########################################################
MVb1 = Pairwise_index(MVFit$B[,1:pdim])
FCb1 = Pairwise_index(FCFit$B[,1:pdim])
FC1b1 = Pairwise_index(FCFit1$B[,1:pdim])
ORb1 = Pairwise_index(ORFit$B[,1:pdim])

i=1
JI = c(MVb1[,i],FCb1[,i],FC1b1[,i],ORb1[,i])
model = rep(c('Multi-view regression','Integrative factor regression',
              'Lasso regression','Latent factor regression'),each = 1485)
data1 = data.frame(model,JI)

p1 = ggplot(data1,aes(x=model,y=JI,group=model,col=model))+
  geom_violin(width=1.4) +
  geom_boxplot(width=0.1, color="grey", alpha=0.5) + 
  scale_fill_viridis(discrete = TRUE)+ggtitle('Jaccard similarity coefficient')+ylab('Values')

i=2
JI = c(MVb1[,i],FCb1[,i],FC1b1[,i],ORb1[,i])
model = rep(c('Multi-view regression','Integrative factor regression',
              'Lasso regression','Latent factor regression'),each = 1485)
data1 = data.frame(model,JI)
p2 = ggplot(data1,aes(x=model,y=JI,group=model,col=model))+
  geom_violin(width=1.4) +
  geom_boxplot(width=0.1, color="grey", alpha=0.5) + 
  scale_fill_viridis(discrete = TRUE)+ggtitle('Otsuka–Ochiai coefficient')+ylab('Values')

i=3
JI = c(MVb1[,i],FCb1[,i],FC1b1[,i],ORb1[,i])
model = rep(c('Multi-view regression','Integrative factor regression',
              'Lasso regression','Latent factor regression'),each = 1485)
data1 = data.frame(model,JI)
p3 = ggplot(data1,aes(x=model,y=JI,group=model,col=model))+
  geom_violin(width=1.4) +
  geom_boxplot(width=0.1, color="grey", alpha=0.5) + 
  scale_fill_viridis(discrete = TRUE)+ggtitle('Sørensen–Dice coefficient')+ylab('Values')

ggpubr::ggarrange(p1,p2,p3,common.legend = T,ncol = 1,nrow = 3)