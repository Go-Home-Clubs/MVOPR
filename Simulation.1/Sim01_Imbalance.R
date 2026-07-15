rm(list = ls.str(mode = 'numeric'))

nlambda = 100
ratio1=0.9
set.seed(354)

SNR_range = seq(from=20,to=30,by=1)

r_range = 1:10
r = 9

n = 200
p = 50
q = 300

SR1 = (p-10)/p
SR2 = (q-10)/q

X = rep(1,n)
M1 = mvtnorm::rmvnorm(n, mean = c(rep(0,p)), sigma = diag(p))
w1 = create_sparsew1(p,q,r=r,ratio = 0.7,size = 1)

###############################################################################
sparse_indices <- sample(1:p, size = SR1*p)
w11 <- runif(p,1,2)*sample(c(-1,1),p,replace = T)
w11[sparse_indices] <- 0
w11 = matrix(w11) 

sparse_indices <- sample(1:q, size = SR2*q)
w31 <- runif(q,1,2)*sample(c(-1,1),q,replace = T)
w31[sparse_indices] <- 0
w31 = matrix(w31)

M2 <- as.matrix(M1) %*% w1 
Noise = mvtnorm::rmvnorm(n, mean = c(rep(0,q)), sigma = diag(q))
ratio = (Frobenius_Norm(M2)/Frobenius_Norm(Noise))
noiseM2 = Noise * ratio/1
M2 <- as.matrix(M1) %*% w1 + noiseM2

Y <- as.matrix(M1) %*% w11 + M2 %*% w31

##
M1w1_FB = Frobenius_Norm(as.matrix(M1) %*% w1)
Noise_FB = Frobenius_Norm(noiseM2)
M2_FB = Frobenius_Norm(M2)

lambda_limit = sqrt(-(M1w1_FB*M2_FB)**2/(sum(as.matrix(M1) %*% w1 * noiseM2)**2 - (M1w1_FB*Noise_FB)**2))

lambda = 1

a = M1w1_FB**2
b = 2*lambda*sum(as.matrix(M1) %*% w1 * noiseM2)
c = lambda**2 * Noise_FB**2 - M2_FB**2
Sigma = b**2 - 4*a*c

C = (-b + sqrt(b**2 - 4*a*c))/(2*a)

a1 = sum(as.matrix(M1) %*% w1 * noiseM2)/M1w1_FB**2
a2 = sum(as.matrix(M1) %*% w1 * noiseM2)**2/M1w1_FB**4
a3 = Noise_FB**2/M1w1_FB**2
##############################################################################
print(Frobenius_Norm(as.matrix(M1) %*% w11)/Frobenius_Norm(M2 %*% w31))
###############################################################################

for(Sigma in c(0.7)){
  ##
  SNR1 = c()
  i = 1
  AUC_Summary = rep(0,18)
  
  for(SNR in SNR_range){
    ## Form M2
    lambda = sqrt(M2_FB**2 / (M1w1_FB**2 * ((SNR + a1)**2 - a2 + a3)))
    
    a = M1w1_FB**2
    b = 2*lambda*sum(as.matrix(M1) %*% w1 * noiseM2)
    c = lambda**2 * Noise_FB**2 - M2_FB**2
    Sigma = b**2 - 4*a*c
    C = (-b + sqrt(b**2 - 4*a*c))/(2*a)
    
    M2 = C*as.matrix(M1) %*% w1 + lambda*noiseM2
    print(C/lambda)
    
    ## Multi-view regression. Find the rank
    #RR = with(list(M1,M2), rrr(M2, M1, maxrank = 30,ic.type = "GIC" ))#Rank_search(M2, M1,r_range = r_range)
    RR = Rank_search(M2, M1,r_range = 1:10)
    r = which(RR[[2]] == min(RR[[2]]))
    r2 = list(); r2$coef = as.matrix(RR[[1]][r][[1]])
    #r = RR$rank
    #r2 = list(); r2$coef = as.matrix(RR$coef)
    
    s3 = svd(M1%*%r2$coef)
    b3 = s3$u[,1:r] %*% t(s3$u[,1:r])
    
    ## Factor decomposition
    lf = Factor_decompose(list(cbind(M1,M2)))
    
    ## Integrative factor decomposition
    lf1 = Factor_decompose(list(M1,M2))
    xf = do.call(cbind,lf1$Fm); xu = do.call(cbind,lf1$U)
    
    ## Cooperative learning
    X_CL_1 = rbind(cbind(M1,M2),cbind(-M1,M2)) # rho = 1
    X_CL_05 = rbind(cbind(M1,M2),cbind(-sqrt(0.5)*M1,sqrt(0.5)*M2)) # rho = 0.5
    #X_AL = rbind(cbind(M1,M2),cbind(0*M1,0*M2)) # rho = 0
    
    for(Rep in 1:100){
      Len = length(c(w11,w31))
      
      Y <- as.matrix(M1) %*% w11 + M2 %*% w31
      Noise = rnorm(n)
      ratio = (Frobenius_Norm(Y)/Frobenius_Norm(Noise))
      noise = Noise * ratio/100
      Y <- as.matrix(M1) %*% w11 + M2 %*% w31 + noise
      
      Y_CL = rbind(Y,0*Y)
      
      Max_lambda = NULL#LassoFit$lambda[1]
      
      ## Adaptive part
      Adalasso = AdapLasso_Penalty_AUC_internal(X = cbind(M1,M2), Y = Y,w11 = w11,w31 = w31,Max_lambda = Max_lambda,
                                                r = 0,nlambda = nlambda,gamma = 1,ratio=ratio1)
      MV = AdapLasso_Penalty_AUC_internal(X = cbind((diag(1,dim(b3)[1]) - b3)%*%M1,M2 - M1%*%r2$coef,as.matrix(s3$u[,1:r])), 
                                          Y = Y,w11 = w11,w31 = w31,Max_lambda = Max_lambda,
                                          r = r,nlambda = nlambda,gamma = 1,ratio=ratio1)
      Mv3 = AdapLasso_Penalty_AUC_internal(X = cbind(lf$U[[1]],lf$Fm[[1]]), 
                                           Y = Y,w11 = w11,w31 = w31,Max_lambda = Max_lambda,
                                           r = dim(lf$Fm[[1]])[2],nlambda = nlambda,gamma = 1,ratio=ratio1)
      Mv4 = AdapLasso_Penalty_AUC_internal(X = cbind(xu,xf), 
                                           Y = Y,w11 = w11,w31 = w31,Max_lambda = Max_lambda,
                                           r = dim(xf)[2],nlambda = nlambda,gamma = 1,ratio=ratio1)
      Mv5 = AdapLasso_Penalty_AUC_internal(X = X_CL_1, Y = Y_CL,w11 = w11,w31 = w31,Max_lambda = Max_lambda,
                                           r = 0,nlambda = nlambda,gamma = 1,ratio=ratio1)
      Mv6 = AdapLasso_Penalty_AUC_internal(X = X_CL_05, Y = Y_CL,w11 = w11,w31 = w31,Max_lambda = Max_lambda,
                                           r = 0,nlambda = nlambda,gamma = 1,ratio=ratio1)
      
      i=1
      print(rbind(Adalasso[[i]],Mv3[[i]],Mv4[[i]],Mv5[[i]],Mv6[[i]],MV[[i]]))
      AUC_Summary = rbind(AUC_Summary,c(Adalasso[[i]],Mv3[[i]],Mv4[[i]],Mv5[[i]],Mv6[[i]],MV[[i]]))
    }
  }
}


Total_AUC = AUC_Summary[-1,c(1,4,7,10,13,16)]
Beta1_AUC = AUC_Summary[-1,c(1,4,7,10,13,16)+1]
Beta2_AUC = AUC_Summary[-1,c(1,4,7,10,13,16)+2]

#m1 = factor(rep(c('AdapLasso','Factor_AdapLasso','IntegFactor_AdapLasso','Cooperative(rho=1)','Cooperative(rho=0.5)','MV_AdapLasso'),each=800), 
#            levels=c('AdapLasso','Factor_AdapLasso','IntegFactor_AdapLasso','Cooperative(rho=1)','Cooperative(rho=0.5)','MV_AdapLasso'))
Total_AUC = AUC_Summary[-1,c(1,10,13,7,4,16)]
Beta1_AUC = AUC_Summary[-1,c(1,10,13,7,4,16)+1]
Beta2_AUC = AUC_Summary[-1,c(1,10,13,7,4,16)+2]

m1 = factor(rep(c('Cooperative(rho=0)','Cooperative(rho=0.5)','Cooperative(rho=1)','IntegFactor','Factor','Multi-view'),each=600), 
            levels=c('Cooperative(rho=0)','Cooperative(rho=0.5)','Cooperative(rho=1)','IntegFactor','Factor','Multi-view'))

Total_AUC = reshape2::melt(Total_AUC); Total_AUC$model = m1
Beta1_AUC = reshape2::melt(Beta1_AUC); Beta1_AUC$model = m1
Beta2_AUC = reshape2::melt(Beta2_AUC); Beta2_AUC$model = m1

Total_AUC$SNR = rep(SNR_range,each=100)
Beta1_AUC$SNR = rep(SNR_range,each=100)
Beta2_AUC$SNR = rep(SNR_range,each=100)

p1 = ggplot(Total_AUC,aes(x=SNR,y=value,fill=model,col=model))+
  geom_violin()+ggtitle('Total AUC')+geom_boxplot()+geom_jitter(alpha=0.5)

p2 = ggplot(Beta1_AUC,aes(x=SNR,y=value,fill=model,col=model))+
  geom_violin()+ggtitle('Beta1 AUC')+geom_boxplot()+geom_jitter(alpha=0.5)

p3 = ggplot(Beta2_AUC,aes(x=SNR,y=value,fill=model,col=model))+
  geom_violin()+ggtitle('Beta2 AUC')+geom_boxplot()+geom_jitter(alpha=0.5)
ggpubr::ggarrange(p1,p2,p3,nrow = 1,ncol = 3,common.legend = T)

########################################################################################
p1 = ggplot(Total_AUC,aes(x=as.factor(SNR),y=value,fill=model))+
  geom_boxplot()+ggtitle('Total AUC')+xlab('SNR')
p2 = ggplot(Beta1_AUC,aes(x=as.factor(SNR),y=value,fill=model))+
  geom_boxplot()+ggtitle('Beta1 AUC')+xlab('SNR')
p3 = ggplot(Beta2_AUC,aes(x=as.factor(SNR),y=value,fill=model))+
  geom_boxplot()+ggtitle('Beta2 AUC')+xlab('SNR')
ggpubr::ggarrange(p1,p2,p3,nrow = 1,ncol = 3,common.legend = T)

###############################################
Total_AUC = AUC_Summary[-1,c(1,10,13,7,4,16)]
Beta1_AUC = AUC_Summary[-1,c(1,10,13,7,4,16)+1]
Beta2_AUC = AUC_Summary[-1,c(1,10,13,7,4,16)+2]

Model_Names = c('Adaptive Lasso','Cooperative(rho=0.5)','Cooperative(rho=1)','IntegFactor','Factor','MVOPR')

i = 1
AUC_total = Total_AUC

colnames(AUC_total) = Model_Names

d1 = reshape2::melt(AUC_total)
d1$snr = rep(SNR_range,each=100)
colnames(d1) = c('Null','Model','Value','SNR')

data1 <- data_summary(d1, varname="Value", 
                      groupnames=c("Model", "SNR"))

i = 2
AUC_total = Beta1_AUC

colnames(AUC_total) = Model_Names

d1 = reshape2::melt(AUC_total)
d1$snr = rep(SNR_range,each=100)
colnames(d1) = c('Null','Model','Value','SNR')

data2 <- data_summary(d1, varname="Value", 
                      groupnames=c("Model", "SNR"))

i = 3
AUC_total = Beta2_AUC

colnames(AUC_total) = Model_Names

d1 = reshape2::melt(AUC_total)
d1$snr = rep(SNR_range,each=100)
colnames(d1) = c('Null','Model','Value','SNR')

data3 <- data_summary(d1, varname="Value", 
                      groupnames=c("Model", "SNR"))

p1 = ggplot(data1, aes(x=SNR, y=Value, group=Model, color=Model)) + 
  geom_line() +
  geom_point()+
  geom_errorbar(aes(ymin=Value-sd, ymax=Value+sd), width=0.5,
                position=position_dodge(0.1))+ggtitle('Total AUC')+ theme(legend.position="none")
p2 = ggplot(data2, aes(x=SNR, y=Value, group=Model, color=Model)) + 
  geom_line() +
  geom_point()+
  geom_errorbar(aes(ymin=Value-sd, ymax=Value+sd), width=0.5,
                position=position_dodge(0.1))+ggtitle('AUC (M1)')+ theme(legend.position="none")
p3 = ggplot(data3, aes(x=SNR, y=Value, group=Model, color=Model)) + 
  geom_line() +
  geom_point()+
  geom_errorbar(aes(ymin=Value-sd, ymax=Value+sd), width=0.5,
                position=position_dodge(0.1))+ggtitle('AUC (M2)')+ theme(legend.position="none")
ggpubr::ggarrange(p1,p2,p3,ncol = 3,nrow = 1,common.legend = T)

###############################################
Total_AUC = AUC_Summary[-1,c(1,10,13,7,4,16)]
Beta1_AUC = AUC_Summary[-1,c(1,10,13,7,4,16)+1]
Beta2_AUC = AUC_Summary[-1,c(1,10,13,7,4,16)+2]