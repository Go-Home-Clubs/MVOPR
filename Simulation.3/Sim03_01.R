set.seed(354)
n = 200
p1 = 100
nlambda = 100
ratio1=0.9

r=1

SNR_range = seq(from = 3,to = 7,by=0.5)

M1 = mvtnorm::rmvnorm(n, mean = c(rep(0,p1)), sigma = diag(p1))

B1 = create_sparsew1(p1,p1,r=r,ratio = 0.5,size = 2)
B2 = create_sparsew1(p1,p1,r=r,ratio = 0.7,size = 1)
B3 = create_sparsew1(p1,p1,r=r,ratio = 0.7,size = 1)

##############################################################################
sparse_indices <- sample(1:p1, size = 0.9*p1)
w11 <- runif(p1,1,2)*sample(c(-1,1),p1,replace = T)
w11[sparse_indices] <- 0
w11 = matrix(w11) 

sparse_indices <- sample(1:p1, size = 0.9*p1)
w21 <- runif(p1,1,2)*sample(c(-1,1),p1,replace = T)
w21[sparse_indices] <- 0
w21 = matrix(w21)

sparse_indices <- sample(1:p1, size = 0.9*p1)
w31 <- runif(p1,1,2)*sample(c(-1,1),p1,replace = T)
w31[sparse_indices] <- 0
w31 = matrix(w31)
###########################################################################

SNR1 = c()

i = 1

AUC1 = c(0,0,0,0)
AUC2 = c(0,0,0,0)
AUC3 = c(0,0,0,0)
AUC4 = c(0,0,0,0)
AUC5 = c(0,0,0,0)
AUC6 = c(0,0,0,0)
Type = c(0,0,0)

M2 = M1 %*% B1 
Noise = mvtnorm::rmvnorm(n, mean = c(rep(0,p1)), sigma = ar1_cor(p1,0.9))
ratio = (Frobenius_Norm(M2)/Frobenius_Norm(Noise))
noiseM2 = Noise * ratio/3
M2 = M1 %*% B1 + noiseM2

M3 = M2 %*% B2 + M1 %*% B3
Noise = mvtnorm::rmvnorm(n, mean = c(rep(0,p1)), sigma = ar1_cor(p1,0.9))
ratio = (Frobenius_Norm(M3)/Frobenius_Norm(Noise))
noiseM2 = Noise * ratio/5
M3 = M2 %*% B2 + M1 %*% B3 + noiseM2

Adap_weight = T

r1 = Rank_search(M2,M1,r_range = 1:5)
r2 = Rank_search(M3,cbind(M2,M1),r_range = 1:5)

rank1 = which(r1[[2]] == min(r1[[2]]))
r1$coef = as.matrix(r1[[1]][rank1][[1]])

rank2 = which(r2[[2]] == min(r2[[2]]))
r2$coef = as.matrix(r2[[1]][rank2][[1]])
#rr = Rank_search(D1,S1,r_range = 1:10)
#rr1 = Rank_search(D2,cbind(D1,S1),r_range = 1:10)
MB = M1 %*% (cbind(r1$coef,r2$coef[(p1+1):(2*p1),]))
s = svd(MB)
P = s$u[,1:qr(MB)$rank] %*% t(s$u[,1:qr(MB)$rank])
M2_res = M2 - M1%*%r1$coef
dnus2 = M2_res %*% r2$coef[1:p1,]; Snus2 = svd(dnus2); nus2 = Snus2$u[,1:qr(M2_res %*% r2$coef[1:p1,])$rank]
P2 = nus2 %*% t(nus2)

Nus = cbind(nus2,s$u[,1:qr(MB)$rank])
Res = cbind((diag(dim(M1)[1]) - P) %*% M1,
            (diag(dim(M2_res)[1]) - P2) %*% M2_res,
            M3 - cbind(M2,M1)%*%r2$coef)
######################################################################
## Factor decomposition
s2 = Factor_decompose(list(M1,M2,M3))
U = do.call(cbind,s2$U)
Fm = do.call(cbind,s2$Fm)

s1 = Factor_decompose(list(cbind(M1,M2,M3)))

## 
X_CL_1 = rbind(cbind(M1,M2,M3),cbind(-M1,M2,0*M3),
               cbind(-M1,0*M2,M3),cbind(0*M1,-M2,M3)) # rho = 1
X_CL_05 = rbind(cbind(M1,M2,M3),cbind(-sqrt(0.5)*M1,sqrt(0.5)*M2,0*M3),
                cbind(-sqrt(0.5)*M1,0*M2,sqrt(0.5)*M3),cbind(0*M1,-sqrt(0.5)*M2,sqrt(0.5)*M3)) # rho = 0.5
for(SNR in SNR_range){
  
  for(Rep in 1:100){
    Len = length(c(w11,w31))
    
    Y <- as.matrix(M1) %*% w11 + M2 %*% w21 + M3 %*% w31
    Noise = rnorm(n)
    ratio = (Frobenius_Norm(Y)/Frobenius_Norm(Noise))
    noise = Noise * ratio/SNR
    Y <- as.matrix(M1) %*% w11 + M2 %*% w21 + M3 %*% w31 + noise
    
    Y_CL = rbind(Y,0*Y,0*Y,0*Y)
    
    ##
    Max_lambda = NULL
    
    ## Adaptive part
    Adalasso = AdapLasso_Penalty_AUC_Multi_Modalities(X = cbind(M1,M2,M3), 
                                                      Y = Y,w11 = w11,w21 = w21,w31 = w31,Max_lambda = Max_lambda,
                                                      r = 0,nlambda = nlambda,gamma = 1,ratio=ratio1,Adap_weight = Adap_weight)
    
    Adalasso1 = AdapLasso_Penalty_AUC_Multi_Modalities(X = cbind(s1$U[[1]],s1$Fm[[1]]), 
                                                       Y = Y,w11 = w11,w21 = w21,w31 = w31,Max_lambda = Max_lambda,
                                                       r = dim(s1$Fm[[1]])[2],nlambda = nlambda,gamma = 1,ratio=ratio1,Adap_weight = Adap_weight)
    
    Mv4 = AdapLasso_Penalty_AUC_Multi_Modalities(X = cbind(U,Fm), 
                                                 Y = Y,w11 = w11,w21 = w21,w31 = w31,Max_lambda = Max_lambda,
                                                 r = dim(Fm)[2],nlambda = nlambda,gamma = 1,ratio=ratio1,Adap_weight = Adap_weight)
    
    CL1 = AdapLasso_Penalty_AUC_Multi_Modalities(X = X_CL_1, 
                                                 Y = Y_CL,w11 = w11,w21 = w21,w31 = w31,Max_lambda = Max_lambda,
                                                 r = 0,nlambda = nlambda,gamma = 1,ratio=ratio1,Adap_weight = Adap_weight)
    CL05 = AdapLasso_Penalty_AUC_Multi_Modalities(X = X_CL_05, 
                                                 Y = Y_CL,w11 = w11,w21 = w21,w31 = w31,Max_lambda = Max_lambda,
                                                 r = 0,nlambda = nlambda,gamma = 1,ratio=ratio1,Adap_weight = Adap_weight)
    
    MV = AdapLasso_Penalty_AUC_Multi_Modalities(X = cbind(Res,Nus), 
                                                Y = Y,w11 = w11,w21 = w21,w31 = w31,Max_lambda = Max_lambda,
                                                r = dim(Nus)[2],nlambda = nlambda,gamma = 1,ratio=ratio1,Adap_weight = Adap_weight)
    ###
    AUC1 = rbind(AUC1,Adalasso[[i]])
    AUC2 = rbind(AUC2,CL1[[i]])
    AUC3 = rbind(AUC3,CL05[[i]])
    AUC4 = rbind(AUC4,Adalasso1[[i]])
    AUC5 = rbind(AUC5,as.numeric(Mv4[[i]]))
    AUC6 = rbind(AUC6,as.numeric(MV[[i]]))
    
    print(rbind(Adalasso[[i]],CL1[[i]],CL05[[i]],Adalasso1[[i]],Mv4[[i]],MV[[i]]))
    
    #print(ggpubr::ggarrange(p1,p2,p3,p4,p5,p6,
    #                        p7,p8,p9,p10,p11,p12,ncol = 3,nrow = 4))
    
    #print(Plot_path(Beta_lasso,c(w11,w31)))
    
    print('------')
  }
}

i = 1
AUC_total = data.frame(AUC1[-1,i],AUC2[-1,i],AUC3[-1,i],AUC4[-1,i],AUC5[-1,i],AUC6[-1,i])

colnames(AUC_total) = c('AdapLasso','Cooperative(rho=1)','Cooperative(rho=0.5)','Factor','IntegFactor','Multi-view')

d1 = reshape2::melt(AUC_total)
d1$snr = rep(SNR_range,each=100)
colnames(d1) = c('Model','Value','SNR')

data1 <- data_summary(d1, varname="Value", 
                      groupnames=c("Model", "SNR"))

i = 2
AUC_total = data.frame(AUC1[-1,i],AUC2[-1,i],AUC3[-1,i],AUC4[-1,i],AUC5[-1,i],AUC6[-1,i])

colnames(AUC_total) = c('AdapLasso','Cooperative(rho=1)','Cooperative(rho=0.5)','Factor','IntegFactor','Multi-view')

d1 = reshape2::melt(AUC_total)
d1$snr = rep(SNR_range,each=100)
colnames(d1) = c('Model','Value','SNR')

data2 <- data_summary(d1, varname="Value", 
                      groupnames=c("Model", "SNR"))

i = 3
AUC_total = data.frame(AUC1[-1,i],AUC2[-1,i],AUC3[-1,i],AUC4[-1,i],AUC5[-1,i],AUC6[-1,i])

colnames(AUC_total) = c('AdapLasso','Cooperative(rho=1)','Cooperative(rho=0.5)','Factor','IntegFactor','Multi-view')

d1 = reshape2::melt(AUC_total)
d1$snr = rep(SNR_range,each=100)
colnames(d1) = c('Model','Value','SNR')

data3 <- data_summary(d1, varname="Value", 
                      groupnames=c("Model", "SNR"))

i = 4
AUC_total = data.frame(AUC1[-1,i],AUC2[-1,i],AUC3[-1,i],AUC4[-1,i],AUC5[-1,i],AUC6[-1,i])

colnames(AUC_total) = c('AdapLasso','Cooperative(rho=1)','Cooperative(rho=0.5)','Factor','IntegFactor','Multi-view')

d1 = reshape2::melt(AUC_total)
d1$snr = rep(SNR_range,each=100)
colnames(d1) = c('Model','Value','SNR')

data4 <- data_summary(d1, varname="Value", 
                      groupnames=c("Model", "SNR"))


p1 = ggplot(data1, aes(x=SNR, y=Value, group=Model, color=Model)) + 
  geom_line() +
  geom_point()+
  geom_errorbar(aes(ymin=Value-sd, ymax=Value+sd), width=0.3,
                position=position_dodge(0.1))+ggtitle('Total AUC')+ theme(legend.position="none")
p2 = ggplot(data2, aes(x=SNR, y=Value, group=Model, color=Model)) + 
  geom_line() +
  geom_point()+
  geom_errorbar(aes(ymin=Value-sd, ymax=Value+sd), width=0.3,
                position=position_dodge(0.1))+ggtitle('AUC (M1)')+ theme(legend.position="none")
p3 = ggplot(data3, aes(x=SNR, y=Value, group=Model, color=Model)) + 
  geom_line() +
  geom_point()+
  geom_errorbar(aes(ymin=Value-sd, ymax=Value+sd), width=0.3,
                position=position_dodge(0.1))+ggtitle('AUC (M2)')+ theme(legend.position="none")
p4 = ggplot(data4, aes(x=SNR, y=Value, group=Model, color=Model)) + 
  geom_line() +
  geom_point()+
  geom_errorbar(aes(ymin=Value-sd, ymax=Value+sd), width=0.3,
                position=position_dodge(0.1))+ggtitle('AUC (M3)')

ggpubr::ggarrange(p1,p2,p3,p4,ncol = 2,nrow = 2,common.legend = T)