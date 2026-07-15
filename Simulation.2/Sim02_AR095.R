rm(list = ls.str(mode = 'numeric'))

nlambda = 100
ratio1=0.9
set.seed(354)

SNRy = 5
SNRx = 5
Sigma = 0.95

r = 1
RN = 30

n = 200
p = 300
q = 300

SR = (p-10)/p

X = rep(1,n)
M1 = mvtnorm::rmvnorm(n, mean = c(rep(0,p)), sigma = diag(p))
w1 = create_sparsew1(p,q,r=r,ratio = 0.5,size = 2)

###############################################################################
sparse_indices <- sample(1:p, size = SR*p)
w11 <- runif(p,1,2)*sample(c(-1,1),p,replace = T)
w11[sparse_indices] <- 0
w11 = matrix(w11)

sparse_indices <- sample(1:q, size = SR*q)
w31 <- runif(q,1,2)*sample(c(-1,1),q,replace = T)
w31[sparse_indices] <- 0
w31 = matrix(w31)

M2 <- as.matrix(M1) %*% w1
Noise = mvtnorm::rmvnorm(n, mean = c(rep(0,q)), sigma = ar1_cor(q,Sigma))
ratio = (Frobenius_Norm(M2)/Frobenius_Norm(Noise))
noiseM2 = Noise * ratio/1
M2 <- as.matrix(M1) %*% w1 + noiseM2

Y <- as.matrix(M1) %*% w11 + M2 %*% w31

###############################################################################
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

SNR1 = c()
i = 1
AUC_Summary = rep(0,18)

for(SNRy in SNRy){
  ## Form M2
  lambda = sqrt(M2_FB**2 / (M1w1_FB**2 * ((SNRx + a1)**2 - a2 + a3)))

  a = M1w1_FB**2
  b = 2*lambda*sum(as.matrix(M1) %*% w1 * noiseM2)
  c = lambda**2 * Noise_FB**2 - M2_FB**2
  Sigma = b**2 - 4*a*c
  C = (-b + sqrt(b**2 - 4*a*c))/(2*a)

  M2 = C*as.matrix(M1) %*% w1 + lambda*noiseM2
  print(C/lambda)

  ## Multi-view regression. Rank selection
  RR = Rank_search(M2, M1)
  r = which(RR[[2]] == min(RR[[2]]))
  r2 = list(); r2$coef = as.matrix(RR[[1]][r][[1]])

  s3 = svd(M1%*%r2$coef)
  b3 = s3$u[,1:r] %*% t(s3$u[,1:r])

  ## Factor decomposition
  lf = Factor_decompose(list(cbind(M1,M2)))
  lf1 = Factor_decompose(list(M1,M2))
  xf = do.call(cbind,lf1$Fm); xu = do.call(cbind,lf1$U)

  ## Cooperative learning
  X_CL_1 = rbind(cbind(M1,M2),cbind(-M1,M2)) # rho = 1
  X_CL_05 = rbind(cbind(M1,M2),cbind(-sqrt(0.5)*M1,sqrt(0.5)*M2)) # rho = 0.5
  X_AL = rbind(cbind(M1,M2),cbind(0*M1,0*M2)) # rho = 0

  for(Rep in 1:100){
    Len = length(c(w11,w31))

    Y <- as.matrix(M1) %*% w11 + M2 %*% w31
    Noise = rnorm(n)
    ratio = (Frobenius_Norm(Y)/Frobenius_Norm(Noise))
    noise = Noise * ratio/SNRy
    Y <- as.matrix(M1) %*% w11 + M2 %*% w31 + noise

    Y_CL = rbind(Y,0*Y)

    Max_lambda = NULL#LassoFit$lambda[1]

    ## Adaptive part
    Adalasso = AdapLasso_Penalty_AUC_internal(X = X_AL, Y = Y_CL,w11 = w11,w31 = w31,Max_lambda = Max_lambda,
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
    Mv5 = AdapLasso_Penalty_AUC_internal(X = X_CL_05, Y = Y_CL,w11 = w11,w31 = w31,Max_lambda = Max_lambda,
                                         r = 0,nlambda = nlambda,gamma = 1,ratio=ratio1)
    Mv6 = AdapLasso_Penalty_AUC_internal(X = X_CL_1, Y = Y_CL,w11 = w11,w31 = w31,Max_lambda = Max_lambda,
                                         r = 0,nlambda = nlambda,gamma = 1,ratio=ratio1)

    i=1
    print(rbind(Adalasso[[i]],Mv3[[i]],Mv4[[i]],Mv5[[i]],Mv6[[i]],MV[[i]]))
    AUC_Summary = rbind(AUC_Summary,c(Adalasso[[i]],Mv3[[i]],Mv4[[i]],Mv5[[i]],Mv6[[i]],MV[[i]]))
  }
}

Total_AUC = AUC_Summary[-1,c(1,10,13,7,4,16)]
Beta1_AUC = AUC_Summary[-1,c(1,10,13,7,4,16)+1]
Beta2_AUC = AUC_Summary[-1,c(1,10,13,7,4,16)+2]

m1 = factor(rep(c('Cooperative(rho=0)','Cooperative(rho=0.5)','Cooperative(rho=1)','Factor','IntegFactor','Multi-view'),each=100),
            levels=c('Cooperative(rho=0)','Cooperative(rho=0.5)','Cooperative(rho=1)','Factor','IntegFactor','Multi-view'))

Total_AUC = reshape2::melt(Total_AUC); Total_AUC$model = m1
Beta1_AUC = reshape2::melt(Beta1_AUC); Beta1_AUC$model = m1
Beta2_AUC = reshape2::melt(Beta2_AUC); Beta2_AUC$model = m1

p1 = ggplot(Total_AUC,aes(x=Var2,y=value,group=model,col=model))+
  geom_violin()+ggtitle('Total AUC')+geom_boxplot()+geom_jitter(alpha=0.5)

p2 = ggplot(Beta1_AUC,aes(x=Var2,y=value,group=model,col=model))+
  geom_violin()+ggtitle('Beta1 AUC')+geom_boxplot()+geom_jitter(alpha=0.5)

p3 = ggplot(Beta2_AUC,aes(x=Var2,y=value,group=model,col=model))+
  geom_violin()+ggtitle('Beta2 AUC')+geom_boxplot()+geom_jitter(alpha=0.5)
ggpubr::ggarrange(p1,p2,p3,nrow = 1,ncol = 3,common.legend = T)

#####################################################################################
AUC_Summary_AR95 = AUC_Summary
