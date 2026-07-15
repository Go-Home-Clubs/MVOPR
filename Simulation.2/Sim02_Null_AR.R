rm(list = ls.str(mode = 'numeric'))

nlambda = 100
ratio1=0.9
set.seed(354)
Balance = 1
Km = 2

r = 1
RN = 30

n = 200
p = 100
q = 100

X = rep(1,n)

M1 = mvtnorm::rmvnorm(n, mean = c(rep(0,p)), sigma = ar1_cor(p,0.7))
M2 = mvtnorm::rmvnorm(n, mean = c(rep(0,q)), sigma = ar1_cor(p,0.7))

###############################################################################
sparse_indices <- sample(1:p, size = 0.9*p)
w11 <- runif(p,1,2)*sample(c(-1,1),p,replace = T)
w11[sparse_indices] <- 0
w11 = matrix(w11) 

sparse_indices <- sample(1:q, size = 0.9*q)
w31 <- runif(q,1,2)*sample(c(-1,1),q,replace = T)
w31[sparse_indices] <- 0
w31 = matrix(w31)

Y <- as.matrix(M1) %*% w11 + M2 %*% w31

##############################################################################
print(Frobenius_Norm(as.matrix(M1) %*% w11)/Frobenius_Norm(M2 %*% w31))
###############################################################################

A1T = c(0,0,0)
A2T = c(0,0,0)
A3T = c(0,0,0)

Mv11T = c(0,0,0)
Mv12T = c(0,0,0)
Mv13T = c(0,0,0)

Mv21T = c(0,0,0)
Mv22T = c(0,0,0)
Mv23T = c(0,0,0)

Mv31T = c(0,0,0)
Mv32T = c(0,0,0)
Mv33T = c(0,0,0)

for(Sigma in c(0.9)){
  ##
  i = 1
  AUC_Summary = rep(0,18)
  
  for(SNR in 1){
    #r2 = rrpack::cv.srrr(M2, M1,nrank = r)
    RR = Rank_search(M2, M1,r_range = 1:5)
    r = which(RR[[2]] == min(RR[[2]]))
    r2 = list(); r2$coef = as.matrix(RR[[1]][r][[1]])
    
    s3 = svd(M1%*%r2$coef)
    b3 = s3$u[,1:r] %*% t(s3$u[,1:r])
    
    lf = Factor_decompose(list(cbind(M1,M2)))
    fn1 = if(length(lf$Fm)==0) 0 else dim(lf$Fm[[1]])[2]
    XF1 = if(length(lf$Fm)==0) lf$U[[1]] else cbind(lf$U[[1]],lf$Fm[[1]])
    
    lf1 = Factor_decompose(list(M1,M2))
    xf = do.call(cbind,lf1$Fm)
    xu = do.call(cbind,lf1$U)
    fn2 = if(length(lf1$Fm)==0) 0 else dim(xf)[2]
    XF2 = if(length(lf1$Fm)==0) xu else cbind(xu,xf)
    
    ## Cooperative learning
    X_CL_1 = rbind(cbind(M1,M2),cbind(-M1,M2)) # rho = 1
    X_CL_05 = rbind(cbind(M1,M2),cbind(-sqrt(0.5)*M1,sqrt(0.5)*M2)) # rho = 0.5
    X_AL = rbind(cbind(M1,M2),cbind(0*M1,0*M2)) 
    
    for(Rep in 1:100){
      Len = length(c(w11,w31))
      
      Y <- as.matrix(M1) %*% w11 + M2 %*% w31
      Noise = rnorm(n)
      ratio = (Frobenius_Norm(Y)/Frobenius_Norm(Noise))
      noise = Noise * ratio/1
      Y <- as.matrix(M1) %*% w11 + M2 %*% w31 + noise
      
      Max_lambda = NULL#LassoFit$lambda[1]
      Y_CL = rbind(Y,0*Y)
      
      ## Adaptive part
      Adalasso = AdapLasso_Penalty_AUC_internal(X = X_AL, Y = Y_CL,w11 = w11,w31 = w31,Max_lambda = Max_lambda,
                                                r = 0,nlambda = nlambda,gamma = 1,ratio=ratio1)
      MV = AdapLasso_Penalty_AUC_internal(X = cbind((diag(1,dim(b3)[1]) - b3)%*%M1,M2 - M1%*%r2$coef,as.matrix(s3$u[,1:r])), 
                                          Y = Y,w11 = w11,w31 = w31,Max_lambda = Max_lambda,
                                          r = r,nlambda = nlambda,gamma = 1,ratio=ratio1)
      Mv3 = AdapLasso_Penalty_AUC_internal(X = XF1, 
                                           Y = Y,w11 = w11,w31 = w31,Max_lambda = Max_lambda,
                                           r = fn1,nlambda = nlambda,gamma = 1,ratio=ratio1)
      Mv4 = AdapLasso_Penalty_AUC_internal(X = XF2, 
                                           Y = Y,w11 = w11,w31 = w31,Max_lambda = Max_lambda,
                                           r = fn2,nlambda = nlambda,gamma = 1,ratio=ratio1)
      Mv5 = AdapLasso_Penalty_AUC_internal(X = X_CL_05, Y = Y_CL,w11 = w11,w31 = w31,Max_lambda = Max_lambda,
                                           r = 0,nlambda = nlambda,gamma = 1,ratio=ratio1)
      Mv6 = AdapLasso_Penalty_AUC_internal(X = X_CL_1, Y = Y_CL,w11 = w11,w31 = w31,Max_lambda = Max_lambda,
                                           r = 0,nlambda = nlambda,gamma = 1,ratio=ratio1)
      
      i=1
      print(rbind(Adalasso[[i]],Mv3[[i]],Mv4[[i]],Mv5[[i]],Mv6[[i]],MV[[i]]))
      AUC_Summary = rbind(AUC_Summary,c(Adalasso[[i]],Mv3[[i]],Mv4[[i]],Mv5[[i]],Mv6[[i]],MV[[i]]))
    }
  }
}