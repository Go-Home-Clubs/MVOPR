rm(list = ls.str(mode = 'numeric'))

nlambda = 100
ratio1=0.9
set.seed(354)

SNR_range = seq(from=3,to=10,by=1)

r_range = 1:10
r = 3

n = 200
p = 300
q = 300

f1 = 7
f2 = 7

SR = (p-10)/p

## Factors
U1 = mvtnorm::rmvnorm(n, mean = c(rep(0,p)), sigma = 4*diag(p))
U2 = mvtnorm::rmvnorm(n, mean = c(rep(0,q)), sigma = 4*diag(q))

F1 = mvtnorm::rmvnorm(n, mean = c(rep(0,f1)), sigma = diag(f1))
F2 = mvtnorm::rmvnorm(n, mean = c(rep(0,f2)), sigma = diag(f2))
Fs = mvtnorm::rmvnorm(n, mean = c(rep(0,r)), sigma = diag(r))

Lambda1 = mvtnorm::rmvnorm(f1, mean = c(rep(0,p)), sigma = diag(p))
Lambda1s = mvtnorm::rmvnorm(r, mean = c(rep(0,p)), sigma = 81*diag(p))
Lambda2 = mvtnorm::rmvnorm(f2, mean = c(rep(0,q)), sigma = diag(q))
Lambda2s = mvtnorm::rmvnorm(r, mean = c(rep(0,q)), sigma = 81*diag(q))

M1 = F1 %*% Lambda1 + Fs %*% Lambda1s + U1
M2 = F2 %*% Lambda2 + Fs %*% Lambda2s + U2

###############################################################################
sparse_indices <- sample(1:p, size = SR*p)
w11 <- runif(p,1,2)*sample(c(-1,1),p,replace = T)
w11[sparse_indices] <- 0
w11 = matrix(w11) 

sparse_indices <- sample(1:q, size = SR*q)
w31 <- runif(q,1,2)*sample(c(-1,1),q,replace = T)
w31[sparse_indices] <- 0
w31 = matrix(w31)

Y <- as.matrix(M1) %*% w11 + M2 %*% w31

##############################################################################
print(Frobenius_Norm(as.matrix(M1) %*% w11)/Frobenius_Norm(M2 %*% w31))
###############################################################################
RR <- with(list(M1,M2), rrr(M2, M1, maxrank = 10,ic.type = "GIC", penaltySVD = "ann"))
r = RR$rank
r2 = list(); r2$coef = as.matrix(RR$coef)

s3 = svd(M1%*%r2$coef)
b3 = s3$u[,1:r] %*% t(s3$u[,1:r])

lf = Factor_decompose(list(cbind(M1,M2)))

lf1 = Factor_decompose(list(M1,M2))
xf = do.call(cbind,lf1$Fm); xu = do.call(cbind,lf1$U)

X_CL_1 = rbind(cbind(M1,M2),cbind(-M1,M2)) # rho = 1
X_CL_05 = rbind(cbind(M1,M2),cbind(-sqrt(0.5)*M1,sqrt(0.5)*M2)) # rho = 0.5
#X_AL = rbind(cbind(M1,M2),cbind(0*M1,0*M2)) # rho = 0

for(Sigma in c(0.7)){
  ##
  SNR1 = c()
  i = 1
  AUC_Summary = rep(0,18)
  
  for(SNR in SNR_range){
    for(Rep in 1:100){
      Len = length(c(w11,w31))
      
      Y <- as.matrix(M1) %*% w11 + M2 %*% w31
      Noise = rnorm(n)
      ratio = (Frobenius_Norm(Y)/Frobenius_Norm(Noise))
      noise = Noise * ratio/SNR
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