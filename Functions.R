Rank_search<-function(Y,X,r_range = 1:10){
  X = as.matrix(X); Y = as.matrix(Y)
  n = dim(X)[1]; p = dim(X)[2]; q = dim(Y)[2]
  Coef = list()
  GIC = c()

  for(r in r_range){
    r1 = rrpack::srrr(Y,X,
                      nrank = r,method = 'glasso',ic.type = 'GIC') #method = 'glasso',ic.type = 'GIC'

    d1 = data.frame(r1$ic.path,r1$lambda); d2 = d1[,5:6]
    GIC = c(GIC,d2[order(d2$X5),][1,1])
    print(d2[order(d2$X5),][1,1])
    Coef = append(Coef,list(as.matrix(r1$coef)))
  }
  re = list(Coef,GIC)
  return(re)
}

create_sparsew1<-function(p,q,r,ratio=0.5,size=0.5){
  sparse_indices <- sample(1:p, size = ratio*p)
  w1Bvector1 = c()
  for(i in 1:r){
    w1Bvector <- runif(p,-size,size)
    w1Bvector[sparse_indices] <- 0
    w1Bvector1 = c(w1Bvector1,w1Bvector)
  }
  w1B = matrix(w1Bvector1,nrow = p,ncol = r)
  w1A = matrix(runif(r*q,-1,1),nrow = r,ncol = q) #runif(r*q,-1,1)
  w1 = w1B%*%w1A

  #w1 = apply(w1, 2, function(col) col / sqrt(sum(col^2)))

  return(w1)
}

CS_cor<-function(n,rho){
  cov = matrix(rho,n,n) + diag(1-rho,n,n)
  return(cov)
}

ar1_cor <- function(n, rho) {
  exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) -
                    (1:n - 1))
  return(rho^exponent)
}

Unstr_cor <- function(n) {
  cor1 = matrix(runif(n*n,-1,1),n,n)
  cor1[lower.tri(cor1)] = t(cor1)[lower.tri(cor1)]
  diag(cor1) = 1

  return(cor1)
}

Toeplitz_cor <- function(n,size=1){
  first_row <- c(1,runif(n-1,-size,size))
  cov_matrix <- Matrix::toeplitz(first_row)

  return(cov_matrix)
}

Frobenius_Norm<-function(A){
  d<-sqrt(sum(diag(t(A)%*%A)))
  return(d)
}

data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

Classification_report<-function(fit,true){
  d = data.frame(fit = c(fit!=0),true = c(true!=0))

  fp<-dim(d[fit!=0&true==0,])[1]
  fn<-dim(d[fit==0&true!=0,])[1]
  tp<-dim(d[fit!=0&true!=0,])[1]
  tn<-dim(d[fit==0&true==0,])[1]

  return(c(fp,fn,tp,tn))
}

AdapLasso_Penalty<-function(X,Y,nlambda=50,r=0,gamma=1,family="gaussian"){
  p = dim(X)[2]
  Index = c(rep(1,p-r),rep(0,r))

  kridge = glmnet::cv.glmnet(x = X, y = Y,
                             alpha=0,intercept = FALSE,
                             standardize = FALSE,nlambda = nlambda,
                             penalty.factor = Index,family = family )
  kridge1 = glmnet::glmnet(x = X, y = Y,
                           lambda = kridge$lambda.min,alpha=0,
                           intercept = FALSE,standardize = FALSE,
                           penalty.factor = Index,family = family )
  best_ridge_coef <- as.numeric((kridge1$beta))^(gamma)

  lassoCV = glmnet::cv.glmnet(x = X, y = Y,alpha=1,
                              intercept = FALSE,standardize = FALSE,
                              nlambda = nlambda,family = family,
                              penalty.factor = c(1 / abs(best_ridge_coef[1:(p-r)]),rep(0,r)))
  LassoFit = glmnet::glmnet(x = X, y = Y,lambda = lassoCV$lambda.min,alpha=1,
                            intercept = FALSE,standardize = FALSE,family = family,
                            penalty.factor = c(1 / abs(best_ridge_coef[1:(p-r)]),rep(0,r)))
  return(LassoFit)
}

Lasso_Penalty<-function(X,Y,nlambda=50,r=1){
  p = dim(X)[2]
  Index = c(rep(1,p-r),rep(0,r))

  best_ridge_coef <- rep(1,p)

  lassoCV = glmnet::cv.glmnet(x = X, y = Y,alpha=1,
                              intercept = FALSE,standardize = FALSE,
                              nlambda = nlambda,
                              penalty.factor = c(1 / abs(best_ridge_coef[1:(p-r)]),rep(0,r)))
  LassoFit = glmnet::glmnet(x = X, y = Y,lambda = lassoCV$lambda.min,alpha=1,
                            intercept = FALSE,standardize = FALSE,
                            penalty.factor = c(1 / abs(best_ridge_coef[1:(p-r)]),rep(0,r)))
  return(LassoFit)
}

Lasso_Penalty_AUC_internal<-function(X,Y,w11,w31,Max_lambda=NULL,r=1,nlambda=50,ratio=0.7){
  p = dim(X)[2]
  Index = c(rep(1,p-r),rep(0,r))

  if(is.null(Max_lambda)){
    LassoFit = glmnet::glmnet(x = X, y = Y,alpha=1,nlambda = nlambda,
                              intercept = FALSE,standardize = FALSE,
                              penalty.factor = Index)

    Max_lambda = LassoFit$lambda[1]
  }

  lambda_range = c(Max_lambda*ratio^c(0,1:(nlambda-2)),1e-4)

  Sen = c(); Sen1 = c(); Sen2 = c()
  Spe = c(); Spe1 = c(); Spe2 = c()

  lambda_len = length(lambda_range)
  p = length(w11); q = length(w31)
  Beta = rep(0,(p+q))

  for(i in 1:lambda_len){
    LassoFit = glmnet::glmnet(x = X, y = Y,alpha=1,lambda = lambda_range[i],
                              intercept = FALSE,standardize = FALSE,
                              penalty.factor = Index)
    beta = as.numeric(LassoFit$beta[1:(p+q)])

    suppressWarnings({tab1 = caret::confusionMatrix(as.factor(c(beta==0)),
                                                    as.factor(c(w11,w31)==0),positive='FALSE')})
    suppressWarnings({tab2 = caret::confusionMatrix(as.factor(c(beta[1:p]==0)),
                                                    as.factor(c(w11)==0),positive='FALSE')})
    suppressWarnings({tab3 = caret::confusionMatrix(as.factor(c(beta[(p+1):(p+q)]==0)),
                                                    as.factor(c(w31)==0),positive='FALSE')})
    Sen = c(Sen,as.numeric(tab1$byClass[1]));Sen1 = c(Sen1,as.numeric(tab2$byClass[1]));Sen2 = c(Sen2,as.numeric(tab3$byClass[1]))
    Spe = c(Spe,as.numeric(tab1$byClass[2]));Spe1 = c(Spe1,as.numeric(tab2$byClass[2]));Spe2 = c(Spe2,as.numeric(tab3$byClass[2]))
    Beta = rbind(Beta,beta)
  }
  d1 = data.frame(FPR = 1-Spe,TPR = Sen)
  d2 = data.frame(FPR = 1-Spe1,TPR = Sen1)
  d3 = data.frame(FPR = 1-Spe2,TPR = Sen2)

  suppressWarnings({auc1 = tryCatch(as.numeric(DescTools::AUC(x = d1$FPR,y=d1$TPR,method = 'trapezoid')),
                                    error = function(e) return(0))})
  suppressWarnings({auc2 = tryCatch(as.numeric(DescTools::AUC(x = d2$FPR,y=d2$TPR,method = 'trapezoid')),
                                    error = function(e) return(0))})
  suppressWarnings({auc3 = tryCatch(as.numeric(DescTools::AUC(x = d3$FPR,y=d3$TPR,method = 'trapezoid')),
                                    error = function(e) return(0))})

  suppressWarnings({AUC_value = c(auc1,auc2,auc3)})

  results = list(AUC_value,d1,d2,d3,Beta)

  return(results)
}

AdapLasso_Penalty_AUC_internal<-function(X,Y,w11,w31,Max_lambda=NULL,
                                         r=1,gamma=1,nlambda=50,ratio=0.7,
                                         family = "gaussian"){
  p = dim(X)[2]
  Index = c(rep(1,p-r),rep(0,r))

  kridge = glmnet::cv.glmnet(x = X, y = Y,
                             alpha=0,intercept = FALSE,
                             standardize = FALSE,nlambda = nlambda,
                             penalty.factor = Index,family = family)
  kridge1 = glmnet::glmnet(x = X, y = Y,
                           lambda = kridge$lambda.min,alpha=0,
                           intercept = FALSE,standardize = FALSE,
                           penalty.factor = Index,family = family)
  best_ridge_coef <- as.numeric((kridge1$beta))^(gamma)

  if(is.null(Max_lambda)){
    LassoFit = glmnet::glmnet(x = X, y = Y,alpha=1,nlambda = nlambda,
                              intercept = FALSE,standardize = FALSE,family = family,
                              penalty.factor = c(1 / abs(best_ridge_coef[1:(p-r)]),rep(0,r)))

    Max_lambda = LassoFit$lambda[1]
  }

  lambda_range = c(Max_lambda*ratio^c(0,1:(nlambda-2)),1e-4)

  Sen = c(); Sen1 = c(); Sen2 = c()
  Spe = c(); Spe1 = c(); Spe2 = c()

  lambda_len = length(lambda_range)
  p = length(w11); q = length(w31)
  Beta = rep(0,(p+q))

  for(i in 1:lambda_len){
    LassoFit = glmnet::glmnet(x = X, y = Y,alpha=1,lambda = lambda_range[i],
                                         intercept = FALSE,standardize = FALSE,family = family,
                                         penalty.factor = c(1 / abs(best_ridge_coef[1:(p+q)]),rep(0,r)))
    beta = as.numeric(LassoFit$beta[1:(p+q)])

    suppressWarnings({tab1 = caret::confusionMatrix(as.factor(c(beta==0)),
                                                    as.factor(c(w11,w31)==0),positive='FALSE')})
    suppressWarnings({tab2 = caret::confusionMatrix(as.factor(c(beta[1:p]==0)),
                                                    as.factor(c(w11)==0),positive='FALSE')})
    suppressWarnings({tab3 = caret::confusionMatrix(as.factor(c(beta[(p+1):(p+q)]==0)),
                                                    as.factor(c(w31)==0),positive='FALSE')})
    Sen = c(Sen,as.numeric(tab1$byClass[1]));Sen1 = c(Sen1,as.numeric(tab2$byClass[1]));Sen2 = c(Sen2,as.numeric(tab3$byClass[1]))
    Spe = c(Spe,as.numeric(tab1$byClass[2]));Spe1 = c(Spe1,as.numeric(tab2$byClass[2]));Spe2 = c(Spe2,as.numeric(tab3$byClass[2]))
    Beta = rbind(Beta,beta)
  }
  d1 = data.frame(FPR = 1-Spe,TPR = Sen)
  d2 = data.frame(FPR = 1-Spe1,TPR = Sen1)
  d3 = data.frame(FPR = 1-Spe2,TPR = Sen2)

  suppressWarnings({auc1 = tryCatch(as.numeric(DescTools::AUC(x = d1$FPR,y=d1$TPR,method = 'trapezoid')),
                                    error = function(e) return(0))})
  suppressWarnings({auc2 = tryCatch(as.numeric(DescTools::AUC(x = d2$FPR,y=d2$TPR,method = 'trapezoid')),
                                    error = function(e) return(0))})
  suppressWarnings({auc3 = tryCatch(as.numeric(DescTools::AUC(x = d3$FPR,y=d3$TPR,method = 'trapezoid')),
                                    error = function(e) return(0))})
  suppressWarnings({AUC_value = c(auc1,auc2,auc3)})

  results = list(AUC_value,d1,d2,d3,Beta)

  return(results)
}

Create_Eigen_Design<-function(n,p,eigenvalue){
  M = mvtnorm::rmvnorm(n, mean = c(rep(0,p)), sigma = diag(p))

  Q <- qr.Q(qr(M))
  #A <- Q %*% diag(eigenvalue) %*% t(Q)
  A <- Q %*% sqrt(diag(eigenvalue))

  return(A)
}

Factor_decompose<-function(X,K=NULL,Km_U=NULL){
  ## Input X should be list. Each element represent one modality
  m = length(X)
  Fm_list = list(); Sigma_list = list(); U_list = list()

  for(i in 1:m){
    S = eigen(X[[i]] %*% t(X[[i]]))
    n = dim(X[[i]])[1]
    if(is.null(K)){
      if(is.null(Km_U)){
        Km_U = min(dim(X[[i]])[1],dim(X[[i]])[2])/2
      }

      loss = Find_Km(X[[i]],S,Km_U = Km_U)
      loss = loss[order(loss$loss),]; Km_est = loss[1,2]
    }else{
      Km_est = K
      loss = 0
    }

    if(Km_est == 0){
      U_est = X[[i]]
      Fm_est = NULL
      Sigma_est = NULL
    }else{
      Fm_est = sqrt(n)*as.matrix(S$vectors[,1:Km_est])
      Sigma_est = (1/n)*t(X[[i]]) %*% Fm_est
      U_est = X[[i]] - Fm_est%*%t(Sigma_est)
    }
    Km_U = NULL
    Fm_list[[i]] = Fm_est; U_list[[i]] = U_est; Sigma_list[[i]] = Sigma_est
  }
  result_list = list(Fm = Fm_list, U = U_list,Sigma = Sigma_list, loss = loss)
  return(result_list)
}

Find_Km<-function(X,S,Km_U = 10){
  EV = S$values; n = dim(X)[1]
  loss = c(Km_loss(X))

  for(j in 1:Km_U){
    Fm = sqrt(n)*as.matrix(S$vectors[,1:j])
    loss = c(loss,Km_loss(X,Fm))
  }
  return(data.frame(loss,Km = 0:Km_U))
}

Km_loss<-function(X,Fm=NULL){
  X = as.matrix(X)

  n = dim(X)[1]
  pm = dim(X)[2]

  if(is.null(Fm)){
    loss = log((1/(n*pm)) * Frobenius_Norm(X)**2)
  }else{
    Fm = as.matrix(Fm); k = dim(Fm)[2]
    loss = log((1/(n*pm)) * Frobenius_Norm(X - (Fm %*% t(Fm))%*%X/n)**2) + k*gPenalty(n,pm)
  }

  return(loss)
}

gPenalty<-function(n,pm){
  return(((n+pm)/(n*pm))*log(n*pm/(n+pm)))
}

gPenalty1<-function(n,pm){
  return(((n+pm)/(n*pm))*log(min(n,pm)))
}

Multi_view_Integrative_Factor<-function(M1,M2,Nusiance,Y){
  M_l = list(M1,M2)

  l = Factor_decompose(M_l)
  Fm = cbind(l$Fm[[1]],l$Fm[[2]])
  U = cbind(l$U[[1]],l$U[[2]])

  Nuisance = cbind(Fm,Nusiance); r = dim(Nuisance)[2]

  fit = AdapLasso_Penalty(cbind(U,Nuisance),Y,r = r,gamma = 1)
  return(fit)
}


Integrative_Factor_AUC_internal<-function(M1,M2,Nusiance,Y,w11,w31,gamma=1,nlambda=50,ratio=0.7){
  #M_l = list(M1,M2)

  #l = Facor_decompose(M_l)
  #Fm = cbind(l$Fm[[1]],l$Fm[[2]])
  #U = cbind(l$U[[1]],l$U[[2]])

  M_l = list(M1)
  l = Factor_decompose(M_l,Km_U = round(dim(M1)[2]/10,0))
  Fm = l$Fm[[1]]
  U = l$U[[1]]

  Nuisance = cbind(Fm,Nusiance); r = dim(Nuisance)[2]

  X = cbind(U,M2,Nuisance)
  p = dim(X)[2]
  Index = c(rep(1,p-r),rep(0,r))

  kridge = glmnet::cv.glmnet(x = X, y = Y,
                             alpha=0,intercept = FALSE,
                             standardize = FALSE,nlambda = nlambda,
                             penalty.factor = Index)
  kridge1 = glmnet::glmnet(x = X, y = Y,
                           lambda = kridge$lambda.min,alpha=0,
                           intercept = FALSE,standardize = FALSE,
                           penalty.factor = Index)
  best_ridge_coef <- as.numeric((kridge1$beta))^(gamma)

  LassoFit = glmnet::glmnet(x = X, y = Y,alpha=1,nlambda = nlambda,
                            intercept = FALSE,standardize = FALSE,
                            penalty.factor = c(1 / abs(best_ridge_coef[1:(p-r)]),rep(0,r)))

  Max_lambda = LassoFit$lambda[1]; lambda_range = c(Max_lambda*ratio^c(0,1:(nlambda-2)),1e-4)

  Sen = c(); Sen1 = c(); Sen2 = c()
  Spe = c(); Spe1 = c(); Spe2 = c()

  lambda_len = length(lambda_range)
  p = length(w11); q = length(w31)

  for(i in 1:lambda_len){
    fit = AdapLasso_Penalty_AUC(X,Y,r = r,gamma = 1,lambda = lambda_range[i])
    beta = as.numeric(fit$beta[1:(p+q)])

    suppressWarnings({tab1 = caret::confusionMatrix(as.factor(c(beta==0)),
                                                    as.factor(c(w11,w31)==0),positive='FALSE')})
    suppressWarnings({tab2 = caret::confusionMatrix(as.factor(c(beta[1:p]==0)),
                                                    as.factor(c(w11)==0),positive='FALSE')})
    suppressWarnings({tab3 = caret::confusionMatrix(as.factor(c(beta[(p+1):(p+q)]==0)),
                                                    as.factor(c(w31)==0),positive='FALSE')})
    Sen = c(Sen,as.numeric(tab1$byClass[1]));Sen1 = c(Sen1,as.numeric(tab2$byClass[1]));Sen2 = c(Sen2,as.numeric(tab3$byClass[1]))
    Spe = c(Spe,as.numeric(tab1$byClass[2]));Spe1 = c(Spe1,as.numeric(tab2$byClass[2]));Spe2 = c(Spe2,as.numeric(tab3$byClass[2]))
  }
  d1 = data.frame(FPR = 1-Spe,TPR = Sen)
  d2 = data.frame(FPR = 1-Spe1,TPR = Sen1)
  d3 = data.frame(FPR = 1-Spe2,TPR = Sen2)

  suppressWarnings({auc1 = tryCatch(as.numeric(DescTools::AUC(x = d1$FPR,y=d1$TPR,method = 'trapezoid')),
                                    error = function(e) return(0))})
  suppressWarnings({auc2 = tryCatch(as.numeric(DescTools::AUC(x = d2$FPR,y=d2$TPR,method = 'trapezoid')),
                                    error = function(e) return(0))})
  suppressWarnings({auc3 = tryCatch(as.numeric(DescTools::AUC(x = d3$FPR,y=d3$TPR,method = 'trapezoid')),
                                    error = function(e) return(0))})

  suppressWarnings({AUC_value = c(auc1,auc2,auc3)})
  results = list(AUC_value,d1,d2,d3)
  return(results)
}

# k3 = Multi_view_Integrative_Factor((diag(1,dim(b2)[1]) - b2)%*%M1,M2 - M1%*%r2$coef,as.matrix(s2$u[,1:r]),Y = Y)
Integ_Factor_Lasso_AUC_internal<-function(M1,M2,Y,w11,w31,Max_lambda=NULL,Nuisance=NULL,Weights=F,gamma=1,nlambda=50,ratio=0.7){
  M_l = list(M1,M2)

  l = Factor_decompose(M_l)
  Fm = if(length(l$Fm) == 0) NULL else if(length(l$Fm) == 1) cbind(l$Fm[[1]]) else cbind(l$Fm[[1]],l$Fm[[2]])
  Fm = cbind(Fm,Nuisance)
  U = cbind(l$U[[1]],l$U[[2]])

  X = cbind(U,Fm)
  p = dim(X)[2];

  r = if(is.null(Fm)) 0 else dim(as.matrix(Fm))[[2]]

  Index = c(rep(1,p-r),rep(0,r))

  if (Weights) {
    Weights = rep(1,(p-r))
  }else {
    kridge = glmnet::cv.glmnet(x = X, y = Y,
                               alpha=0,intercept = FALSE,
                               standardize = FALSE,nlambda = nlambda,
                               penalty.factor = Index)
    kridge1 = glmnet::glmnet(x = X, y = Y,
                             lambda = kridge$lambda.min,alpha=0,
                             intercept = FALSE,standardize = FALSE,
                             penalty.factor = Index)
    best_ridge_coef <- as.numeric((kridge1$beta))^(gamma)
    Weights = 1/abs(best_ridge_coef[1:(p-r)])
  }

  if(is.null(Max_lambda)){
    LassoFit = glmnet::glmnet(x = X, y = Y,alpha=1,nlambda = nlambda,
                              intercept = FALSE,standardize = FALSE,
                              penalty.factor = c(Weights,rep(0,r)))

    Max_lambda = LassoFit$lambda[1]
  }

  lambda_range = c(Max_lambda*ratio^c(0,1:(nlambda-2)),1e-4)

  Sen = c(); Sen1 = c(); Sen2 = c()
  Spe = c(); Spe1 = c(); Spe2 = c()

  lambda_len = length(lambda_range)
  p = length(w11); q = length(w31)
  Beta = rep(0,(p+q))

  for(i in 1:lambda_len){
    fit = glmnet::glmnet(x = X, y = Y,alpha=1,lambda = lambda_range[i],
                                    intercept = FALSE,standardize = FALSE,
                                    penalty.factor = c(Weights,rep(0,r)))
    beta = as.numeric(fit$beta[1:(p+q)])

    suppressWarnings({tab1 = caret::confusionMatrix(as.factor(c(beta==0)),
                                                    as.factor(c(w11,w31)==0),positive='FALSE')})
    suppressWarnings({tab2 = caret::confusionMatrix(as.factor(c(beta[1:p]==0)),
                                                    as.factor(c(w11)==0),positive='FALSE')})
    suppressWarnings({tab3 = caret::confusionMatrix(as.factor(c(beta[(p+1):(p+q)]==0)),
                                                    as.factor(c(w31)==0),positive='FALSE')})
    Sen = c(Sen,as.numeric(tab1$byClass[1]));Sen1 = c(Sen1,as.numeric(tab2$byClass[1]));Sen2 = c(Sen2,as.numeric(tab3$byClass[1]))
    Spe = c(Spe,as.numeric(tab1$byClass[2]));Spe1 = c(Spe1,as.numeric(tab2$byClass[2]));Spe2 = c(Spe2,as.numeric(tab3$byClass[2]))
    Beta = rbind(Beta,beta)
  }
  d1 = data.frame(FPR = 1-Spe,TPR = Sen)
  d2 = data.frame(FPR = 1-Spe1,TPR = Sen1)
  d3 = data.frame(FPR = 1-Spe2,TPR = Sen2)

  suppressWarnings({auc1 = tryCatch(as.numeric(DescTools::AUC(x = d1$FPR,y=d1$TPR,method = 'trapezoid')),
                                    error = function(e) return(0))})
  suppressWarnings({auc2 = tryCatch(as.numeric(DescTools::AUC(x = d2$FPR,y=d2$TPR,method = 'trapezoid')),
                                    error = function(e) return(0))})
  suppressWarnings({auc3 = tryCatch(as.numeric(DescTools::AUC(x = d3$FPR,y=d3$TPR,method = 'trapezoid')),
                                    error = function(e) return(0))})

  suppressWarnings({AUC_value = c(auc1,auc2,auc3)})

  results = list(AUC_value,d1,d2,d3,Beta)
  return(results)
}

Factor_Lasso_AUC_internal<-function(M1,M2,Y,w11,w31,Max_lambda=NULL,Nuisance=NULL,Weights=NULL,gamma=1,nlambda=50,ratio=0.7){
  M_l = list(M1)

  l = Factor_decompose(M_l)
  Fm = if(length(l$Fm) == 0) NULL else cbind(l$Fm[[1]])
  U = cbind(l$U[[1]]);Fm = cbind(Fm,Nuisance)

  X = cbind(U,M2,Fm);
  p = dim(X)[2]; r = if(is.null(Fm)) 0 else dim(as.matrix(Fm))[[2]]
  Index = c(rep(1,p-r),rep(0,r))

  if (is.null(Weights)) {
    Weights = rep(1,(p-r))
  }else {
    kridge = glmnet::cv.glmnet(x = X, y = Y,
                               alpha=0,intercept = FALSE,
                               standardize = FALSE,nlambda = nlambda,
                               penalty.factor = Index)
    kridge1 = glmnet::glmnet(x = X, y = Y,
                             lambda = kridge$lambda.min,alpha=0,
                             intercept = FALSE,standardize = FALSE,
                             penalty.factor = Index)
    best_ridge_coef <- as.numeric((kridge1$beta))^(gamma)
    Weights = 1/abs(best_ridge_coef[1:(p-r)])
  }

  if(is.null(Max_lambda)){
    LassoFit = glmnet::glmnet(x = X, y = Y,alpha=1,nlambda = nlambda,
                              intercept = FALSE,standardize = FALSE,
                              penalty.factor = c(Weights,rep(0,r)))

    Max_lambda = LassoFit$lambda[1]
  }

  lambda_range = c(Max_lambda*ratio^c(0,1:(nlambda-2)),1e-4)

  Sen = c(); Sen1 = c(); Sen2 = c()
  Spe = c(); Spe1 = c(); Spe2 = c()

  lambda_len = length(lambda_range)
  p = length(w11); q = length(w31)
  Beta = rep(0,(p+q))

  for(i in 1:lambda_len){
    fit = glmnet::glmnet(x = X, y = Y,alpha=1,lambda = lambda_range[i],
                         intercept = FALSE,standardize = FALSE,
                         penalty.factor = c(Weights,rep(0,r)))
    beta = as.numeric(fit$beta[1:(p+q)])

    suppressWarnings({tab1 = caret::confusionMatrix(as.factor(c(beta==0)),
                                                    as.factor(c(w11,w31)==0),positive='FALSE')})
    suppressWarnings({tab2 = caret::confusionMatrix(as.factor(c(beta[1:p]==0)),
                                                    as.factor(c(w11)==0),positive='FALSE')})
    suppressWarnings({tab3 = caret::confusionMatrix(as.factor(c(beta[(p+1):(p+q)]==0)),
                                                    as.factor(c(w31)==0),positive='FALSE')})
    Sen = c(Sen,as.numeric(tab1$byClass[1]));Sen1 = c(Sen1,as.numeric(tab2$byClass[1]));Sen2 = c(Sen2,as.numeric(tab3$byClass[1]))
    Spe = c(Spe,as.numeric(tab1$byClass[2]));Spe1 = c(Spe1,as.numeric(tab2$byClass[2]));Spe2 = c(Spe2,as.numeric(tab3$byClass[2]))
    Beta = rbind(Beta,beta)
  }
  d1 = data.frame(FPR = 1-Spe,TPR = Sen)
  d2 = data.frame(FPR = 1-Spe1,TPR = Sen1)
  d3 = data.frame(FPR = 1-Spe2,TPR = Sen2)

  suppressWarnings({auc1 = tryCatch(as.numeric(DescTools::AUC(x = d1$FPR,y=d1$TPR,method = 'trapezoid')),
                                    error = function(e) return(0))})
  suppressWarnings({auc2 = tryCatch(as.numeric(DescTools::AUC(x = d2$FPR,y=d2$TPR,method = 'trapezoid')),
                                    error = function(e) return(0))})
  suppressWarnings({auc3 = tryCatch(as.numeric(DescTools::AUC(x = d3$FPR,y=d3$TPR,method = 'trapezoid')),
                                    error = function(e) return(0))})

  suppressWarnings({AUC_value = c(auc1,auc2,auc3)})

  results = list(AUC_value,d1,d2,d3,Beta)
  return(results)
}


Plot_path<-function(beta,w11,w31){
  W = c(w11,w31); p1=length(w11)
  n = dim(beta)[1]; p = dim(beta)[2]

  d = reshape2::melt(beta); d$lambda = 1:n
  d$W = as.factor(rep(W!=0,each=n))
  d$part = c(rep('w11',p1*n),rep('w31',(p-p1)*n))

  p1 = ggplot(d,aes(x=lambda,y=value,group=Var2, color=W ))+
    geom_line()
  p2 = ggplot(d[d$part=='w11',],aes(x=lambda,y=value,group=Var2, color=W  ))+
    geom_line()
  p3 = ggplot(d[d$part=='w31',],aes(x=lambda,y=value,group=Var2, color=W  ))+
    geom_line()

  P = ggpubr::ggarrange(p1,p2,p3,nrow = 3,ncol = 1)

  return(P)
}

spline_ROC<-function(X,nlambda=50,newX = seq(from = 0,to=1,by=0.01)){
  l = dim(X)[1];
  nD <- data.frame(FPR = newX)
  fSpline = c(0,0)

  for(j in 1:(l/nlambda)){
    data = X[(1+(j-1)*nlambda):(nlambda*j),]
    spline_model <- lm(TPR ~ bs(FPR, df = 5),data)
    fTPR = predict(spline_model,newdata = nD)
    df = data.frame(FPR = newX,TPR = fTPR)
    fSpline = rbind(fSpline,df)
  }
  fSpline = fSpline[-1,]
  return(fSpline)
}

CI_Spline<-function(X,nlambda=50,newX = seq(from = 0,to=1,by=0.01)){
  S1 = spline_ROC(X,nlambda,newX)
  low = c(); up = c(); Mu = c()

  for(j in newX){
    S2 = S1[S1$FPR==j,]
    sd = sd(S2$TPR); mu = mean(S2$TPR)

    low = c(low,mu-1.96*sd); up = c(up,mu+1.96*sd)
    Mu = c(Mu,mu)
  }
  d1 = data.frame(newX,Mu,low,up)
  return(d1)
}

MSE<-function(b1,b2){
  return(sum((b1-b2)**2)/length(b1))
}

AdapLasso_Penalty_AUC_Multi_Modalities<-function(X,Y,w11,w21,w31,Max_lambda=NULL,Adap_weight=F,r=1,gamma=1,nlambda=50,ratio=0.7){
  p = dim(X)[2]
  Index = c(rep(1,p-r),rep(0,r))

  if(Adap_weight == T){
    kridge = glmnet::cv.glmnet(x = X, y = Y,
                               alpha=0,intercept = FALSE,
                               standardize = FALSE,nlambda = nlambda,
                               penalty.factor = Index)
    kridge1 = glmnet::glmnet(x = X, y = Y,
                             lambda = kridge$lambda.min,alpha=0,
                             intercept = FALSE,standardize = FALSE,
                             penalty.factor = Index)
    best_ridge_coef <- as.numeric((kridge1$beta))^(gamma)
  }else{
    best_ridge_coef = rep(1,p-r)
  }

  if(is.null(Max_lambda)){
    LassoFit = glmnet::glmnet(x = X, y = Y,alpha=1,nlambda = nlambda,
                              intercept = FALSE,standardize = FALSE,
                              penalty.factor = c(1 / abs(best_ridge_coef[1:(p-r)]),rep(0,r)))

    Max_lambda = LassoFit$lambda[1]
  }

  lambda_range = c(Max_lambda*ratio^c(0,1:(nlambda-2)),1e-4)

  Sen = c(); Sen1 = c(); Sen2 = c(); Sen3 = c()
  Spe = c(); Spe1 = c(); Spe2 = c(); Spe3 = c()

  lambda_len = length(lambda_range)
  p1 = length(w11); p2 = length(w21); p3 = length(w31)
  Pp = (p1 + p2 + p3)
  Beta = rep(0,(p1 + p2 + p3))

  for(i in 1:lambda_len){
    LassoFit = glmnet::glmnet(x = X, y = Y,alpha=1,lambda = lambda_range[i],
                              intercept = FALSE,standardize = FALSE,
                              penalty.factor = c(1 / abs(best_ridge_coef[1:(Pp)]),rep(0,r)))
    beta = as.numeric(LassoFit$beta[1:Pp])

    suppressWarnings({tab1 = caret::confusionMatrix(as.factor(c(beta==0)),
                                                    as.factor(c(w11,w21,w31)==0),positive='FALSE')})
    suppressWarnings({tab2 = caret::confusionMatrix(as.factor(c(beta[1:p1]==0)),
                                                    as.factor(c(w11)==0),positive='FALSE')})
    suppressWarnings({tab3 = caret::confusionMatrix(as.factor(c(beta[(p1+1):(p1+p2)]==0)),
                                                    as.factor(c(w21)==0),positive='FALSE')})
    suppressWarnings({tab4 = caret::confusionMatrix(as.factor(c(beta[(p1+p2+1):Pp]==0)),
                                                    as.factor(c(w31)==0),positive='FALSE')})
    Sen = c(Sen,as.numeric(tab1$byClass[1]));Sen1 = c(Sen1,as.numeric(tab2$byClass[1]));Sen2 = c(Sen2,as.numeric(tab3$byClass[1]));Sen3 = c(Sen3,as.numeric(tab4$byClass[1]))
    Spe = c(Spe,as.numeric(tab1$byClass[2]));Spe1 = c(Spe1,as.numeric(tab2$byClass[2]));Spe2 = c(Spe2,as.numeric(tab3$byClass[2]));Spe3 = c(Spe3,as.numeric(tab4$byClass[2]))
    Beta = rbind(Beta,beta)
  }
  d1 = data.frame(FPR = 1-Spe,TPR = Sen)
  d2 = data.frame(FPR = 1-Spe1,TPR = Sen1)
  d3 = data.frame(FPR = 1-Spe2,TPR = Sen2)
  d4 = data.frame(FPR = 1-Spe3,TPR = Sen3)

  suppressWarnings({auc1 = tryCatch(as.numeric(DescTools::AUC(x = d1$FPR,y=d1$TPR,method = 'trapezoid')),
                                    error = function(e) return(0))})
  suppressWarnings({auc2 = tryCatch(as.numeric(DescTools::AUC(x = d2$FPR,y=d2$TPR,method = 'trapezoid')),
                                    error = function(e) return(0))})
  suppressWarnings({auc3 = tryCatch(as.numeric(DescTools::AUC(x = d3$FPR,y=d3$TPR,method = 'trapezoid')),
                                    error = function(e) return(0))})
  suppressWarnings({auc4 = tryCatch(as.numeric(DescTools::AUC(x = d4$FPR,y=d4$TPR,method = 'trapezoid')),
                                    error = function(e) return(0))})

  suppressWarnings({AUC_value = c(auc1,auc2,auc3,auc4)})

  results = list(AUC_value,d1,d2,d3,d4,Beta)

  return(results)
}
