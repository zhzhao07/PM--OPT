TOTAL_combo <- function(
    x,                # 数据框
    resp, trt,        # 响应、处理列名
    u, std.u,         # 协变量及其标准化列名
    model.function,   # e.g. hierarchical.models4
    nmodels,
    nTree = 1000, nodeSize = 10
) {
  # 1. 标准化协变量（如尚未）

  # 2. FACT matching 得到 grp
  Y  <- x[[resp]]
  W  <- x[[trt]]
  D  <- proximityScore(as.matrix(x[, std.u, drop=FALSE]), Y, W,
                       nTree = nTree, nodeSize = nodeSize)
  pm <- factMatch(D, 3, 3)$solution
  pm <- factMatchPrune(pm, D,
                       indexTreatment = which(W == 1),
                       indexControl   = which(W == 0),
                       Y = Y)
  grp <- matchToGroup(
    matchMatrix    = pm$matchMatrixPrune,
    indexControl   = which(W == 0),
    indexTreatment = which(W == 1),
    Y               = Y
  )
  
  # 3. 正确提取伪对索引
  idx_pairs <- which(!is.infinite(grp$tauMatrix), arr.ind = TRUE)
  t_idx <- mapply(function(r,c) grp$treatmentMatrix[r,c],
                  idx_pairs[,1], idx_pairs[,2])
  c_idx <- mapply(function(r,c) grp$controlMatrix[r,c],
                  idx_pairs[,1], idx_pairs[,2])
  tau_combo <- mapply(function(r,c) grp$tauMatrix[r,c],
                      idx_pairs[,1], idx_pairs[,2])
  
  # 4. 训练所有基础模型
  base_models <- model.function(x[,c(u,resp,trt)])
  
  # 5. 构造 residual matrix ehatt
  ehatt <- matrix(NA, nrow = length(tau_combo), ncol = nmodels)
  for (j in seq_len(nmodels)) {
    df_t <- data.frame(x[t_idx, u, drop = FALSE], Tt = 1)
    df_c <- data.frame(x[c_idx, u, drop = FALSE], Tt = 0)
    pred_t <- predict(base_models[[j]], newdata = df_t)
    pred_c <- predict(base_models[[j]], newdata = df_c)
    ehatt[, j] <- tau_combo - (pred_t - pred_c)
  }
  
  # 6. 求 OPT 权重 via solve.QP
  Dmat <- crossprod(ehatt) + diag(1e-6, nmodels)
  # estimate sigma^2 from first model residuals
  resid1 <- Y - predict(base_models[[1]], newdata = x)
  sigma2 <- sum(resid1^2) / (nrow(x) - length(coef(base_models[[nmodels]])))
  dvec  <- matrix(c(-sigma2*lct(x[,c(u,resp,trt)])),nmodels,1) #
  
  Amat <- cbind(
    rep(1, nmodels),    # sum(w)=1
    diag(nmodels),      # w >= 0
    -diag(nmodels)      # w <= 1
  )
  bvec <- c(1, rep(0, nmodels), rep(-1, nmodels))
  
  sol <- quadprog::solve.QP(Dmat, dvec, Amat, bvec, meq = 1)
  w_opt <- pmax(sol$solution, 0)
  w_opt <- w_opt / sum(w_opt)
  
  return(list(w_OPT = w_opt))
}




library(Biobase)
library(quadprog)
library(mgcv)
library(mvtnorm)

#rm(list=ls())

p=6#2
M=p;#变量加截距项个数
MM=2^p

hierarchical.models64 <- function(x) {
  # 确保数据为 data.frame 格式
  x <- as.data.frame(x)
  
  # 提取响应变量和处理变量
  Y <- x[, 1]
  Tt <- x[, 2]
  
  # 提取协变量部分并重命名为 U1, U2, ...
  original_covariates <- x[, 3:ncol(x), drop = FALSE]
  num_covs <- ncol(original_covariates)
  new_cov_names <- paste0("U", seq_len(num_covs))
  colnames(original_covariates) <- new_cov_names
  
  # 用新协变量名替代
  covariates <- new_cov_names
  
  # 辅助函数：生成所有协变量子集
  power_set <- function(set) {
    n <- length(set)
    subsets <- unlist(lapply(0:n, function(k) {
      combn(set, k, simplify = FALSE)
    }), recursive = FALSE)
    return(subsets)
  }
  
  models <- list()
  model_id <- 1
  main_subsets <- power_set(covariates)
  
  for (main_vars in main_subsets) {
    formula_terms <- c("Tt")  # 必含 Tt
    
    if (length(main_vars) > 0) {
      # 添加主效应
      formula_terms <- c(formula_terms, main_vars)
      # 添加交互项
      interaction_terms <- paste0("Tt:", main_vars)
      formula_terms <- c(formula_terms, interaction_terms)
    }
    
    # 构造公式
    model_formula <- as.formula(paste("Y ~", paste(formula_terms, collapse = " + ")))
    
    # 构造数据框供 lm 使用
    model_df <- data.frame(Y = Y, Tt = Tt, original_covariates)
    
    # 拟合模型
    models[[model_id]] <- lm(model_formula, data = model_df)
    model_id <- model_id + 1
  }
  
  return(models)
}

hierarchical.models4=function(x)
{
  Y=x[,1];Tt=x[,2];u1=x[,3];u2=x[,4];
  sol1=lm(formula = Y ~ Tt+u1+u2+Tt:u1+Tt:u2)
  sol2=lm(formula = Y ~ Tt+u1+Tt:u1)
  sol3=lm(formula = Y ~ Tt+u2+Tt:u2)
  sol4=lm(formula = Y ~ Tt)
  list(sol1,sol2,sol3,sol4)
}

lct=function(x)
{lc=1:nmodels
for(i in 1:nmodels){lc[i]=length(hierarchical.models64(x)[[i]]$coefficients)}
lc
}



varf=function(x){
  # xx=rmvnorm(x,rep(0,P),Sigma)
  # Tt=rbinom(x,1,0.5) #
  # #Tte =   1/(1 + exp(xx[,1] + xx[,2]))
  # #Tt=rbinom(x,1,Tte)
  # #e=rnorm(x,0,1);#e2=rnorm(x,0,1);e=Tt*e1+(1-Tt)*e2;
  # e1=rnorm(x,0,1);#rskewlap(x)
  # e2=rskewlap(x);e=Tt*e1+(1-Tt)*e2;
  # e=Tt*e1+(1-Tt)*e2;#rnorm(x,0,1)#
  # Y=Tt*(h1(gam1,xx))+(1-Tt)*(h2(gam2,xx))+e#
  
  xx = matrix(runif(P * x, -1, 1), ncol = P)
  prop =  exp(xx[,1])/(1 + exp(xx[,1]))
  Tt = rbinom(x, 1, prop)
  tau = delt(gam1,xx[,1:P])
  mu0 = h2(gam2,xx[,1:P])
  noise = rnorm(x, 0, 1)
  Y = Tt * tau + mu0 + noise
  r= seq(0.1,0.9,0.2)
  c=sqrt(c(r*var(e))/c(c(1-r)*c(var(Y-e))) )
  c
}

TOTAL <- function(
    x,                ###data frame on which to apply the data splitting
    resp,             ###name of the response variable in x 
    trt,              ###name of the treatment variable in x.  Should =1 if treated, 0 otherwise.
    u,                ###names of the covariates used in the model
    std.u,            ###names of the standardized covariates used to compute distances 
    nsplits=100,      ###number of data splittings over which to average
    dev.fraction=0.5, ###fraction of the data to be used for model development 
    model.function,   ###function should take a dataset x, construct the candidate procedures
    ###being evaluated on x, and return a list of models on which the 
    ###predict() function can be used
    nmodels           ###length of the model list that will be created by model.function
)
{
  ehatt=matrix(NA,ncol=nmodels,nrow=nrow(x)); Aic=Bic=rep(NA,nmodels)
  trn.dat.trt<-subset(x, get(trt)==1);  trn.dat.ctrl<-subset(x, get(trt)==0)
  nearest.neighbor<-matchpt(as.matrix(trn.dat.trt[,std.u]),as.matrix(trn.dat.ctrl[,std.u]))[,1]
  trn.dat.trt.paired<-trn.dat.trt[,c(resp,trt,u)]#被配对的样本   是处理元
  trn.dat.sim.ctrl<-trn.dat.ctrl[nearest.neighbor,c(resp,trt,u)]; #trn.dat.sim.ctrl[,trt]<-0#和处理组匹配的对照组
  nearest.neighbor1<-matchpt(as.matrix(trn.dat.ctrl[,std.u]), as.matrix(trn.dat.trt[,std.u]))[,1]
  trn.dat.ctrl.paired<-trn.dat.ctrl[,c(resp,trt,u)]  #被配对的样本   是对照元
  trn.dat.sim.trt<-trn.dat.trt[nearest.neighbor1,c(resp,trt,u)]; #trn.dat.sim.trt[,trt]<-1#和对照组匹配的处理组
  
  ###Append results 
  trn.dat.trta<-rbind(trn.dat.trt.paired, trn.dat.sim.trt)#处理组  处理元和配对处理元
  trn.dat.ctrla<-rbind(trn.dat.sim.ctrl,trn.dat.ctrl.paired)#对照组  配对对照元和对照元
  trn.nn.delta.tilde=trn.dat.trta[,resp]-trn.dat.ctrla[,resp]
  trn.true=rbind(trn.dat.trt.paired[,c(u,trt)],trn.dat.ctrl.paired[,c(u,trt)])
  trn.truet=trn.truec=trn.true; trn.truet$Tt=1; trn.truec$Tt=0
  lc=lct(x[,c(u,resp,trt)])#c(6,5,5,4,4,3,4,3,2)
  
  for (j in 1:nmodels) { 
    ehatt[,j]=trn.nn.delta.tilde-(predict(model.function(x[,c(resp,trt,u)])[[j]],trn.truet)-predict(model.function(x[,c(resp,trt,u)])[[j]], trn.truec))
    Aic[j]=(length(x$Y)*log(mean((x$Y-predict(model.function(x[,c(resp,trt,u)])[[j]],x))^2))+(lc[j])*2)
    Bic[j]=(length(x$Y)*log(mean((x$Y-predict(model.function(x[,c(resp,trt,u)])[[j]],x))^2))+(lc[j])*log(length(x$Y)))
  }
  aic=which.min(Aic); bic=which.min(Bic)
  wa=exp(0.5*min(Aic)-0.5*Aic)/sum(exp(0.5*min(Aic)-0.5*Aic)) #S-AIC权
  wb=exp(0.5*min(Bic)-0.5*Bic)/sum(exp(0.5*min(Bic)-0.5*Bic)) #S-BIC权重
  
  sigmatc=sum((x$Y-predict(model.function(x[,c(resp,trt,u)])[[nmodels]],x))^2)/(length(x$Y)-lc[nmodels])#summary(hierarchical.models4(x)[[1]])$sigma^2
  a1 <- t(ehatt[,1:(nmodels)]) %*% ehatt[,1:(nmodels)]+diag(0.000001,nmodels)
  a2<-matrix(c(-sigmatc*lc),nmodels,1) #
  a3 <- t(rbind(matrix(1,nrow=1,ncol=nmodels),diag(nmodels),-diag(nmodels)))
  a4 <- rbind(1,matrix(0,nrow=nmodels,ncol=1),matrix(-1,nrow=nmodels,ncol=1))
  QP <- solve.QP(a1,a2,a3,a4,1)
  mw <- QP$solution
  mw <- as.matrix(mw)
  mw <- mw*(mw>0)
  mw <- mw/sum(mw)  
  
  
  ntrain=nrow(x)
  cate.se.cv<-matrix(NA,nrow=nsplits,ncol=nmodels)#
  w=matrix(NA,ncol=nmodels,nrow=nsplits)
  
  for (i in 1:nsplits) {
    ###Create training and validation data
    val.inds<-sample(1:ntrain,round(ntrain*(1-dev.fraction)),replace=F)
    trn.data<-x[-val.inds,]   
    trn.data.trt<-subset(trn.data, get(trt)==1)　　
    trn.data.ctrl<-subset(trn.data, get(trt)==0)
    val.data<-x[val.inds,]  #should contain u, std.u, resp and trt    
    val.data.trt<-subset(val.data, get(trt)==1)
    val.data.ctrl<-subset(val.data, get(trt)==0)
    nval<-length(val.inds)
    
    ###Match each obs in the treatment group to its nearest control
    nearest.neighbors<-matchpt(as.matrix(val.data.trt[,std.u]), 
                               as.matrix(val.data.ctrl[,std.u]))
    nn.indices<-nearest.neighbors[,1]
    
    val.data.trt.paired1<-val.data.trt[,c(u,resp,trt)]
    val.nn.delta.tilde<-val.data.trt.paired1[,resp] - val.data.ctrl[nn.indices,resp]
    val.nn.yt.tilde=val.data.trt.paired1[,resp]
    val.nn.yc.tilde=val.data.ctrl[nn.indices,resp]
    val.data.trt.paire=val.data.ctrl[nn.indices,c(u,resp,trt)]
    val.data.sim.ctrl<-val.data.trt.paired1[,c(u,trt)]; val.data.sim.ctrl[,trt]<-0
    
    ###Now match each obs in the control group to its nearest treated obs
    nearest.neighbors<-matchpt(as.matrix(val.data.ctrl[,std.u]), 
                               as.matrix(val.data.trt[,std.u]))
    nn.indices<-nearest.neighbors[,1]
    
    val.data.ctrl.paired1<-val.data.ctrl[,c(u,resp,trt)]  
    val.nn.delta.tilde2<-val.data.trt[nn.indices,resp] - val.data.ctrl.paired1[,resp]
    val.nn.yt.tilde2=val.data.trt[nn.indices,resp]
    val.nn.yc.tilde2=val.data.ctrl.paired1[,resp]
    val.data.ctrl.paire=val.data.trt[nn.indices,c(u,resp,trt)]
    
    val.data.sim.trt<-val.data.ctrl.paired1[,c(u,trt)]; val.data.sim.trt[,trt]<-1
    
    ###Append results 
    val.nn.delta.tilde<-c(val.nn.delta.tilde, val.nn.delta.tilde2)
    val.nn.yta.tilde=c(val.nn.yt.tilde, val.nn.yt.tilde2)
    val.nn.yca.tilde=c(val.nn.yc.tilde, val.nn.yc.tilde2)
    
    val.data.trt.paired<-rbind(val.data.trt.paired1[,c(u,trt)], val.data.sim.trt)
    val.data.sim.ctrl<-rbind(val.data.sim.ctrl, val.data.ctrl.paired1[,c(u,trt)])
    
    ###Build models on development data
    models.train<-model.function(trn.data[,c(u,resp,trt)] )
    
    nn.sspe<-rep(NA,nmodels)      
    sigmaht=matrix(NA,nrow=length(val.inds),ncol=nmodels)
    
    for (j in 1:nmodels) { 
      sigmaht[,j]=sqrt((val.nn.yta.tilde-predict(models.train[[j]], rbind(val.data.trt.paired1,val.data.ctrl.paire)))^2+
                         (val.nn.yca.tilde-predict(models.train[[j]], rbind(val.data.trt.paire,val.data.ctrl.paired1)))^2)
      val.delta.hats<-predict(models.train[[j]], val.data.trt.paired)-predict(models.train[[j]], val.data.sim.ctrl)
      nn.sspe[j]<-sum((val.delta.hats - val.nn.delta.tilde)^2)
      w[i,j]=prod((10*dnorm((val.nn.delta.tilde-val.delta.hats)/(sigmaht[,j])))/(sigmaht[,j]))
    }
    w[i,]=w[i,]/sum(w[i,])
    cate.se.cv[i,]<-nn.sspe/length(val.nn.delta.tilde)
  }
  wn=colMeans(w)#
  TE.CV<-apply(cate.se.cv,2,mean)
  return(list(which.min(TE.CV),wn,mw,aic,bic,wa,wb))  
}


#rho=c(0,0.5,0.8)
nd=c(100,200,400)
for (tt in 1:3)
{
  set.seed(1)
  n=nd[tt];P=6;nmodels=64
  # Sigma=diag(P)
  # rrho=rho[tt];    for(u in 1:(P-1)){for(v in 1:(P-1)){Sigma[u,v] = rrho^abs(u-v)}}
  # Sigma[P,P]=1
   #gam1=gam2=c(0.5,1.5)
   gam1=gam2 = rep(0, P); #beta[(P+1)] = c(1); 
   gam1[1:5] = c(2,1.5,0,0,0); gam2[1:5] = -c(2,1.5,-1,-0.8,-0.5); 
  # h1<-function(ce,X){
  #   as.matrix(ce[1]*((X[,1])^2)+ce[2]*((X[,2]))+ce[1]*((X[,1]))+ce[2]*(X[,2])^2)
  # }
  # h2<-function(ce,X){
  #   as.matrix(ce[1]*((X[,1])^2)+ce[2]*((X[,2])) )
  # }
  # delt<-function(ce,X){
  #   as.matrix(ce[1]*((X[,1]))+ce[2]*(X[,2])^2)
  # }
  h1<-function(ce1,ce2,X){
    as.matrix( cbind(X) %*% ce1 + 2 * (X[,1])^2)+ as.matrix( cbind(X) %*% ce2)
  }
  h2<-function(ce,X){
    as.matrix( cbind(X) %*% ce + 2 * (X[,1])^2)
    #as.matrix(ce[1]*((X[,1])^2)+ce[2]*((X[,2])) ) mu0
  }
  delt<-function(ce,X){
    as.matrix( cbind(X) %*% ce)
    }
  
  K=100;N=1000000
  mse=matrix(0,K,nmodels)
  optmse=saicmse=sbicmse=teemmse=optNmse=tecvmse=aicmse=bicmse=rep(0,K);A=NULL
  c=varf(10000000)
  for (r in 1:length(c)){
    gamm1=c[r]*gam1
    gamm2=c[r]*gam2
    for(k in 1:K)#循环次数
    {  
      # xx=rmvnorm(n,rep(0,P),Sigma); 
      # #Tte =   1/(1 + exp(xx[,1] + xx[,2]))
      # Tt=rbinom(n,1,0.5)
      # e1=rnorm(n,0,1);#rskewlap(n)
      # e2=rskewlap(n)#rnorm(n,0,.5);#e2=rnorm(n,0,1);
      # e=Tt*e1+(1-Tt)*e2;#rnorm(n,0,1)#Tt*e1+(1-Tt)*e2;
      # Y=Tt*(h1(gamm1,xx[,1:2]))+(1-Tt)*(h2(gamm2,xx[,1:2]))+e
      
      xx = matrix(runif(P * n, -1, 1), ncol = P)
      prop =  exp(xx[,1])/(1 + exp(xx[,1]))
      Tt = rbinom(n, 1, prop)
      tau = delt(gamm1,xx[,1:P])
      mu0 = h2(gamm2,xx[,1:P])
      noise = rnorm(n, 0, 1)
      Y = Tt * tau + mu0 + noise
      
      xxy=cbind(Y,Tt,xx[,1:P]);  std.uu=scale(xxy[,3:(P+2)])
      X=as.data.frame(cbind(xxy,std.uu));  colnames(X)=c("Y","Tt",paste0("u", seq_len(P)),paste0("u", seq_len(P),".u"))
      wn=TOTAL(x=X, resp="Y", trt="Tt", u=paste0("u", seq_len(P)),
               std.u=paste0("u", seq_len(P),".u"), nsplits=100, 
               dev.fraction=0.5,
               model.function=function(x){hierarchical.models64(x)}, nmodels=nmodels)
      
      res <- TOTAL_combo(
        x              = X,
        resp="Y", trt="Tt", u=paste0("u", seq_len(P)),
        std.u=paste0("u", seq_len(P),".u"),
        model.function=function(x){hierarchical.models64(x)},
        nmodels        = 64
      )
      wnN=res$w_OPT
      xnew= matrix(runif(P * N, -1, 1), ncol = P)
        #rmvnorm(N,rep(0,P-1),Sigma[-P,-P])
      xnewt=cbind(1,xnew);  xnewc=cbind(0,xnew)
      colnames(xnewt)=colnames(xnewc)=c("Tt",paste0("u", seq_len(P)))
      deltahat=matrix(0,N,nmodels);deltatrue=delt(gamm2,xnew)
      for(j in 1:nmodels){
        deltahat[,j]=(predict(hierarchical.models64(X[,1:(P+2)])[[j]],as.data.frame(xnewt))
                      -predict(hierarchical.models64(X[,1:(P+2)])[[j]],as.data.frame(xnewc)))
      }
      mse[k,]=colMeans((deltatrue%*%rep(1,nmodels)-deltahat)^2)
      tecvmse[k]=mean((deltatrue-deltahat[,wn[[1]]])^2) #mse[k,wn[[1]]]
      teemmse[k]=mean((deltatrue-deltahat%*%wn[[2]])^2) #mse[k,]%*%wn[[2]] 
      optmse[k]=mean((deltatrue-deltahat%*%wn[[3]])^2) #mse[k,]%*%wn[[3]] 
      optNmse[k]=mean((deltatrue-deltahat%*%wnN)^2) #mse[k,]%*%wn[[3]] 
      
      aicmse[k]=mean((deltatrue-deltahat[,wn[[4]]])^2) #mse[k,wn[[4]]]
      bicmse[k]=mean((deltatrue-deltahat[,wn[[5]]])^2) #mse[k,wn[[5]]]
      saicmse[k]=mean((deltatrue-deltahat%*%wn[[6]])^2) # mse[k,]%*%wn[[6]] 
      sbicmse[k]=mean((deltatrue-deltahat%*%wn[[7]])^2) #mse[k,]%*%wn[[7]] 
      print(k)
      
    }
    ARisk=round(cbind(mean(aicmse),mean(bicmse),mean(saicmse),mean(sbicmse),mean(tecvmse),mean(teemmse),mean(optmse),mean(optNmse)),digits = 3)
    colnames(ARisk)=c("AIC","BIC","cAIC","BMA","TECV","TEEM","OPT","OPTN")
    Models=round(matrix(colMeans(mse),1,64),digits = 3);colnames(Models)=1:nmodels
    B=cbind(Models,ARisk)
    A=rbind(A,B)
    print(r)
  }
  write.table(A,paste('design',1,'k',K,'n',n,'.txt',sep=''))
}
