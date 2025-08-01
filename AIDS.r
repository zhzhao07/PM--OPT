
TOTAL_combo1 <- function(
    x,                # 数据框
    resp, trt,        # 响应、处理列名
    u, std.u,         # 协变量及其标准化列名
    model.function,   # e.g. hierarchical.models4
    nmodels,
    nTree = 2000, nodeSize = 20
) {
  
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
  base_models <- model.function(x)
  
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
  resid1 <- Y - predict(base_models[[nmodels]], newdata = x)
  sigma2 <- sum(resid1^2) / (nrow(x) - length(coef(base_models[[nmodels]])))
  dvec  <- matrix(c(-sigma2*lct(x)),nmodels,1) #
  
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


rtt=as.data.frame(read.table("first_data.txt", header=TRUE))
# attach(rtt)
#xx=rtt
xy=as.matrix(rtt)
#as.numeric(x[2][which(rt[2]=="t"),])
xy[which(rtt[,2]=="t"),2]=1
xy[which(rtt[,2]=="c"),2]=0
xx=as.data.frame(matrix(as.numeric(xy),ncol=5,byrow=FALSE))
colnames(xx)=colnames(rtt)
z=as.data.frame(xx)
colnames(z)=c("Y","Tt","CD4_0","RNA_0","Age")

u=c("CD4_0","RNA_0","Age")
std.u=c("CD4_0.u","RNA_0.u","Age.u")
x.u=z[,u]

std.uu=scale(x.u)
colnames(std.uu)=c("CD4_0.u","RNA_0.u","Age.u")


X1=cbind(z,std.uu)
Y1=X1$Y
#ot=1:5#c(5,6,14,15,22)
MM=22#length(ot)

hierarchical.models3=function(x)
{
  Y=x[,1];Tt=x[,2];CD4_0=x[,3];RNA_0=x[,4]; Age=x[,5];
  sol1=lm(formula =Y ~ CD4_0)
  sol2=lm(formula =Y ~ Age + CD4_0)
  sol3=lm(formula =Y ~ RNA_0 + CD4_0)
  sol4=lm(formula =Y ~ Age + RNA_0 + CD4_0)
  
  sol5=lm(formula =Y ~ Tt + CD4_0)
  sol6=lm(formula =Y ~ Tt+ CD4_0 + Tt:CD4_0)
  sol7=lm(formula =Y ~ Tt+ CD4_0 +Age)
  sol8=lm(formula =Y ~ Tt+ CD4_0 + RNA_0)
  
  
  sol9=lm(formula =Y ~ Tt + Age + Tt:CD4_0 + CD4_0)
  sol10=lm(formula =Y ~ Tt + Age + RNA_0 + CD4_0)
  sol11=lm(formula =Y ~ Tt + RNA_0 + Tt:CD4_0 + CD4_0)
  sol12=lm(formula =Y ~ Tt + Age + Tt:Age+ CD4_0)
  sol13=lm(formula =Y ~ Tt + RNA_0 + Tt:RNA_0+ CD4_0)
  
  sol14=lm(formula =Y ~ Tt + Age + Tt:CD4_0 + Tt:Age + CD4_0)
  sol15=lm(formula =Y ~ Tt + RNA_0 + Tt:CD4_0 + Tt:RNA_0 + CD4_0)
  sol16=lm(formula =Y ~ Tt + RNA_0 + Tt:CD4_0 + Age + CD4_0)
  sol17=lm(formula =Y ~ Tt + RNA_0 + Tt:Age + Age + CD4_0)
  sol18=lm(formula =Y ~ Tt + RNA_0 + Tt:RNA_0 + Age + CD4_0)
  
  sol19=lm(formula =Y ~ Tt+ RNA_0  + Age + Tt:CD4_0 + Tt:Age + CD4_0)
  sol20=lm(formula =Y ~ Tt+ RNA_0  + Age + Tt:CD4_0 + Tt:RNA_0 + CD4_0)
  sol21=lm(formula =Y ~ Tt+ RNA_0  + Age + Tt:Age + Tt:RNA_0 + CD4_0)
  
  sol22=lm(formula =Y ~ Tt+ RNA_0  + Age + Tt:Age + Tt:RNA_0+ Tt:CD4_0 + CD4_0)
  list(sol1,sol2,sol3,sol4,sol5,sol6,sol7,sol8,sol9,sol10,sol11,sol12,sol13,
       sol14,sol15,sol16,sol17,sol18,sol19,sol20,sol21,sol22)
}


TOTAL1 <- function(
    x,                ###data frame on which to apply the data splitting
    resp,             ###name of the response variable in x 
    trt,              ###name of the treatment variable in x.  Should =1 if treated, 0 otherwise.
    u,                ###names of the covariates used in the model
    std.u,            ###names of the standardized covariates used to compute distances 
    nsplits,      ###number of data splittings over which to average
    dev.fraction=0.5, ###fraction of the data to be used for model development 
    model.function,   ###function should take a dataset x, construct the candidate procedures
    ###being evaluated on x, and return a list of models on which the 
    ###predict() function can be used
    nmodels,           ###length of the model list that will be created by model.function
    deltatrue,theta)
{
  deltahatt=ehatt=matrix(NA,ncol=nmodels,nrow=nrow(x)); Aic=Bic=rep(NA,nmodels)
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
  trn.truet=trn.truec=as.data.frame(trn.true); trn.truet$Tt=1; trn.truec$Tt=0
  lc=lct(x)#c(6,5,5,4,4,3,4,3,2)
  #deltahatt=ehatt
  for (j in 1:nmodels) { 
    deltahatt[,j]=(predict(model.function(x)[[j]],trn.truet)-predict(model.function(x)[[j]], trn.truec))
    
    ehatt[,j]=trn.nn.delta.tilde-deltahatt[,j]#(predict(model.function(x)[[j]],trn.truet)-predict(model.function(x)[[j]], trn.truec))
    Aic[j]=(length(x$Y)*log(mean((x$Y-predict(model.function(x)[[j]],x))^2))+(lc[j])*2)
    Bic[j]=(length(x$Y)*log(mean((x$Y-predict(model.function(x)[[j]],x))^2))+(lc[j])*log(length(x$Y)))
  }
  aic=which.min(Aic); bic=which.min(Bic)
  wa=exp(0.5*min(Aic)-0.5*Aic)/sum(exp(0.5*min(Aic)-0.5*Aic)) #S-AIC权
  wb=exp(0.5*min(Bic)-0.5*Bic)/sum(exp(0.5*min(Bic)-0.5*Bic)) #S-BIC权重
  
  sigmatc=sum((x$Y-predict(model.function(x)[[nmodels]],x))^2)/(length(x$Y)-lc[nmodels])#summary(hierarchical.models3(x)[[62]])$sigma^2
  a1 <- t(ehatt[,1:(nmodels)]) %*% ehatt[,1:(nmodels)]+diag(0.000001,nmodels)
  a2<-matrix(c(-sigmatc*lc),nmodels,1) #
  a3 <- t(rbind(matrix(1,nrow=1,ncol=nmodels),diag(nmodels),-diag(nmodels)))
  a4 <- rbind(1,matrix(0,nrow=nmodels,ncol=1),matrix(-1,nrow=nmodels,ncol=1))
  QP <- solve.QP(a1,a2,a3,a4,1)
  mw <- QP$solution
  mw <- as.matrix(mw)
  mw <- mw*(mw>0)
  mw <- mw/sum(mw)  
  
  w=matrix(NA,ncol=nmodels,nrow=nsplits)
  ntrain<-nrow(x)
  cate.se.cv<-matrix(NA,nrow=nsplits,ncol=nmodels)
  
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
    
    val.data.trt.paired<-val.data.trt[,c(u,resp,trt)]
    val.nn.delta.tilde<-val.data.trt.paired[,resp] - val.data.ctrl[nn.indices,resp]
    val.data.trt1=val.data.trt.paired[,resp];val.data.ctrl1=val.data.ctrl[nn.indices,resp]
    val.data.sim.ctrl<-val.data.trt.paired[,c(u,trt)]; val.data.sim.ctrl[,trt]<-0
    
    ###Now match each obs in the control group to its nearest treated obs
    nearest.neighbors<-matchpt(as.matrix(val.data.ctrl[,std.u]), 
                               as.matrix(val.data.trt[,std.u]))
    nn.indices<-nearest.neighbors[,1]
    
    val.data.ctrl.paired<-val.data.ctrl[,c(u,resp,trt)]  
    val.nn.delta.tilde2<-val.data.trt[nn.indices,resp] - val.data.ctrl.paired[,resp] 
    val.data.trt2=val.data.trt[nn.indices,resp]; val.data.ctrl2=val.data.ctrl.paired[,resp] 
    val.data.sim.trt<-val.data.ctrl.paired[,c(u,trt)]; val.data.sim.trt[,trt]<-1
    
    ###Append results 
    val.data.trta=c(val.data.trt1,val.data.trt2); val.data.ctrla=c(val.data.ctrl1,val.data.ctrl2)
    val.nn.delta.tilde<-c(val.nn.delta.tilde, val.nn.delta.tilde2)
    val.data.trt.paired<-rbind(val.data.trt.paired[,c(u,trt)], val.data.sim.trt)
    val.data.sim.ctrl<-rbind(val.data.sim.ctrl, val.data.ctrl.paired[,c(u,trt)])
    
    ###Build models on development data
    models.train<-model.function(trn.data)
    
    nn.sspe<-rep(NA,nmodels)      
    sigmaht=matrix(NA,nrow=length(val.inds),ncol=nmodels)
    
    for (j in 1:nmodels) { 
      val.delta.hats<-predict(models.train[[j]], val.data.trt.paired) - 
        predict(models.train[[j]], val.data.sim.ctrl)
      nn.sspe[j]<-sum((val.delta.hats - val.nn.delta.tilde)^2)
      sigmaht[,j]=sqrt((val.data.trta-predict(models.train[[j]], val.data.trt.paired))^2+
                         (val.data.ctrla-predict(models.train[[j]], val.data.sim.ctrl))^2)
      w[i,j]=prod((20*dnorm((val.nn.delta.tilde-val.delta.hats)/(sigmaht[,j])))/(sigmaht[,j]))
    }
    w[i,]=w[i,]/sum(w[i,])
    cate.se.cv[i,]<-nn.sspe/length(val.nn.delta.tilde)
  }
  wn=colMeans(w)#
  TE.CV<-apply(cate.se.cv,2,mean)
  wx=list(which.min(TE.CV),wn,mw,aic,bic,wa,wb)
  
  # tecvdeltahat=deltahatt[,wx[[1]]]
  # teemdeltahat=deltahatt%*%wx[[2]]
  # optdeltahat=deltahatt%*%wx[[3]]
  # aicdeltahat=deltahatt[,wx[[4]]]
  # bicdeltahat=deltahatt[,wx[[5]]]
  # saicdeltahat=deltahatt%*%wx[[6]]
  # sbicdeltahat=deltahatt%*%wx[[7]]
  
  #mse=colMeans((deltatrue%*%rep(1,nmodels)-deltahatt)^2)
  #tecvmse=mean((deltatrue-deltahatt[,wx[[1]]])^2) #mse[k,wn[[1]]]
  #teemmse=mean((deltatrue-deltahatt%*%wx[[2]])^2) #mse[k,]%*%wn[[2]] 
  #optmse=mean((deltatrue-deltahatt%*%wx[[3]])^2) #mse[k,]%*%wn[[3]] 
  #aicmse=mean((deltatrue-deltahatt[,wx[[4]]])^2) #mse[k,wn[[4]]]
  #bicmse=mean((deltatrue-deltahatt[,wx[[5]]])^2) #mse[k,wn[[5]]]
  #saicmse=mean((deltatrue-deltahatt%*%wx[[6]])^2) # mse[k,]%*%wn[[6]] 
  #sbicmse=mean((deltatrue-deltahatt%*%wx[[7]])^2) #mse[k,]%*%wn[[7]] 
  #ARisk=round(cbind((aicmse),(bicmse),(saicmse),(sbicmse),(tecvmse),(teemmse),(optmse)),digits = 4)
  #colnames(ARisk)=c("AIC","BIC","cAIC","BMA","TECV","TEEM","OPT")
  #Models=round(matrix((mse),1,nmodels),digits = 4);colnames(Models)=1:nmodels
  #B=cbind(Models,ARisk)
  # if(theta==1)  return((wx))   
  #else  return((B))   
  #  return(list(tecvdeltahat,teemdeltahat,optdeltahat,aicdeltahat,bicdeltahat,saicdeltahat,sbicdeltahat))
  return((wx))  }




TOTAL <- function(
    x,                ###data frame on which to apply the data splitting
    resp,             ###name of the response variable in x 
    trt,              ###name of the treatment variable in x.  Should =1 if treated, 0 otherwise.
    u,                ###names of the covariates used in the model
    std.u,            ###names of the standardized covariates used to compute distances 
    nsplits,      ###number of data splittings over which to average
    dev.fraction=0.5, ###fraction of the data to be used for model development 
    model.function,   ###function should take a dataset x, construct the candidate procedures
    ###being evaluated on x, and return a list of models on which the 
    ###predict() function can be used
    nmodels,           ###length of the model list that will be created by model.function
    deltatrue,theta)
{
  deltahatt=ehatt=matrix(NA,ncol=nmodels,nrow=nrow(x)); Aic=Bic=rep(NA,nmodels)
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
  trn.truet=trn.truec=as.data.frame(trn.true); trn.truet$Tt=1; trn.truec$Tt=0
  lc=lct(x)#c(6,5,5,4,4,3,4,3,2)
  #deltahatt=ehatt
  for (j in 1:nmodels) { 
    deltahatt[,j]=(predict(model.function(x)[[j]],trn.truet)-predict(model.function(x)[[j]], trn.truec))
    
    ehatt[,j]=trn.nn.delta.tilde-deltahatt[,j]#(predict(model.function(x)[[j]],trn.truet)-predict(model.function(x)[[j]], trn.truec))
    Aic[j]=(length(x$Y)*log(mean((x$Y-predict(model.function(x)[[j]],x))^2))+(lc[j])*2)
    Bic[j]=(length(x$Y)*log(mean((x$Y-predict(model.function(x)[[j]],x))^2))+(lc[j])*log(length(x$Y)))
  }
  aic=which.min(Aic); bic=which.min(Bic)
  wa=exp(0.5*min(Aic)-0.5*Aic)/sum(exp(0.5*min(Aic)-0.5*Aic)) #S-AIC权
  wb=exp(0.5*min(Bic)-0.5*Bic)/sum(exp(0.5*min(Bic)-0.5*Bic)) #S-BIC权重
  
  sigmatc=sum((x$Y-predict(model.function(x)[[nmodels]],x))^2)/(length(x$Y)-lc[nmodels])#summary(hierarchical.models3(x)[[62]])$sigma^2
  a1 <- t(ehatt) %*% ehatt+diag(0.00001,MM)
  a2<-matrix(c(-sigmatc*lc),MM,1) #
  a3 <- t(rbind(matrix(1,nrow=1,ncol=MM),diag(MM),-diag(MM)))
  a4 <- rbind(1,matrix(0,nrow=MM,ncol=1),matrix(-1,nrow=MM,ncol=1))
  QP <- solve.QP(a1,a2,a3,a4,1)
  mw1 <- QP$solution
  mw1 <- as.matrix(mw1)
  mw1 <- mw1*(mw1>0)
  mw1 <- mw1/sum(mw1)  
  mw=mw1#rep(0,nmodels)
  #mw[ot]=mw1
  res <- TOTAL_combo1(
    x=X1, resp="Y", trt="Tt", u=c("CD4_0","RNA_0","Age"),
    std.u=c("CD4_0.u","RNA_0.u","Age.u"), 
    model.function=function(x){hierarchical.models3(x)}, nmodels=22
  )
  wnN=res$w_OPT
  
  w=matrix(NA,ncol=nmodels,nrow=nsplits)
  ntrain<-nrow(x)
  cate.se.cv<-matrix(NA,nrow=nsplits,ncol=nmodels)
  
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
    
    val.data.trt.paired<-val.data.trt[,c(u,resp,trt)]
    val.nn.delta.tilde<-val.data.trt.paired[,resp] - val.data.ctrl[nn.indices,resp]
    val.data.trt1=val.data.trt.paired[,resp];val.data.ctrl1=val.data.ctrl[nn.indices,resp]
    val.data.sim.ctrl<-val.data.trt.paired[,c(u,trt)]; val.data.sim.ctrl[,trt]<-0
    
    ###Now match each obs in the control group to its nearest treated obs
    nearest.neighbors<-matchpt(as.matrix(val.data.ctrl[,std.u]), 
                               as.matrix(val.data.trt[,std.u]))
    nn.indices<-nearest.neighbors[,1]
    
    val.data.ctrl.paired<-val.data.ctrl[,c(u,resp,trt)]  
    val.nn.delta.tilde2<-val.data.trt[nn.indices,resp] - val.data.ctrl.paired[,resp] 
    val.data.trt2=val.data.trt[nn.indices,resp]; val.data.ctrl2=val.data.ctrl.paired[,resp] 
    val.data.sim.trt<-val.data.ctrl.paired[,c(u,trt)]; val.data.sim.trt[,trt]<-1
    
    ###Append results 
    val.data.trta=c(val.data.trt1,val.data.trt2); val.data.ctrla=c(val.data.ctrl1,val.data.ctrl2)
    val.nn.delta.tilde<-c(val.nn.delta.tilde, val.nn.delta.tilde2)
    val.data.trt.paired<-rbind(val.data.trt.paired[,c(u,trt)], val.data.sim.trt)
    val.data.sim.ctrl<-rbind(val.data.sim.ctrl, val.data.ctrl.paired[,c(u,trt)])
    
    ###Build models on development data
    models.train<-model.function(trn.data)
    
    nn.sspe<-rep(NA,nmodels)      
    sigmaht=matrix(NA,nrow=length(val.inds),ncol=nmodels)
    
    for (j in 1:nmodels) { 
      val.delta.hats<-predict(models.train[[j]], val.data.trt.paired) - 
        predict(models.train[[j]], val.data.sim.ctrl)
      nn.sspe[j]<-sum((val.delta.hats - val.nn.delta.tilde)^2)
      sigmaht[,j]=sqrt((val.data.trta-predict(models.train[[j]], val.data.trt.paired))^2+
                         (val.data.ctrla-predict(models.train[[j]], val.data.sim.ctrl))^2)
      w[i,j]=prod((20*dnorm((val.nn.delta.tilde-val.delta.hats)/(sigmaht[,j])))/(sigmaht[,j]))
    }
    w[i,]=w[i,]/sum(w[i,])
    cate.se.cv[i,]<-nn.sspe/length(val.nn.delta.tilde)
  }
  wn=colMeans(w)#
  TE.CV<-apply(cate.se.cv,2,mean)
  wx=list(which.min(TE.CV),wn,mw,aic,bic,wa,wb)
  
  
  mse=colMeans((deltatrue%*%rep(1,nmodels)-deltahatt)^2)
  tecvmse=mean((deltatrue-deltahatt[,wx[[1]]])^2) #mse[k,wn[[1]]]
  teemmse=mean((deltatrue-deltahatt%*%wx[[2]])^2) #mse[k,]%*%wn[[2]] 
  optmse=mean((deltatrue-deltahatt%*%wx[[3]])^2) #mse[k,]%*%wn[[3]] 
  optNmse=mean((deltatrue-deltahatt%*%wnN)^2) #mse[k,]%*%wn[[3]] 
  aicmse=mean((deltatrue-deltahatt[,wx[[4]]])^2) #mse[k,wn[[4]]]
  bicmse=mean((deltatrue-deltahatt[,wx[[5]]])^2) #mse[k,wn[[5]]]
  saicmse=mean((deltatrue-deltahatt%*%wx[[6]])^2) # mse[k,]%*%wn[[6]] 
  sbicmse=mean((deltatrue-deltahatt%*%wx[[7]])^2) #mse[k,]%*%wn[[7]] 
  ARisk=round(cbind((aicmse),(bicmse),(saicmse),(sbicmse),(tecvmse),(teemmse),(optmse),optNmse),digits = 4)
  colnames(ARisk)=c("AIC","BIC","cAIC","BMA","TECV","TEEM","OPT","OPTN")
  Models=round(matrix((mse),1,nmodels),digits = 4);colnames(Models)=1:nmodels
  B=cbind(Models,ARisk)
  if(theta==1)  return((wx))   
  else  return((B))   
}

lct=function(x)
{lc=1:nmodels
for(i in 1:nmodels){lc[i]=length(hierarchical.models3(x)[[i]]$coefficients)}
lc
}



set.seed(3)
n=length(Y1);P=3;nmodels=22
K=100;
#mse=matrix(0,K,nmodels)
#optmse=saicmse=sbicmse=teemmse=tecvmse=aicmse=bicmse=rep(0,K);
A=matrix(0,K,8+nmodels)
for(k in 1:K)#循环次数
  
{  
  trt="Tt"
  x.trt<-subset(X1, get(trt)==1); 
  x.ctrl<-subset(X1, get(trt)==0)
  trn.true=rbind(x.trt,x.ctrl)[,c(u,trt)]
  trn.truet=trn.truec=as.data.frame(trn.true); trn.truet$Tt=1; trn.truec$Tt=0
  # gamm=hierarchical.models2(X)[[Z[[1]]]]$coef
  mu=predict(hierarchical.models3(X1)[[19]],X1);sig=sqrt(mean((Y1-predict(hierarchical.models3(X1)[[19]],X1))^2))
  #summary(hierarchical.models3(X1)[[19]])$sigma
  e=rnorm(n,0,sig);#e2=rnorm(n,0,1);e=Tt*e1+(1-Tt)*e2;
  deltatrue=predict(hierarchical.models3(X1)[[19]],trn.truet)-predict(hierarchical.models3(X1)[[19]],trn.truec)
  Y=mu+e
  
  u=c("CD4_0","RNA_0","Age")
  Xy=cbind(Y,X1[,c("Tt",u,std.u)])
  colnames(Xy)=c("Y","Tt","CD4_0","RNA_0","Age","CD4_0.u","RNA_0.u","Age.u")
  X=as.data.frame(Xy)
  
  rest=TOTAL(x=X, resp="Y", trt="Tt", u=c("CD4_0","RNA_0","Age"),
             std.u=c("CD4_0.u","RNA_0.u","Age.u"), nsplits=100, 
             dev.fraction=0.5,
             model.function=function(x){hierarchical.models3(x)}, nmodels=22,deltatrue=as.matrix(deltatrue),theta=0)
  A[k,]=as.matrix(rest)
  print(k)
}
