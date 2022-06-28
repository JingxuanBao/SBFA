dwl.new <- function(X=NULL,y=NULL,XX=NULL,Xy=NULL,yy=NULL)
{
  dwl <- new.env(parent=emptyenv())
  
  if ( !is.null(X) & !is.null(y) )
  {
    dwl$X <- X
    dwl$y <- as.vector(y)
    
    dwl$n <- nrow(X)
    dwl$p <- ncol(X)
    dwl$m <- min(dwl$n,dwl$p)
    
    dwl$XxMax <- min(3*dwl$n,dwl$p)
    dwl$Xx <- matrix(0,dwl$p,dwl$XxMax)
    dwl$XxIdx <- rep(0,dwl$p)
    dwl$XxCnt <- 0
    
    dwl$Xy <- as.vector(t(X) %*% y)
    dwl$yy <- sum(y^2)
  }
  else if ( !is.null(XX) & !is.null(Xy) )
  {
    dwl$p <- ncol(XX)
    dwl$m <- dwl$p
    
    dwl$XxMax <- dwl$p
    dwl$Xx <- as.matrix(XX)
    dwl$XxIdx <- 1:dwl$p
    dwl$XxCnt <- dwl$p
    
    dwl$Xy <- as.vector(Xy)
    dwl$yy <- as.numeric(yy)
  }
  else
  {
    message("X and y or XX and Xy must be provided.")
    stop()
  }
  
  dwl$lam <- rep(max(abs(dwl$Xy))*1.2,dwl$p)
  
  dwl$Idx <- rep(0,dwl$p)
  dwl$A <- numeric()
  dwl$nA <- 0
  
  dwl$B <- numeric()
  dwl$S <- numeric()
  
  dwl$C <- dwl$Xy
  dwl$iXXa <- matrix(0,dwl$m,dwl$m)
  
  dwl$coef = rep(0,dwl$p)
  dwl$sign = rep(0,dwl$p)
  
  dwl
}


dwl.lookahead <- function(dwl,tlam)
{
  ######################
  # Find the direction #
  ######################
  
  dwl$dlam <- tlam - dwl$lam
  
  # calculate dwl$dB/dalpha and dC/dalpha
  if ( dwl$nA > 0 )
  {
    dwl$dB <- as.vector(-dwl$iXXa[1:dwl$nA,1:dwl$nA] %*% (dwl$S * dwl$dlam[dwl$A]))
    dC <- as.vector(-dwl.getXXI(dwl,dwl$A) %*% dwl$dB)
  }
  else
    dC <- rep(0,dwl$p)
  
  
  ######################
  # How far can we go? #
  ######################
  
  # find breakpoint
  dwl$alpha <- 1
  dwl$type <- 0
  dwl$idx <- 0
  
  if ( dwl$nA > 0 )
  {
    pbp0 = -dwl$B/dwl$dB
    for ( l in 1:dwl$nA )
      if ( (dwl$B[l]+dwl$dB[l])*dwl$S[l] < 0 & pbp0[l] < dwl$alpha )
      {
        dwl$alpha <- pbp0[l]
        dwl$type <- 0
        dwl$idx <- dwl$A[l]
      }
  }
  
  pbp1 = (dwl$lam-dwl$C)/(dC-dwl$dlam)
  pbp2 = -(dwl$lam+dwl$C)/(dC+dwl$dlam)
  for ( k in 1:dwl$p )
    if ( dwl$Idx[k] == 0 )
    {
      if ( dwl$C[k]+dC[k] > dwl$lam[k]+dwl$dlam[k] & pbp1[k] < dwl$alpha )
      {
        dwl$alpha <- pbp1[k]
        dwl$type <- 1
        dwl$idx <- k
      }
      if ( dwl$C[k]+dC[k] < -dwl$lam[k]-dwl$dlam[k] & pbp2[k] < dwl$alpha )
      {
        dwl$alpha <- pbp2[k]
        dwl$type <- -1
        dwl$idx <- k
      }
    }
  
}


dwl.updateBC <- function(dwl)
{
  dwl$coef <- rep(0,dwl$p)
  if ( dwl$nA > 0 )
  {
    dwl$B <- as.vector(dwl$iXXa[1:dwl$nA,1:dwl$nA] %*% ( dwl$Xy[dwl$A] - dwl$S*dwl$lam[dwl$A] ))
    dwl$C <- dwl$Xy - as.vector(dwl.getXXI(dwl,dwl$A) %*% dwl$B)
    dwl$coef[dwl$A] <- dwl$B
  }
  else
  {
    dwl$B <- numeric()
    dwl$C <- dwl$Xy
  }
  dwl$sign <- sign(dwl$coef)
}


dwl.fit <- function(dwl,tlam)
{
  # Eat free lunch
  for ( i in 1:dwl$p )
  {
    if ( dwl$Idx[i] == 0 & dwl$lam[i] < tlam[i] )
      dwl$lam[i] <- tlam[i]
  }
  
  niter = 0
  repeat
  {
    niter = niter + 1
    
    ##############
    # Look ahead #
    ##############
    
    lh = dwl.lookahead(dwl,tlam)
    
    
    ############################
    # Add or remove a variable #
    ############################
    
    if ( dwl$alpha < 1 )
      if ( dwl$type == 0 )
        dwl.remove(dwl,dwl$idx)
    else
    {
      dwl$B <- dwl$B + dwl$alpha*dwl$dB
      dwl.add(dwl,dwl$idx,dwl$type)
    }
    
    
    ##########
    # Update #
    ##########
    
    dwl$lam <- dwl$lam + dwl$dlam*dwl$alpha
    dwl.updateBC(dwl)
    
    
    if ( dwl$alpha ==  1 )
      break
  }
  
  
  list(coef=dwl$coef,sign=dwl$sign,niter=niter)
}



dwl.add = function(dwl,k,sgn)
{
  b = dwl.getXXI(dwl,k)
  
  if ( dwl$nA > 0 )
  {
    a = dwl$iXXa[1:dwl$nA,1:dwl$nA]
    del = drop(a %*% b[dwl$A])
    d = drop(b[k] - crossprod(del,b[dwl$A]))
    
    if ( d < 1e-8 )
    {
      # message("Warning: numerical instability")
      
      pos = which.max(del*sgn/dwl$B)
      dwl.remove(dwl,dwl$A[pos])
      
      if ( dwl$nA > 0 )
      {
        a = dwl$iXXa[1:dwl$nA,1:dwl$nA]
        del = drop(a %*% b[dwl$A])
        d = drop(b[k] - crossprod(del,b[dwl$A]))
      }
    }
  }
  
  # Now add k
  if ( dwl$nA > 0 )
  {
    dwl$iXXa[1:dwl$nA,1:dwl$nA] <- a + del %*% t(del) / d
    dwl$iXXa[1:dwl$nA,dwl$nA+1] <- -del / d
    dwl$iXXa[dwl$nA+1,1:dwl$nA] <- -del / d
    dwl$iXXa[dwl$nA+1,dwl$nA+1] <- 1/d
  }
  else
  {
    dwl$iXXa[1] <- 1/b[k]
  }
  
  dwl$nA <- dwl$nA+1
  dwl$Idx[k] <- dwl$nA
  dwl$A <- c(dwl$A,k)
  dwl$S <- c(dwl$S,sgn)
}



dwl.remove = function(dwl,k)
{
  l = dwl$Idx[k]
  dwl$Idx[k] <- 0
  if ( l<dwl$nA )
    dwl$Idx[dwl$A[(l+1):dwl$nA]] <- dwl$Idx[dwl$A[(l+1):dwl$nA]] - 1
  dwl$A <- dwl$A[-l]
  dwl$S <- dwl$S[-l]
  
  if ( dwl$nA>1 )
  {
    a = dwl$iXXa[1:dwl$nA,1:dwl$nA]
    b = a[,l]
    dwl$iXXa[1:(dwl$nA-1),1:(dwl$nA-1)] <- a[-l,-l] - b[-l] %*% t(b[-l]) / b[l]
  }
  dwl$iXXa[,dwl$nA] <- 0
  dwl$iXXa[dwl$nA,] <- 0
  
  dwl$nA <- dwl$nA-1
}



dwl.getXXI <- function(dwl,I)
{
  for ( k in I )
    if ( dwl$XxIdx[k] == 0 )
    {
      dwl$XxCnt <- dwl$XxCnt + 1
      if ( dwl$XxCnt > dwl$XxMax )
      {
        oldmax = dwl$XxMax
        oldXx = dwl$Xx
        dwl$XxMax <- min(oldmax*2,dwl$p)
        dwl$Xx <- matrix(0,dwl$p,dwl$XxMax)
        dwl$Xx[,1:oldmax] <- oldXx
      }
      dwl$XxIdx[k] <- dwl$XxCnt
      dwl$Xx[,dwl$XxCnt] <- t(dwl$X) %*% dwl$X[,k]
    }
  
  as.matrix(dwl$Xx[,dwl$XxIdx[I]])
}

init <- function(X,type,param)
{
  p = nrow(X)
  n = ncol(X)
  
  Y = matrix(0,p,n)
  for ( j in 1:p )
  {
    if ( type[j] == 0 )
    {
      Y[j,] = X[j,]
    }
    if ( type[j] == 1 )
    {
      pbar = (X[j,]+1)/(param[j]+2)
      Y[j,] = log(pbar/(1-pbar))
    }
    if ( type[j] == 2 )
    {
      pbar = (X[j,]+1)/(param[j]+X[j,]+2)
      Y[j,] = log(pbar/(1-pbar))
    }
  }
  return(Y)
}

llk <- function(X,mu,type,param)
{
  p = nrow(X)
  n = ncol(X)
  
  l = 0
  for ( j in 1:p )
  {
    if ( type[j] == 0 )
      l = l - n*log(mean((X[j,]-mu[j,])^2))/2
    else if ( type[j] == 1 )
      l = l + sum(dbinom(X[j,],param[j],1/(1+exp(-mu[j,])),TRUE))
    else if ( type[j] == 2 )
      l = l + sum(dnbinom(X[j,],param[j],1/(1+exp(mu[j,])),log=TRUE))
  }
  l
}

GenLambert <- function(t,r,a)
{
  x = t+a/r
  while (TRUE)
  {
    f = (x-t)*(exp(x/2)+r)-a
    g = (x/2-t/2+1)*exp(x/2)+r
    x = x - f/g
    if ( abs(f/g) < 1e-10 ) break
  }
  x
}

SBFA_EM <- function(X,type,param,E,L,v1,v2,a_omega,b_omega,m.init=1,scale=T,W.init=NULL,eps=1e-3,maxIter=500)
{
  # Obtain dimension for X
  p = nrow(X)
  n = ncol(X)
  print(paste0("The input data matrix has ",p," rows (features) and ",n," columns (samples)."))
  print(paste0("The number of factors are assumed to be ",L))
  
  # get adjency matrix for omega
  E = unique(rbind(E,E[,c(2,1)]))
  adj = matrix(0,nrow = p,ncol = p)
  adj[E] = 1
  
  # initialization (see unified likelihood function)
  psi = matrix(0,p,n)
  kappa = matrix(0,p,n)
  bji = matrix(4,p,n) 
  for ( j in 1:p ){
    if ( type[j] == 0 ){
      psi[j,] = X[j,]
    }
    if ( type[j] == 1 ){
      kappa[j,] = X[j,]-param[j]/2
      bji[j,] = param[j]
    }
    if ( type[j] == 2 ){
      kappa[j,] = (X[j,]-param[j])/2
      bji[j,] = param[j]+X[j,]
    }
  }
  
  
  # initialize Y 
  Y = init(X,type,param) 
  
  # initialize m
  if ( length(m.init) == p ){
    m = m.init
  }else if ( m.init == 0 ){
    m = rep(0,p)
  }else if ( m.init == 1 ){
    m = apply(Y,1,function(x) mean(x,0.2))
  }else if ( m.init == 2 ){
    m = apply(Y,1,median)
  }else if ( m.init == 3 ){
    m = apply(Y,1,meam)
  }else{
    m = rep(0,p)
  }
  
  # Find scale
  if ( scale ){
    g = sqrt(apply((Y-m)^2,1,mean)) 
    g[which(g==0)]<-1
  }else{
    g = rep(1,p)
  }
  
  # W and Z init 
  WZinit = (Y-m)/g
  if ( !is.null(W.init) ){
    W = W.init
    WW = t(W)%*%W
    Z = chol2inv(chol(WW)) %*% t(W) %*% WZinit
  }else{
    init_svd = svd(WZinit,L,L)
    W = init_svd$u %*% diag(init_svd$d[1:L],L) / sqrt(n)
    Z = t(init_svd$v) * sqrt(n)
  }
  Sigma_Z = array(diag(1,L),c(L,L,n))
  precision_Z = Sigma_Z
  
  # Initialization for omega
  omega = matrix(0,p,p)
  omega[E] = -rgamma(nrow(E),a_omega,b_omega)
  
  omega[lower.tri(omega)] = 0
  omega = omega+t(omega)
  diag(omega) = 1-rowSums(omega)
  
  # Initialization for Sigma_alpha
  Sigma_alpha = replicate(L,v2*chol2inv(chol(omega))) 
  
  # Initialization for mu_alpha
  mu_alpha = matrix(v1, nrow = p,ncol = L)  
  
  # Initialization for mu_rho (Use prior mean)
  mu_rho = bji/4
  
  # Initialization for varphi (values will not be used)
  varphi = matrix(0,nrow = p,ncol = n)
  
  # Initialization for mu_Z (values will not be used)
  mu_Z = Z
  
  iter = 0
  while ( iter < maxIter ){
    iter = iter + 1
    W_old=W
    if(iter%%100==0||iter==1){
      print(iter)
    }
    
    # E-step for omega
    A = matrix(a_omega,nrow = p,ncol = p)
    B = matrix(b_omega,nrow = p,ncol = p)
    for (j in 1:p) {
      for (k in 1:p) {
        B[j,k] = B[j,k] + (sum((mu_alpha[j,] - mu_alpha[k,])^2) + sum(Sigma_alpha[k,k,]) + sum(Sigma_alpha[j,j,]))/(2*v2)
      }
    }
    omega<- -A/B #calculate mean
    diag(omega) = 0
    omega = omega * adj
    
    diag(omega) = 1-rowSums(omega)
    
    # E-step for mu_alpha
    for (j in 1:p) {
      for (l in 1:L) {
        a = abs(W[j,l]) * exp( 0.5 * Sigma_alpha[j,j,l] )
        b = 1/v2 * omega[j,j]
        c = 1/v2 * omega[j,-j]%*%mu_alpha[-j,l] - (v1+v2)/v2
        #cat(a,b,c,"\n")
        #cat(j,a/b,c/b,"\n")
        mu_alpha[j,l] = - lambertW0(a/b * exp(-c/b)) - c/b
      }
    }
    
    # E-step for Sigma_alpha
    for (l in 1:L) {
      Delta = Sigma_alpha[,,l]%*%omega
      for (j in 1:p) {
        Delta[-j,] = Delta[-j,] - outer(Sigma_alpha[-j,j,l],omega[j,])
        if (sum(omega[j,-j]==0)==(p-1)){
          oSo = 0
          c = 0
          Sigma_alpha[-j,j,l] = 0
          Sigma_alpha[j,-j,l] = 0
        }else{
          oSo = crossprod(omega[j,-j],Delta[-j,j])
          numerator = v2 - sqrt(v2^2 + 4 * oSo * Sigma_alpha[j,j,l])
          denominator = 2 * oSo
          c = numerator/denominator
          Sigma_alpha[-j,j,l] = as.numeric(c) * Delta[-j,j]
          Sigma_alpha[j,-j,l] = Sigma_alpha[-j,j,l]
        }
        
        t.alpha = c^2*oSo
        if (W[j,l]!=0){
          a = 1/(abs(W[j,l]) * exp(mu_alpha[j,l]))
          r = omega[j,j]/v2 * a
          Sigma_alpha[j,j,l] = GenLambert(t.alpha,r,a)
        }else{
          Sigma_alpha[j,j,l] = t.alpha + v2/omega[j,j]
        }
        Delta[-j,] = Delta[-j,] + outer(Sigma_alpha[-j,j,l],omega[j,])
        Delta[j,] = Sigma_alpha[j,,l]%*%omega
      }
    }
    
    # E-step for precision_Z
    for (i in 1:n) {
      precZ = diag(L) + t(W)%*%(W*mu_rho[,i])
      Sigma_Z[,,i] = chol2inv(chol(precZ))
    }
    temp = t(W) %*% ((psi-m)*mu_rho+kappa)
    for (i in 1:n) {
      mu_Z[,i] =  Sigma_Z[,,i] %*% temp[,i] 
    }
    WSW = apply(array((W%*%matrix(Sigma_Z,L))*matrix(W,p,n*L),c(p,L,n)),c(1,3),sum)
    varphi2 = WSW + (W%*%mu_Z+m - psi)^2
    varphi = sqrt(varphi2)
    
    for (j in 1:p) {
      if (type[j]==0){
        mu_rho[j,] = (param[j]+n)/(param[j]+sum(varphi2[j,]))
      }else{
        for (i in 1:n) {
          if(varphi[j,i]<=100){
            if(varphi[j,i]<=0.0001){
              mu_rho[j,i] = bji[j,i]*(1+varphi[j,i]/2)/2/(2+varphi[j,i]+varphi[j,i]^2/2)
            }else{
              mu_rho[j,i] = ( bji[j,i] * (exp(varphi[j,i])-1) )/( 2 * varphi[j,i] * (exp(varphi[j,i]) + 1) )
            }
          }else{
            mu_rho[j,i] = ( bji[j,i] )/( 2 * varphi[j,i] )
          }
        }
      }
    }
    
    # M-step for m
    for (j in 1:p) {
      temp = 0
      for (i in 1:n) {
        temp = temp + mu_rho[j,i] * (t(W[j,]) %*% mu_Z[,i])
      }
      m[j] = 1/sum(mu_rho[j,]) * (sum(kappa[j,]) - temp + sum(mu_rho[j,] * psi[j,]))
    }
    
    for (j in 1:p) {
      A_DWL = apply(Sigma_Z*rep(mu_rho[j,],each=L^2),1:2,sum) + mu_Z%*%(t(mu_Z)*mu_rho[j,])
      B_DWL = drop(mu_Z%*%((psi[j,]-m[j])*mu_rho[j,]+kappa[j,]))
      dwl = dwl.new(XX=A_DWL,Xy=B_DWL)
      dwl.fit(dwl,exp(mu_alpha[j,]+0.5*Sigma_alpha[j,j,]))
      W[j,] = dwl$coef
      rm(dwl)
    }
    
    mu = m + W%*%mu_Z
    df = sum(W!=0)
    LLK = llk(X,mu,type,param)
    BIC = -2*LLK + log(n)*df
    
    if (max(abs(W-W_old))<eps){
      break
    }
    
  }
  
  return(list(p=p,n=n,L=L,E=E,m=m,W=W,mu_alpha=mu_alpha,mu_rho=mu_rho,mu_Z=mu_Z,Sigma_alpha=Sigma_alpha,Sigma_Z=Sigma_Z,Omega=omega,BIC=BIC,iter=iter))
}