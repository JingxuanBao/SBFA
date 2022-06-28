rm(list = ls())

L = c(2,3,4) 
v1 = c(4,5) 
v2 = log(2) 
a_omega = c(1,3,5)
b_omega = 1
maxIter=500
edge_prop<-1

n_fea<-20 # number of features in each modality

set.seed(23)
# First dataset - 20 features, 3 latent variables
W1 = matrix(0,nrow = n_fea,ncol = 3)
W1[1:(n_fea/2),1] = runif(n_fea/2,min = 1.5,max = 2.5)
W1[(n_fea/2+1):n_fea,2] = runif(n_fea/2,min = 1.5,max = 2.5)
# Second Dataset - 20 features, 3 latent variables
W2 = matrix(0,nrow = n_fea,ncol = 3)
W2[1:(n_fea/2),1] = runif(n_fea/2,min = 1.5,max = 2.5)
W2[(n_fea/2+1):n_fea,3] = runif(n_fea/2,min = 1.5,max = 2.5)
# Third Dataset - 20 features, 3 latent variables
W3 = matrix(0,nrow = n_fea,ncol = 3)
W3[1:(n_fea/2),2] = runif(n_fea/2,min = 1.5,max = 2.5)
W3[(n_fea/2+1):n_fea,3] = runif(n_fea/2,min = 1.5,max = 2.5)

# Ground Truth W - combining W1,W2,W3 - 90 features
W_true=rbind(W1,W2,W3)
# Randomly assign negative values to the ground truth W
W_true=W_true*matrix(sign(runif(n_fea*3*3,min = -1,max = 1)),nrow = nrow(W_true),ncol = ncol(W_true))

# Generate factor matrix Z
Z_true = matrix(rnorm(3*160),nrow = 3)

# Generate the mean of data distribution
mu_true = W_true%*%Z_true

prob<- 1/(1+exp(-mu_true))
trail<-matrix(sample(1:15,size=(nrow(prob)),replace = T),nrow = nrow(prob),ncol = ncol(prob),byrow = F)

X1 = matrix(NA,nrow = nrow(mu_true)/3,ncol = ncol(mu_true))
for (i in 1:nrow(X1)) {
  for (j in 1:ncol(X1)) {
    X1[i,j] = rbinom(1,size = 1,prob = prob[i,j])
  }
}
X2 = matrix(NA,nrow = nrow(mu_true)/3,ncol = ncol(mu_true))
for (i in 1:nrow(X2)) {
  for (j in 1:ncol(X2)) {
    X2[i,j] = rbinom(1,size = trail[(nrow(X1)+i),j],prob = prob[(nrow(X1)+i),j])
  }
}
X3 = matrix(NA,nrow = nrow(mu_true)/3,ncol = ncol(mu_true))
for (i in 1:nrow(X3)) {
  for (j in 1:ncol(X3)) {
    X3[i,j] = rnorm(1,mean = mu_true[(nrow(X1)+nrow(X2)+i),j],sd=2)
    #X3[i,j] = rnorm(1,mean = mu_true[(nrow(X1)+nrow(X2)+i),j],sd=2)
  }
}

# Generate the data according to the assumed distribution
X_simu=rbind(X1,X2,X3)

# Set parameters
type = c(rep(1,n_fea*2),rep(0,n_fea))
param = c(rep(1,n_fea),trail[n_fea+1:(2*n_fea),1],rep(2,n_fea))

# Generate pathway information
resGraph<-vector(mode = "list",length = 6)
j=0
firstIdxMat<-which(W1[1:(n_fea/2),]!=0,arr.ind = T)
for (i in unique(firstIdxMat[,'col'])) {
  j=j+1
  tempIdx<-firstIdxMat[,1][firstIdxMat[,'col']==i]
  resGraph[[j]]<-combn(tempIdx,2)
}
firstIdxMat<-which(W1[(n_fea/2+1):(n_fea),]!=0,arr.ind = T)
for (i in unique(firstIdxMat[,'col'])) {
  j=j+1
  tempIdx<-firstIdxMat[,1][firstIdxMat[,'col']==i]
  resGraph[[j]]<-combn(tempIdx+n_fea/2,2)
}
secondIdxMat<-which(W2[1:(n_fea/2),]!=0,arr.ind = T)
for (i in unique(secondIdxMat[,'col'])) {
  j=j+1
  tempIdx<-secondIdxMat[,1][secondIdxMat[,'col']==i]
  resGraph[[j]]<-combn(tempIdx+n_fea,2)
}
secondIdxMat<-which(W2[(n_fea/2+1):n_fea,]!=0,arr.ind = T)
for (i in unique(secondIdxMat[,'col'])) {
  j=j+1
  tempIdx<-secondIdxMat[,1][secondIdxMat[,'col']==i]
  resGraph[[j]]<-combn(tempIdx+3/2*n_fea,2)
}
thirdIdxMat<-which(W3[1:(n_fea/2),]!=0,arr.ind = T)
for (i in unique(thirdIdxMat[,'col'])) {
  j=j+1
  tempIdx<-thirdIdxMat[,1][thirdIdxMat[,'col']==i]
  resGraph[[j]]<-combn(tempIdx+n_fea*2,2)
}
thirdIdxMat<-which(W3[(n_fea/2+1):n_fea,]!=0,arr.ind = T)
for (i in unique(thirdIdxMat[,'col'])) {
  j=j+1
  tempIdx<-thirdIdxMat[,1][thirdIdxMat[,'col']==i]
  resGraph[[j]]<-combn(tempIdx+n_fea*2.5,2)
}
resGraph<-resGraph[!sapply(resGraph,is.null)]
num_comb<-unique(sapply(resGraph, ncol))

if (edge_prop=="full"){
  E = t(bind_cols(resGraph))
  row.names(E)<-NULL
  E = rbind(E,E[,c(2,1)])
}else if (n_fea/2*as.numeric(edge_prop)>=num_comb) {
  warning("Number of edges is larger than fully connected graph, set teh graph as fully connected graph!")
  E = t(bind_cols(resGraph))
  row.names(E)<-NULL
  E = rbind(E,E[,c(2,1)])
}else if(n_fea/2*as.numeric(edge_prop)>0){
  resGraph = lapply(resGraph, function(x,num_comb_=num_comb,
                                       n_fea_=n_fea,
                                       edge_prop_=edge_prop){
    x[,sample(x=1:num_comb_,size = round(n_fea_/2*as.numeric(edge_prop_)),replace = F)]
  })
  E = t(bind_cols(resGraph))
  row.names(E)<-NULL
  E = rbind(E,E[,c(2,1)])
}else if (n_fea/2*as.numeric(edge_prop)==0){
  edge_prop<-as.numeric(edge_prop)
  E<-matrix(NA,nrow = edge_prop*nrow(W_true),ncol = 2) 
}
E<-E[order(E[,1],E[,2]),]

# This may take a while to run
autotuneRes<-SBFA_EM_AUTOT(X=X_simu,
                           type=type,
                           param=param,
                           E=E,
                           L=L,
                           v1=v1,
                           v2=v2,
                           a_omega=a_omega,
                           b_omega=b_omega,
                           maxIter = maxIter,
                           W.init = NULL)
