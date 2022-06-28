SBFA_EM_AUTOT<-function(X,type,param,E,L,v1,v2,a_omega,b_omega,m.init=1,scale=T,W.init=NULL,eps=1e-6,maxIter=500){
  len_v1<-length(v1)
  len_v2<-length(v2)
  len_L<-length(L)
  len_a_omega<-length(a_omega)
  len_b_omega<-length(b_omega)
  totalIter<-len_v1*len_v2*len_L*len_a_omega*len_b_omega
  print(paste0("Total number of iteration is ",totalIter))
  
  result_Mat<-data.frame(matrix(NA,nrow = 6,ncol = totalIter))
  row.names(result_Mat)<-c("L","a_omega","b_omega","v1","v2","BICs")
  colnames(result_Mat)<-1:totalIter
  
  curBIC<-Inf
  cur_iter<-0
  for (idx_L in 1:len_L) {
    for (idx_a_omega in 1:len_a_omega) {
      for (idx_b_omega in 1:len_b_omega) {
        for (idx_v1 in 1:len_v1) {
          for (idx_v2 in 1:len_v2) {
            cur_iter<-cur_iter+1
            print(paste0("In iteration: ",cur_iter))
            tryCatch({
              temp_GBFA<-SBFA_EM(X,
                                 type,
                                 param,
                                 E,
                                 L[idx_L],
                                 v1[idx_v1],
                                 v2[idx_v2],
                                 a_omega[idx_a_omega],
                                 b_omega[idx_b_omega],
                                 m.init=m.init,
                                 scale=scale,
                                 W.init=W.init,
                                 eps=eps,
                                 maxIter=maxIter)
              result_Mat[1,cur_iter]<-L[idx_L]
              result_Mat[2,cur_iter]<-a_omega[idx_a_omega]
              result_Mat[3,cur_iter]<-b_omega[idx_b_omega]
              result_Mat[4,cur_iter]<-v1[idx_v1]
              result_Mat[5,cur_iter]<-v2[idx_v2]
              result_Mat[6,cur_iter]<-temp_GBFA$BIC
              
              #saveRDS(temp_GBFA,file = paste0(savePath,"_L",L[idx_L],"_v1",v1[idx_v1],"_v2",v2[idx_v2],".rds") )
              
              if(temp_GBFA$BIC<=curBIC){
                curBIC<-temp_GBFA$BIC
                bestModel<-temp_GBFA
                bestL<-L[idx_L]
                bestv1<-v1[idx_v1]
                bestv2<-v2[idx_v2]
                best_a_omega<-a_omega[idx_a_omega]
                best_b_omega<-b_omega[idx_b_omega]
                bestIter<-cur_iter
                print(paste0("Best model (BIC) is detected at ",cur_iter,"/",totalIter))
              }
            },error=function(msg){
              print("Iteration abort!")
              print(msg)
            })
            
          }
        }
      }
    }
  }
  return(list(tuningProcess=result_Mat,bestIter=bestIter,bestBIC=curBIC,bestModel=bestModel,bestL=bestL,bestv1=bestv1,bestv2=bestv2,best_a_omega=best_a_omega,best_b_omega=best_b_omega))
}