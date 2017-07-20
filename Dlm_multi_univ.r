
dyn.load("~/Documentos/Analisis_grupo/multivariate_DLM/multiva_dlm.so")
MultiDlm <- function(Vm0, Vc0, Npar, Vbeta, Y, coord, Nneigh, Stop, X, n01, Cbeta, VecProb, conditional){
  if(Stop<=dim(X)[1] & Stop>0){
    tep = .C("multiDlm", as.double(Vm0), as.double(Vc0), as.integer(Npar), as.double(Vbeta), as.double(t(Y)), 
             as.double(t(coord)), dim(Y), dim(coord), as.integer(Nneigh), as.integer(Stop), 
             as.double(t(X)), dim(X), as.double(n01), as.double(Cbeta), res1=as.double(t(VecProb)), 
             as.integer(conditional))
    res11 = matrix(tep$res1, byrow = TRUE, ncol=dim(Y)[2])
    #res11 = matrix(tep$res12, byrow = TRUE, ncol=1)
    return(significance=res11)}else{return(c("Error in Stop: subscript out of bounds"))}
}




Input.Data <- function(mask1, data1){
  
  dime.data <- dim(data1)
  pos1 <- which(mask1==1, arr.ind = TRUE)
  vol.input.c <- sapply(1:dim(pos1)[1], function(i){data1[pos1[i,1], pos1[i,2], pos1[i,3],1:dime.data[4]]})
  
  return(list(vol.imput.c = vol.input.c, posi = pos1))
  
}
