library(oro.nifti)
source("~/Documentos/Analisis_grupo/multivariate_DLM/Dlm_multi_univ.r")
#MultiDlm(): multivariate dlm
#Input.Data(): a function that allows to organize the input data in a proper way
##########
# COTROL #
##########
#Covariates
covaria <- read.csv("~/Documentos/C_codes/covariables.csv", header=FALSE, sep="")
str(covaria)

#DATA FROM CONTROL 1
ffd.c <- readNIfTI("~/Documentos/Analisis_grupo/Controles/C000393_BOLD_fMRS_SENSE_FSL_10_1.feat/std.nii.gz")@.Data
mask.c1 <- ifelse(ffd.c[,,,2]!=0, 1, 0)
res <- sapply(1:dim(mask.c1)[3], function(i)sum(mask.c1[,,i]))
posi.no.signal <- which(res!=0)

new.vol <- ffd.c[,,posi.no.signal,]
dim(new.vol)
new.mask <- mask.c1[,,posi.no.signal]
dim(new.mask)


ffd.dlm.f = array(0, c(2,dim(new.vol)[1:3]))
dim(ffd.dlm.f)

#ORGANIZE THE TIME SERIES AND POSITIONS (INSIDE THE VOL) IN MATRICES
##########
#CONTROL##
##########
system.time(
for(i in 1:dim(new.vol)[3]){
  
  if(i!=dim(new.vol)[3]){
    posi <- (i-1):(i+1)
    posi1 <- posi[posi>0]
  }else{posi1 <- (dim(new.vol)[3]-1):dim(new.vol)[3]}
  
  
  
  Aux.new.vol <- new.vol[,,posi1,]
  Aux.new.mask <- new.mask[,,posi1]
  
  
  data.input.c1 <- Input.Data(mask1 = Aux.new.mask, data1 = Aux.new.vol)
  y.data.input.c1 <- data.input.c1$vol.imput.c
  y.data.input.c11 <- apply(y.data.input.c1, 2, function(x)(x-mean(x))/sd(x))
  
  
  posi.data.input.c1 <- data.input.c1$posi
  dim(posi.data.input.c1)
  #rm(list=c("ffd.c1", "ffd.p1"))
  
  #Vm0: prior mean value
  #Vc0: prior variance value
  #Npar: number of regressors
  #Vbeta: Value of the descount factor
  #Y: matrix containing the voxel's series
  #coord: matrix containing the voxel position into the brain array
  #cluster: The minimum number of voxels into the neighborhood
  #Nneigh: The number of neighbors into the prefixed neighboring geometry: in the case of cube shape, the number will be 27
  #Stop: Stop point or last period of time analyzed with the DLM
  #X: Desing matrix
  #n01: prior hyperparameter
  #Cbeta: prior hyperparameter
  #VecProb: output vector containing the final result associated with the FBST test
  #Evidence: prefixed evidence value for the FBST test
  
  
  
  system.time( Mult_dlm  <- MultiDlm(Vm0 = 0, Vc0 = 100, Npar = 2, Vbeta = sqrt(1/0.95), Y = y.data.input.c11, 
                                    coord = t(posi.data.input.c1), corr1 = 0.75, cluster = 7, Nneigh = 27, 
                                    Stop = 90, X = as.matrix(covaria), n01 = 1, Cbeta = 0.90, 
                                    VecProb = matrix(rep(0,2*dim(y.data.input.c11)[2]), ncol=dim(y.data.input.c11)[2]), Evidencia = 0.01)
  )
  
  

  ffd.dlm.c = array(0, c(2,dim(Aux.new.vol)[1:3]))
  #dim(ffd.dlm.c)
  for(j in 1:dim(Mult_dlm)[2]){
    ffd.dlm.c[1:2, posi.data.input.c1[j, 1], posi.data.input.c1[j, 2], posi.data.input.c1[j, 3]]= Mult_dlm[1:2,j]
  }
  
  
  
  ffd.dlm.f[,,,i]<- ffd.dlm.c[,,,which(posi1==i)]
  
  
  
})
 


dim(ffd.dlm.f)
ffd.dlm.f2 = array(0, c(2,dim(ffd.c)[1:3]))
ffd.dlm.f2[ , , , posi.no.signal] <- ffd.dlm.f
dim(ffd.dlm.f2)



test.final <- function(x){
  if(x[2]>0&x[2]!=10000){
    res <- pnorm(0, x[1], sqrt(x[2]), lower.tail = FALSE)
    res2 <- ifelse(res>0.99, 4,1)
  }else res2 <- 1
  return(res2)
}


final.vol <- apply(ffd.dlm.f2, 2:4, test.final)
dim(final.vol)




ffd.ref <- readNIfTI("~/Documentos/Analisis_grupo/Controles/C000393_BOLD_fMRS_SENSE_FSL_10_1.feat/reg/standard.nii.gz")
mask.ref <- ifelse(ffd.ref!=0, 1, 0)
Z.visual.c <- nifti(final.vol*mask.ref, datatype=16)
#Z.visual <- test11
yrange.c <- c(3.9, max(Z.visual.c, na.rm=TRUE))
#setEPS()
#postscript(paste("Animation/plot",j,".eps", sep=""))
overlay(ffd.ref, ifelse(Z.visual.c>3.9, Z.visual.c, NA), zlim.x=range(ffd.ref)*0.95, zlim.y=yrange.c, xlab="1", axes=TRUE, col.y = "red")







##################
#univariate case##
##################
ffd.uni <- readNIfTI("Documentos/Analisis_grupo/grupos_out/control1.nii.gz")@.Data
dim(ffd.uni)

test.final <- function(x){
  if(x[2]>0){
    res <- pnorm(0, x[1], sqrt(x[2]), lower.tail = FALSE)
    res2 <- ifelse(res>0.99, 4,1)
  }else res2 <- 1
  return(res2)
}


Map.uni <- apply(ffd.uni, 2:4, test.final)
dim(Map.uni)

ffd.ref <- readNIfTI("~/Documentos/Analisis_grupo/Controles/C000393_BOLD_fMRS_SENSE_FSL_10_1.feat/reg/standard.nii.gz")
mask.uni <- ifelse(ffd.ref!=0, 1, 0)
Z.visual.c <- nifti(Map.uni*mask.uni, datatype=16)
#Z.visual <- test11
yrange.c <- c(3.9, max(Z.visual.c, na.rm=TRUE))
#setEPS()
#postscript(paste("Animation/plot",j,".eps", sep=""))
overlay(ffd.ref, ifelse(Z.visual.c>3.9, Z.visual.c, NA), zlim.x=range(ffd.ref)*0.95, zlim.y=yrange.c, xlab="1", axes=TRUE, col.y = "red")


