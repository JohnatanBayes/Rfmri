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




data.input.c1 <- Input.Data(mask1 = mask.c1, data1 = ffd.c)
y.data.input.c1 <- data.input.c1$vol.imput.c
y.data.input.c11 <- apply(y.data.input.c1, 2, function(x)(x-mean(x))/sd(x))
dim(y.data.input.c11)

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



system.time(Mult_dlm  <- MultiDlm(Vm0 = 0, Vc0 = 100, Npar = 2, Vbeta = sqrt(1/0.95), Y = y.data.input.c11, 
                      coord = t(posi.data.input.c1), corr1 = 0.80, cluster = 5, Nneigh = 27, 
                      Stop = 90, X = as.matrix(covaria), n01 = 1, Cbeta = 0.90, 
                      VecProb = matrix(rep(0,2*dim(y.data.input.c11)[2]), ncol=dim(y.data.input.c11)[2]), Evidencia = 0.01))




ffd.dlm.c = array(0, c(2,dim(ffd.c)[1:3]))
#dim(ffd.dlm.c)
for(j in 1:dim(Mult_dlm)[2]){
  ffd.dlm.c[1:2, posi.data.input.c1[j, 1], posi.data.input.c1[j, 2], posi.data.input.c1[j, 3]]= Mult_dlm[1:2,j]
}

dim(ffd.dlm.c)
max(ffd.dlm.c)
min(ffd.dlm.c)


test.final <- function(x){
  if(x[2]>0&x[2]!=10000){
    res <- pnorm(0, x[1], sqrt(x[2]), lower.tail = FALSE)
    res2 <- ifelse(res>0.99, 4,1)
  }else res2 <- 1
  return(res2)
}


final.vol <- apply(ffd.dlm.c, 2:4, test.final)
dim(final.vol)




ffd.ref <- readNIfTI("~/Documentos/Analisis_grupo/Controles/C000393_BOLD_fMRS_SENSE_FSL_10_1.feat/reg/standard.nii.gz")
mask.ref <- ifelse(ffd.ref!=0, 1, 0)
Z.visual.c <- nifti(final.vol*mask.ref, datatype=16)
#Z.visual <- test11
yrange.c <- c(3.9, max(Z.visual.c, na.rm=TRUE))
#setEPS()
#postscript(paste("Animation/plot",j,".eps", sep=""))
overlay(ffd.ref, ifelse(Z.visual.c>3.9, Z.visual.c, NA), zlim.x=range(ffd.ref)*0.95, zlim.y=yrange.c, xlab="1", axes=TRUE, col.y = "red")

