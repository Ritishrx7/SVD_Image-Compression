# Le Code de la décomposition SVD d'une image: Ritish et Andrea


# Libraries used
library(png)
library(base)
library(BBmisc)
library(stats)
library(SpatialPack)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Exemple de chargement d'un fichier PNG
img_mont<-readPNG("mont.png", info = TRUE)
img_mont


#On affiche l'image en plot
plot(img_mont, xlab="i", ylab = "j")

#L'image transformé a une matrice
img_matrix<-as.matrix(img_mont)
img_matrix


#On centre et réduit les valeurs 
norm_image <-BBmisc::normalize(img_matrix, method="standardize", range = c(0,1))
norm_image



#Plot de l'image normalisé
plot(norm_image,xlab="i", ylab = "j")





# On montre le conditionnement de la matrice




#SVD using the standardized matrix


# On decompose la matrice de l'image avec la SVD sous la forme S * SIGMA * T^T
svd_decomp<-svd(norm_image)
svd_decomp

#La Matrice S
svd_s<-svd_decomp$u
svd_s

#La Matrice diagonale Sigma des valeurs singulières
svd_sigma<-diag(svd_decomp$d)
svd_sigma

# La matrice T transpose
svd_t<-svd_decomp$v
svd_tp<-t(svd_t)
svd_tp

#La SVD multiplication
svd_matrix<-svd_s %*% svd_sigma %*% svd_tp
svd_matrix






#Equal test to check if svd matrix = norm_img matrix
Equal_test<-all.equal(norm_image,svd_matrix)
Equal_test





#Reconstruction de l'image avec la SVD: La fonction qui extrait le niveau de details de l'image avec standardization

svd_image<-function(n)
{
#LA decomposition SVD
image_final<-svd_decomp$u[,1:n] %*% diag(svd_decomp$d[1:n]) %*% t(svd_decomp$v[,1:n])

#L'image avec la matrice de SVD. On fait une rotation de 90`

image_final_rot<-t(image_final)[,nrow(image_final):1]

#Image en gray pour avoir les valeurs entre 0 et 1
return(image_sortie<-image(image_final_rot,col=grey(seq(0,1,length=180)),axes=FALSE))
}


# On affiche l'evolution de la decomposition SVD avec une grille de 3x3 avec 9 état de itérations

par(mfrow=c(3,3),mar=rep(4,4))
img_final_n <- sapply(c(2,20,30,50,90),svd_image)


#On change la grille d'affichage de graphique à défaut
dev.off()
par(mfrow=c(1,2))




#Reconstruction de l'image avec la SVD avec la standardization

svd_image(2)  # 110.3 kb
svd_image(20) # 147.3 kb
svd_image(30) # 150.7 kb
svd_image(50) # 150   kb
svd_image(90) # 151.3 kb



#Graphics

# On visualize les valeurs singulieres au carées sur la somme pour voir la variance expliquée de chaque valeur singulière

plot(svd_decomp$d^2/sum(svd_decomp$d^2),type="o",xlab="Sum of Square Sigma", ylab ="Contribution of the singular value", col="dark red")
axis(2,col="dark red")

par(mfrow=c(1,2))

# On créer la séconde graphique
plot(x=c(1:100),y=lapply(c(1:100),FUN = function(n) {sum(svd_decomp$d[1:n]^2)/sum(svd_decomp$d^2)}),type="o",col="dark blue",ylab ="Singular value contribution to image accuracy",xlab="Sum of Sigma squared",bty="o")




#Fonctions MSE et PNSR

#La somme des erreurs carrées de notre decomposition SVD
MSE<-function(f,m)
{

g<-svd_decomp$u[,1:m] %*% diag(svd_decomp$d[1:m]) %*% t(svd_decomp$v[,1:m])
mse<-(1/(nrow(f)*m))*sum(sum(f[nrow(f):m,nrow(f):m]-g[nrow(f):m,nrow(f):m])^2)
return(mse)
}

#Peak Signal to noise ratio
PNSR<-function(f,m)
{
  g<-svd_decomp$u[,1:m] %*% diag(svd_decomp$d[1:m]) %*% t(svd_decomp$v[,1:m])
  mse<-(1/(nrow(f)*m))*sum(sum(f[nrow(f):m,nrow(f):m]-g[nrow(f):m,nrow(f):m])^2)
  pnsr<-20*(log10(max(f[1:m,1:m])/sqrt(mse)))
  return(pnsr)
}

#MSE ET PNSR for rank n

MSE(norm_image,2)
MSE(norm_image,20)
MSE(norm_image,30)
MSE(norm_image,50)
MSE(norm_image,90)

PNSR(norm_image,2)
PNSR(norm_image,20)
PNSR(norm_image,30)
PNSR(norm_image,50)
PNSR(norm_image,90)



#Graphics of MSE and PNSR

#rang de la matrice SVD
Rank_M<-c(2:100)

#Graphique du Mean square error
MSE_norm <-sapply(Rank_M,MSE,f=norm_image)
plot(x=Rank_M,y=MSE_norm,type="l",col="dark blue",ylab ="Mean Square Error (MSE)",xlab="Rank nxp of reduced SVD matrix",bty="o")

#Graphique de Peak Signal to Noise Ratio
PNSR_norm <- sapply(Rank_M,PNSR,f=norm_image)
plot(x=Rank_M,y=PNSR_norm,type="l",col="red",ylab ="Peak Signal to Noise ratio (PSNR)",xlab="Rank nxp of reduced SVD matrix",bty="o",ylim =c(-45,85))






#L'image de la matrice original A*A^T et A^T*A (Image)

original1 <-(img_matrix %*% t(img_matrix))
image(original1,col=grey(seq(0,1,length=180)),axes=FALSE)

original2 <- (t(img_matrix) %*% img_matrix)
image(original,col=grey(seq(0,1,length=180)),axes=FALSE)





#IMAGE DECOMPOSITION without standardization


#SVD without standardization
svd_nodenoise <-svd(img_matrix)
svd_matrix_nodenoise<-svd_nodenoise$u %*% diag(svd_nodenoise$d) %*% t(svd_nodenoise$v)
sigma_nostd<-diag(svd_nodenoise$d)







#SVD function with no standardization

svd_image2<-function(n)
{
  #LA decomposition SVD
  image_final2<-svd_nodenoise$u[,1:n] %*% diag(svd_nodenoise$d[1:n]) %*% t(svd_nodenoise$v[, 1:n])
  
  #L'image avec la matrice de SVD. On fait une rotation de 90`
  
  image_final_rot2<-t(image_final2)[,nrow(image_final2):1]
  
  #Image en gray pour avoir les valeurs entre 0 et 1
  image_sortie2<-image(image_final_rot2,col=grey(seq(0,1,length=180)),axes=FALSE)
}


svd_image2(2)  # 116.2 kb
svd_image2(20) # 159.8 kb
svd_image2(30) # 161.3 kb
svd_image2(50) # 168.4 kb
svd_image2(90) # 170.4 kb




MSE2<-function(f,m)
{
  
  g<-svd_nodenoise$u[,1:m] %*% diag(svd_nodenoise$d[1:m]) %*% t(svd_nodenoise$v[,1:m])
  mse2<-(1/(nrow(f)*m))*sum(sum(f[nrow(f):m,nrow(f):m]-g[nrow(f):m,nrow(f):m])^2)
  return(mse2)
}

PNSR2<-function(f,m)
{
  g<-svd_nodenoise$u[,1:m] %*% diag(svd_nodenoise$d[1:m]) %*% t(svd_nodenoise$v[,1:m])
  mse2<-(1/(nrow(f)*m))*sum(sum(f[nrow(f):m,nrow(f):m]-g[nrow(f):m,nrow(f):m])^2)
  pnsr2<-20*(log10(max(f[1:m,1:m])/sqrt(mse2)))
  return(pnsr2)
}


Rank_M2<-c(2:100)
#Graphique du Mean square error
MSE_norm2 <-sapply(Rank_M2,MSE2,f=img_matrix)
plot(x=Rank_M2,y=MSE_norm2,type="l",col="dark blue",ylab ="Mean Square Error (MSE)",xlab="Rank nxp of reduced SVD matrix",bty="o")

#Graphique de Peak Signal to Noise Ratio
PNSR_norm2 <- sapply(Rank_M2,PNSR2,f=img_matrix)
plot(x=Rank_M2,y=PNSR_norm2,type="l",col="red",ylab ="Peak Signal to Noise ratio (PSNR)",xlab="Rank nxp of reduced SVD matrix",bty="o", ylim =c(-45,85))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#


#CONDITION NUMBER OF THE MATRIX

#Avant la normalisation
cond_num<-kappa(img_mont, exact = TRUE)
cond_num

#Apres la normalisation
cond_num_norm<-kappa(norm_image, exact = TRUE)
cond_num_norm


dev.off()
#Condition number function of matrix with rank 2:n
condition_number <- function(f,k)
{
  return(kappa(f[1:k,1:k],exact = TRUE))
  
}
cond_range <- c(1:100)




#Condition number of matrix with rank 2:n with image without normalised matrix

kappa_range <- sapply(cond_range,condition_number,f=sigma_nostd)
plot(x=cond_range,y=kappa_range,type="l",col="dark green",ylab ="Condition number SVD without Standardization",xlab="Condition number of SVD rank nxp",bty="o")

#Condition number of matrix with rank 2:n with image of the SVD decomposition matrix with normalised matrix

kappa_range <- sapply(cond_range,condition_number,f=svd_sigma)
plot(x=cond_range,y=kappa_range,type="l",col="dark green",ylab ="Condition number SVD with Standardization",xlab="Condition number of SVD rank nxp",bty="o")

#Scenario analysis of the SVD standardized matrix of the condition number

#Condition number of n singular values (Standardized)

kappa(svd_sigma[1:2,1:2], exact = TRUE)
kappa(svd_sigma[1:20,1:20], exact = TRUE)
kappa(svd_sigma[1:30,1:30], exact = TRUE)
kappa(svd_sigma[1:50,1:50], exact = TRUE)
kappa(svd_sigma[1:90,1:90], exact = TRUE)

#Condition number of n singular values (Not Standardized)

kappa(sigma_nostd[1:2,1:2], exact = TRUE)
kappa(sigma_nostd[1:20,1:20], exact = TRUE)
kappa(sigma_nostd[1:30,1:30], exact = TRUE)
kappa(sigma_nostd[1:50,1:50], exact = TRUE)
kappa(sigma_nostd[1:90,1:90], exact = TRUE)



#Condition number of AT*A using svd singular values
kappa(t(svd_sigma[1:2,1:2])%*%svd_sigma[1:2,1:2], exact = TRUE)
kappa(t(svd_sigma[1:20,1:20])%*%svd_sigma[1:20,1:20], exact = TRUE)
kappa(t(svd_sigma[1:30,1:30])%*%svd_sigma[1:30,1:30], exact = TRUE)
kappa(t(svd_sigma[1:50,1:50])%*%svd_sigma[1:50,1:50], exact = TRUE)
kappa(t(svd_sigma[1:90,1:90])%*%svd_sigma[1:90,1:90], exact = TRUE)

#Scenario Analysis
#Standardized Image matrix with SVD

MSE(norm_image,2)
MSE(norm_image,20)
MSE(norm_image,30)
MSE(norm_image,50)
MSE(norm_image,90)


PNSR(norm_image,2)
PNSR(norm_image,20)
PNSR(norm_image,30)
PNSR(norm_image,50)
PNSR(norm_image,90)


#Original Image matrix with SVD

MSE2(img_matrix,2)
MSE2(img_matrix,20)
MSE2(img_matrix,30)
MSE2(img_matrix,50)
MSE2(img_matrix,90)


PNSR2(img_matrix,2)
PNSR2(img_matrix,20)
PNSR2(img_matrix,30)
PNSR2(img_matrix,50)
PNSR2(img_matrix,90)



#%%%%%%%%%%%%%%%%%%%%%%%%%%COMPUTATION TIME OF SVD%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#Not standardized SVD

system.time(svd_image2(2))
system.time(svd_image2(20))
system.time(svd_image2(30))
system.time(svd_image2(50))
system.time(svd_image2(90))

#Not standardized SVD

system.time(svd_image(2))
system.time(svd_image(20))
system.time(svd_image(30))
system.time(svd_image(50))
system.time(svd_image(90))


