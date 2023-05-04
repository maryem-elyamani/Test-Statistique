library(ggplot2)
library(rmutil)
library(matlib)
library(splines)


# Définir les paramètres

n = 50
num_x = 50
num_y = 50

d=1
k=2

alpha_x <- 1
beta_x <- 1
alpha_y <- 1
beta_y <- 1

x <- seq(from=0, to=1, length.out=num_x)
y <- seq(from=0, to=1, length.out=num_y)
num_splines = k+d+1
X <- sort(rbeta(n, alpha_x, beta_x))
Y <- sort(rbeta(n, alpha_y, beta_y))


Base_Spline<-function(X){
  N_X <- bs(X, df = num_splines, degree = d, intercept = TRUE)
  N_X
}
print(dim(Base_Spline(X)))

M_hat<-function(X,num_splines,n){
 
  M_x = matrix(0,nrow = num_splines, ncol = num_splines) 
  for (i in (1:n)){
    M_x = M_x + Base_Spline(X)[i,] %*% t(Base_Spline(X)[i,]) 
  }
  M_x<-M_x/n
}




D<- function(X,Y,n) {
  
  # la matrice M
  M = matrix(0,nrow = num_splines, ncol = num_splines)
  for (i in (1:num_splines)){
    for (j in (1:num_splines)){
      M[i,j] = sum(Base_Spline(X)[,i]*Base_Spline(X)[,j])*(max(Base_Spline(X))-min(Base_Spline(X)))
    }
    M <- M/n
  }

  
  Vecteur_1<-matrix(1,1,num_splines)
  vect_n<-n*Vecteur_1
  vect_n
  Mh_X<-M_hat(X,num_splines,n)
  Mh_Y<-M_hat(Y,num_splines,n)
  D<-vect_n%*%(Mh_X-Mh_Y)%*%inv(Mh_Y)%*%M%*%inv(Mh_Y)%*%(Mh_X-Mh_Y)%*%t(Vecteur_1)
  D
}

D(X,Y,n)



set.seed(123) # pour la reproductibilité

M = 1000
n = 50
num_x = num_y = n
alpha_x = alpha_y = 1
beta_x = beta_y = 1
k = 2
d = 1

Dobs <- numeric(M) # initialiser le vecteur Dobs

for (j in 1:M) {
  # générer les données X et Y
  X <- runif(num_x, 0, 1)
  Y <- runif(num_y, 0, 1)
  # calculer la statistique du test D pour ces données
  Dobs[j] <-D(X, Y, n)
}

# afficher les premières valeurs de Dobs
Dobs
