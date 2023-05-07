library(ggplot2)
library(rmutil)
library(matlib)
library(splines)




# returns a vector equidistal distributed numbers from 0 to 1
# corresponds to x in script
equidist_points<-function(n){
  seq(from=0, to=1, length.out=n)
}


num_splines<-function(d,k){
  d+k+1
}

# calculates the base spline functions given the random data X and
# returns b-spline-functions given the parameters d and k
B_Splines<-function(X,d,k){
  N_X <- bs(X, df = num_splines(d,k), degree = d, intercept = TRUE)
  N_X
}


# calculates M hat as described in the script, for given samples X and
# parameters d and k
M_hat<-function(X,d,k,n){
  b_spls = B_Splines(X,d,k)
  M_x = matrix(0,nrow = num_splines(d,k), ncol = num_splines(d,k)) 
  #n = length(X)
  
  for (i in (1:n)){
    M_x = M_x + b_spls[i,] %*% t(b_spls[i,]) 
  }
  M_x<-M_x/n
  M_x
}


# Calculates the test satistic D as described in the script for the two sample
# vectors X and Y. num_u refers to the number of equidistant points in vector
# x (see script)
D<- function(X,Y,num_e=1000,d=1,k=2) {
  
  n = length(X)
  
  # vector just containing ones
  vec_1<-matrix(1,1,num_splines(d,k))
  
  # equidistal points
  x <- equidist_points(num_e)
  
  # approximations of the given M matrixes
  M <- M_hat(x,d,k,n)
  Mh_X<-M_hat(X,d,k,n)
  Mh_Y<-M_hat(Y,d,k,n)
  
  # calculation of the statistic
  D_XY<-n*vec_1%*%(Mh_X-Mh_Y)%*%inv(Mh_Y)%*%M%*%inv(Mh_Y)%*%(Mh_X-Mh_Y)%*%t(vec_1)
  D_XY
}

# creates M random samples of of D of size n, in the case, where H_0 is true 
# returns samples, and scth quantiles for each iteration step 
D_dist_H0 <-function(M,n, q=.95){
  if (M >= 1000){print("Generating data, this might take some time ...")}
  
  Dobs <- numeric(M) # initialiser le vecteur Dobs
  sc_s <- numeric(M)
  
  for (j in 1:M) 
  {
    # générer les données X et Y
    X <- runif(n, 0, 1)
    Y <- runif(n, 0, 1)
    # calculer la statistique du test D pour ces données
    Dobs[j] <-D(X, Y)
    sc_s[j] <- quantile(Dobs[1:j], q)
    #if (j %% 1000 == 0){print(paste("M:",j,"sc:",sc_s[j]))}
  }
  lst = list(sc_s=sc_s,Dobs=Dobs)
  lst
}


