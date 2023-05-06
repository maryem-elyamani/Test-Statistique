library(ggplot2)
library(rmutil)
library(matlib)
library(splines)




#n = 50
#num_x = 1000

#d=1
#k=2

#num_splines = k+d+1


# x <- seq(from=0, to=1, length.out=num_x)
# 
# X=c(0.645612, 0.66022, 0.555223, 0.368246, 0.801863, 0.188038, 0.837897, 
#     0.872886, 0.179029, 0.312767, 0.184966, 0.782871, 0.196528, 0.111166, 
#     0.950258, 0.0409111, 0.0230874, 0.737858, 0.381488, 0.693835, 
#     0.271724, 0.779804, 0.624836, 0.0160502, 0.831894, 0.850674, 
#     0.624769, 0.214307, 0.590515, 0.273547, 0.716782, 0.510014, 0.608093, 
#     0.342415, 0.458901, 0.969161, 0.0731037, 0.445583, 0.939654, 
#     0.918209, 0.186333, 0.29723, 0.741534, 0.802549, 0.71817, 0.848047, 
#     0.430392, 0.26375, 0.697366, 0.165481)
# 
# Y=c(0.65977, 0.0316515, 0.872726, 0.944487, 0.595143, 0.770168, 
#     0.923727, 0.256481, 0.372617, 0.430102, 0.721855, 0.244513, 0.635611, 
#     0.651448, 0.52978, 0.516194, 0.351899, 0.360599, 0.800771, 0.533902, 
#     0.342824, 0.14494, 0.442621, 0.903258, 0.508918, 0.70013, 0.941215, 
#     0.0815881, 0.971915, 0.918777, 0.43683, 0.814787, 0.711265, 
#     0.0630177, 0.209283, 0.212617, 0.532894, 0.734362, 0.777764, 
#     0.162219, 0.428761, 0.430543, 0.571107, 0.568565, 0.666056, 0.96149, 
#     0.886778, 0.00261489, 0.541221, 0.717969)




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


