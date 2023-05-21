library(ggplot2)
library(rmutil)
library(matlib)
library(splines)


# Définir les paramètres

n = 50
num_x = 1000

d=1
k=2

alpha_x <- 1
beta_x <- 1
alpha_y <- 1
beta_y <- 1

x <- seq(from=0, to=1, length.out=num_x)
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
  #M = matrix(0,nrow = num_splines, ncol = num_splines)
  #for (i in (1:num_splines)){
  #  for (j in (1:num_splines)){
  #    M[i,j] = sum(Base_Spline(X)[,i]*Base_Spline(X)[,j])*(max(Base_Spline(X))-min(Base_Spline(X)))
  #  }
  #  M <- M/n
  #}

  
  Vecteur_1<-matrix(1,1,num_splines)
  vect_n<-n*Vecteur_1
  vect_n
  M <- M_hat(x,num_splines,n)
  Mh_X<-M_hat(X,num_splines,n)
  Mh_Y<-M_hat(Y,num_splines,n)
  D<-vect_n%*%(Mh_X-Mh_Y)%*%inv(Mh_Y)%*%M%*%inv(Mh_Y)%*%(Mh_X-Mh_Y)%*%t(Vecteur_1)
  D
}

D(X,Y,n)





M = 10000
n = 50

Dobs <- numeric(M) # initialiser le vecteur Dobs
sc_s <- numeric(M)

for (j in 1:M) 
  {
  # générer les données X et Y
  X <- runif(n, 0, 1)
  Y <- runif(n, 0, 1)
  # calculer la statistique du test D pour ces données
  Dobs[j] <-D(X, Y, n)
  sc_s[j] <- quantile(Dobs[1:j], 0.95)
  if (j %% 50 == 0){print(paste("M:",j,"sc:",sc_s[j]))}
}
Dobs
print("finished")
plot(sc_s, type="l")
print("finished: ")
print(paste(" quantile M=1000 ", sc_s[1000]))
print(paste(" quantile M=5000 ", sc_s[5000]))
print(paste(" quantile M=10000 ", sc_s[10000]))
ggplot(data = data.frame(sc_s[100:M])) +
  geom_line(aes(x = seq_along(sc_s[100:M]), y = sc_s[100:M], color = "Threshold value")) +
  xlab("Iteration") +
  ylab("Threshold value (sc)") +
  ggtitle("Threshold value (sc) over iterations") +
  scale_color_manual(values = c("Threshold value" = "magenta"))



