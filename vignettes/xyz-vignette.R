## ---- echo=FALSE---------------------------------------------------------
source('vignettecode.R')
set.seed(1)

## ---- echo=TRUE----------------------------------------------------------
#set dimensions
n<-100;p<-1000;
#create data
X<-matrix(2*rbinom(n*p,1,0.5)-1,n,p)
Y<-X[,1]*X[,2]

#find top 10 interactions
result<-xyz_search(X,Y,L=5,N=10,binary=TRUE,negative=TRUE)
#the first element contains the interaction pairs the second element contains their strength
print(result)

## ---- echo=TRUE,fig.show='hold',fig.width=6, fig.height=6----------------
#set dimensions
n<-100;p<-1000;
#create data
X<-matrix(rnorm(n*p),n,p)

#build model
Y<-3*X[,5]+2*X[,1]*X[,2]-3*X[,7]*X[,4]+rnorm(n)

#find top 10 interactions
result<-xyz_regression(X,Y,L=10,n_lambda=10,alpha=0.9)
#the first element contains the main effects and the third the interaction effects, we look at the fifth lambda
print(result)
plot(result)
#predict
predict(result,matrix(rnorm(10*p),10,p))

