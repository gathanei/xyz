# xyz algorithm

R package of the xyz algorithm for fast interaction search in high-dimensional data.

## what is it?

Finding interactions in high-dimensional data can be computationally very expensive. If the data set
has p variables then naive search incurs a quadratic cost in p.

The xyz algorithm finds strong interaction in almost linear time.

### example

A simple example in R:
```
library(xyz)
n<-300
p<-1000
#construct a binary matrix
X<-matrix(sample(c(-1,1),replace=TRUE,n*p),n,p)
#set an interaction of the pair (1,2)
Y<-X[,1]*X[,2]
#run the interaction search
result<-xyz_search(X,Y,L=10,N=10,binary=TRUE,negative=TRUE)
#print the result
print(result)
```

### how to install

This package is only on github and not yet on CRAN. To install run the following code:

```
library(devtools)
install_github("gathanei/xyz")
```

### reference

[xyz-paper](https://arxiv.org/abs/1610.05108)
