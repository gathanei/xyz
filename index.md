---
layout: dis_template
---

The xyz package implements the [xyz](https://arxiv.org/abs/1610.05108) algorithm for fast interaction search in high-dimensional data.

# A simple problem

Given a data matrix \\( X \in \mathbb{R}^{n\times p} \\) and response vector \\( Y \in \mathbb{R}^n \\). As a simple example we want to fit the following model:

$$Y_i = \beta_{jk} X_{ij}X_{ik}+\varepsilon_i$$

The problem is we don't know the exact interaction pair \\( (j,k) \\) . If we would loop through all possible pairs, we get a quadratic runtime \\(\mathcal{O}(np^2)\\). The xyz algorithm provably returns the correct interaction pair in subquadratic time, that is \\( \mathcal{O}(np^{\alpha}) \\) with \\( \alpha < 2 \\).

# How to install

You can install the package either by compiling it yourself or by using the devtools R library.
{% highlight R %}
library(devtools)
install_github("gathanei/xyz")
{% endhighlight %}

# How to use

We generate the above model.

{% highlight R %}
library(xyz)
n<-300
p<-1000
#construct a binary matrix
X<-matrix(sample(c(-1,1),replace=TRUE,n*p),n,p)
#set an interaction of the pair (1,2)
Y<-X[,1]*X[,2]+rnorm(n)
#run the interaction search
result<-xyz_search(X,Y,L=10,N=10,binary=TRUE,negative=TRUE)
#print the result
print(result)
{% endhighlight %}
