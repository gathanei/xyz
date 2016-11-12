---
layout: dis_template
---

The xyz package implements the [xyz](https://arxiv.org/abs/1610.05108) algorithm for fast interaction search in high-dimensional data.

# A simple problem

Given a data matrix \\( X \in \mathbb{R}^{n\times p} \\) and response vector \\( Y \in \mathbb{R}^n \\). As a simple example we want to fit the following model:

$$Y_i = \beta_{jk} X_{ij}X_{ik}+\varepsilon_i.$$

The problem is we don't know the exact interaction pair \\( (j,k) \\) . If we would loop through all possible pairs, we get a quadratic runtime \\(\mathcal{O}(np^2)\\). The xyz algorithm provably returns the correct interaction pair in subquadratic time, that is \\( \mathcal{O}(np^{\alpha}) \\) with \\( \alpha < 2 \\). More elaborate models can be considered, see package vignette for details.

# How to install

You can either install the package directly from CRAN (using ```install.packages("xyz") ```) or you can get
the most up to date version by installing from github:
{% highlight R %}
library(devtools)
install_github("gathanei/xyz")
{% endhighlight %}

# How to use

We generate the model \\( Y_i=2 X_{i1}X_{i2}+\varepsilon_i \\), where \\(X \in \\{-1,1\\}^{n \times p} \\) and \\( \varepsilon_i \\) is iid Gaussian noise.

{% highlight R %}
library(xyz)
n<-300
p<-1000
#construct a data matrix X
X<-matrix(sample(c(-1,1),replace=T,n*p),n,p)
#interaction of the pair (1,2)
Y<-2*X[,1]*X[,2]+rnorm(n)
#run the interaction search
result<-xyz_search(X,Y,L=10,N=10,binary=T,negative=T)
#print the result
print(result)
{% endhighlight %}

For a practical verification of the subquadratic runtime see the plot below.

{: .center}
![asymptotic_runtime]({{ site.baseurl }}/assets/asymptotic_runtime.png)

*It depicts three runs of xyz on a data set with different interaction strengths. Green is the weakest interaction
and purple the strongest. Recorded is the time taking to discover the correct interaction.*

The package also implements interaction search for high-dimensional regression (see the function `xyz_regression`).

# References

[The xyz algorithm for fast interaction search in high-dimensional data. (Gian-Andrea Thanei, Nicolai Meinshausen, Rajen D. Shah)](https://arxiv.org/abs/1610.05108)
