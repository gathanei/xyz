#' Interaction search
#' @param X A matrix.
#' @param Y A vector.
#' @param L An integer indicating how many projection steps are performed.
#' @param N A integer, controlling the number of pairs that will be returned in the end.
#' @param binary  A logical indicating if X is binary or continuous.
#' @param negative A logical indicating if also negative interactions should be searched for.
#' @return \code{N} strongest interactions between \code{X} and \code{Y} after \code{L} projections.
#' @examples
#' n<-300
#' p<-1000
#' #construct a binary matrix
#' X<-matrix(sample(c(-1,1),replace=TRUE,n*p),n,p)
#' #set an interaction of the pair (1,2)
#' Y<-X[,1]*X[,2]
#' #run the interaction search
#' result<-xyz_search(X,Y,L=10,N=10,binary=TRUE,negative=TRUE)
#' #print the result
#' print(result)
#' @export
xyz_search<-function(X,Y,L=10,N=100,binary=TRUE,negative=TRUE) {
  L<-round(L)
  if(L < 1) {
    stop("Number of runs has to be at least 1.")
  }
  s<-30
  #do translating checks here
  if(!is.matrix(X)) {
    X<-as.matrix(X)
  }
  n<-dim(X)[1]
  p<-dim(X)[2]
  if(!is.vector(Y)) {
    stop("Y has to be a vector.")
  }
  if(length(Y) != n) {
    stop("Y and X have to have the same number of rows")
  }
  if(n < 10) {
    stop(paste("You have ",n," samples. The number of samples should at least be 10.",sep=""))
  }
  result<-list()
  if(binary) {
    if(!((sum(X==1)+sum(X==-1))==prod(dim(X)))) {
      stop("X is not binary.")
    }
  }
  result<-interaction_search(X, Y, L, N, negative,binary)
  result[[1]]<-result[[1]]+1
  class(result)<-"xyz_search_result"
  return(result)

  stop("You reached the end of the function and it won't return anything. This is not good.")
}

#' @export
summary.xyz_search_result<-function(object,...,maxsum=10) {
  if(length(object[[2]])< 1) {
    output<-"no interactions discovered"
    if(length(object[[2]])==1) {
      cat("intereaction pair: (",object[[1]][1,1],",",object[[1]][2,1],") strength: ",object[[2]][1],sep="")
    }
  }
  l<-min(maxsum,length(object[[2]]))
  if(l > 1) {
    for(i in 1:l) {
      cat("intereaction pair: (",object[[1]][1,i],",",object[[1]][2,i],") strength: ",object[[2]][i],"\n",sep="")
    }
  }
}

#' @export
print.xyz_search_result<-function(x,...) {
  return(summary(x,maxsum=length(x[[2]])))
}

#' @importFrom graphics plot
#' @export
plot.xyz_search_result<-function(x,...) {
  return("no plot functionality implemented")
}

