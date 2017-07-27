context("Gets previous results")

test_that("Gives results from vignette", {
  #####
  # Example 1
  set.seed(18990833)
  n <- 100;p<-1000
  X <- matrix(2*rbinom(n*p,1,0.5)-1,n,p)
  Y <- X[,1]*X[,2]
  result <- xyz_search(X,Y,L=5,N=10,binary=TRUE,negative=TRUE)

  # save_to_test(result, "p1e3_ex")
  expect_equal(result, read_to_test("p1e3_ex"), tolerance = 1.49e-08)

  #####
  # Example 2
  set.seed(1337876)

  n <- 100;p<-1000;
  X <- matrix(rnorm(n*p),n,p)
  Y <- 3*X[,5]+2*X[,1]*X[,2]-3*X[,7]*X[,4]+rnorm(n)

  result <- xyz_regression(X,Y,L=10,n_lambda=10,alpha=0.9)

  # save_to_test(result, "p1e4_ex")
  expect_equal(result, read_to_test("p1e4_ex"), tolerance = 1.49e-08)
})
