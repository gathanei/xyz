#include <Rcpp.h>
#include <numeric>
#include <utility>
//[[Rcpp::plugins(cpp11)]]

using namespace Rcpp;

typedef std::pair<double, int> paired_doubleint;

inline bool cmp_double_smaller(paired_doubleint left, paired_doubleint right) {
  return left.first < right.first;
}

inline bool cmp_double_bigger(paired_doubleint left, paired_doubleint right) {
  return left.first > right.first;
}

inline int round_to_int(double x) {
  int result = 0;
  if(x > 0) x += 0.5;
  if(x < 0) x -= 0.5;
  result = (int)x;
  return result;
}

/* order (ascending) function for IntegerVector */
inline IntegerVector order_vector(NumericVector x, bool decreasing) {
  IntegerVector result(x.size());
  std::iota(result.begin(), result.end(), 0);
  if(decreasing){
    auto comparator = [&x](int a, int b){ return x[a] > x[b]; };
    std::sort(result.begin(), result.end(), comparator);
  } else{
    auto comparator = [&x](int a, int b){ return x[a] < x[b]; };
    std::sort(result.begin(), result.end(), comparator);
  }

  return result;
}

/* sort PairMatrix using an already computed order */
inline void sort_using_order_intmat(IntegerMatrix &x, const IntegerVector &x_o) {
  // Consider using std::swap as in https://stackoverflow.com/a/838789/5861244
  for(int i = 0; i < x.nrow(); ++i){
    IntegerVector temp = x(i, _);
    auto x_o_it = x_o.begin();
    for(int j = 0; j < x.ncol(); ++j, ++x_o_it)
      x(i, j) = temp[*x_o_it];
  }
}

//[[Rcpp::export]]
IntegerVector sample_uniform(int range, int n) {
    if(range < 1)
      stop("range is zero or negative");

    if(n < 1)
      stop("number of samples is zero or negative");
    NumericVector temp = runif(n,-0.499,range-0.501);
    IntegerVector results(n);

    for (int i = 0; i < n; ++i) results[i] = round_to_int(temp[i]);

    return results;
}

//[[Rcpp::export]]
IntegerVector sample_int_replace(NumericVector probabilities, int n) {

    if(n < 1) stop("number of samples is zero or negative");

    int m = probabilities.size();


    IntegerVector alias(m);

    NumericVector probabilities1 = clone(probabilities);

    double avg = sum(probabilities1);

    probabilities1 = probabilities1*(m/avg);

    if(var(probabilities1) < 0.01) {
      IntegerVector results(n);
      NumericVector temp = runif(n,-0.499,m-0.501);
      for(int i = 0; i < n; ++i) {
        results[i] = round_to_int(temp[i]);
      }
      return results;
    }

    NumericVector probs(m);

    int index_small = 0;
    int index_large = 0;

    IntegerVector small(m);
    IntegerVector large(m);

    for (int i = 0; i < m; ++i) {
        if(probabilities1[i] > 1.0) {
            ++index_large;
            large[index_large] = i;

        } else {
            ++index_small;
            small[index_small] = i;

        }
    }

    int less = 0;
    int more = 0;
    while (index_small > 0 && index_large > 0) {
        less = small[index_small];
        more = large[index_large];
        --index_small;
        --index_large;
        probs[less] = probabilities1[less];
        alias[less] = more;
        probabilities1[more] = probabilities1[more]+probabilities1[less]-1;

        if (probabilities1[more] > 1.0) {
            ++index_large;
            large[index_large] = more;
        } else {
            ++index_small;
            small[index_small] = more;
        }
    }
    while (index_small > 0) {
        probs[small[index_small]]=1;
        --index_small;
    }
    while (index_large > 0) {
        probs[large[index_large]]=1;
        --index_large;
    }

    IntegerVector results(n);
    NumericVector temp = runif(n,-0.499,m-0.501);
    NumericVector temp1 = runif(n,0,1);
    int col = 0;
    for (int i = 0; i < n; ++i) {
        col = round_to_int(temp[i]);
        if (temp1[i] < probs[col]) {
            results[i] = col;
        } else {
            results[i] = alias[col];
        }
    }
    return results;
}

/* entrywise multiplication of X and Y */
//[[Rcpp::export]]
NumericMatrix prod_matrix_vector(IntegerMatrix X, NumericVector r) {
    int n = X.nrow(); int p = X.ncol(); NumericMatrix result(n,p);

    for(int j = 0; j < p; ++j) {
      for(int i = 0; i < n; ++i) {
        result(i,j)=X(i,j)*r[i];
        }
    }

    return result;
}

/* sum over columns of matrix for specific indexes */
//[[Rcpp::export]]
NumericVector colsum_index(NumericMatrix X,IntegerVector indexes) {
    int p = X.ncol();
    NumericVector result(p);
    double temp = 0.0;
    for(int i = 0; i < p; ++i) {
        temp = 0.0;
        for(int j = 0; j < indexes.size(); ++j) {
            temp += X(indexes[j],i);
        }
        result[i] = temp;
    }
    return result;
}

//
void copy_vector_to_column(NumericMatrix X, NumericVector Y, int k) {
    for (int i = 0; i < X.nrow(); ++i) {
        X(i,k)=Y[i];
    }
}

/* sum over columns of matrix for specific indexes */
//[[Rcpp::export]]
NumericVector absolute_covariates(NumericMatrix X, NumericVector Y) {
    int p = X.ncol();
    int n = X.nrow();
    NumericVector abs_cov(p);

    double temp = 0.0;
    for (int l = 0; l < p; ++l) {
        temp = 0.0;
        for (int i = 0; i < n; ++i) {
            temp += X(i,l)*Y[i];
        }
        abs_cov[l]=std::abs(temp)/n;
    }
    return abs_cov;
}

//[[Rcpp::export]]
NumericVector absolute_covariates_pairs(IntegerMatrix pairs, NumericMatrix X, NumericVector Y) {
    int n = X.nrow();

    int nr_pairs = pairs.ncol();
    NumericVector abs_cov(nr_pairs);

    int l = 0;
    int k = 0;
    double temp = 0.0;
    for(int j = 0; j < nr_pairs; ++j) {
      l = pairs(0,j);
      k = pairs(1,j);
      temp = 0.0;
      for(int i = 0; i < n; ++i) {
          temp += Y[i]*X(i,l)*X(i,k);
      }
      abs_cov[j]=std::abs(temp)/n;
    }
    return abs_cov;
}


typedef std::pair<int, int> paired_int;

bool cmp_paired_int(paired_int left, paired_int right) {
  if(left.first < right.first) return true;

  if(left.first == right.first) {
    if(left.second < right.second) {
      return true;
    }
  }
  return false;
}

/* order (ascending) function for IntegerVector */
//[[Rcpp::export]]
IntegerMatrix clean_pairs(IntegerMatrix pairs) {
    int n_pairs = pairs.ncol();

    if(n_pairs == 1) return pairs;

    std::vector<paired_int> pairs1;
    pairs1.reserve(n_pairs);

    for(int i = 0; i < n_pairs; ++i) {
        if(pairs(0,i) <= pairs(1,i)) {
          pairs1.push_back(std::make_pair(pairs(0,i), pairs(1,i)));
        } else {
          pairs1.push_back(std::make_pair(pairs(1,i), pairs(0,i)));
        }
    }
    std::sort(pairs1.begin(), pairs1.end(), cmp_paired_int);

    pairs1.erase(std::unique(pairs1.begin(), pairs1.end()), pairs1.end());

    IntegerMatrix result(2,pairs1.size());
    for(int i = 0; i < result.ncol(); ++i) {
        result(0,i) = pairs1[i].first;
        result(1,i) = pairs1[i].second;
    }
    return result;
}

/* translates a continuous matrix to binary finding the closest binary representation */
//[[Rcpp::export]]
IntegerMatrix translate_to_binary(NumericMatrix X, int max_number_of_iterations) {
    int p = X.ncol(); int n = X.nrow();

    double convergence_ratio = 0.025;

    NumericVector colSum(p); LogicalVector cluster(n); IntegerMatrix X_bin(n,p);

    for (int i = 0; i < n; ++i) colSum = colSum+X(i,_);

    colSum = colSum/n;

    double center1 = 0.0; double center2 = 0.0;

    double temp_center1 = 0.0; double temp_center2 = 0.0;

    double barrier = 0.0;

    int counter = 0; int ratio_changed = 1; int count1 = 0; int count2 = 0;

    for(int j = 0; j < p; ++j) {
        barrier = colSum[j];
        for (int i = 0; i < n; ++i) {
            if (X(i,j) < barrier) {
                cluster(i)=false; ++count1; center1 += X(i,j);
            } else {
                cluster(i)=true; ++count2; center2 += X(i,j);
            }
        }

        center1 = center1/count1; center2 = center2/count2;

        counter = 0;
        while (counter < max_number_of_iterations) {

            temp_center1 = center1*count1; temp_center2 = center2*count2;

            ratio_changed = 0;

            for(int i = 0; i < n; ++i) {
                if(std::abs(X(i,j)-center1) < std::abs(X(i,j)-center2)) {
                    if (cluster(i)) {
                        cluster(i)=false;
                        temp_center1 += X(i,j); temp_center2 -= X(i,j);
                        count1 += 1; count2 -= 1; ++ratio_changed;
                    }
                } else {
                    if (!cluster(i)) {
                        cluster(i)=true;
                        temp_center1 -= X(i,j); temp_center2 += X(i,j);
                        count1 -= 1; count2 += 1; ++ratio_changed;
                    }
                }
            }
            if (count1 > 0) {
                center1 = temp_center1/count1;
            } else {
                break;
            }
            if (count2 > 0) {
                center2 = temp_center2/count2;
            } else {
                break;
            }
            ++counter;
            if (ratio_changed/n < convergence_ratio) {
                break;
            }
        }
        for (int i = 0; i < n; ++i) {
            if (cluster(i)) {
                X_bin(i,j) =1;
            } else {
                X_bin(i,j) = -1;
            }
        }
    }
    return X_bin;
}

/* estimate interaction frequency of non interacting pairs */
//[[Rcpp::export]]
NumericVector estimate_background_interaction_frequency(IntegerMatrix X, IntegerVector Y, int number_of_samples) {
    int n = X.nrow(); int p = X.ncol();
    double temp = 0.0; int cord1 = 0; int cord2 = 0;

    if (number_of_samples < 2) {
        number_of_samples = 2;
    }

    NumericVector frequencies(number_of_samples);
    IntegerVector samples1 = sample_uniform(p,number_of_samples); IntegerVector samples2 = sample_uniform(p,number_of_samples);

    for (int j = 0; j < number_of_samples; ++j) {
        temp = 0.0; cord1 = samples1[j]; cord2 = samples2[j];
        for(int i = 0; i < n; ++i) {
          if(Y[i]==X(i,cord1)*X(i,cord2)) ++temp;
        }
        frequencies[j]= temp/n;
    }
    return frequencies;
}

//[[Rcpp::export]]
List find_strongest_pairs(List pairs, NumericMatrix X, NumericVector Y, int max_number_of_pairs) {
  int size_of_list = pairs.size();

  if(max_number_of_pairs > 500) stop("You consider too many pairs, usually one would only consider 10");

  int total_number = size_of_list*max_number_of_pairs;
  IntegerMatrix temp_pairs(2,total_number);
  NumericVector temp_strength(total_number);

  NumericVector temps;
  IntegerVector otemps;
  int counter = 0;
  for(int m = 0; m < size_of_list; ++m) {
    IntegerMatrix tempp = pairs(m);
    temps = absolute_covariates_pairs(tempp,X,Y);
    otemps = order_vector(temps,true);
    sort_using_order_intmat(tempp,otemps);
    for(int j = 0; j <  std::min(max_number_of_pairs,tempp.ncol()); ++j) {
      temp_pairs(0,counter) = tempp(0,j);
      temp_pairs(1,counter) = tempp(1,j);
      ++counter;
    }
  }
  //fix for the bug that made (1,1) an intercation pair most of the time
  IntegerMatrix tpairs(2,counter);

  for(int i = 0; i < counter; ++i) {
    tpairs(0,i) = temp_pairs(0,i);
    tpairs(1,i) = temp_pairs(1,i);
  }
  temp_pairs = clean_pairs(tpairs);

  temp_strength = absolute_covariates_pairs(temp_pairs,X,Y);
  otemps = order_vector(temp_strength,true);
  sort_using_order_intmat(temp_pairs,otemps);

  int length = std::min(max_number_of_pairs,(int) temp_strength.size());
  IntegerMatrix tempp(2,length);
  temps =  NumericVector(length);

  for(int j =0; j < length; ++j) {
    tempp(0,j)=temp_pairs(0,j);
    tempp(1,j)=temp_pairs(1,j);
    temps[j]=temp_strength[otemps[j]];
  }
  List result(2);
  result(0)=tempp;
  result(1)=temps;
  return result;
}

//[[Rcpp::export]]
IntegerMatrix equalpairs(NumericVector u, NumericVector v, IntegerVector ou, IntegerVector ov, int max_number_of_pairs) {
    //set sizes of array
    int nu = u.size();
    int nv = v.size();

    //init two lists to store pairs
    std::list<int> pairs_u;
    std::list<int> pairs_v;

    //set pointers
    int start = 0;
    int j = 0;

    //set counter
    int count = 0;

    //set precision epsilon
    double eps = 0.0000001;

    //start looping through u vector
    for(int i = 0; i < nu; ++i) {

        //increase if too small
        while(v[start]<u[i]-eps && start < nv-1) {
            ++start;
        }

        //if close consider the pairs that might be close
        if(std::abs(v[start]-u[i]) < eps) {
            j = start;
            while(std::abs(v[j]-u[i]) < eps) {
                //add pairs that are clsoe
                pairs_u.push_front(ou[i]);
                pairs_v.push_front(ov[j]);

                ++j;
                ++count;
                if(j >= nv) {
                    break;
                }
            }
        }
        //if there are too many pairs kill the search
        if(count > max_number_of_pairs) {
            break;
        }
    }
    int n = 0;
    //fill pairs in a 2x|pairs| matrix
    if(pairs_u.size() > 0) {
      IntegerMatrix pairs(2,pairs_u.size());
      while(!pairs_u.empty()) {
          pairs(0,n)=pairs_u.back();
          pairs_u.pop_back();
          pairs(1,n)=pairs_v.back();
          pairs_v.pop_back();
          ++n;
      }
      return pairs;
    }
    IntegerMatrix pairs(2,1);
    return pairs;
}

//[[Rcpp::export]]
List projected_equal_pairs(IntegerMatrix X, NumericVector Y, int number_of_runs, int max_number_of_collisions, bool negative) {

    int n = X.nrow();
    int p = X.ncol();

    int length_of_pairs = number_of_runs;
    //now start with the actual task
    NumericVector abs_Y = abs(Y);
    IntegerVector sign_Y = sign(Y);

    NumericVector frequencies = estimate_background_interaction_frequency(X,sign_Y,std::min(p,2000));

    double s2 = mean(frequencies);
    //double s2 = 0.5;
    if(s2 > 0.9) {Rcout << "s2=" << s2 <<" "; stop("background interaction strength seems to be unusually high");}
    if(s2 < 0.1) {Rcout << "s2=" << s2 <<" "; stop("background interaction strength seems to be unusually low");}

    float pp = p;
    int size_of_subsample = round_to_int(-log(pp)/log(s2));

    if(size_of_subsample < 1) stop("calculated sub sample size is below 1");
    if(size_of_subsample > n) stop("calculated sub sample size is above n");

    NumericVector r = runif(n);
    NumericVector ry(n);
    for(int i = 0; i < n; ++i) {
        ry[i]= r[i]*sign_Y[i];
    }
    NumericMatrix Xr = prod_matrix_vector(X,r);

    NumericMatrix Zr = prod_matrix_vector(X,ry);
    if(negative) {
      length_of_pairs = 2*number_of_runs;
    }
    List pairs(length_of_pairs);
    IntegerVector indexes;
    NumericVector x(p);
    NumericVector z(p);
    NumericVector z_minus(p);
    IntegerVector order_x(p);
    IntegerVector order_z(p);
    IntegerVector order_z_minus(p);

    int counter = 0;
    for (int l = 0; l < number_of_runs; ++l) {
        indexes = sample_int_replace(abs_Y,size_of_subsample);
        x = colsum_index(Zr,indexes);
        z = colsum_index(Xr,indexes);


        order_x = order_vector(x,false);
        order_z = order_vector(z,false);

        x.sort();
        z.sort();
        pairs(counter) = clean_pairs(equalpairs(x,z,order_x,order_z,max_number_of_collisions));
        ++counter;
        if(negative) {
          int j = 0;
          for(int i = p-1; i >= 0; --i) {
              z_minus[i] = -z[j];
              order_z_minus[i] = order_z[j];
              ++j;
          }
          pairs(counter) = clean_pairs(equalpairs(x,z_minus,order_x,order_z_minus,max_number_of_collisions));
          ++counter;
        }
    }
    return pairs;
}

//[[Rcpp::export]]
List naive_interaction_search(NumericMatrix X, NumericVector Y, int max_number_of_pairs) {
  int n = X.nrow();
  int p = X.ncol();
  IntegerMatrix pairs(2,(p*(p+1))/2);
  NumericVector strength((p*(p+1))/2);

  int count = 0;
  for(int k = 0; k < p; ++k) {
    for(int l = k; l < p; ++l) {
      double temp = 0;
      for(int i = 0; i < n; ++i) {
          temp += Y[i]*X(i,l)*X(i,k);
      }
      pairs(0,count)=k;
      pairs(1,count)=l;
      strength(count)=std::abs(temp);
      ++count;
    }
  }
  IntegerVector order_pairs = order_vector(strength,true);
  sort_using_order_intmat(pairs,order_pairs);
  strength = strength[order_pairs];
  List result(2);

  if(max_number_of_pairs < strength.size()) {
    IntegerMatrix temp_pairs(2,max_number_of_pairs);
    NumericVector temp_strength(max_number_of_pairs);

    for(int i = 0; i < max_number_of_pairs; ++i) {
      temp_pairs(0,i)=pairs(0,i);
      temp_pairs(1,i)=pairs(1,i);
      temp_strength[i]=strength[i];

    }
    result(0) = temp_pairs;
    result(1) = temp_strength;
    return result;
  }

  result(0) = pairs;
  result(1) = strength;
  return result;
}

// [[Rcpp::export]]
List interaction_search(NumericMatrix X, NumericVector Y, int number_of_runs, int max_number_of_pairs,bool negative,bool binary) {
  int n = X.nrow();
  int p = X.ncol();
  List result(2);
  if(p < 20) {
    result = naive_interaction_search(X,Y,max_number_of_pairs);
  } else {
    IntegerMatrix X_binary(n,p);
    if(binary) {
      X_binary = (IntegerMatrix) X;
    } else {
      X_binary = translate_to_binary(X,10);
    }
    int max_number_of_collisions = 2*p;
    List pairs = projected_equal_pairs(X_binary, Y,number_of_runs, max_number_of_collisions, negative);
    result = find_strongest_pairs(pairs,X,Y,max_number_of_pairs);
  }
  return result;
}

//[[Rcpp::export]]
List interaction_search_low_level(IntegerMatrix X_binary,NumericMatrix X, NumericVector Y, int number_of_runs,int max_number_of_pairs) {
  int p = X.ncol();
  List result(2);
  if(p < 20) {
    result = naive_interaction_search(X,Y,max_number_of_pairs);
  } else {
    int max_number_of_collisions = 2*p;
    List pairs = projected_equal_pairs(X_binary, Y,number_of_runs, max_number_of_collisions, true);
    result = find_strongest_pairs(pairs,X,Y,max_number_of_pairs);
  }
  return result;
}

//[[Rcpp::export]]
double soft_threshold(double beta_tilde, double normalization, double lambda, double alpha) {
    double denom = normalization+lambda*(1.0-alpha);
    double s = ((beta_tilde > 0)-(beta_tilde < 0))*std::max(std::abs(beta_tilde)-lambda*alpha,0.0);
    return s/denom;
}

//[[Rcpp::export]]
NumericVector create_lambda_sequence(double lambda_max, int n_lambda, double factor_eps_inv) {
    double eps_inv = factor_eps_inv; NumericVector lambdas(n_lambda);
    if (n_lambda < 2) stop("n_lambda has to be at least two");
    lambdas[0] = lambda_max;

    double b = std::log(eps_inv)/(n_lambda-1);
    for (int i = 1; i < n_lambda; ++i) lambdas[i] = lambda_max*std::exp(-i*b);

    return lambdas;
}

/*scans main effects and includes new potential candidates*/
//[[Rcpp::export]]
bool scan_main_effects(const NumericMatrix &X, const NumericVector &Y, const NumericVector &residuals,List main_effects, List beta_main,
                       const NumericVector &lambdas, double alpha, int r, int add_max, bool strong) {
    /*check if sequential strong rule is applied at very beginning of lambda sequence*/
    int p = X.ncol();

    if (strong && r == 0) strong = false; //stop("strong rule cannot be applied at lambda max");

    IntegerVector temp_main_effects = main_effects[r];
    NumericVector temp_beta_main = beta_main[r];

    /*find effects already included*/
    LogicalVector already_contained(p);
    for (int l = 0; l < temp_main_effects.size(); ++l) {
        already_contained[temp_main_effects[l]] = true;
    }

    NumericVector covs = absolute_covariates(X,residuals);
    IntegerVector order_covs = order_vector(covs,true);

    double threshold = alpha*lambdas[r];

    if (strong) threshold = alpha*(2*lambdas[r]-lambdas[r-1]);

    std::vector<int> coefs;
    coefs.reserve(p);
    int ptr = 0;
    int count = 0;
    /*pick vectors that are violating threshold condition*/
    while (covs[order_covs[ptr]] >= threshold && ptr < p && count < add_max) {
        if (!already_contained[order_covs[ptr]]) {
            coefs.push_back(order_covs[ptr]);
            ++count;
        }
        ++ptr;
    }

    /*if no vectors are added then communicate that nothing has changed*/
    if (coefs.size() == 0) {
        return false;
    }

    if(temp_main_effects.size() == 0) {
        /*no main effects have yet been collected*/
        temp_main_effects = IntegerVector(coefs.begin(),coefs.end());
        temp_beta_main = NumericVector(temp_main_effects.size());
        main_effects[r] = temp_main_effects;
        beta_main[r] = temp_beta_main;
        return true;
    }

    IntegerVector new_main_effects(temp_main_effects.size()+coefs.size());
    IntegerVector new_beta_main(temp_beta_main.size()+coefs.size());

    main_effects[r] = IntegerVector(temp_main_effects.size()+coefs.size());
    beta_main[r] = IntegerVector(temp_beta_main.size()+coefs.size());

    for(int l = 0; l < temp_main_effects.size(); ++l) {
      new_main_effects[l] = temp_main_effects[l];
      new_beta_main[l] = temp_beta_main[l];
    }

    ptr = 0;
    int sum_size = temp_main_effects.size()+coefs.size();
    for (int l = temp_main_effects.size(); l < sum_size; ++l) {
        new_main_effects[l] = coefs.at(ptr);
        ++ptr;
    }
    main_effects[r] = new_main_effects;
    beta_main[r] = new_beta_main;

    return true;
}

//[[Rcpp::export]]
NumericVector scale_intr(NumericMatrix X, int pair_x, int pair_y) {
    int n = X.nrow();
    NumericVector intr(n);
    double mean = 0.0;
    double var = 0.0;
    for (int i = 0; i < n; ++i) {
        intr[i]=X(i,pair_x)*X(i,pair_y);
        mean += intr[i];
    }
    mean = mean/n;
    for (int i = 0; i < n; ++i) {
        intr[i]=intr[i]-mean;
        var += intr[i]*intr[i];
    }
    if (std::abs(var) > 0 && n > 1) {
            intr = intr/std::sqrt(var/(n-1));
    }
    return intr;
}

//[[Rcpp::export]]
bool scan_intr_effects(const NumericMatrix &X, const NumericVector &Y, const IntegerMatrix &X_bin, const NumericVector &residuals,
                       List intr_effects, List beta_intr, NumericMatrix &intr_vars,
                       const NumericVector &lambdas, double alpha, int r, int projections, bool strong) {

    int n = X.nrow();

    double threshold = alpha*lambdas[r];
    if (strong) threshold = alpha*(2*lambdas[r]-lambdas[r-1]);

    List result_interaction_search(2);

    result_interaction_search = interaction_search_low_level(X_bin,X,residuals,projections,20);
    IntegerMatrix pairs_is = result_interaction_search(0);

    int number_considered = std::min(20,pairs_is.ncol());

    NumericVector covariates(number_considered);
    NumericMatrix intr_vars_new(n,number_considered);

    IntegerMatrix temp_intr_effects = intr_effects[r];
    NumericVector temp_beta_intr = beta_intr[r];
    NumericVector temp_itr(n);

    double sum1 = 0.0;

    int count = 0;

    LogicalVector consider(number_considered);

    std::vector<int> indexes;
    indexes.reserve(number_considered);


    for (int l = 0; l < number_considered; ++l) {
        sum1 = 0.0;
        temp_itr = scale_intr(X,pairs_is(0,l),pairs_is(1,l));
                for (int i = 0; i < n; ++i) {
            sum1 += temp_itr[i]*residuals[i];
            intr_vars_new(i,l) = temp_itr[i];
        }

        covariates[l] = sum1;
        if (std::abs(sum1) > threshold) {

            bool duplicate = false;

            for (int k = 0; k < temp_intr_effects.ncol(); ++k) {
                if (temp_intr_effects(0,k) == pairs_is(0,l) && temp_intr_effects(1,k) == pairs_is(1,l)) {
                    duplicate = true;
                    break;
                }
                if (temp_intr_effects(1,k) == pairs_is(0,l) && temp_intr_effects(0,k) == pairs_is(1,l)) {
                    duplicate = true;
                    break;
                }
            }

            if(!duplicate) {
                indexes.push_back(l);
                ++count;
            }
        }

    }
    if (count == 0) {
        return false;
    }
    int nr_intr = temp_intr_effects.ncol()+count;

    IntegerMatrix new_intr_effects(2,nr_intr);
    NumericVector new_beta_intr(nr_intr);

    for (int l = 0; l < temp_intr_effects.ncol(); ++l) {
        new_intr_effects(0,l) = temp_intr_effects(0,l);
        new_intr_effects(1,l) = temp_intr_effects(1,l);
        new_beta_intr[l] = temp_beta_intr[l];
    }

    count = 0;
    for (int l = temp_intr_effects.ncol(); l < nr_intr; ++l) {
        new_intr_effects(0,l) = pairs_is(0,indexes.at(count));
        new_intr_effects(1,l) = pairs_is(1,indexes.at(count));
        new_beta_intr[l] = 0.0;
        ++count;
    }
    intr_effects[r] = new_intr_effects;
    beta_intr[r] = new_beta_intr;
    return true;
}

//[[Rcpp::export]]
void update_intr_final(List intr_effects, List beta_intr) {
    int r = intr_effects.size()-1;

    IntegerMatrix temp_eff = intr_effects[r];
    NumericVector temp_bet(beta_intr[r]);
    intr_effects(0) = temp_eff;
    beta_intr(0) = temp_bet;
}

//[[Rcpp::export]]
NumericVector calculate_xbeta(const NumericMatrix &X, const NumericVector &Y,
                                  const NumericVector &intercept,
                                  const List main_effects, const List beta_main,
                                  const List intr_effects, const List beta_intr, const NumericMatrix &intr_vars,
                                  int r) {

    int n = X.nrow();
    NumericVector xbeta(n);

    /* add intercept */
    for (int i = 0; i < n; ++i) {
        xbeta[i] = intercept[r];
    }

    /* deduct main effects */
    IntegerVector temp_main_effects(main_effects[r]);
    NumericVector temp_beta_main(beta_main[r]);
    if (temp_main_effects.size() > 0) {
        for (int l = 0; l < temp_main_effects.size(); ++l) {
            for (int i = 0; i < n; ++i) {
                xbeta[i] = xbeta[i]+X(i,temp_main_effects[l])*temp_beta_main[l];
            }
        }
    }
    /* deduct intr effects */
    NumericVector temp_beta_intr(beta_intr[r]);
    if (temp_beta_intr.size() > 0) {
        for (int l = 0; l < temp_beta_intr.size(); ++l) {
            for (int i = 0; i < n; ++i) {
                xbeta[i] = xbeta[i]+intr_vars(i,l)*temp_beta_intr[l];
            }
        }
    }
    return xbeta;
}

/* calculate residuals */
//[[Rcpp::export]]
NumericVector calculate_residuals(const NumericMatrix &X, const NumericVector &Y,
                               const NumericVector &intercept,
                               const List main_effects, const List beta_main,
                               const List intr_effects, const List beta_intr, const NumericMatrix &intr_vars,
                               int r) {
    int n = X.nrow();
    NumericVector res(n);
    res = Y-calculate_xbeta(X, Y, intercept, main_effects, beta_main, intr_effects, beta_intr, intr_vars, r);
    return res;
}

//[[Rcpp::export]]
int iterate(const NumericMatrix &X, const NumericVector &Y, NumericVector &residuals,
            const NumericVector &intercept,
            const List main_effects, List beta_main,
            const List intr_effects, List beta_intr,const NumericMatrix &intr_vars,
            const NumericVector &weights, const NumericVector &lambdas, double alpha, int r, int maxiter_inner) {

    int n = X.nrow();

    /*count number of iterations*/
    int count_iterations = 0;

    int steps_to_update_residuals = 1000;

    double update_precision = 0.0001;

    /* update residuals before starting to iterate */
    residuals = calculate_residuals(X,Y,intercept,main_effects,beta_main,intr_effects, beta_intr, intr_vars,r);

    IntegerVector temp_main_effects(main_effects[r]);
    IntegerMatrix temp_intr_effects = intr_effects[r];
    /* update: note that we decide three situations*/
    int nr_main = temp_main_effects.size();
    int nr_intr = temp_intr_effects.ncol();


    /* update main effects */
    NumericVector temp_beta_main = beta_main[r];

    NumericVector temp_beta_intr = beta_intr[r];

    bool converged = false;

    double sum1 = 0.0;
    double sum2 = 0.0;
    double temp = 0.0;
    double beta_tilde = 0.0;
    int variable = 0;


    while (!converged && count_iterations < maxiter_inner) {
        converged = true;

        /* update main effects*/
        if (nr_main > 0) {
            for (int l = 0; l < nr_main; ++l) {
                variable = temp_main_effects[l];
                sum1 = 0.0;
                sum2 = 0.0;
                for (int i = 0; i < n; ++i) {
                    temp = weights[i]*X(i,variable);
                    sum1 += temp*residuals[i];
                    sum2 += temp*X(i,variable);
                }
                beta_tilde = soft_threshold(sum1+temp_beta_main[l],sum2,lambdas[r]*std::sqrt(2.0),alpha);


                if (std::abs(beta_tilde - temp_beta_main[l]) > update_precision) {
                    converged = false;

                    temp = temp_beta_main[l]-beta_tilde;
                    /* update residuals */
                    for (int i = 0; i < n; ++i) {
                        residuals[i] = residuals[i]+X(i,variable)*temp;
                    }
                    temp_beta_main[l] = beta_tilde;
                    //residuals = calculate_residuals(X,Y,intercept,main_effects,beta_main,intr_effects, beta_intr, intr_vars,r);

                }
            }
        }

        /* update intr effects */
        if (nr_intr > 0) {
            for (int l = 0; l < nr_intr; ++l) {
                sum1 = 0.0;
                sum2 = 0.0;
                for (int i = 0; i < n; ++i) {
                    temp = weights[i]*intr_vars(i,l);
                    sum1 += temp*residuals[i];
                    sum2 += temp*intr_vars(i,l);
                }
                beta_tilde = soft_threshold(sum1+temp_beta_intr[l],sum2,lambdas[r],alpha);
                if (std::abs(beta_tilde - temp_beta_intr[l]) > update_precision) {
                    converged = false;

                    temp = temp_beta_intr[l]-beta_tilde;
                    /* update residuals */
                    for (int i = 0; i < n; ++i) {
                        residuals[i] = residuals[i]+intr_vars(i,l)*temp;
                    }
                    temp_beta_intr[l] = beta_tilde;
                    //residuals = calculate_residuals(X,Y,intercept,main_effects,beta_main,intr_effects, beta_intr, intr_vars,r);
                }
            }
        }
        ++count_iterations;
        if(count_iterations % steps_to_update_residuals == 0) {
            residuals = calculate_residuals(X,Y,intercept,main_effects,beta_main,intr_effects, beta_intr, intr_vars,r);
        }
    }

    /* propagate results */
    if(nr_main > 0) {
        beta_main[r] = temp_beta_main;
    }

    if(nr_intr > 0) {
        beta_intr[r] = temp_beta_intr;
    }
    return count_iterations;
}

/* build interaction variables*/
//[[Rcpp::export]]
NumericMatrix update_intr_vars(const NumericMatrix &X, List intr_effects, bool standardize, int r) {
    int n = X.nrow();
    IntegerMatrix temp_intr_effects = intr_effects[r];
    NumericMatrix intr_vars(n,temp_intr_effects.ncol());

    for(int l = 0; l < intr_vars.ncol(); ++l) {
        int pair_x = temp_intr_effects(0,l);
        int pair_y = temp_intr_effects(1,l);

        NumericVector intr(n);
        double mean = 0.0;
        double var = 0.0;
        for (int i = 0; i < n; ++i) {
            intr[i]=X(i,pair_x)*X(i,pair_y);
            mean += intr[i];
        }
        if(standardize) {
          mean = mean/n;
          for (int i = 0; i < n; ++i) {
              intr[i]=intr[i]-mean;
              var += intr[i]*intr[i];
          }
          if (std::abs(var) > 0 && n > 1) {
              intr = intr/std::sqrt(var/(n-1));
              intr_vars(_,l) = intr;
          }
        }
    }
    return intr_vars;
}

/* clean effects from those with zero entries */
//[[Rcpp::export]]
void clean_all_effects(List main_effects, List beta_main,
                       List intr_effects, List beta_intr,
                       int r) {

    double precision = 0.001;

    int count = 0;

    IntegerVector temp_main_effects(main_effects[r]);
    NumericVector temp_beta_main(beta_main[r]);

    if (temp_beta_main.size() > 0) {
        for (int l = 0; l < temp_beta_main.size(); ++l) {
            if (std::abs(temp_beta_main[l])> precision) {
                ++count;
            }
        }
        if(count == 0) {
            count = 1;
        }
        IntegerVector new_main_effects(count);
        NumericVector new_beta_main(count);

        count = 0;
        for (int l = 0; l < temp_beta_main.size(); ++l) {
            if (std::abs(temp_beta_main[l])> precision) {
                new_main_effects[count] = temp_main_effects[l];
                new_beta_main[count] = temp_beta_main[l];
                ++count;
            }
        }
        main_effects[r] =new_main_effects;
        beta_main[r] = new_beta_main;
    }

    count = 0;

    IntegerMatrix temp_intr_effects = intr_effects[r];
    NumericVector temp_beta_intr(beta_intr[r]);

    if (temp_beta_intr.size() > 0) {
        for (int l = 0; l < temp_beta_intr.size(); ++l) {
            if (std::abs(temp_beta_intr[l])> precision) {
                ++count;
            }
        }
        /* if there are no effects just take the zero effect as dummy effect*/
        if (count == 0) {
            count = 1;
        }
        IntegerMatrix new_intr_effects(2,count);
        NumericVector new_beta_intr(count);
        count = 0;
        for (int l = 0; l < temp_beta_intr.size(); ++l) {
            if (std::abs(temp_beta_intr[l])> precision) {
                new_intr_effects(0,count) = temp_intr_effects(0,l);
                new_intr_effects(1,count) = temp_intr_effects(1,l);
                new_beta_intr[count] = temp_beta_intr[l];
                ++count;
            }
        }
        intr_effects[r] =new_intr_effects;
        beta_intr[r] = new_beta_intr;
    }
}

//[[Rcpp::export]]
void warm_start(List main_effects, List beta_main,
                List intr_effects, List beta_intr,
                int r) {
    main_effects[r+1] = main_effects[r];
    beta_main[r+1]= beta_main[r];
    intr_effects[r+1] = intr_effects[r];
    beta_intr[r+1] = beta_intr[r];
}

// [[Rcpp::export]]
List gaussiglmnet(NumericMatrix X, NumericVector Y, NumericVector lambdas, double alpha, bool standardize, int max_main_effects, int max_interaction_effects, int max_outer, int number_of_nnis_runs) {
    int n = X.nrow(); int n_lambda = lambdas.size();

    int maxiter_inner=1000; int add_max_main = 100;

    IntegerMatrix X_bin = translate_to_binary(X,10);

    List main_effects(n_lambda); List intr_effects(n_lambda); List beta_main(n_lambda); List beta_intr(n_lambda);

    for (int r = 0; r < n_lambda; ++r) { main_effects[r] = IntegerVector(1); intr_effects[r] = IntegerMatrix(2,1); beta_main[r] = NumericVector(1); beta_intr[r] = NumericVector(1); }

    NumericVector intercept(n_lambda);  NumericMatrix intr_vars(n,1); NumericVector weights(n); NumericVector residuals(Y);

    intercept[0]=mean(Y);

    for(int i = 0; i < n; ++i)  weights[i]=1.0/n;

    NumericVector covs = absolute_covariates(X,residuals);

    if(lambdas[0] < 0) lambdas = create_lambda_sequence(max(covs)-0.001,n_lambda,10.0);

    bool changed = false;

    for (int r = 0; r < n_lambda; ++r) {
        intercept[r] = mean(Y);
        for (int iter = 0; iter < max_outer; ++iter) {

            changed = true;

            changed = scan_main_effects(X,Y,residuals,main_effects,beta_main,lambdas,alpha,r,add_max_main,true);

            iterate(X,Y,residuals,intercept,main_effects,beta_main,intr_effects,beta_intr,intr_vars,weights,lambdas,alpha,r,maxiter_inner);

            residuals = calculate_residuals(X,Y,intercept,main_effects,beta_main,intr_effects,beta_intr,intr_vars,r);

            changed = changed & scan_intr_effects(X,Y,X_bin,residuals,intr_effects,beta_intr,intr_vars,lambdas,alpha,r,number_of_nnis_runs,true);

            intr_vars = update_intr_vars(X,intr_effects,standardize,r);

            iterate(X,Y,residuals,intercept,main_effects,beta_main,intr_effects,beta_intr,intr_vars,weights,lambdas,alpha,r,maxiter_inner);

            residuals = calculate_residuals(X,Y,intercept,main_effects,beta_main,intr_effects,beta_intr,intr_vars,r);

            if (!changed) {
                break;
            }
        }

        clean_all_effects(main_effects,beta_main,intr_effects,beta_intr,r);

        intr_vars = update_intr_vars(X,intr_effects,standardize,r);

        iterate(X,Y,residuals,intercept,main_effects,beta_main,intr_effects,beta_intr,intr_vars,weights,lambdas,alpha,r,maxiter_inner);

        if(r < n_lambda-1) warm_start(main_effects,beta_main,intr_effects,beta_intr,r);
    }

    int last_entry = intr_effects.size()-1;
    /*start the second loop*/
    for (int r = 0; r < n_lambda; ++r) {
            intr_effects[r] = intr_effects[last_entry];

            beta_intr[r] = beta_intr[last_entry];

            scan_main_effects(X,Y,residuals,main_effects,beta_main,lambdas,alpha,r,add_max_main,true);

            iterate(X,Y,residuals,intercept,main_effects,beta_main,intr_effects,beta_intr,intr_vars,weights,lambdas,alpha,r,maxiter_inner);

            residuals = calculate_residuals(X,Y,intercept,main_effects,beta_main,intr_effects,beta_intr,intr_vars,r);

            iterate(X,Y,residuals,intercept,main_effects,beta_main,intr_effects,beta_intr,intr_vars,weights,lambdas,alpha,r,maxiter_inner);

            residuals = calculate_residuals(X,Y,intercept,main_effects,beta_main,intr_effects,beta_intr,intr_vars,r);

            clean_all_effects(main_effects,beta_main,intr_effects,beta_intr,r);
    }
    return List::create(main_effects,beta_main,intr_effects,beta_intr,lambdas,intercept);
}
