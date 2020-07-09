#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector rdirichlet_cpp(NumericVector alphas){
  NumericVector y(alphas.size());
  for(int i=0; i<alphas.size(); i++){
    y[i] = rgamma(1,alphas[i],1)[0];
  }
  return(y/sum(y));
}

// [[Rcpp::export]]
List cpp_gibbs_collapsed_pem_type(IntegerVector arm, // 0 or 1
                             NumericMatrix time,
                             NumericMatrix event,
                             IntegerVector tumor_type, // 0,1,2,3,...
                             int T, // number of tumor types
                             int k, // number of PD-L1 classes
                             int H, // number of piecewise exponential cut-points
                             IntegerVector classes,
                             IntegerMatrix X,
                             NumericVector Nclass,
                             NumericVector alpha,
                             NumericMatrix shape,
                             NumericMatrix rate,
                             int NITER,
                             bool verbose){
  // Vectors to save the generated values
  int NARMS = 2;
  NumericMatrix pi(NITER,k);  
  NumericMatrix CL(NITER,classes.size());
  NumericVector lambda(NITER*NARMS*k*H*T);
  // Initialization of the person-time and event indicators
  double pt[NARMS][k][H][T];
  double ev[NARMS][k][H][T];
  double pt_temp;
  double ev_temp;
  int index;
  memset(pt, 0, sizeof(double)*NARMS*k*H*T);
  memset(ev, 0, sizeof(double)*NARMS*k*H*T);
  for(int i=0; i<classes.size(); i++){
    for(int j=0; j<H; j++){
      pt[ arm[i] ][ classes[i] ][ j ][ tumor_type[i] ] += time(i,j);
      ev[ arm[i] ][ classes[i] ][ j ][ tumor_type[i] ] += event(i,j);
    }
  }
  
  // Runs the collapsed data-augmentation Gibbs sampler
  IntegerVector frame = seq_len(k)-1;
  for(int iter=0; iter<NITER; iter++){
    if(verbose) Rcout << "Iteration " << iter+1 << " of " << NITER << "\n";
    // impute the membership indicators
    for(int i=0; i<classes.size(); i++){
      // removes i-th subject from sufficient statistics
      Nclass[classes[i]] -= 1;
      for(int h=0; h<H; h++){
        pt[ arm[i] ][ classes[i] ][ h ][ tumor_type[i] ] -= time(i,h);
        ev[ arm[i] ][ classes[i] ][ h ][ tumor_type[i] ] -= event(i,h);
      }
      // Computes the predictive probabilities
      NumericVector lpr(k); 
      NumericVector lmult(k);
      for(int j=0; j<k; j++){
        lpr[j] = X(i,j)*(alpha[j]+Nclass[j]);
        lmult[j] = 0;
        for(int h=0; h<H; h++){
          ev_temp = ev[ arm[i] ][ j ][ h ][ tumor_type[i] ];
          pt_temp = pt[ arm[i] ][ j ][ h ][ tumor_type[i] ];
          lmult[j] += event(i,h)*log(shape(j,h)+ ev_temp);
          lmult[j] -= (shape(j,h)+ev_temp+event(i,h))*log(rate(j,h)+pt_temp+time(i,h));
          lmult[j] += (shape(j,h)+ev_temp)*log(rate(j,h)+pt_temp);
        }
      }
      lpr = log(lpr) + lmult;
      lpr = lpr - max(lpr) - log(sum(exp(lpr-max(lpr))));
      // assigns new memebership indicator
      IntegerVector samp;
      NumericVector probs;
      probs = exp(lpr);
      samp = RcppArmadillo::sample(frame,1,TRUE,probs);
      classes[i] = samp[0]; 
      // updates sufficient statistics
      Nclass[classes[i]] += 1;
      for(int h=0; h<H; h++){
        pt[ arm[i] ][ classes[i] ][ h ][ tumor_type[i] ] += time(i,h);
        ev[ arm[i] ][ classes[i] ][ h ][ tumor_type[i] ] += event(i,h);
      }
    }
    // Saves the imputed classes
    CL(iter,_) = classes;
    // updates the class probabilities
    pi(iter,_) = rdirichlet_cpp(alpha + Nclass); 
    // Updates the survival distribution parameters
    // FIX THIS
    for(int a=0; a<NARMS; a++){
      for(int j=0; j<k; j++){
        for(int h=0; h<H; h++){
          for(int t=0; t<T; t++){
            index = iter + NITER*a + NITER*NARMS*j + NITER*NARMS*k*h + NITER*NARMS*k*H*t;
            lambda( index ) = rgamma(1,
                    shape(j,h) + ev[ a ][ j ][ h ][ t ],
                    1/(rate(j,h) + pt[ a ][ j ][ h ][ t ]))[0];
          }
        }
      }
    }
  }  
  // Output the generated values
  List out;
  out["pi"] = pi;
  out["lambda"] = lambda;
  out["classes"] = CL;
  return out;
}


