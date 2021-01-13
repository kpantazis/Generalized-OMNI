#'
#' Computation of coefficients of covariance matrices for Weighted OMNI 
#' @param w  the weight vector
#' @param m  the number of graphs
#' @param s1  the s1-th block-row of the OMNI matrix
#' @param s2  the s2-th block-row of the OMNI matrix
#' @param alpha cumulative weight vector 
#' 
#' @return the coefficient of the covariance matrix among the estimated 
#'         latent positions of s1-th and s2-th graphs
#' 
#' @references The Importance of Being Correlated:Implications of
#'  Dependence in Joint Spectral Inference across Multiple Networks
#'
#' @author Konstantinos Pantazis <kpantazi@umd.edu>
#' @export 

    #' Compute alpha
    #' @param w the weight vector
    #' @param m the number of the graphs
    #' @param s the block-row of the OMNI matrix
    #'
    #' @return a vector of length \eqn{m}
    #'
    alpha <- function(w,m,s){
      summa <- 0
      alpha <- rep(0,m)
      for (k in 1:m){
        summa <- w[s]/(w[s]+w[k]) + summa
      }
      summa <- summa + 1/2
      for (l in 1:m){
        if (s==l){
          alpha[l] <- summa
          }
        else{
          alpha[l] <- w[l]/(w[l]+w[s])
          }
        }
    return(alpha)
    }
    
#' Coefficients 
coeff_weighted <- function(s1,s2,w,m,alpha){
    return(sum((alpha(w,m,s1)-alpha(w,m,s2))^2)/m^2)
}

#Classical OMNI example
# w= c(1,1,1,1)
# m = 4
# coeff_weighted(1,4,w,m,alpha)
