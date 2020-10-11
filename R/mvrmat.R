#' Design Matrices for Modeling Elements of Regressograms
#'
#' @description mvrmat returns design matrices for polynomial models used for modeling elements of modified Cholesky block decomposition.
#'
#' @param time a vector with T equally or unequally spaced time points.
#' @param j a positive integer for the number of outcomes.
#' @param r a positive integer for the order of time-lag polynomial in constructing design matrix for regressogram.
#' @param s a positive integer for the order of time polynomial in constructing design matrix for innovariogram.
#'
#' @return design matrices for models selected for regressogram and innovariogram. Specifically,
#' @return
#' \itemize{
#'   \item U is the design matrix for time-lag polynomials selected for regressogram.
#'   \item V is the design matrix for time polynomials selected for log innovariogram.
#' }
#'
#' @usage mvrmat(time,j,r,s)
#'
#' @references Kohli, P. Garcia, T. and Pourahmadi, M. 2016 Modeling the Cholesky Factors of Covariance Matrices of Multivariate Longitudinal Data, Journal of Multivariate Analysis, 145, 87-100.
#'
#' @export
#'
#' @examples data(Tcells)
#' time <- c(0, 2, 4, 6, 8, 18, 24, 32, 48, 72)
#' t <- length(time)
#' j <- 4
#' n <- 44
#' r <- 3
#' s <- 2
#'design_mat <- mvrmat(time,j,r,s)

mvrmat <- function(time,j,r,s){
  t2 <- 0
  t <- length(time)
  for (i in 2:t){
    for(k in 1:(i-1)){
      t2 <- c(t2,(time[i]-time[k]))
    }
  }
  t2 <- t2[-1]
  if(r>1){
    U1_mat <- cbind(rep(1,length(t2)),poly(t2,deg=r))
  }else{
    U1_mat <- cbind(rep(1,length(t2)))
  }
  BigU_mat <- rep(1,(ncol(U1_mat)*j))
  for(i in 1:nrow(U1_mat)){
    BigU_mat <- rbind(BigU_mat,t(kronecker(diag(1,j),U1_mat[i,],FUN="*")))
  }
  BigU <- BigU_mat[-1,]
  V <- poly(time,degree=s)
  V.mat <- cbind(rep(1,t),V[,c(1:s)])
  list(U=BigU,V=V.mat)
}
