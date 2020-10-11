#' Multivariate Cholesky Block Decomposition
#'
#' @description mvchol returns modified Cholesky block decomposition (MCBD) of the sample covariance matrix for multivariate longitudinal data with j outcomes.
#'
#' @param S Sample covariance matrix of the multivariate longitudinal data.
#' @param t positive integer indicating the number of repeated measurements.
#' @param j positive integer representing the number of outcomes.
#'
#' @return Elements from MCBD of the covariance matrix:
#' \itemize{
#'   \item T.mat is the regression coefificients matrix \eqn{T} matrix in the MCBD.
#'   \item Phit is \eqn{\Phi_{t}} matrix with regression coefficients stacked as \eqn{\phi_{t1}, \phi_{t2}, ..., \phi_{tk}}.
#'   \item bigPhi is \eqn{\Phi} matrix with regression coefficients stacked as \eqn{\Phi_{2}, \Phi_{3}, ..., \Phi_{T}}.
#'   \item D is the innovation variance matrix with block diagonal elements \eqn{D_1, D_2, \ldots , D_t}.
#'   \item Dt are the j*j innovation variance matrices.
#'}
#'
#' @usage mvchol(S,t,j)
#'
#'@references Kohli, P. Garcia, T. and Pourahmadi, M. 2016 Modeling the Cholesky Factors of Covariance Matrices of Multivariate Longitudinal Data, Journal of Multivariate Analysis, 145, 87-100.
#'
#' @export
#'
#' @examples data(Tcells)
#' S <- cov(Tcells)
#' t <- 10
#' j <- 4
#' MCBD <- mvchol(S,t,j)
#' names(MCBD)

mvchol <- function(S,t,j){
  sigma.sample <- S
  phi.step1 <-  phi.step2  <- array(0,dim=c(t,j*(t-1),j))
  phi.step2[1,,] <- phi.step2[1,,] <- matrix(0,nrow=j*(t-1))
  bigphi <- matrix(0,nrow=t*(t-1)*j/2,ncol=j)
  Phit <- array(0,dim=c(t,(t-1)*j,j))
  for (i in 2:t){
    phi.step1[i,1:(j*(i-1)),] <- solve(S[1:(j*(i-1)),1:(j*(i-1))])%*%S[1:(j*(i-1)),(j*(i-1)+1):(j*i)]
    for (k in 1:(i-1)){
      phi.step2[i,(j*(k-1)+1):(j*k),] <- t(phi.step1[i,(j*(k-1)+1):(j*k),])
    }
    bigphi[((i-2)*(i-1)*j/2+1):(j*i*(i-1)/2),] <- phi.step2[i,1:((i-1)*j),]
    Phit[i,1:((i-1)*j),] <- bigphi[(j/2*(i-1)*(i-2)+1):(j/2*i*(i-1)),]
  }
  T.mat <- diag(1,(t*j))
  for(i in 2:t){
    for (k in 1:(i-1)){
      T.mat[((i-1)*j+1):(i*j),((k-1)*j+1):(k*j)] <- -Phit[i,((k-1)*j+1):(k*j),]
    }
  }
  Lbar <- t(chol(S))
  Dbar <- (diag(diag(Lbar)))^2
  L <-  (diag(diag(Lbar))) %*% (solve(Lbar))
  Tbar <- matrix(0,nrow=nrow(S),ncol=ncol(S))
  for (i in 1:t){
    Tbar[(j*(i-1)+1):(j*i),(j*(i-1)+1):(j*i)] <- L[(j*(i-1)+1):(j*i),(j*(i-1)+1):(j*i)]
  }
  D <- (solve(Tbar)) %*% Dbar %*% (solve(t(Tbar)))
  Dt <- array(0,dim=c(t,j,j))
  for (i in 1:t){
    Dt[i,,] <- D[(j*(i-1)+1):(j*i),(j*(i-1)+1):(j*i)]
  }
  list(T.mat=T.mat,Phit=Phit,bigphi=bigphi,D=D,Dt=Dt)
}
