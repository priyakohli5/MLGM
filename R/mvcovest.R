#' Estimated Covariance Matrix for Multivariate Longitudinal Data
#'
#' @description mvcovest returns estimated covariance matrix estimate using the estimates for the elements of modified Cholesky block decomposition (MCBD).
#'
#' @param time a vector with T equally or unequally spaced time points.
#' @param j a positive integer for the number of outcomes.
#' @param gamma.mat matrix with rj rows and j columns for estimates for parameters in regressogram models. Here r is the order of time-lag polynomial used for regressogram plots.
#' @param alpha.mat Vector with sj(j+1)/2 estimates for parameters in log innovariogram models. Here s is the order of the time polynomial used for log innovariogram plots.
#'
#' @return
#' \itemize{
#'   \item bigT is the estimated regression coefficients matrix in the MCBD.
#'   \item D.mat is the estimated D matrix in the MCBD.
#'   \item Sigma is the estimated covariance matrix.
#'   \item logDt.list block of log innovation matrix in a list format.
#'   \item logDt elements of log innovation matrix in a matrix format.
#' }
#'
#' @usage mvcovest(time,j,gamma.mat,alpha.mat,U,V)
#'
#' @references Kohli, P. Garcia, T. and Pourahmadi, M. 2016 Modeling the Cholesky Factors of Covariance Matrices of Multivariate Longitudinal Data, Journal of Multivariate Analysis, 145, 87-100.
#'
#' @references Kohli, P. Du, X. and Shen, H. 2020+ Multivariate Longitudinal Graphical Models (MLGM): An R Package for Visualizing and Modeling Mean \& Dependence Patterns in Multivariate Longitudinal Data, submitted.
#'
#' @importFrom Matrix bdiag
#'
#' @import expm
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
#' design_mat <- mvrmat(time,j,r,s)
#' gene.names <- c("FYB", "CD69", "IL2RG", "CDC2")
#' MVR <- mvr(Tcells,time,j,n,inno=FALSE,inverse=FALSE,loginno=TRUE,plot=FALSE,pch.plot=19,par1.r = 4,par2.r = 4,par1.d=4,par2.d=4)
#' j2 <- j^2
#' lags <- 0
#' for(i in 1:(t-1)){
#' time.lags <- diff(time,i)
#' lags <- c(lags,rep(time.lags))
#' }
#' lags <- lags[-1]
#' gamma.init <- matrix(0,nrow=j2,ncol=r+1)
#' # par(mfrow=c(4,4))
#' for(i in 1:j2){
#'  data.phi <- MVR$Phi.plot[[i]]
#'  best.model <- polyfit(data.phi,lags,degree=r,plot=FALSE,title=gene.names[i],xlabel="Time points",ylabel="Expression Response",indicator=FALSE,index=1,intercept=TRUE)
#'  gamma.init[i,] <- best.model$output$coefficients[1:(r+1)]
#'  # if wanted to see initial fits
#'  # best.model <- polyfit(data.phi,lags,degree=r,plot=TRUE,title=gene.names[i],xlabel="Time points",ylabel="Expression Response",indicator=FALSE,index=1,intercept=TRUE)
#'  # legend("bottomright",legend=c("observed","fitted"),col=c("black","red"),lty=c(1,2),lwd=c(2,2),cex=0.7)
#' }
#' gamma.mat <- matrix(0,((r+1)*j),j)
#' for(i in 1:j){
#' row <- ((i-1)*j)+1
#' col <- ((i-1)*(r+1)) +1
#' gamma.mat[col:(i*(r+1)),1:j] <- t(gamma.init[row:(i*j),(1:(r+1))])
#' }
#' # par(mfrow=c(4,3))
#' j3 <- j*(j+1)/2
#' alpha.init <- 0
#' for(i in 1:j3){
#'  data.logD <- MVR$logD.elements[i,]
#'  model.best <- polyfit(data.logD,time,degree=s,plot=FALSE,title="",xlabel="Time points",ylabel="Expression Response",indicator=FALSE,index=4,intercept=TRUE)
#'  alpha.init <- c(alpha.init,model.best$output$coefficients[1:(s+1)])
#'  # if wanted to see initial fits
#'  # model.best <- polyfit(data.logD,time,degree=s,plot=TRUE,title="",xlabel="Time points",ylabel="Expression Response",indicator=FALSE,index=4,intercept=TRUE)
#'  # legend("topright",legend=c("observed","fitted"),col=c("black","red"),lty=c(1,2),lwd=c(2,2),cex=0.8)
#' }
#' alpha.init <- alpha.init[-1]
#' alpha.mat <- matrix(alpha.init,nrow=(s+1),ncol=j3)
#' Cov_est <- mvcovest(time,j,gamma.mat,alpha.mat,design_mat$U,design_mat$V)$Sigma

mvcovest <- function(time,j,gamma.mat,alpha.mat,U,V){
  T <- time
  t <- length(time)
  bigPhi <- U %*% gamma.mat
  Phit <- array(0,dim=c(t,(t-1)*j,j)) 	# Phi(2),Phi(3)...where Phi(t=2) is a j*j matrix, Phi(t=3) is a (2j)*j matrix...
  for (m in 2:t){
    Phit[m,1:((m-1)*j),] <- bigPhi[(j/2*(m-1)*(m-2)+1):(j/2*m*(m-1)),]
  }
  bigT <- diag(1,(t*j)) #T*J#
  for(m in 2:t){
    for (k in 1:(m-1)){
      bigT[((m-1)*j+1):(m*j),((k-1)*j+1):(k*j)] <- -Phit[m,((k-1)*j+1):(k*j),]
    }
  }

  n <- j*(j+1)/2
  Delta <- t(V %*% alpha.mat)
  Dpoly.hat <- matrix(0,nrow=n,ncol=t)
  for(i in 1:n){
    Dpoly.hat[i,] <- V %*% alpha.mat[,i]
  }
  logDt<- array(0,dim=c(t,j,j))
  row.dim <- 0
  for(i in 1:j){
    j.l1 <- j-i+1
    row.dim <- c(row.dim,rep(i,j.l1))
  }
  row.dim <- row.dim[-1]
  col.dim <- 0
  for(i in 1:j){
    col.dim <- c(col.dim,i:j)
  }
  col.dim <- col.dim[-1]

  for(k in 1:n){
    r <- row.dim[k]
    c <- col.dim[k]
    logDt[1:t,r,c] <- Dpoly.hat[k,]
  }
  logDt.list <- list(0)
  Dt.list <- list(0)
  for(i in 2:j){
    for(k in 1:(i-1)){
      logDt[1:t,i,k] <- logDt[1:t,k,i]
    }
  }
  for(i in 1:t){
    logDt.list[[i]] <- logDt[i,,]
    Dt.list[[i]] <- expm(logDt.list[[i]])
  }
  D.mat <- as.matrix(bdiag((Dt.list)))
  Sigma <- solve(bigT) %*% (D.mat) %*% t(solve(bigT))
  Sigma <- round(Sigma,7)
  list(bigT=bigT,D.mat=D.mat,Sigma=Sigma,logDt.list=logDt.list,logDt=logDt)
}
