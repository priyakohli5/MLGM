#' Least Squares Estimates for Mean and Covariance in Multivariate Longitudinal Data
#'
#' @description mvlsfit returns least squares estimates for the mean and covariance using an iterative process.
#'
#' @param data a data frame (or matrix) with n rows for subjects and T columns for the repeated measurements.
#' @param time a vector with T equally or unequally spaced time points.
#' @param j a positive integer for the number of outcomes.
#' @param p a positive integer for the order of time polynomial used for average or the number of covariates.
#' @param r a positive integer for the order of time-lag polynomial in constructing design matrix for regressogram.
#' @param s a positive integer for the order of time polynomial in constructing design matrix for innovariogram.
#' @param lags a vector with lags corresponding to all time points.
#' @param X.mat design matrix for modeling the mean.
#' @param U design matrix for modeling the regressogram elements.
#' @param V design matrix for modeling the log innovariogram elements.
#' @param beta.init a vector with initial values for the mean parameter.
#' @param gamma.init a matrix with rj rows and j columns with initial estimates for parameters in the regressogram models.
#' @param alpha.init a vector with initial estimates for parameters in the log innovariogram models.
#' @param beta a vector with initial values for the mean parameters. To initial the algorithm the default is a vector of zeros.
#' @param tol a double or real number indicating the tolerance to decide convergence of mean and covariance parameters. The default is 0.001.
#' @param Nmax a positive integers indicating the maximum number of iterations for the algorithm if it does not converge before. The default is 500.
#'
#' @return estimates for mean and covariance parameters in matrix form, estimated mean vector and covariance matrix, and number of iterations. Specifically,
#' \itemize{
#'   \item output is the matrix with least square estimates for the mean and covariance parameters and their standard errors.
#'   \item mean is the estimated mean vector for the multivariate longitudinal data.
#'   \item sigma is the estimated covariance matrix for the multivariate longitudinal data.
#'   \item gamma.mat has the estimated regressogram parameters in a matrix form.
#'   \item alpha.mat has the estimated log innovariogram oarameters in a matrix format.
#'   \item Nmax is the number of iterations.
#' }
#'
#' @usage mvlsfit(data,time,j,p,r,s,lags,X.mat,U,V,beta.init,gamma.init,alpha.init,beta=c(rep(10,length(beta.init))),tol=0.001,Nmax=100)
#'
#' @references Kohli, P. Garcia, T. and Pourahmadi, M. 2016 Modeling the Cholesky Factors of Covariance Matrices of Multivariate Longitudinal Data, Journal of Multivariate Analysis, 145, 87-100.
#'
#' @references Kohli, P. Du, X. and Shen, H. 2020+ Multivariate Longitudinal Graphical Models (MLGM): An R Package for Visualizing and Modeling Mean \& Dependence Patterns in Multivariate Longitudinal Data, submitted.
#'
#' @export
#'
#' @import gdata
#'
#' @examples data(Tcells)
#' time <- c(0, 2, 4, 6, 8, 18, 24, 32, 48, 72)
#' t <- length(time)
#' j <- 4
#' n <- 44
#' p <- 4
#' r <- 3
#' s <- 2
#' X <- poly(time,degree=p)
#' X <- cbind(rep(1,t),X)
#' X.mat <- matrix(0,nrow=1,ncol=(p+1)*j)
#' for(i in 1:t){
#'  X.mat <- rbind(X.mat,kronecker(diag(j),t(X[i,])))
#' }
#' X.mat <- X.mat[-1,]
##removing 15th and 20th columns of X.mat
#' X.mat <- X.mat[,-c(15,20)]
#' design_mat <- mvrmat(time,j,r,s)
#' gene.names <- c("FYB", "CD69", "IL2RG", "CDC2")
#' degree.mean <- c(4,4,3,3)
#' beta.init <- 0
#' for(i in 1:j){
#'  data.gene <- Tcells[,seq(i, ncol(Tcells), j)]
#'  mean.gene <- apply(data.gene,2,mean)
#'  model <-  polyfit(mean.gene,time,degree=degree.mean[i],plot=FALSE,title="",xlabel="Time points",ylabel="Expression Response",indicator=FALSE, intercept=TRUE,pch.plot=19,lwd.fit=2,lty.fit=2)
#'  nc <- degree.mean[i] + 1
#'  beta.init <- c(beta.init,model$output$coefficients[1:nc])
#' }
#' beta.init <- beta.init[-1]
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
#' results <- mvlsfit(Tcells,time,j,p,r,s,lags,X.mat,design_mat$U,design_mat$V,beta.init,gamma.init,alpha.init,beta=matrix(0,nrow=length(beta.init),ncol=1),tol=0.001,Nmax=500)
#' par.est <- results$output
#' mean.est <- results$mean
#' cov.est <- results$sigma

mvlsfit <- function(data,time,j,p,r,s,lags,X.mat,U,V,beta.init,gamma.init,alpha.init,beta=matrix(0,nrow=length(beta.init),ncol=1),tol=0.001,Nmax=500){
  data <- as.matrix(data)
  Nsub <- nrow(data)
  n <- ncol(data)
  t <- length(time)
  beta.new <- matrix(beta.init,ncol=1)
  mean.init <- apply(data,2,mean)
  sigma.init <- cov(data)
  sigma <- sigma.init
  gamma <- unmatrix(gamma.init,byrow=TRUE)
  alpha <- alpha.init

  N <- 0
  beta0 <- 1e6
  gamma0 <- 1e6
  alpha0 <- 1e6

  while(max(abs(beta0-beta.new), abs(alpha0-alpha), abs(gamma0-gamma))>tol & N<Nmax){
    alpha0 <- matrix(alpha,ncol=1)
    gamma0 <- matrix(gamma,ncol=1)
    beta0 <- matrix(beta,ncol=1)
    CholElements <- mvchol(sigma,t,j)
    Phi <- CholElements$bigphi
    Phit <- CholElements$Phit
    D <- CholElements$D
    Dt <- CholElements$Dt
    term1 <- 0
    term2 <- 0
    Nsub <- nrow(data)
    for (i in 1:Nsub){
      term1 <- term1 + t(X.mat) %*% X.mat
      term2 <- term2 + t(X.mat) %*% data[i,]
    }
    beta <- solve(term1) %*% term2
    data.mean <- apply(data,2,mean)
    error.beta <- data.mean-(X.mat%*%beta)
    sigma.est <- sum(error.beta^2)/(nrow(data)-length(beta.init))
    if(round(sigma.est,3)==0){
      sigma.est <- 0.0001
    }
    beta.var.mat <- sigma.est * solve(t(X.mat) %*% X.mat)
    beta.var <- diag(beta.var.mat)
    beta.se <- sqrt(round(beta.var,7))
    gamma.mat <- gamma.mat <- (solve(t(U)%*%U))%*%(t(U)%*% Phi)
    gamma <- unmatrix(gamma.mat,byrow=TRUE)
    error.gamma <- Phi-(U%*%gamma.mat)
    sigma.est <- sum(error.gamma^2)/(nrow(data)-length(gamma.init))
    if(round(sigma.est,3)==0){
      sigma.est <- 0.0001
    }
    gamma.var.mat <- sigma.est * solve(t(U)%*%U)
    gamma.var <- unmatrix(gamma.var.mat,byrow=TRUE)
    gamma.var <- abs(round(gamma.var,7))
    gamma.var[gamma.var==0] <- 0.01
    gamma.se <- sqrt(round(gamma.var,7))
    t1 <- 1:t
    n <- (j*(j+1))/2
    n <- j*(j+1)/2
    logD <- matrix(0,nrow=n,ncol=t)
    log.Dt <- array(0,dim=dim(Dt))
    row.dim <-  rep(1:j,times=c(j:1))
    col.dim <- 0
    for(i in 1:j){
      col.dim <- c(col.dim,i:j)
    }
    col.dim <- col.dim[-1]

    for(i in 1:t){
      Dt.mat <- Dt[i,,]
      ev <- eigen(Dt.mat)
      index <- order(ev$values,decreasing=TRUE)
      log.Nt <- diag(log(ev$values[index]))
      Pt <- ev$vectors[,index]
      log.Dt[i,,] <- t(Pt) %*% log.Nt %*% Pt
    }
    for(k in 1:n){
      r <- row.dim[k]
      c <- col.dim[k]
      logD[k,] <-  log.Dt[1:t,r,c]
    }
    logD.mat <- t(logD)
    alpha.mat <- (solve(t(V)%*%V))%*%(t(V)%*% logD.mat)
    alpha <- unmatrix(alpha.mat,byrow=TRUE)
    error.alpha <- logD.mat-(V%*%alpha.mat)
    sigma.est <- sum(error.alpha^2)/(nrow(data)-length(alpha.init))
    if(round(sigma.est,3)==0){
      sigma.est <- 0.0001
    }
    alpha.var.mat <- sigma.est * solve(t(V)%*%V)
    alpha.var <- unmatrix(alpha.var.mat,byrow=TRUE)
    alpha.var <- abs(round(alpha.var,7))
    alpha.mat[alpha.mat==0] <- 0.01
    alpha.se <- sqrt(round(alpha.var,7))
    alpha.se <- rep(alpha.se,n)
    EstCov <- mvcovest(time,j,gamma.mat,alpha.mat,U,V)
    sigma <- EstCov$Sigma
    N <- N+1
  }
  mean <- X.mat %*% beta
  output <- matrix(0,nrow=2,ncol=sum(c(length(beta),length(gamma),length(alpha))))
  colnames(output) <- c(rep("beta",length(beta)),rep("gamma",length(gamma)),
                        rep("alpha",length(alpha)))
  rownames(output) <- c("estimates","est.se")
  output["estimates",] <- c(beta,gamma,alpha)
  output["est.se",] <- c(beta.se, gamma.se[1:length(gamma)], alpha.se[1:length(alpha)])
  list(output=round(output,digits=4),mean=mean,sigma=sigma,gamma.mat=gamma.mat,alpha.mat=alpha.mat,Nmax=Nmax)
}
