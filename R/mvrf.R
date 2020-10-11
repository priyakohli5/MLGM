#' Multivariate Regressograms with Estimated Values
#'
#' @description mvrf returns multivariate regressograms for multivariate longitudinal data with j outcomes. It includes the estimated values for the multivariate regressograms.
#'
#' @param data a data frame (or matrix) with n rows for subjects and T columns for the repeated measurements.
#' @param time vector with T equally or unequally spaced time points.
#' @param j a positive integer for the number of outcomes.
#' @param N a positive integer for the number of subjects.
#' @param inno a logical indicating if elements of innovation variance matrices be used innovariogram plot. The default is FALSE.
#' @param inverse a logical indicating if elements of inverse innovation variance matrices be used innovariogram plot. The default is FALSE.
#' @param loginno a logical indicating if elements of log innovation variance matrices be used innovariogram plot. The default is TRUE.
#' @param sigma.fitted  Tj X Tj, positive definite estimate of the covariance matix for multivariate longitudinal data with j outcomes.
#' @param mcol color option for the average response. The default is black.
#' @param fcol color option for the estimated average response. The default is red.
#' @param pch.plot a integer indicating type of symbols to be used in  multivariate regressograms. The default is 19 for a solid dot.
#' @param par1.r a positive integer indicating number of rows in multiple regressogram plots. The default is 2.
#' @param par2.r a positive integer indicating number of columns in multiple regressogram plots. The default is 2.
#' @param par1.d a positive integer indicating number of rows in multiple innovariogram plots. The default is 2.
#' @param par2.d a positive integer indicating number of columns in multiple innovariogram plots. The default is 2.
#' @param lwd.fit integer for line width of the estimated values. The default is 2.
#' @param lty.fit integer for line width of the estimated values. The default is 2.
#'
#' @usage mvrf(data,time,j,N,inno=FALSE,inverse=FALSE,loginno=TRUE,sigma.fitted,mcol="black",fcol="red",pch.plot=19,par1.r=2,par2.r=2,par1.d=2,par2.d=2,lwd.fit=2,lty.fit=1)
#'
#'@references Kohli, P. Garcia, T. and Pourahmadi, M. 2016 Modeling the Cholesky Factors of Covariance Matrices of Multivariate Longitudinal Data, Journal of Multivariate Analysis, 145, 87-100.
#' @references Kohli, P. Du, X. and Shen, H. 2020+ Multivariate Longitudinal Graphical Models (MLGM): An R Package for Visualizing and Modeling Mean \& Dependence Patterns in Multivariate Longitudinal Data, submitted.
#'
#' @return Multivariate regressograms with fitted values are returned and following elements of modified Cholesky block decomposition:
#' \itemize{
#'  \item Phit are the correlation coefficient matrices obtained from sample covariance matrix.
#'  \item Phi.plot the elements from jXj Phit matrices obtained from sample covariance matrix. This is in a list format where each element of the list represents elements for each regressogram plot.
#'  \item Phi.fitted.plot the elements from jXj Phit matrices obtained from estimated covariance matrix. This is in a list format where each element of the list represents the estimated elements for each regressogram plot.
#'  \item D.elements are the innovation variance matrices obtained from sample covariance matrix. This is in matrix format where each row represents the elements for each innovariogram plot. These elements are not plotted by default.
#'  \item D.fitted.elements are the innovation variance matrices obtained from estimated covariance matrix. This is in matrix format where each row represents the estimated elements for each innovariogram plot. These elements are not plotted by default.
#'  \item Dinv.elements are the inverse innovation variance matrices obtained from sample covariance matrix. This is in matrix format where each row represents the elements for each innovariogram plot. These elements are not plotted by default.
#'  \item Dinv.fitted.elements are the inverse innovation variance matrices obtained from estimated covariance matrix. This is in matrix format where each row represents the estimated elements for each innovariogram plot. These elements are not plotted by default.
#'  \item logD.elements are the log innovation variance matrices obtained from sample covariance matrix. This is in matrix format where each row represents the elements for each innovariogram plot. These elements are plotted by default.
#'  \item logD.fitted.elements are the log innovation variance matrices obtained from estimated covariance matrix. This is in matrix format where each row represents the estimated elements for each innovariogram plot. These elements are plotted by default.
#' }
#'
#' @export
#'
#' @examples data(Tcells)
#' time <- c(0, 2, 4, 6, 8, 18, 24, 32, 48, 72)
#' j <- 4
#' N <- 44
#' data(Cov_est)
#' sigma.fitted <- Cov_est
#' MVRF <- mvrf(Tcells,time,j,N,inno=FALSE,inverse=FALSE,loginno=TRUE,as.matrix(Cov_est),mcol="black",fcol="red",pch.plot=19,par1.r=4,par2.r=4,par1.d=4,par2.d=3,lwd.fit=2,lty.fit=1)

mvrf <- function(data,time,j,N,inno=FALSE,inverse=FALSE,loginno=TRUE,sigma.fitted,mcol="black",fcol="red",pch.plot=19,par1.r=2,par2.r=2,par1.d=2,par2.d=2,lwd.fit=2,lty.fit=1){
  sigma.sample <- cov(data)
  t <- length(time)
  CholElements <- mvchol(sigma.sample,t,j)
  CholElements.fitted <- mvchol(sigma.fitted,t,j)
  Phit <- CholElements$Phit
  Phit.fitted <- CholElements.fitted$Phit
  Phi.lags.all <-   Phi.fitted.lags.all <- list(0)
  Lag.all <- Lag.fitted.all <- list(0)
  row.dim <- rep(1:j,each=j)
  col.dim <- rep(1:j,j)
  for(i in 1:(t-1)){
    phi.lags <- array(0,dim=c((t-i),j,j))
    phi.fitted.lags <- array(0,dim=c((t-i),j,j))
    count <- 0
    for(b in (i+1):t){
      count <- count + 1
      m <- (b-i)
      phi.lags[count,,] <- Phit[b,(j*(m-1)+1):(j*m),]
      phi.fitted.lags[count,,] <- Phit.fitted[b,(j*(m-1)+1):(j*m),]
    }
    Phi.lags.all[[i]] <- phi.lags
    Phi.fitted.lags.all[[i]] <- phi.fitted.lags
    lag.elements <- lag.fitted.elements <- matrix(0,nrow=(j^2),ncol=count)

    for(k in 1:(j^2)){
      r <- row.dim[k]
      c <-	 col.dim[k]
      lag.elements[k,] <-  Phi.lags.all[[i]][1:count,r,c]
      lag.fitted.elements[k,] <- Phi.fitted.lags.all[[i]][1:count,r,c]
    }
    Lag.all[[i]] <- lag.elements
    Lag.fitted.all[[i]] <- lag.fitted.elements
  }

  combine <- function(y,x,t){
    form <- 0
    for(i in 1:(t-1)){
      form <- c(form,x[[i]][y,])
    }
    form <- form[-1]
    return(form)
  }

  Phi.plot <- lapply(1:(j^2),combine,Lag.all,t)
  Phi.fitted.plot <- lapply(1:(j^2),combine,Lag.fitted.all,t)

  ylab.all <- 0
  count2 <- 0
  main.all <- 0
  for(i in 1:j){
    for(k in 1:j){
      count2 <- count2 + 1
      num <- c(i,k)
      num_i <- c(i)
      num_k <- c(k)
      ylab.all <- c(ylab.all,bquote(phi[paste(.(i),.(k))]))
      main.all <- c(main.all,letters[count2])
    }
  }
  ylab.all <- ylab.all[-1]
  main.all <- main.all[-1]

  lags <- 0
  for(i in 1:(t-1)){
    time.lags <- diff(time,i)
    lags <- c(lags,rep(time.lags))
  }
  lags <- lags[-1]

  par(mfrow=c(par1.r,par2.r))
  for(i in 1:(j^2)) {
    plot(lags,Phi.plot[[i]],xlab="Time Lags",main=main.all[i],ylab=ylab.all[i],ylim=c(min(Phi.plot[[i]])-1,max(Phi.plot[[i]])+1),
         pch=pch.plot,col=mcol)
    par(new=TRUE)
    plot(lags,Phi.fitted.plot[[i]],lwd=lwd.fit,type="l",lty=lty.fit,col=fcol,xlab="",ylab="",ylim=c(min(Phi.plot[[i]])-1,max(Phi.plot[[i]])+1),main="")
  }

  D <- CholElements$D
  Dt <- CholElements$Dt
  D.fitted <- CholElements.fitted$D
  Dt.fitted <- CholElements.fitted$Dt
  Dt.inv <- array(0,dim=c(t,j,j))
  Dt.fitted.inv <- array(0,dim=c(t,j,j))

  for (i in 1:t){
    Dt.inv[i,,] <- solve(Dt[i,,])
    Dt.fitted.inv[i,,] <- solve(Dt.fitted[i,,])
  }
  t1 <- 1:t
  n <- j*(j+1)/2
  D.elements <- matrix(0,nrow=n,ncol=t)
  Dinv.elements <- matrix(0,nrow=n,ncol=t)
  D.fitted.elements <- matrix(0,nrow=n,ncol=t)
  Dinv.fitted.elements <- matrix(0,nrow=n,ncol=t)
  row.dim <-  rep(1:j,times=c(j:1))
  col.dim <- 0
  for(i in 1:j){
    col.dim <- c(col.dim,i:j)
  }
  col.dim <- col.dim[-1]
  for(k in 1:n){
    r <- row.dim[k]
    c <- col.dim[k]
    D.elements[k,] <-  Dt[1:t,r,c]
    Dinv.elements[k,] <-  Dt.inv[1:t,r,c]
    D.fitted.elements[k,] <-  Dt.fitted[1:t,r,c]
    Dinv.fitted.elements[k,] <-  Dt.fitted.inv[1:t,r,c]
  }

  if(inno==TRUE){
    count2 <- 0
    main.all <- 0
    ylab.all <- 0
    for(i in 1:j){
      for(k in 1:j){
        count2 <- count2 + 1
        num_i <- c(i)
        num_k <- c(k)
        main.all <- c(main.all,letters[count2])
        ylab.all <- c(ylab.all,bquote(d[paste(.(i),.(k))]))
      }
    }
    main.all <- main.all[-1]
    ylab.all <- ylab.all[-1]
    lags <- 0
    innoplot.list <- list(0)
    par(mfrow=c(par1.d,par2.d))
    for(i in 1:n){
      plot(time,D.elements[i,],xlab="Time Points",ylab=ylab.all[i],ylim=c(min(D.elements[i,])-1,max(D.elements[i,])+1),
           main=main.all[i],pch=pch.plot,col=mcol)
      par(new=TRUE)
      plot(time,D.fitted.elements[i,],lwd=lwd.fit,type="l",lty=lty.fit,col=fcol,xlab="",ylab="",main="",ylim=c(min(D.elements[i,])-1,max(D.elements[i,])+1))
    }
  }

  if(inverse==TRUE){
    count2 <- 0
    ylab.all <- 0
    main.all <- 0
    for(i in 1:j){
      for(k in 1:j){
        count2 <- count2 + 1
        num_i <- c(i)
        num_k <- c(k)
        main.all <- c(main.all,letters[count2])
        ylab.all <- c(ylab.all,bquote(m[paste(.(i),.(k))]))
      }
    }
    main.all <- main.all[-1]
    ylab.all <- ylab.all[-1]
    par(mfrow=c(par1.d,par2.d))
    for(i in 1:n){
      plot(time,Dinv.elements[i,],xlab="Time Points",ylab=ylab.all[i],ylim=c(min(Dinv.elements[i,])-1,max(Dinv.elements[i,])+1),
           main=main.all[i],pch=pch.plot,col=mcol)
      par(new=TRUE)
      plot(time,Dinv.fitted.elements[i,],lwd=lwd.fit,type="l",lty=lty.fit,col=fcol,xlab="",ylab="",main="",ylim=c(min(Dinv.elements[i,])-1,max(Dinv.elements[i,])+1))
    }
    list(Phi.plot=Phi.plot,Phi.fitted.plot=Phi.fitted.plot,Dinv.elements=Dinv.elements,Dinv.fitted.elements=Dinv.fitted.elements)

  }

  if(loginno==TRUE){
    logD.elements <- matrix(0,nrow=n,ncol=t)
    log.Dt <- array(0,dim=dim(Dt))
    logD.fitted.elements <- matrix(0,nrow=n,ncol=t)
    log.Dt.fitted <- array(0,dim=dim(Dt))
    for(i in 1:t){
      Dt.mat <- Dt[i,,]
      Dt.fitted.mat <- Dt.fitted[i,,]
      ev <- eigen(Dt.mat)
      ev.fitted <- eigen(Dt.fitted.mat)
      index <- order(ev$values,decreasing=TRUE)
      index.fitted <- order(ev.fitted$values,decreasing=TRUE)
      log.Nt <- diag(log(ev$values[index]))
      log.Nt.fitted <- diag(log(ev.fitted$values[index]))
      Pt <- ev$vectors[,index]
      Pt.fitted <- ev.fitted$vectors[,index.fitted]
      log.Dt[i,,] <- t(Pt) %*% log.Nt %*% Pt
      log.Dt.fitted[i,,] <- t(Pt.fitted) %*% log.Nt.fitted %*% Pt.fitted
    }
    for(k in 1:n){
      r <- row.dim[k]
      c <- col.dim[k]
      logD.elements[k,] <-  log.Dt[1:t,r,c]
      logD.fitted.elements[k,] <- log.Dt.fitted[1:t,r,c]
    }
    count2 <- 0
    main.all <- 0
    ylab.all <- 0
    for(i in 1:j){
      for(k in 1:j){
        count2 <- count2 + 1
        num_i <- c(i)
        num_k <- c(k)
        main.all <- c(main.all,letters[count2])
        ylab.all <- c(ylab.all,bquote(d[paste(.(i),.(k))]))
      }
    }
    main.all <- main.all[-1]
    ylab.all <- ylab.all[-1]
    lags <- 0
    par(mfrow=c(par1.d,par2.d))
    for(i in 1:n){
      loginnoplot.list<- plot(time,logD.elements[i,],xlab="Time Points",ylab=ylab.all[i],ylim=c(min(logD.elements[i,])-1,max(logD.elements[i,])+1),
                              main=main.all[i],pch=pch.plot,col=mcol)
      par(new=TRUE)
      plot(time,logD.fitted.elements[i,],lwd=lwd.fit,type="l",lty=lty.fit,col=fcol,xlab="",ylab="",main="",ylim=c(min(logD.elements[i,])-1,max(logD.elements[i,])+1))
    }
    list(Phi.plot=Phi.plot,Phi.fitted.plot=Phi.fitted.plot,logD.elements=logD.elements,logD.fitted.elements=logD.fitted.elements)
  }else{
    list(Phi.plot=Phi.plot,Phi.fitted.plot=Phi.fitted.plot,D.elements=D.elements,D.fitted.elements=D.fitted.elements)
  }
}
