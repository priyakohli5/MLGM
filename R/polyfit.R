#' Fit Polynomial Functions of Time-lag and Time
#'
#' @description polyfit returns results of fitting polynomial or indicator functions of time or time-lag function.
#'
#' @param data a data frame (or matrix) with n rows for subjects and T columns for the repeated measurements.
#' @param time a vector with T equally or unequally spaced time points.
#' @param degree positive integer indicating the order of the polynomial function. The default is 2.
#' @param plot a logical indicating whether observed and fitted values be plotted or not. The default is TRUE.
#' @param title a character string indicating title of the profile plot. The default is blank.
#' @param xlabel a character string indicating label for the x-axis. The default is blank.
#' @param ylabel a character string indicating label for the y-axis. The default is blank.
#' @param indicator a logical indicating whether an indicator function of time will be used or not. The default is FALSE.
#' @param index a positive integer indicating the time point(s) or lag at which indicator should be 1. It is 0 for all other time points. The default is 1.
#' @param intercept a logical indicating whether intercept term should be included in the indicator model or not. The default is TRUE.
#' @param timepch.plot a integer indicating type of symbols to be used in  plot. The default is 19 for a solid dot.
#' @param lwd.fit integer for line width of the estimated average. The default is 2.
#' @param lty.fit integer for line width of the estimated average. The default is 2.
#'
#' @return Plot of observed and fitted model id plot is TRUE and summary of the fitted linear model. More specifically,
#' \itemize{
#'   \item output is the summary of the fitted linear model.
#'   \item AIC is the Akaike's information criterion for the fitted model.
#'   \item BIC is the Baye's information criterion for the fitted model.
#'   \item fitted is a vector of estimated values.
#' }
#'
#' @usage polyfit(data,time,degree=2,plot=TRUE,title="",xlabel="",ylabel="",indicator,index=1,intercept="TRUE",pch.plot=19,lwd.fit=2,lty.fit=2)
#'
#' @export
#'
#' @examples data(Tcells)
#' time <- c(0, 2, 4, 6, 8, 18, 24, 32, 48, 72)
#' j <- 4
#' n <- 44
#' gene.names <- c("FYB", "CD69", "IL2RG", "CDC2")
#' degree.order <- c(4,4,3,3)
#' par(mfrow=c(2,2))
#' for(i in 1:j){
#'  data.gene <- Tcells[,seq(i, ncol(Tcells), j)]
#'  mean.gene <- apply(data.gene,2,mean)
#'  polyfit(mean.gene,time,degree=degree.order[i],plot=TRUE,title=gene.names[i],xlabel="Time points",ylabel="Expression Response",indicator=FALSE,index=1,intercept=TRUE,pch.plot=19,lwd.fit=2,lty.fit=2)
#'  legend("bottomright",legend=c("observed","fitted"),col=c("black","red"),lty=c(1,2),lwd=c(2,2),cex=0.75)
#'}


polyfit <- function(data,time,degree,plot=TRUE,title="",xlabel="",ylabel="",indicator=FALSE,index=1,intercept=TRUE,pch.plot=19,lwd.fit=2,lty.fit=2){
  if(indicator==FALSE){
    poly.fit <- lm(data~poly(time,degree=degree))
  }else{
    n <- length(index)
    if(n==1){
      X.Indi <- rep(0,length(time))
      index <- which(abs(data)==max(abs(data)))
      X.Indi[index] <- 1
    }
    if(n>1){
      X.Indi <- rep(0,length(time))
      for(i in 1:n){
        X.Indi[which(time==index[i])]<-1
      }
    }
    if(intercept=="TRUE"){
      poly.fit <- lm(data~X.Indi)
    }else{
      poly.fit <- lm(data~X.Indi-1)
    }
  }
  output <- summary(poly.fit)
  AIC <- AIC(poly.fit)
  BIC <- BIC(poly.fit)
  fitted <- poly.fit$fitted.values
  if(plot==TRUE){
    plot(time,data,ylab=ylabel,xlab=xlabel,pch=pch.plot,main=title,ylim=c(min(data)-1,max(data)+1))
    par(new=TRUE)
    plot(time,fitted,type="l",lwd=lwd.fit,col="red",lty=lty.fit,xaxt="n",yaxt="n",main="",xlab="",ylab="",ylim=c(min(data)-1,max(data)+1))
  }
  list(output=output,AIC=AIC,BIC=BIC,fitted=fitted)
}

