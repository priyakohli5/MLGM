#' Profile Plot for a Longitudinal Outcome with it's Estimated Mean
#'
#' @description mvpf returns profile plot for longitudinal outcome. It includes the observed and estimated mean values.
#'
#' @param data a data frame (or matrix) with n rows for subjects and T columns for the repeated measurements.
#' @param time vector with T equally or unequally spaced time points.
#' @param mean.fitted vector with estimated average values.
#' @param title a character string indicating title of the profile plot. The default is blank.
#' @param xlabel a character string indicating label for the x-axis. The default is blank.
#' @param ylabel a character string indicating label for the y-axis. The default is blank.
#' @param scol color option for lines representing subjects. The default is gray.
#' @param mcol color option for the average response. The default is black.
#' @param fcol color option for the estimated average response. The default is red.
#' @param lwd.mean integer for line width of the average. The default is 2.
#' @param lwd.fit integer for line width of the estimated average. The default is 2.
#' @param lty.fit integer for line width of the estimated average. The default is 2.
#'
#' \itemize{
#'   \item profile plot with observed and estimated mean values.
#' }
#'
#' @usage mvpf(data,time,mean.fitted,title="",xlabel="",ylabel="",scol="gray", mcol="black",fcol="red",lwd.mean=2,lwd.fit=2,lty.fit=2)
#'
#' @export
#'
#' @examples data(Tcells)
#' time <- c(0, 2, 4, 6, 8, 18, 24, 32, 48, 72)
#' j <- 4
#' n <- 44
#' gene.names <- c("FYB", "CD69", "IL2RG", "CDC2")
#' par(mfrow=c(2,2))
#' for(i in 1:j){
#' data.gene <- Tcells[,seq(i, ncol(Tcells), j)]
#' mean.gene <- apply(data.gene,2,mean)
#'  mvpf(data.gene,time,mean.gene,title=gene.names[i],xlabel="Hours",ylabel="Expression Response",scol="gray", mcol="black",fcol="red",lwd.mean=2,lwd.fit = 3,lty.fit =2)
#' }

mvpf <- function(data,time,mean.fitted,title="",xlabel="",ylabel="",scol="gray", mcol="black",fcol="red",f2col="blue",lwd.mean=2,lwd.fit=2,lty.fit=2){
  plot(time,data[1,],col="gray",type="l",ylab=ylabel,xlab=xlabel,main=title,ylim=range(data))
  for(i in 2:nrow(data)){
    lines(time,data[i,],col=scol)
  }
  lines(time,apply(data,2,mean),col=mcol,lwd=lwd.mean)
  lines(time,mean.fitted,col=fcol,lwd=lwd.fit,lty=lty.fit)
}

