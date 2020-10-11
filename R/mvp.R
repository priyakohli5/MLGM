#' Profile Plot for a Longitudinal Outcome
#'
#' @description mvp returns profile plot for a longitudinal outcome with repeated measurements. This plot shows the change in variable over the repeated measurements for each subject and the average of all subjects.
#'
#' @param data a data frame (or matrix) with n rows for subjects and T columns for the repeated measurements.
#' @param time a vector with T equally or unequally spaced time points.
#' @param mean a logical indicating whether average response should be plotted. The default is TRUE.
#' @param title a character string indicating title of the profile plot. The default is blank.
#' @param xlabel a character string indicating label for the x-axis. The default is blank.
#' @param ylabel a character string indicating label for the y-axis. The default is blank.
#' @param scol color option for lines representing subjects. The default is gray.
#' @param mcol color option for the average response. The default is black.
#' @param plot logical indicating whether profile plot is returned or not. The default is TRUE.
#' @param lwd.mean integer for line width of the average. The default is 2.
#'
#' @return
#' \itemize{
#'   \item mean.data is the average response.
#'   \item profile plot if plot is TRUE.
#' }
#'
#' @usage mvp(data,time,mean=TRUE,title="",xlabel="",ylabel="",scol="gray", mcol="black",plot=TRUE,lwd.mean=2)
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
#'  mvp(Tcells[,seq(i, ncol(Tcells), j)],time,mean=TRUE,title=gene.names[i],xlabel="Time points",ylabel="Expression Response",scol="gray",mcol="black",plot=TRUE,lwd.mean=2)
#' }

mvp <- function(data,time,mean=TRUE,title="",xlabel="",ylabel="",scol="gray", mcol="black",plot=TRUE,lwd.mean=2){
  plot(time,data[1,],col=scol,type="l",ylab=ylabel,xlab=xlabel,main=title,ylim=range(data))
  mean.data <- apply(data,2,mean)
  if(plot==TRUE){
    for(i in 2:nrow(data)){
      lines(time,data[i,],col=scol)
    }
    if (mean==TRUE){
      lines(time,mean.data,col=mcol,lwd=lwd.mean)
    }
  }
  return(mean.data)
}
