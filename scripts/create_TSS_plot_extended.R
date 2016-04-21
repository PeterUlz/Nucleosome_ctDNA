args <- commandArgs()
data<-read.csv(args[4],header=TRUE,sep="\t")

#define sliding window smoothing function
slideFunct <- function(data, window, step){
  total <- length(data)
  spots <- seq(from=1, to=(total-window), by=step)
  result <- vector(length = length(spots))
  for(i in 1:length(spots)){
    result[i] <- mean(data[spots[i]:(spots[i]+window)])
  }
  return(result)
}

png(filename = paste(args[4],".extended.png",sep=""), width = 6500, height = 1180, units = "px", pointsize = 20, bg = "white", res = NA)
nf <- layout(matrix(c(1,2,3,4,5), 1, 5, byrow=TRUE))

#First plot zoomed in version
plot(data$Position,data$Mean.Cov,type="l",xlab="Position relative to TSS",ylab="Read depth",lwd=5,cex.lab=3,cex.axis=2)
polygon(c(data$Position,rev(data$Position)),c(data$LowerBound,rev(data$UpperBound)),col = "grey75", border = FALSE)
abline(v=0,lty=2)
lines(data$Position,data$Mean.Cov)

#Second plot sliding window (window-size: 10)
data_smooth_10<-slideFunct(data$Mean.Cov,10,1)
plot(data$Position[5:1995],data_smooth_10,type="l",xlab="Position relative to TSS",ylab="Read depth",col="red",lwd=5,cex.lab=3,cex.axis=2)
abline(v=0,lty=2)

#Third Plot sliding window (window-size: 100)
data_smooth_100<-slideFunct(data$Mean.Cov,100,1)
plot(data$Position[50:1950],data_smooth_100,xlab="Position relative to TSS",ylab="Read depth",type="l",col="green",lwd=5,cex.lab=3,cex.axis=2)
abline(v=0,lty=2)

#Fourth Plot sliding window (window-size: 400)
data_smooth_100<-slideFunct(data$Mean.Cov,400,1)
plot(data$Position[200:1800],data_smooth_100,xlab="Position relative to TSS",ylab="Read depth",type="l",col="blue",lwd=5,cex.lab=3,cex.axis=2)
abline(v=0,lty=2)

#Fifth plot using zero as baseline
mean_cov<-mean(data$Mean.Cov)
plot(data$Position,data$Mean.Cov,type="l",ylim=c(0,max(data$Mean.Cov)+5),xlab="Position relative to TSS",ylab="Read depth",lwd=3,cex.lab=3,cex.axis=2)
polygon(c(data$Position,rev(data$Position)),c(data$LowerBound,rev(data$UpperBound)),col = "grey75", border = FALSE)
lines(data$Position,rep(mean_cov,length(data$Position)),ylim=c(0,max(data$Mean.Cov)+5),lwd=7,col="blue",xlab="Position relative to TSS",ylab="Read depth (mean)")
lines(data$Position,data$Mean.Cov)
text(mean(data$Position),mean_cov+1,labels=round(mean_cov,digits=2),cex=3)
abline(v=0,lty=2)



dev.off()


