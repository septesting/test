
Bandi5 <- function(x0=x0,dx,nx,DT,bw,na,avec)  {
	require(KernSmooth)
# Set up constants and useful preliminaries
DT <- DT
	SF <- 1/(bw*sqrt(2*pi))  # scale factor for kernel calculation
	x02 <- x0*x0 # second power of x
	x04 <- x0 * x0 * x0 * x0
	dx2 <- dx*dx # second power of dx
	dx3 <- dx2*dx
	dx4 <- dx2*dx2  # fourth power of dx
	dx5 <- dx3*dx2
	dx6 <- dx2*dx4  # sixth power of dx
# Compute matrix of kernel values
	Kmat <- matrix(0,nrow=na,ncol=nx)
	for(i in 1:(nx)) {  # loop over columns (x0 values)
		Kmat[,i] <- SF*exp(-0.5*(x0[i]-avec)*(x0[i]-avec)/(bw*bw))
	}
# Compute M1, M2, M4, moment ratio and components of variance for each value of a
	M1.a <- rep(0,na)
	M2.a <- rep(0,na)
	M3.a <- rep(0,na)
	M4.a <- rep(0,na)
	M5.a <- rep(0,na)
	M6.a <- rep(0,na)
	M6M4r <- rep(0,na)  # vector to hold column kernel-weighted moment ratio
	mean.a <- rep(0,na) # centering of conditional variance
	SS.a <- rep(0,na)  # sum of squares
	for(i in 1:na) {  # loop over rows (a values)
		Ksum <- sum(Kmat[i,])  # sum of weights
		M1.a[i] <- (1/DT)*sum(Kmat[i,]*dx)/Ksum
		M2.a[i] <- 0.5*(1/DT)*sum(Kmat[i,]*dx2)/Ksum
		M3.a[i] <- (1/DT)*sum(Kmat[i,]*dx3)/Ksum
		M4.a[i] <- (1/DT)*sum(Kmat[i,]*dx4)/Ksum
		M5.a[i] <- (1/DT)*sum(Kmat[i,]*dx5)/Ksum
		M6.a[i] <- (1/DT)*sum(Kmat[i,]*dx6)/Ksum
		M6M4r[i] <- M6.a/M4.a[i]
		mean.a[i] <- sum(Kmat[i,]*x0[2:(nx+1)])/Ksum 
#SS.a[i] <- sum(Kmat[i,]*x02[2:(nx+1)])/Ksum 
		SS.a[i] <- sum(Kmat[i,]*x04[2:(nx+1)])/Ksum
#Y4.a[i] <- sum(Kmat[i,]*x04[2:(length(x0))])/Ksum 
#Y2.a[i] <- sum(Kmat[i,]*x02[2:(length(x0))])/Ksum 
		
	}
# Compute conditional variance
#S2.x <- SS.a - (mean.a*mean.a) # sum of squares minus squared mean
	S2.x <- 1/4 * log((SS.a/((sd(x0)^4))*3)) # sum of squares minus squared mean
# Compute jump frequency, diffusion and drift
	sigma2.Z <-(((M6M4r)/(5))) # average the column moment ratios
	lamda.Z <- (M4.a/(3*sigma2.Z*sigma2.Z)) #* sqrt(sigma2.Z)
	sigma2.dx <- (M2.a -(lamda.Z*sigma2.Z))
# set negative diffusion estimates to zero
	diff.a <- ifelse(sigma2.dx>0,sigma2.dx,0)
	sigma2.dx <- M2.a   # total variance of dx
	sigma3.dx <- M3.a
	sigma4.dx <- M4.a
	sigma5.dx <- M5.a
	sigma6.dx <- M6.a
	mu.a <- M1.a
	D3.a <- M3.a 
	outlist <- list(mu.a,sigma2.dx,sigma3.dx,sigma4.dx,sigma5.dx,sigma6.dx,diff.a,sigma2.Z,lamda.Z,S2.x,D3.a)

	return(outlist)
} 




Mehrddj<-function(timeseries,bandwidth=0.3,na=50,DT,logtransform=FALSE,interpolate=FALSE){

	

timeseries<-ts(timeseries) 

d <- as.data.frame(timeseries)
timeseries[is.na(d)] <- 0

if (dim(timeseries)[2]==1 ){


Y=timeseries[,1]
print(sd(Y))
g=sd(Y)
#Y=(Y - mean(Y))/sd(Y)
#Y=(Y - mean(Y))
                                         #for (japan), removed Diff from timeseries
Y=Y
timeindex=1:length(Y)#dim(timeseries)[1]
}else if(dim(timeseries)[2]==2){
Y<-timeseries[,2]
print(sd(Y))
g=sd(Y)
#Y=(Y - mean(Y))/sd(Y)
Y=Y
#Y=(Y - mean(Y))
timeindex<-timeseries[,1]
}else{
warning("not right format of timeseries input")
}

		
	# Interpolation
	if (interpolate){
		YY<-approx(timeindex,Y,n=lenfth(Y),method="linear")
		Y<-YY$y
		}else{
		Y<-Y}
			
	# Log-transformation
	if (logtransform){
		Y<-log(Y+1)}
		
	# Preliminaries
	Xvec1<-Y
	Tvec1<-timeindex
	dXvec1 <- diff(Y)
	DT <- DT
	#DT <- Tvec1[2]-Tvec1[1]
	bw <- bandwidth*sd(Xvec1) # bandwidth 
	alow <- min(Xvec1)
	ahigh <- max(Xvec1)
	na <- na
	avec <- seq(alow,ahigh,length.out=na)
	nx <- length(dXvec1)

	# Bandi-type estimates
	ParEst <- Bandi5(Xvec1,dXvec1,nx,DT,bw,na,avec)
	Drift.vec <- ParEst[[1]]
	TotVar.dx.vec <- ParEst[[2]]
	M3.vec <- ParEst[[3]]
	M4.vec <- ParEst[[4]]
	M5.vec <- ParEst[[5]]
	M6.vec <- ParEst[[6]]
	Diff2.vec <- ParEst[[7]]
	Sigma2Z <- ParEst[[8]]
 	LamdaZ.vec <- DT * ParEst[[9]]
	S2.vec <- ParEst[[10]]
	D3.vec <- ParEst[[11]]


    
	TotVar.i <- approx(x=avec,y=TotVar.dx.vec,xout=Xvec1)
	TotVar.t <- TotVar.i$y
    Drift.i <- approx(x=avec,y=Drift.vec,xout=Xvec1)
    Drift.t <- Drift.i$y     
    #M2.i <- approx(x=avec,y=M2.vec,xout=Xvec1)
	#M2.t <- M2.i$y
	M3.i <- approx(x=avec,y=M3.vec,xout=Xvec1)
	M3.t <- M3.i$y
	M4.i <- approx(x=avec,y=M4.vec,xout=Xvec1)
	M4.t <- M4.i$y
	M5.i <- approx(x=avec,y=M5.vec,xout=Xvec1)
	M5.t <- M5.i$y
	M6.i <- approx(x=avec,y=M6.vec,xout=Xvec1)
	M6.t <- M6.i$y
	Diff2.i <- approx(x=avec,y=Diff2.vec,xout=Xvec1)
	Diff2.t <- Diff2.i$y
		#M2.i <- approx(x=avec,y=M2.vec,xout=Xvec1)
		M3.i <- approx(x=avec,y=M3.vec,xout=Xvec1)
		M4.i <- approx(x=avec,y=M4.vec,xout=Xvec1)
		M5.i <- approx(x=avec,y=M5.vec,xout=Xvec1)
		M6.i <- approx(x=avec,y=M6.vec,xout=Xvec1)
        D3.i <- approx(x=avec,y=D3.vec,xout=Xvec1)
	D3.t <- D3.i$y
	Lamda.i <-  approx(x=avec,y=LamdaZ.vec,xout=Xvec1)
	Lamda.t <- Lamda.i$y
	Lamda.t <- ifelse(Lamda.t < (1/DT), Lamda.t,1)
	LamdaZ.vec <- ifelse(LamdaZ.vec < (1/DT), LamdaZ.vec,1)
    S2.i <- approx(x=avec,y=S2.vec,xout=Xvec1)
	S2.t <- S2.i$y
    

ss <- matrix(NA,nrow=na, ncol=9)
ss[,1] <- avec
ss[,2] <- Drift.vec
ss[,3] <- TotVar.dx.vec
ss[,4] <- M3.vec
ss[,5] <- M4.vec
ss[,6] <- M5.vec
ss[,7] <- M6.vec
ss[,8] <- Diff2.vec
#ss[,10] <- TotVar.dx.vec
ss[,9] <- LamdaZ.vec

write.table(ss, file = "D:/U/Codes/Non-stationary/ss.txt", row.names = FALSE , col.names = FALSE)
#write.table(x,Drift.vec,file = "D:/U/Codes/Non-stationary/5.txt", row.names = FALSE , col.names = FALSE)
#write.table(x,Diff2.vec,file = "D:/U/Codes/Non-stationary/6.txt", row.names = FALSE , col.names = FALSE)
#write.table(x,TotVar.dx.vec,file = "D:/U/Codes/Non-stationary/7.txt", row.names = FALSE , col.names = FALSE)
#write.table(x,LamdaZ.vec,file = "D:/U/Codes/Non-stationary/8.txt", row.names = FALSE , col.names = FALSE)
#write.table(timeseries,file = "D:/U/Codes/Non-stationary/9.txt", row.names = FALSE , col.names = FALSE)
write.table(Xvec1,file = "D:/U/Codes/Non-stationary/data.txt", row.names = FALSE , col.names = FALSE)
write.table(Lamda.t, file = "D:/U/Codes/Non-stationary/jump-rate.txt", row.names = FALSE , col.names = FALSE)
#write.table(S2.t, file = "D:/U/Codes/Non-stationary/13.txt", row.names = FALSE , col.names = FALSE)
write.table(Drift.t, file = "D:/U/Codes/Non-stationary/D1.txt", row.names = FALSE , col.names = FALSE)
write.table(TotVar.t, file = "D:/U/Codes/Non-stationary/D2.txt", row.names = FALSE , col.names = FALSE)
write.table(Sigma2Z, file = "D:/U/Codes/Non-stationary/Sigma2.txt", row.names = FALSE , col.names = FALSE)




    

# Plot the data
	dev.new()
	par(mfrow=c(5,1),mar=c(3, 3, 2, 2),mgp=c(1.5,0.5,0),oma=c(1,1,1,1))
	#plot(Tvec1,Xvec1,type='l',col='black',lwd=2,xlab='',ylab='original data')
	#grid()
	#plot(Tvec1[1:length(Tvec1)-1],dXvec1,type='l',col='black',lwd=2,xlab='time',ylab='first-diff data')
	#grid()
	#plot(Tvec1,Drift.t,type='l',lwd=1,col='red',xlab='time',ylab='local drift')
	#plot(Tvec1,Diff2.t,type='l',lwd=1,col='blue',xlab='time',ylab='local diffusion')
	#plot(Tvec1,Lamda.t,type='l',lwd=1,col='red',xlab='time',ylab='jump rate')
	#plot(Tvec1,S2.t,type='l',lwd=1,col='black',xlab='time',ylab='Non-Gaussian Parameter')

	plot(Tvec1,Xvec1,type='l',col='black',lwd=2,xlab='',ylab='original data')
	plot(Xvec1,TotVar.t,type='l',lwd=1,col='blue',xlab='X',ylab='Diffusion')
	plot(Xvec1,Drift.t,type='l',lwd=1,col='red',xlab='X',ylab='Drift')
	plot(avec,LamdaZ.vec ,type='l',lwd=1,col='red',xlab='x',ylab='jump rate')
	plot(avec,Sigma2Z,type='l',lwd=1,col='green',xlab='x',ylab='Jump intensity')



	
	
	# Plot indicators versus a
	dev.new()
	par(mfrow=c(4,1),mar=c(3, 3, 2, 2) ,cex.axis=1,cex.lab=1,mgp=c(2,1,0),oma=c(1,1,2,1))
	plot(Tvec1,Xvec1,type='l',col='black',lwd=2,xlab='',ylab='original data')
	plot(Tvec1,TotVar.t,type='l',lwd=1,col='blue',xlab='t',ylab='Diffusion')
	plot(Tvec1,Drift.t,type='l',lwd=1,col='red',xlab='t',ylab='Drift')
	plot(Tvec1,Lamda.t ,type='l',lwd=1,col='red',xlab='t',ylab='jump rate')


	#plot(avec,S2.vec,type='l',lwd=1,col='black',xlab='a',ylab='conditional variance')
	#plot(avec,Sigma2Z,type='l',lwd=1,col='green',xlab='a',ylab='Jump intensity 1')
    #plot(avec,TotVar.dx.vec/(g),type='l',lwd=1,col='blue',xlab='a',ylab='diffusion without jump')
    #plot(avec,D3.vec,type='l',lwd=1,col='green',xlab='a',ylab='Skewness')
	#plot(avec,Drift.vec,type='l',lwd=1,col='red',xlab='a',ylab='Drift')
    #plot(avec,Diff2.vec,type='l',lwd=1,col='red',xlab='a',ylab='diffusion')
	#plot(avec,LamdaZ.vec ,type='l',lwd=1,col='red',xlab='a',ylab='jump rate')
	#mtext("Drift and Diffusion Coefficients ",side=3,line=0.1,outer=TRUE)
	 
	 
	 
	 		
		 
	 
	 
	# Plot indicators versus time
	#dev.new()
	#par(mfrow=c(1,1),mar=c(3, 3, 2, 2),cex.axis=1,cex.lab=1,mgp=c(1.5,0.5,0),oma=c(1,1,2,1))
	#plot(Tvec1,S2.t,type='l',lwd=1,col='black',xlab='time',ylab='conditional variance')
	#plot(Tvec1,D3.vec,type='l',lwd=1,col='blue',xlab='time',ylab='Skewness')
	#plot(Tvec1,TotVar.t,type='l',lwd=1,col='blue',xlab='time',ylab='local diffusion')
   # plot(Sigma2Z.t,type='l',lwd=1,col='green',xlab='a',ylab='Jump intensity')
	#plot(Tvec1,Diff2.t,type='l',lwd=1,col='green',xlab='time',ylab='diffusion')
	#plot(Tvec1,Lamda.t,type='l',lwd=1,col='red',xlab='time',ylab='jump intensity')
	#mtext("DDJ nonparametrics versus time",side=3,line=0.1,outer=TRUE)

	# Output
	nonpar_x<-data.frame(avec,S2.vec,TotVar.dx.vec,Diff2.vec,LamdaZ.vec)
	nonpar_t<-data.frame(Tvec1,S2.t,TotVar.t,Diff2.t,Lamda.t)
	#return(c(nonpar_x,nonpar_t))
    
    Sigma2Z=Sigma2Z
    print(g)
    return(Sigma2Z)

	}