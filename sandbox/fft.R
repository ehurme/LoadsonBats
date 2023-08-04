f<-seq(-50,50,by=1e-2)
fc<-0.3
BW<-0.75
par(mfrow=c(2,1))
curve(BP(x,fc=fc,BW=BW,type="p"),-2,2,ylim=c(-0.2,1)
      ,main="Filterweights" ,xlab="fx",ylab="w" )
curve(BP(x,fc=fc,BW=BW,type="s"),add=TRUE,lty=2)
curve(BP(x,fc=fc,BW=BW,type="b"),add=TRUE,lty=3)
curve(BP(x,fc=fc,BW=BW,type="g"),add=TRUE,lty=4)
abline(v=c(fc,fc+BW,fc-BW),lty=3,col="grey")
#thecorrespondingFourier-Transforms
ty<-c("p","s","b","g")
A0<-integrate(BP,fc=fc,BW=BW,type="s",lower=-2,upper=2)$value
plot(NA,NA,xlab="x",ylab="|A|" ,
     main="correspondingconvolutionkernels" ,
     xlim=2*c(-1,1),ylim=c(0,sqrt(2)*A0/(length(f)*BW*min(diff(f)))) )
for(i in 1:length(ty)) {
  FT<-spec.fft(y=BP(f,fc,BW,type=ty[i]))
lines(FT$fx*length(FT$fx)/diff(range(f)),Mod(FT$A),lty=i)
}

WS <- deconvolve(x = df$time, y = df$accZ_mg)
plot(abs(WS$f), WS$S)



# envelope
##noisysignalwithamplitudemodulation
x<-seq(0,1,length.out=2e2)
#originaldata
y<-(abs(x-0.5))*sin(20*2*pi*x)
ye<-base::Re(envelope(y))
#plotresults
plot(x,y,type="l",lwd=1,col="darkgrey",lty=2,ylab="y",main="Spectralfiltering")
lines(x,ye)
legend("bottomright",c("modulated","envelope"),col=c("grey","black"),lty=c(2,1))


##noisysignalwithamplitudemodulation
x<-seq(0,1,length.out=500)
#originaldata
y_org<-(1+sin(2*2*pi*x))*sin(20*2*pi*x)
#overlaysomenoise
y_noise<-y_org+rnorm(length(x),sd=0.2)
#filterthenoisydata
y_filt<-filter.fft(df$accZ_mg,df$time,fc=12,BW=4,n=50)
#plotresults
layout(1)
plot(df$time, df$accZ_mg,type="l",lwd=1,col="darkgrey",lty=2,ylab="y",main="Spectralfiltering")
lines(df$time, df$accZ_mg,lwd=5,col="grey")
lines(df$time,y_filt)
legend("topright",c("org","noisy","filtered"),col=c("grey","darkgrey","black") ,lty=c(1,2,1),lwd=c(5,1,1))

sf <- spec.fft(y = y_filt, x = df$time)

####noisysignalwithamplitudemodulation####
x<-seq(0,3,length.out=1000)
#originaldata #extendedexamplefromenvelopefunction
y<-1*(abs(x-1.5))*sin(10*2*pi*x)+
  ifelse(x>1.5,sin(15*(1+0.25*(x-1.5))*2*pi*x),0)
ye<-base::Re(envelope(y))
par(mfrow=c(2,1),mar=c(1,3.5,3,3),mgp=c(2.5,1,0))
#plotresults
plot(x,y,type="l",lwd=1,col="darkgrey",lty=2,ylab="y",main="OriginalData",xaxt="n",xlab="")
lines(x,ye)
legend("bottomright",c("modulated","envelope"),col=c("grey","black"),lty=c(2,1))
par(mar=c(3.5,3.5,2,0))
wf<-waterfall(y,x,nf=3)
#rasterImage2(x=wf$x,y=wf$fx,z=wf$A
#,ylim=c(0,60))
plot(wf,ylim=c(0,40),main="Waterfall")

wf <- waterfall(y = y_filt, x = df$time, nf = 3)
plot(wf)


x <- 1:4
fft(x)
fft(fft(x), inverse = TRUE)/length(x)

## Slow Discrete Fourier Transform (DFT) - e.g., for checking the formula
fft0 <- function(z, inverse=FALSE) {
  n <- length(z)
  if(n == 0) return(z)
  k <- 0:(n-1)
  ff <- (if(inverse) 1 else -1) * 2*pi * 1i * k/n
  vapply(1:n, function(h) sum(z * exp(ff*(h-1))), complex(1))
}

relD <- function(x,y) 2* abs(x - y) / abs(x + y)
n <- 2^8
z <- complex(n, rnorm(n), rnorm(n))
## relative differences in the order of 4*10^{-14} :
summary(relD(fft(z), fft0(z)))
summary(relD(fft(z, inverse=TRUE), fft0(z, inverse=TRUE)))


samplingFrequency <- 128
timeInterval <- 1/samplingFrequency
x<-c(0:1048575)*timeInterval
inputSignal <- 60000*sin(2*pi*(0.042)*x+4.5)+100000*sin(2*pi*(0.043)*x+3.3)
ffts <- fft(inputSignal)
freq <- c(0:(length(ffts)-1))*samplingFrequency/length(ffts)
plot(freq,abs(ffts)/(length(ffts)/2),type="h",xlim=c(0.042,0.043))

N <- length(ffts)/2
ind <- 1:ceiling(N)
af <- abs(ffts[ind])/N  ## or Mod(ffts[ind])
i <- which(af %in% rev(sort(af))[1:2])   ## find two max frequencies

ff <- ffts[i]/length(ind-1)
Mod(ff)
atan2(Im(ff), Re(ff))
freq[i]
