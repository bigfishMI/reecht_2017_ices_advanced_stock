##################################
# do simple assessment using TMB #
##################################
require(TMB)

## make sure you create a working directory with the assessment data and the "simple.cpp" file
## then set that as your working directory
setwd("y:/Meeting-Seminar/BIGFISH/")

## ##################################################
## Data/params  prep.:

## read data objects
catches    <- t( as.matrix( read.table("catch_num.dat")))
wts        <- t( as.matrix( read.table("weights.dat")))
tuning     <- t( as.matrix( read.table("tuning.dat")))
mat        <- t( as.matrix( read.table("maturity.dat")))
catch_ton  <- as.matrix(read.table("catch_ton.dat"))
M          <- rep(0.1,times=nrow(catches))

## set the correct dimension names for data objects
dimnames(catches)   <-  list("age"=1:nrow(catches), "year" = 1995:2015)
dimnames(wts)       <- list("age"=1:nrow(wts), "year" = 1995:2015)
dimnames(tuning)    <- list("age"=1:nrow(tuning), "year" = 1995:2015)
dimnames(mat)       <- list("age"=1:nrow(mat), "year" = 1995:2015)
dimnames(catch_ton) <-  list("age"="all", "year" = 1995:2015)

catches[1:5, 1:5]
wts[1:10, 1:5]
tuning[1:10, 1:5]
mat[1:10, 1:10]
M
catch_ton

nages  <- dim(catches)[1]
nyears <- dim(catches)[2]

data <- list(dms=c(nyears,nages),
             catches=catches,
             wts=wts,
             I=tuning,
             mat=mat,
             M=M)

## Create the parameter object:
Fsel <- c(-2,5)
logFyrs   <- seq(-1.1, -0.3, length.out=nyears)
logRec    <- log(rep(800,nyears))
logStartN <- log(rep(400,nages-1))
logq      <- rep(-1, 10)
logsigma  <- rep(1,2)

pars <- list(Fsel=Fsel,
             logFyrs=logFyrs,
             logRec=logRec,
             logStartN=logStartN,
             logq=logq,
             logsigma=logsigma)
## ##################################################
## Linking the C++ model:

## compile it; load it

compile("./simple2.cpp")

dyn.load(dynlib("simple2"))

## create the AD object:
obj <- MakeADFun(data = data, parameters = pars,
                 DLL = "simple2")

## ##################################################
## run it!
system.time(res <- nlminb(obj$par, obj$fn, obj$gr, control=list("iter.max"=200)))


res$par

(AIC <- 2 * res$objective + 2 * length(res$par))

#######################################################
## OK we ran it, but how do we get stuff back out! ####
#######################################################

## generate SD errors using sdreport() method (ADREPORT() variables)
sdrep <- sdreport(obj)

MLE.vec <- sdrep$value
SD.vec <- sdrep$sd
names(SD.vec) <- names(MLE.vec)

str(SD.vec)

head(MLE.vec)
head(SD.vec)

## SSB(y) and F(y) MLE and approximate 95% confidence intervals

ssb.mle <- MLE.vec[names(MLE.vec)=='SSB']
ssb.lq <- ssb.mle-1.96*SD.vec[names(SD.vec)=='SSB']
ssb.uq <- ssb.mle+1.96*SD.vec[names(SD.vec)=='SSB']
Fy.mle <- MLE.vec[names(MLE.vec)=='Fy']
Fy.lq <- Fy.mle-1.96*SD.vec[names(SD.vec)=='Fy']
Fy.uq <- Fy.mle+1.96*SD.vec[names(SD.vec)=='Fy']
Cton.mle <- MLE.vec[names(MLE.vec)=='Cton']
Cton.lq <- Cton.mle-1.96*SD.vec[names(SD.vec)=='Cton']
Cton.uq <- Cton.mle+1.96*SD.vec[names(SD.vec)=='Cton']
Fbar26.mle <- MLE.vec[names(MLE.vec)=='Fbar26']
Fbar26.lq <- Fbar26.mle-1.96*SD.vec[names(SD.vec)=='Fbar26']
Fbar26.uq <- Fbar26.mle+1.96*SD.vec[names(SD.vec)=='Fbar26']

## polygon plots for both

smin <- min(ssb.lq)
smax <- max(ssb.uq)
fmin <- min(Fy.lq)
fmax <- max(Fy.uq)
cmin <- min(c(Cton.lq, catch_ton))
cmax <- max(c(Cton.uq, catch_ton))
fbmin <- min(Fbar26.lq)
fbmax <- max(Fbar26.uq)
yrs <- 1:nyears

X11(width = 20)
par(mfrow=c(1,4))
## SSB:
plot(yrs,rep(NA,nyears),xlab='year',ylab='SSB',ylim=c(smin,smax),main='Spawning stock biomass',type='p')
polygon(c(yrs,rev(yrs)),c(ssb.lq,rev(ssb.uq)),col='cyan',lty=0)
lines(yrs,ssb.mle,lty=1,lwd=1.5,col='red')
## Fmax:
plot(yrs,rep(NA,nyears),xlab='year',ylab='F',ylim=c(fmin,fmax),main='Maximum fishing mortality',type='p')
polygon(c(yrs,rev(yrs)),c(Fy.lq,rev(Fy.uq)),col='cyan',lty=0)
lines(yrs,Fy.mle,lty=1,lwd=1.5,col='red')
## Total catch:
plot(yrs,rep(NA,nyears),xlab='year',ylab='Catch (t)',ylim=c(cmin,cmax),main='Expected catch (tons)',type='p')
polygon(c(yrs,rev(yrs)),c(Cton.lq,rev(Cton.uq)),col='cyan',lty=0)
lines(yrs,Cton.mle,lty=1,lwd=1.5,col='red')
points(yrs, catch_ton, pch = 1, cex = 2)
points(yrs, catch_ton, pch = 16, cex = 0.1)
## Fbar26
plot(yrs,rep(NA,nyears),xlab='year',ylab='Fbar',ylim=c(fbmin,fbmax),
     main='Mean fishing mortality age 2 to 6',type='p')
polygon(c(yrs,rev(yrs)),c(Fbar26.lq,rev(Fbar26.uq)),col='cyan',lty=0)
lines(yrs,Fbar26.mle,lty=1,lwd=1.5,col='red')


## now take a look at the fits and other stuff in the REPORT() variables

## this obtains the report variables from the AD object "obj"

rep <- obj$report()

## catches

library(lattice)

## create a data frame for the results is always useful for plotting

Chat <- rep$Chat
Csd <- rep[['exp(logsigma)']][1]
stdres <- log(catches/Chat)/Csd

c.df <- expand.grid(year=1:nyears,age=1:nages,obs=NA,pred=NA,stdres=NA)
for(y in 1:nyears) {

  c.df[c.df$year == y,'obs'] <- catches[,y]
  c.df[c.df$year == y,'pred'] <- Chat[,y]
  c.df[c.df$year == y,'stdres'] <- stdres[,y]

}

## predicted vs. observed:
key.plot <- list(space='top',col=c('blue','magenta'),text=list(c("Observed","Predicted")))

xyplot(obs+pred~age|as.factor(year),
       c.df,type=c('p','l'),lty=c(0,1),
       pch=c(19,NA),col=c('blue','magenta'),
       lwd=c(0,2),ylab='Catches',key=key.plot)

## bubble plot for absolute magnitude of the standardised residuals:
bub.scale <- 0.5
xyplot(age~as.factor(year),data=c.df,xlab='year',ylab='age',
       grid=TRUE,
       cex=c.df$stdres,
       panel = function(x,y,...,cex) {
         panel.xyplot(x,y,cex=abs(cex[])/bub.scale,pch=21,...)
       })

# alternative: see if SD(stdres) > 1 for each age along years
# ages where SD(stdres) > 1 could be "problematic"

plot(1:nages,apply(stdres,1,sd),xlab='age',ylab='SD(stdres)',type='p',pch=19,col='blue')
abline(h=1,col='magenta',lty=2)

## survey:
Ihat <- rep$U1
Isd <- rep[['exp(logsigma)']][2]
stdres <- log(tuning[1:15,]/Ihat)/Isd

i.df <- expand.grid(year=1:nyears,age=1:nages,obs=NA,pred=NA,stdres=NA)
for(y in 1:nyears) {

  i.df[i.df$year == y,'obs'] <- tuning[,y]
  i.df[i.df$year == y,'pred'] <- Ihat[,y]
  i.df[i.df$year == y,'stdres'] <- stdres[,y]

}

## predicted vs. observed:
key.plot <- list(space='top',col=c('blue','magenta'),text=list(c("Observed","Predicted")))

xyplot(obs+pred~age|as.factor(year),i.df,type=c('p','l'),lty=c(0,1),pch=c(19,NA),col=c('blue','magenta'),lwd=c(0,2),ylab='Index',key=key.plot)

## bubble plot for absolute magnitude of the standardised residuals:
bub.scale <- 0.5
xyplot(age~as.factor(year),data=i.df,xlab='year',ylab='age',
       grid=TRUE,
       cex=i.df$stdres,
       panel = function(x,y,...,cex) {
         panel.xyplot(x,y,cex=abs(cex[])/bub.scale,pch=21,...)
       })

# alternative: see if SD(stdres) > 1 for each age along years
# ages where SD(stdres) > 1 could be "problematic"

plot(1:nages,apply(stdres,1,sd),xlab='age',ylab='SD(stdres)',type='p',pch=19,col='blue')
abline(h=1,col='magenta',lty=2)

## selectivity:
self1 <- rep$self1
plot(1:nages,self1,xlab='ages',ylab='selectivity',main='',type='l',lwd=2,col='magenta')

## recruitment:
R <- rep$N[1,]
plot(1:nyears,R,type='l',xlab='year',ylab='recruitment',main=paste('SD(log(R)) = ',round(sd(log(R)),2),sep=""),lwd=2,col='magenta')

save.image("simple_assessment_TMB.RData")

