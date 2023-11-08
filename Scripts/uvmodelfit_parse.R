library('colorRamps')
#-------- Parse arguments
parseArg <- function( args ){
    argNum <- length(args)
    fileNum <- argNum
    return( list(fileList = args[1:argNum]))
}
#-------- Color Assignment
colIndex <- function(VLSR, VelRange){
    VelSpan  <- diff(VelRange)
    index <- ifelse( VLSR <= VelRange[1], 1,
                ifelse(VLSR >= VelRange[2], VelSpan, floor(VLSR - VelRange[1]) + 1))
    return(index)
}
#-------- Specific parameters
mapSpan <- 0.002
Vsys <- 1492
VelRange <- c(Vsys-150,Vsys+150)
VelSpan  <- diff(VelRange)
SpecRange <- c(1300,1700)
fileName <- 'X109d26e.maser.uvmodelfit.log'
specfile <- 'X109d26e_NGC_1052.SPW3.contsub.NatClean.veloc.txt'
BPA <- -7.351   # Synthesized beam PA (deg) for X109d26e   : Estimated beam: bmin=11.2 mas, bmaj=17.44 mas, bpa=-7.351 degrees
eFlux <- 2.3e-3
Title = '2023-07-10'
#fileName <- "X10a7a20.maser.uvmodelfit.log"
#specfile <- 'X10a7a20_NGC_1052.SPW3.contsub.NatClean.veloc.txt'
#BPA <- 41.53    # Synthesized beam PA (deg) for X10a7a20    : Estimated beam: bmin=10.44 mas, bmaj=11.85 mas, bpa=41.53 degrees
#eFlux <- 2.0e-3
#Title = '2023-07-24'
#-------- Start Program
Arguments <- commandArgs(trailingOnly = T)
if(length(Arguments) > 0){  fileName <- Arguments[1] }
fileLines <- readLines( fileName )
#---- Parse intensity
Flux_index <- grep('^I =', fileLines)
RA_index   <- grep('^x =', fileLines)
DEC_index  <- grep('^y =', fileLines)
CH_index   <- grep('Writing ', fileLines)
#---- Read spectral data
maserSpec <- read.table(specfile, skip=8)
names(maserSpec) <- c('VLSR','Flux')
#---- Read data
chNum <- length(Flux_index)
Flux <- Fluxerr <- RA <- RAerr <- DEC <- DECerr <- CH <- numeric(chNum)
for(index in 1:chNum){
    lineElement <- strsplit(fileLines[Flux_index[index]], ' ')[[1]]
    Flux[index] <- as.numeric(lineElement[3]); Fluxerr[index] <- as.numeric(lineElement[5])
    lineElement <- strsplit(fileLines[RA_index[index]], ' ')[[1]]
    RA[index]   <- as.numeric(lineElement[3]); RAerr[index] <- as.numeric(lineElement[5])
    lineElement <- strsplit(fileLines[DEC_index[index]], ' ')[[1]]
    DEC[index]  <- as.numeric(lineElement[3]); DECerr[index] <- as.numeric(lineElement[5])
    CH[index]   <- as.numeric(strsplit(fileLines[CH_index[index]], 'cl.')[[1]][2])
}
DF <- data.frame(Flux=Flux, CH=CH, RA=RA, DEC=DEC, eFlux=Fluxerr, eRA=RAerr, eDEC=DECerr)
#---- Adding VLSR info
DF$VLSR <- DF$CH
for( ch in DF$CH ){
    DF[DF$CH == ch,]$VLSR = maserSpec$VLSR[ch]
    DF[DF$CH == ch,]$Flux = mean(maserSpec$Flux[(ch-2):(ch+2)])
}
#---- Adding position along and perpendicular to the jet (PA=68.35 deg, along the continuum jet)
PA <- 68.35* pi / 180.0
CS <- cos(PA); SN <- sin(PA) 
DF$J <- SN* DF$RA  + CS* DF$DEC
DF$D <- SN* DF$DEC - CS* DF$RA
#---- Filter by 3 sigma
FiltDF <- DF[DF$Flux > 3.0*eFlux,]	# 3 sigma for X10a7a20
SpecDF <- DF[((DF$VLSR > min(SpecRange)) & (DF$VLSR < max(SpecRange))),]
pdf(paste(fileName, '.pdf', sep=''))
plot(FiltDF$RA, FiltDF$DEC, xlim=c(mapSpan, -mapSpan), ylim=c(-mapSpan, mapSpan), xlab='Relative RA [arcsec]', ylab='Relative DEC [arcsec]', main=Title, type='n', asp=1)
#---- Error bars
CS <- cos(pi* BPA/ 180.0)
SN <- sin(pi* BPA/ 180.0)
ifelse(mean(FiltDF$eRA) < mean(FiltDF$eDEC), P <- matrix(c(SN^2, CS^2, CS^2, SN^2), nrow=2), P <- matrix(c(CS^2, SN^2, SN^2, CS^2), nrow=2))
Pinv <- solve(P)
FiltDF$eMaj <- 1.0 / sqrt(Pinv[1,1]/FiltDF$eRA^2 + Pinv[1,2]/FiltDF$eDEC^2)
FiltDF$eMin <- 1.0 / sqrt(Pinv[2,1]/FiltDF$eRA^2 + Pinv[2,2]/FiltDF$eDEC^2)
arrows( FiltDF$RA  - SN* FiltDF$eMaj,
        FiltDF$DEC - CS* FiltDF$eMaj,
        FiltDF$RA  + SN* FiltDF$eMaj,
        FiltDF$DEC + CS* FiltDF$eMaj, length=0, col='#7F7F7F3F', lwd=0.3)
arrows( FiltDF$RA  + CS* FiltDF$eMin,
        FiltDF$DEC - SN* FiltDF$eMin,
        FiltDF$RA  - CS* FiltDF$eMin,
        FiltDF$DEC + SN* FiltDF$eMin, length=0, col='#7F7F7F3F', lwd=0.3)
points( FiltDF$RA, FiltDF$DEC, pch=20, cex=50*FiltDF$Flux, col=paste(matlab.like(VelSpan)[colIndex(FiltDF$VLSR, VelRange)], '7F', sep=''), lwd=0 )
dev.off()
#---- Plot Spectrum with Color
pdf(paste(fileName, '.spec.pdf', sep=''))
chWid <- abs(median(diff(SpecDF$VLSR)))
plot(SpecDF$VLSR, SpecDF$Flux, type='n', xlim=SpecRange, ylim=c(-eFlux, max(SpecDF$Flux)), xlab='VLSR [km/s]', ylab='Flux Density [Jy]', main=Title)
grid(lty=2, col='gray', lwd=0.5)
rect( SpecDF$VLSR - 0.5*chWid, 0.0, SpecDF$VLSR + 0.5*chWid, SpecDF$Flux,  col=paste(matlab.like(VelSpan)[colIndex(SpecDF$VLSR, VelRange)], '7F', sep=''), lwd=0, border=NA)
dev.off()
#---- Jet-axis transform
if(0){
plot(FiltDF$J, FiltDF$D, xlim=c(mapSpan,-mapSpan), ylim=c(-mapSpan, mapSpan), xlab='J [arcsec]', ylab='D [arcsec]', main=fileName, type='n', asp=1)
grid(lty=2, col='gray', lwd=0.5)
arrows( FiltDF$J+FiltDF$eRA, FiltDF$D, FiltDF$J-FiltDF$eRA, FiltDF$D, length=0, col='#7F7F7F3F')
arrows( FiltDF$J, FiltDF$D-FiltDF$eDEC, FiltDF$J, FiltDF$D+FiltDF$eDEC, length=0, col='#7F7F7F3F')
points( FiltDF$J, FiltDF$D, pch=20, cex=25*FiltDF$Flux, col=paste(matlab.like(VelSpan)[colIndex(FiltDF$VLSR, VelRange)], '7F', sep='') )
}
#---- Position velocity analysis
fit <- lm(formula=VLSR-Vsys ~ J + D, data=FiltDF, weights = Flux/(eRA*eDEC)) 
pdf(paste(fileName, '.J.pdf', sep=''))
plot(FiltDF$J, FiltDF$VLSR, type='n', xlim=c(0.0015, -0.0015), ylim=SpecRange, xlab='Offset along jet [arcsec]', ylab='VLSR [km/s]', main=Title)
grid(lty=2, col='gray', lwd=0.5)
arrows(FiltDF$J-FiltDF$eRA, FiltDF$VLSR, FiltDF$J+FiltDF$eRA, FiltDF$VLSR, length=0)
points(FiltDF$J, FiltDF$VLSR, pch=20, cex=50*FiltDF$Flux, col=paste(matlab.like(VelSpan)[colIndex(FiltDF$VLSR, VelRange)], 'FF', sep=''))
dev.off()
pdf(paste(fileName, '.D.pdf', sep=''))
plot(FiltDF$D, FiltDF$VLSR, type='n', xlim=c(-0.0015, 0.0015), ylim=SpecRange, xlab='Offset perpendicular to jet [arcsec]', ylab='VLSR [km/s]', main=Title)
grid(lty=2, col='gray', lwd=0.5)
arrows(FiltDF$D-FiltDF$eDEC, FiltDF$VLSR, FiltDF$D+FiltDF$eDEC, FiltDF$VLSR, length=0)
points(FiltDF$D, FiltDF$VLSR, pch=20, cex=50*FiltDF$Flux, col=paste(matlab.like(VelSpan)[colIndex(FiltDF$VLSR, VelRange)], 'FF', sep=''))
dev.off()
