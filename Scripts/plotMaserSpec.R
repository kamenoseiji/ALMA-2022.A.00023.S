library(ggplot2)
FWG <- 2*sqrt(2*log(2))	# Gaussian FWHM/sigma
#-------- for step fill
fillStep <- function(df){
    step <- 0.5* abs(median(diff(df$VLSR)))
    df$lead <- df$VLSR - step
    df$lag  <- df$VLSR + step
    return(df)
}
#-------- Integrated flux density
integFlux <- function(df, velRange){
    DFR <- df[((df$VLSR > velRange[1] ) & (df$VLSR < velRange[2] )),]
    velStep <- abs(median(diff(DFR$VLSR)))
    return( list(sum(DFR$Flux)* velStep, velStep* sqrt(sum(DFR$errFlux^2))))
}
#-------- Weighted mean velocity
weightedMeanVeloc <- function(df, velRange){
    DFR <- df[((df$VLSR > velRange[1] ) & (df$VLSR < velRange[2] )),]
    sumFlux  <- sum(DFR$Flux)
    meanVel <- as.numeric(DFR$Flux %*% DFR$VLSR)/sumFlux
    relVel  <- DFR$VLSR - meanVel
    errVel  <- sqrt(sum(relVel^2 * DFR$errFlux^2)) / sumFlux
    return( list(meanVel, errVel))
}
#-------- Weight (to avoid ozone contamination)
contamWeight <- function(df, maskRange){
	df$Weight <- 1.0
	mask_index <- which( (df$VLSR > min(maskRange)) & (df$VLSR < max(maskRange)))
	df$Weight[mask_index] <- 0.0
	return(df)
}
GaussFn   <- function(x, a, b, c){ return(a* exp(-0.5*((x - b)/c)^2)) }	# a:peak, b:ceter, c:width
LorentzFn <- function(x, a, b, c){ return(a* c^2 / ((x - b)^2 + c^2)) } # a:peak, b:ceter, c:width
#-------- GaussianFit
gaussFit <- function(df, lineRange, compCenter, compFlux, compWidth){
    fitDF <- df[((df$VLSR > min(lineRange)) & (df$VLSR < max(lineRange))),]
    compNum <- length(compCenter)
	formulaText <- 'Flux ~ flux1 * exp(-0.5*((VLSR - cent1) / widt1)^2)' 	# first Gaussian component
	initParam <- list(compFlux[1], compCenter[1], compWidth[1])
	initParamName <- c('flux1', 'cent1', 'widt1')
	for(comp_index in 2:compNum){	#---- Generate fitting parameters of flux, center, and width
		formulaText <- paste(formulaText, sprintf(' + flux%d * exp(-0.5*((VLSR - cent%d) /  widt%d)^2)', comp_index, comp_index, comp_index))
		initParam <- append(initParam, list(compFlux[comp_index], compCenter[comp_index], compWidth[comp_index]))
		initParamName <- append(initParamName, c(sprintf('flux%d', comp_index), sprintf('cent%d', comp_index), sprintf('widt%d', comp_index)))
	}
	names(initParam) <- initParamName
	fit <- nls(formula = as.formula(formulaText), data=fitDF, start=initParam)
	return(fit)
}
#-------- First Gaussian + multiple Lorents 
gaussLorentzFit <- function(df, lineRange, compList, compCenter, compFlux, compWidth){
    fitDF <- df[((df$VLSR > min(lineRange)) & (df$VLSR < max(lineRange))),]
    compNum <- length(compList)
	formulaText <- 'Flux ~ flux1 * exp(-0.5*((VLSR - cent1) / widt1)^2)' 	# first Gaussian component
	if(compList[1] == 'L'){	formulaText <- 'Flux ~ flux1 * widt1^2 / ((VLSR - cent1)^2 + widt1^2)'} 	# first Lorentzian component
	initParam <- list(compFlux[1], compCenter[1], compWidth[1])
	initParamName <- c('flux1', 'cent1', 'widt1')
	for(comp_index in 2:compNum){	#---- Generate fitting parameters of flux, center, and width
		if(compList[comp_index] == 'G'){ formulaText <- paste(formulaText, sprintf(' + flux%d * exp(-0.5*((VLSR - cent%d) /  widt%d)^2)', comp_index, comp_index, comp_index)) }
		if(compList[comp_index] == 'L'){ formulaText <- paste(formulaText, sprintf(' + flux%d * widt%d^2 / ((VLSR - cent%d)^2 + widt%d^2)', comp_index, comp_index, comp_index, comp_index)) }
		initParam <- append(initParam, list(compFlux[comp_index], compCenter[comp_index], compWidth[comp_index]))
		initParamName <- append(initParamName, c(sprintf('flux%d', comp_index), sprintf('cent%d', comp_index), sprintf('widt%d', comp_index)))
	}
	names(initParam) <- initParamName
	fit <- nls(formula = as.formula(formulaText), data=fitDF, weights=1/errFlux^2, start=initParam, control = list(maxiter = 500))
	return(fit)
}
#-------- plot Data and model
plotDataModel <- function(df, FIT, lineRange, compList, Title){
	fitDF <- df[((df$VLSR > min(lineRange)) & (df$VLSR < max(lineRange))),]
	plot(fitDF$VLSR, fitDF$Flux, type='n', xlim=lineRange, ylim=c(-0.02, 0.08), xlab='VLSR [km/s]', ylab='Flux Density [Jy]', main=Title)
	grid(lty=2, col='gray', lwd=0.5)
	arrows(fitDF$VLSR, fitDF$Flux - fitDF$errFlux, fitDF$VLSR, fitDF$Flux + fitDF$errFlux, length=0, lwd=0.5, col='gray')
	points(fitDF$VLSR, fitDF$Flux, pch=20, cex=0.5)
	lines(fitDF$VLSR, predict(FIT, fitDF), col='blue')
	for(comp_index in seq_along(compList)){
		if(compList[comp_index] == 'G'){
			compFlux <- GaussFn(fitDF$VLSR, coef(FIT)[[3*(comp_index-1)+1]], coef(FIT)[[3*(comp_index-1)+2]], coef(FIT)[[3*(comp_index-1)+3]] )
			text_sd <- sprintf('%d : $%.1f \\pm %.1f$ & $%.1f \\pm %.1f$ & $%.1f \\pm %.1f$  & $%.3f \\pm %.3f$ \n', comp_index,
				coef(FIT)[[3*(comp_index-1)+2]],     summary(FIT)$coefficients[,'Std. Error'][[3*(comp_index-1)+2]],		# center velocity
				coef(FIT)[[3*(comp_index-1)+1]]*1e3, summary(FIT)$coefficients[,'Std. Error'][[3*(comp_index-1)+1]]*1e3,	# flux (mJy)
				coef(FIT)[[3*(comp_index-1)+3]]*FWG, summary(FIT)$coefficients[,'Std. Error'][[3*(comp_index-1)+3]]*FWG,   # FWHM (km/s)
				coef(FIT)[[3*(comp_index-1)+3]]* coef(FIT)[[3*(comp_index-1)+1]]* sqrt(2.0*pi),                             # Integrated flux density
                sqrt(2.0*pi)* sqrt( (summary(FIT)$coefficients[,'Std. Error'][[3*(comp_index-1)+3]]/coef(FIT)[[3*(comp_index-1)+3]])^2 + (summary(FIT)$coefficients[,'Std. Error'][[3*(comp_index-1)+1]]/coef(FIT)[[3*(comp_index-1)+1]])^2)* coef(FIT)[[3*(comp_index-1)+3]]* coef(FIT)[[3*(comp_index-1)+1]])

		}
		if(compList[comp_index] == 'L'){
			compFlux <- LorentzFn(fitDF$VLSR, coef(FIT)[[3*(comp_index-1)+1]], coef(FIT)[[3*(comp_index-1)+2]], coef(FIT)[[3*(comp_index-1)+3]] )
			text_sd <- sprintf('%d : $%.1f \\pm %.1f$ & $%.1f \\pm %.1f$ & $%.1f \\pm %.1f$  & $%.3f \\pm %.3f$ \n', comp_index,
				coef(FIT)[[3*(comp_index-1)+2]],     summary(FIT)$coefficients[,'Std. Error'][[3*(comp_index-1)+2]],		# center velocity
				coef(FIT)[[3*(comp_index-1)+1]]*1e3, summary(FIT)$coefficients[,'Std. Error'][[3*(comp_index-1)+1]]*1e3,	# flux (mJy)
				coef(FIT)[[3*(comp_index-1)+3]]*2,   summary(FIT)$coefficients[,'Std. Error'][[3*(comp_index-1)+3]]*2,      # FWHM (km/s)
				coef(FIT)[[3*(comp_index-1)+3]]* coef(FIT)[[3*(comp_index-1)+1]]* pi,                                       # Integrated flux density
                pi* sqrt( (summary(FIT)$coefficients[,'Std. Error'][[3*(comp_index-1)+3]]/coef(FIT)[[3*(comp_index-1)+3]])^2 + (summary(FIT)$coefficients[,'Std. Error'][[3*(comp_index-1)+1]]/coef(FIT)[[3*(comp_index-1)+1]])^2)* coef(FIT)[[3*(comp_index-1)+3]]* coef(FIT)[[3*(comp_index-1)+1]])
		}
		lines(fitDF$VLSR, compFlux, col=comp_index, lty=2)
		cat(text_sd)
	}
	abline(h=-0.01, col='gray')
	arrows(fitDF$VLSR, fitDF$Flux - predict(FIT, fitDF) - fitDF$errFlux - 0.01, fitDF$VLSR, fitDF$Flux - predict(FIT, fitDF) + fitDF$errFlux - 0.01, length=0, lwd=0.5, col='black')
	chiSQ <- sum( ((fitDF$Flux - predict(FIT, fitDF))/fitDF$errFlux)^2 )
	text( min(lineRange), -0.02, pos=4, sprintf('residual chi-squared / d.o.f. = %.1f / %d', chiSQ, nrow(fitDF) - 3*length(compList)))
	#weightCenter <- fitDF$VLSR %*% fitDF$Flux / sum(fitDF$Flux)
	#cat(sprintf('Velocity mean = %.1f\n', weightCenter))
}
#-------- Start
maskRange <- c(1150, 1180)	# Velocity range for ozone contamination
May12 <- read.table('NGC1052_SPW3.veloc.txt')
Jul10 <- read.table('2022.A.00023.S.X109d26e_X12c34.veloc.txt', skip=9, header=FALSE)
Jul24 <- read.table('2022.A.00023.S.X10a7a20_X65a8.veloc.txt', skip=9, header=FALSE)
names(May12) <- c('VLSR', 'Flux')
names(Jul10) <- c('VLSR', 'Flux')
names(Jul24) <- c('VLSR', 'Flux')
#
May12 <- contamWeight(fillStep(May12), maskRange)
Jul10 <- contamWeight(fillStep(Jul10), maskRange)
Jul24 <- contamWeight(fillStep(Jul24), maskRange)
#-------- Spectral noise
May12$errFlux <- 3.3e-3
Jul10$errFlux <- 2.3e-3
Jul24$errFlux <- 2.0e-3
#-------- Component Model (Gaussian or Lorentzian)
modelMay12 <- c('G', 'G')
modelJul10 <- c('G', 'L', 'G', 'G', 'G', 'G', 'L')
modelJul24 <- c('G', 'L', 'G', 'G', 'G', 'G', 'L')
#-------- May12 baseline subtraction
ContFit <- lm(formula=Flux~VLSR, data=May12[((May12$VLSR < 1150) | (May12$VLSR > 1700)),])
baseLine <- predict(ContFit, data.frame(VLSR=May12$VLSR))
May12$Flux <- May12$Flux - baseLine
May12$baseLine <- 0.0
#-------- Integrated flux density
lineRange <- c(1200, 1700)
May12integ <- integFlux(May12, lineRange); May12meanVel <- weightedMeanVeloc(May12, lineRange)
Jul10integ <- integFlux(Jul10, lineRange); Jul10meanVel <- weightedMeanVeloc(Jul10, lineRange)
Jul24integ <- integFlux(Jul24, lineRange); Jul24meanVel <- weightedMeanVeloc(Jul24, lineRange)
cat(sprintf('May12: integFlux = %f +- %f Jy km/s  /  mean = %f +- %f km/s\n', May12integ[1], May12integ[2], May12meanVel[1], May12meanVel[2]))
cat(sprintf('Jul10: integFlux = %f +- %f Jy km/s  /  mean = %f +- %f km/s\n', Jul10integ[1], Jul10integ[2], Jul10meanVel[1], Jul10meanVel[2]))
cat(sprintf('Jul24: integFlux = %f +- %f Jy km/s  /  mean = %f +- %f km/s\n', Jul24integ[1], Jul24integ[2], Jul24meanVel[1], Jul24meanVel[2]))
#cat(sprintf('May12: integFlux = %f +- %f Jy km/s\n', integFlux(May12, lineRange)[1], integFlux(May12, lineRange)[2] ))
#cat(sprintf('Jul24: integFlux = %f +- %f Jy km/s\n', integFlux(Jul10, lineRange)[1], integFlux(Jul10, lineRange)[2] ))
#cat(sprintf('Jul24: integFlux = %f +- %f Jy km/s\n', integFlux(Jul24, lineRange)[1], integFlux(Jul24, lineRange)[2] ))
#-------- Gaussian/Lorenzian components
fitMay12 <- gaussLorentzFit(May12, lineRange, modelMay12, c(1467, 1333), c(0.03, 0.06), c(70, 20))
fitJul10 <- gaussLorentzFit(Jul10, lineRange, modelJul10, c(1540, 1373, 1486, 1525, 1552, 1595, 1607), c(0.040, 0.042, 0.008, 0.014, 0.012, 0.010, 0.03), c(50, 25, 1, 6, 3, 2, 20))
fitJul24 <- gaussLorentzFit(Jul24, lineRange, modelJul24, c(1539, 1380, 1491, 1547, 1567, 1603, 1611), c(0.027, 0.056, 0.007, 0.034, 0.007, 0.046, 0.007), c(75, 26, 10,28, 2, 20, 2))
# fitJul24 <- gaussLorentzFit(Jul24, lineRange, modelJul24, c(1539, 1380, 1491, 1547, 1567, 1611, 1603), c(0.027, 0.056, 0.007, 0.034, 0.007, 0.006, 0.046), c(75, 26, 10,28, 2, 2, 21))
# fitJul24 <- gaussLorentzFit(Jul24, lineRange, modelJul24, c(1539, 1380, 1491, 1567, 1611, 1603, 1547), c(0.027, 0.058, 0.007, 0.007, 0.006, 0.046, 0.034), c(75, 26, 10, 2, 2, 21,28))
pdf('plotSpecComp.May12.pdf'); plotDataModel(May12, fitMay12, lineRange, modelMay12, '2022-05-12'); dev.off(); 
pdf('plotSpecComp.Jul10.pdf'); plotDataModel(Jul10, fitJul10, lineRange, modelJul10, '2023-07-10'); dev.off()
pdf('plotSpecComp.Jul24.pdf'); plotDataModel(Jul24, fitJul24, lineRange, modelJul24, '2023-07-24'); dev.off()
#-------- Red profile by subtracting blue model
RedJul10 <- Jul10; RedJul24 <- Jul24
comp_index <- 2 # subtraction of blue component
RedJul10$Flux <- RedJul10$Flux - LorentzFn(RedJul10$VLSR, coef(fitJul10)[[3*(comp_index-1)+1]], coef(fitJul10)[[3*(comp_index-1)+2]], coef(fitJul10)[[3*(comp_index-1)+3]] )
RedJul24$Flux <- RedJul24$Flux - LorentzFn(RedJul24$VLSR, coef(fitJul24)[[3*(comp_index-1)+1]], coef(fitJul24)[[3*(comp_index-1)+2]], coef(fitJul24)[[3*(comp_index-1)+3]] )
RedJul10integ <- integFlux(RedJul10, lineRange); RedJul10meanVel <- weightedMeanVeloc(RedJul10, lineRange)
RedJul24integ <- integFlux(RedJul24, lineRange); RedJul24meanVel <- weightedMeanVeloc(RedJul24, lineRange)
cat(sprintf('Jul10: Red = %f +- %f Jy km/s   /  mean = %f +- %f km/s\n', RedJul10integ[1], RedJul10integ[2], RedJul10meanVel[1], RedJul10meanVel[2]))
cat(sprintf('Jul24: Red = %f +- %f Jy km/s   /  mean = %f +- %f km/s\n', RedJul24integ[1], RedJul24integ[2], RedJul24meanVel[1], RedJul24meanVel[2]))
if(0){
#-------- Flux scaling
Jul10$baseLine <- 0.075
Jul24$baseLine <- 0.1
Jul10$Flux <- Jul10$Flux + Jul10$baseLine
Jul24$Flux <- Jul24$Flux + Jul24$baseLine
#-------- Date category
May12$Date <- rep('2022-05-12', nrow(May12))
Jul10$Date <- rep('2023-07-10 (+ 0.075 Jy)', nrow(Jul10))
Jul24$Date <- rep('2023-07-24 (+ 0.1 Jy)', nrow(Jul24))
#-------- Plot Maser profiles
Maser <- rbind(May12, Jul10, Jul24)
ggp <- ggplot(data=Maser, mapping=aes(x=VLSR, y=Flux)) + geom_vline(xintercept = 1492.0, linetype = 'dotted') + geom_rect(xmin=1150, xmax=1180, ymin=-0.01, ymax=0.15, fill='gray', alpha=0.4) + geom_step(aes(color=Date), direction='mid') + geom_rect(aes(xmin=lead, xmax=lag, ymin=baseLine, ymax=Flux, fill=Date), alpha=0.5) + xlim(600, 2350) + xlab('VLSR [km/s]') + ylab('Flux Density [Jy]')
ggsave('NGC1052_Maser.pdf', device='pdf', width=11, height=8, units='in')
}
if(0){
#-------- Velocity Drift
velDF <- data.frame( date=as.Date(c('2022-05-12', '2023-07-10', '2023-07-24')),
    VLSR1=c(summary(fitMay12)$coefficients[,'Estimate']['b0'], summary(fitJul10)$coefficients[,'Estimate']['b0'], summary(fitJul24)$coefficients[,'Estimate']['b0']), 
    VLSR2=c(summary(fitMay12)$coefficients[,'Estimate']['b1'], summary(fitJul10)$coefficients[,'Estimate']['b1'], summary(fitJul24)$coefficients[,'Estimate']['b1']), 
    eVLSR1=c(summary(fitMay12)$coefficients[,'Std. Error']['b0'], summary(fitJul10)$coefficients[,'Std. Error']['b0'], summary(fitJul24)$coefficients[,'Std. Error']['b0']), 
    eVLSR2=c(summary(fitMay12)$coefficients[,'Std. Error']['b1'], summary(fitJul10)$coefficients[,'Std. Error']['b1'], summary(fitJul24)$coefficients[,'Std. Error']['b1'])) 
png('velocityShift.png')
plot(velDF$date, velDF$VLSR1, type='n', xlab='Date', ylab='VLSR [km/s]',ylim=c(1300,1600),  main='NGC 1052 H2O maser at 321 GHz')
grid(lty=2, col='gray', lwd=0.5)
arrows( velDF$date, velDF$VLSR1 - velDF$eVLSR1, velDF$date, velDF$VLSR1 + velDF$eVLSR1, length=0)
arrows( velDF$date, velDF$VLSR2 - velDF$eVLSR2, velDF$date, velDF$VLSR2 + velDF$eVLSR2, length=0)
points( velDF$date, velDF$VLSR1, pch=20, col='blue' )
points( velDF$date, velDF$VLSR2, pch=20, col='red' )
fit1 <- lm(VLSR1 ~ date, data=velDF, weights = 1.0 / eVLSR1^2)
abline(fit1)
fit2 <- lm(VLSR2 ~ date, data=velDF, weights = 1.0 / eVLSR2^2)
abline(fit2)
dev.off()
}
