library(colorRamps)
num_points <- 4096
disp3D <- function(df){
    scalePlot <- 0.25
    plot(df$X, df$Y, pch=20, cex=df$radius, lwd=0, col=df$colors, asp=1, xlim=c(-scalePlot, scalePlot), ylim=c(-scalePlot, scalePlot), axes=F, xlab="", ylab="")
}
incline3D <- function(df, angle){
    CS <- cospi(angle/180)
    SN <- sinpi(angle/180)
    df$X <- df$clumpY
    df$Y <- df$clumpZ* CS - df$clumpX* SN
    df$Z <- df$clumpX
    return(df[order(df$Z),])
}
#---- clumpy torus (r, theta, phi)
r <- 0.072 + rexp(num_points, rate=15)
phi <- runif(num_points, min=0, max=2)
x <- r* cospi(phi)
y <- r* sinpi(phi)
z <- rnorm(num_points, mean=0, sd=(r-0.054)*exp(-20*(r-0.072)^2))
cos_theta <- abs(z)/sqrt(r^2 + z^2)
L <- 20*exp( -2*(r - 0.054) )
ClumpPos <- data.frame(clumpX=x, clumpY=y, clumpZ=z, BH = rep(0,num_points), temperature = L * (cos_theta + 0.5) / (r^2 + z^2), radius = rnorm(num_points, mean=0.2, sd=0.02) / (z^2 + r^2 + 0.1))
#-------- Coloring by temperature
ClumpPos$colors <- paste(blue2green2red(128)[floor(128* ClumpPos$temperature / (max(ClumpPos$temperature) + 5)) + 1], "40", sep="")    # alpha = 40
ClumpPos[1,]$colors <- "#000000"  # Core component is black
ClumpPos[1,]$clumpX <- 0.0; ClumpPos[1,]$clumpY <- 0.0; ClumpPos[1,]$clumpZ <- 0.0; ClumpPos[1,]$radius <- 2
#-------- inclination and projection
ClumpPos <- incline3D(ClumpPos, 20)
#-------- Front and Back separation
ClumpBack  <- ClumpPos[ClumpPos$Z <= 0.0,]
ClumpFront <- ClumpPos[ClumpPos$Z >= 0.0,]
rownames(ClumpBack) <- 1:nrow(ClumpBack)
rownames(ClumpFront) <- 1:nrow(ClumpFront)
#-------- Plot Front and Back
pdf("back.pdf")
disp3D(ClumpBack)
dev.off()
pdf("front.pdf")
disp3D(ClumpFront)
dev.off()
#-------- Top view
ClumpPos <- incline3D(ClumpPos, 90)
pdf("top.pdf")
disp3D(ClumpPos)
dev.off()
#-------- Cross cut
ClumpPos <- incline3D(ClumpPos, 0)
ClumpBack  <- ClumpPos[ClumpPos$Z <= 0.0,]
pdf("CrossCut.pdf")
disp3D(ClumpBack)
dev.off()

