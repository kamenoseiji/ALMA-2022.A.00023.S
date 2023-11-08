# library(plotly)
library(scatterplot3d)
library(colorRamps)
num_points <- 8192

#---- clumpy torus (r, theta, phi)
r <- 0.2 + rexp(num_points, rate=2)
# x <- c(r, -r)
phi <- runif(num_points, min=0, max=2)
x <- r* cospi(phi)
y <- r* sinpi(phi)
z <- rnorm(num_points, mean=0, sd=sqrt(0.5* r^2* exp(-0.5*r^2)))
# radius <- abs(rnorm(num_points, mean=1, sd=0.25)) / (z^2 + r^2 + 1)
cos_theta <- abs(z)/sqrt(r^2 + z^2)
L <- 250* exp( -0.1* r )
ClumpPos <- data.frame(clumpX=x, clumpY=y, clumpZ=z, BH = rep(0,num_points), temperature = L * cos_theta / (r^2 + z^2), radius = rnorm(num_points, mean=1, sd=0.1) / (z^2 + r^2 + 10))
ClumpPos$colors <- paste(blue2green2red(128)[floor(128* ClumpPos$temperature / 3000) + 1], "40", sep="")
# axisSet <- list(showgrid = F, showline = F, showaxeslabels = F, showbackground = F, nticks = 0, title = '')
ClumpPos[1,]$clumpX <- 0.0; ClumpPos[1,]$clumpY <- 0.0; ClumpPos[1,]$clumpZ <- 0.0; ClumpPos[1,]$radius <- 0.5; ClumpPos[1,]$temperature <- 2990; ClumpPos[1,]$colors <- "#FF0000FF"
ClumpBack  <- ClumpPos[((ClumpPos$clumpX <= 0.0) & (abs(ClumpPos$clumpY) < 3) & (abs(ClumpPos$clumpZ) < 3)),]
ClumpFront <- ClumpPos[((ClumpPos$clumpX >= 0.0) & (abs(ClumpPos$clumpY) < 3) & (abs(ClumpPos$clumpZ) < 3)),]

pdf("back.pdf")
scatterplot3d(ClumpBack$clumpX, ClumpBack$clumpY, ClumpBack$clumpZ, box=F, axis=F, grid=F, angle=20, pch=20, cex.symbols=2* ClumpBack$radius, color=ClumpBack$colors)
dev.off()
pdf("front.pdf")
scatterplot3d(ClumpFront$clumpX, ClumpFront$clumpY, ClumpFront$clumpZ, box=F, axis=F, grid=F, angle=20, pch=20, cex.symbols=2* ClumpFront$radius, color=ClumpFront$colors)
dev.off()

#figBack <- plot_ly(ClumpBack, type="scatter3d", mode="markers", x = ~clumpX, y = ~clumpY, z = ~clumpZ, color = ~temperature, size=~radius,
#    marker=list(symbol='circle', sizemode='diameter', line=list(width=0)), sizes = c(0, 10))
#figBack <- figBack %>% layout(scene = list(xaxis=list(range=c(-3,3)), yaxis=list(range=c(-3,3)), zaxis=list(range=c(-3,3)), bgcolor="#FFFFFF00", #camera=list(center=list(x=0, y=0, z=0), eye=list(x=1, y=0, z=0.2))))

#orca(figFront, "back.svg")
#
#figFront <- plot_ly(ClumpFront, type="scatter3d", mode="markers", x = ~clumpX, y = ~clumpY, z = ~clumpZ, color = ~temperature, size=~radius,
#    marker=list(symbol='circle', sizemode='diameter', line=list(width=0)), sizes = c(0, 10))
#figFront <- figFront %>% layout(scene = list(xaxis=list(range=c(-3,3)), yaxis=list(range=c(-3,3)), zaxis=list(range=c(-3,3)), bgcolor="#FFFFFF00", camera=list(center=list(x=0, y=0, z=0), eye=list(x=1, y=0, z=0.2))))

#orca(figFront, "front.svg")


# ggsave('hidoi.pdf', device='pdf', width=11, height=11, units='in')
# add_trace(BH, type="scatter3d", mode="markers", x = ~clumpX, y = ~clumpY, z = ~clumpZ, color = 'black', size=~radius,
#    marker=list(symbol='circle', sizemode='diameter', line=list(width=0)), sizes = c(0, 10)) %>%
#layout(xaxis = axisSet)

#fig <- fig %>% layout(ClumpPos[1],  type="scatter3d", mode="markers", x = ~clumpX, y = ~clumpY, z = ~clumpZ, color='red', size=~radius, #marker=list(symbol='circle', sizemode='diameter', line=list(width=0)), sizes = c(0, 10))
#fig <- layout(fig, xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE), yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))
#fig <- layout(fig, yaxis = list(scaleanchor = "x"))
#fig <- layout( fig, shapes = list(
#				list( type = 'circle',
#						xref = 'x', x0=-0.1, x1=0.1,
#						yref = 'y', y0=-0.1, y1=0.1,
#						opacity = 0.7)))
#