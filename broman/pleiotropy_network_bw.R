source("colors_bw.R")
library(broman)
iArrows <- igraph:::igraph.Arrows

pdf("../Figs/pleiotropy_network_bw.pdf", height=6, width=5, pointsize=32)
par(mar=rep(0.1,4),las=1,fg="black",col="black",col.axis="black",col.lab="black",
    bg=bgcolor,bty="n")

plot(0,0, xaxt="n", yaxt="n", xlab="", ylab="", type="n",
     xlim=c(-16-2/3, 150), ylim=c(0, 100), xaxs="i", yaxs="i")
x <- c(25, 75)
y <- seq(10, 90, len=5)
text(x[1], y[2], expression(q[2]))
text(x[1], y[4], expression(q[1]))

text(x[2], y[1], expression(y[3]))
text(x[2], y[3], expression(y[2]))
text(x[2], y[5], expression(y[1]))

arrowcol <- "black"
xd <- 8
yd <- 3.5
arrowlwd <- 5
arrowlen <- 0.3
arrows(x[1]+xd, y[2]-yd, x[2]-xd, y[1]+yd, lwd=arrowlwd,
      col=arrowcol, len=arrowlen)
arrows(x[1]+xd, y[2]+yd, x[2]-xd, y[3]-yd, lwd=arrowlwd,
      col=arrowcol, len=arrowlen)
arrows(x[1]+xd, y[4]-yd, x[2]-xd, y[3]+yd, lwd=arrowlwd,
      col=arrowcol, len=arrowlen)
arrows(x[1]+xd, y[4]+yd, x[2]-xd, y[5]-yd, lwd=arrowlwd,
      col=arrowcol, len=arrowlen)

arrowcol <- "black"
yd2 <- 4
xd2 <- 8
iArrows(x[1]-xd2, y[2]+yd2, x[1]-xd2, y[4]-yd2,
       h.lwd=arrowlwd, sh.lwd=arrowlwd,  sh.col=arrowcol,
       curve=1, width=arrowlwd, size=0.3, code=0)
# code=3 makes it double-headed
# code=0 makes no heads

# adding the heads "manually"
xd3 <- 0.05
yd3 <- 0.05
arrows(x[1]-xd2-xd3, y[4]-yd2-yd3, x[1]-xd2, y[4]-yd2,
       lwd=arrowlwd, col=arrowcol, len=arrowlen)
arrows(x[1]-xd2-xd3, y[2]+yd2+yd3, x[1]-xd2, y[2]+yd2,
       lwd=arrowlwd, col=arrowcol, len=arrowlen)

dev.off()
