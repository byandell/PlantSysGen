source("func.R")
if(bw) source("colors_bw.R") else source("colors.R")

##############################
# DH_fig.R: Illustration of a doubled haploid
##############################

if(bw) {
  pdf(file="../Figs/DH_bw.pdf", width=9.75, height=6.5, pointsize=12, onefile=TRUE)
} else {
  pdf(file="../Figs/DH.pdf", width=9.75, height=6.5, pointsize=12, onefile=TRUE)
}

if(bw) {
  par(mar=rep(0.1,4),las=1,fg="black",col="black",col.axis="black",col.lab="black",
      bg=bgcolor,bty="n")
} else {
  par(mar=rep(0.1,4),las=1,fg="white",col="white",col.axis="white",col.lab="white",
      bg=bgcolor,bty="n")
}
plot(0,0,xlim=c(0,864),ylim=c(0,480),xaxt="n",yaxt="n",xlab="",ylab="",type="n")

rect_size <- list(x=10,y=95)
offset <- 18

chr_rect <- function(x, y, rect_size, offset, paired = TRUE) {
  if(paired) {
    x <- x + c(0, rect_size$x + offset)
    y <- rep(y, 2)
  }
  nodes <- list(xleft = x, ybottom = y)
  nodes$xright <- nodes$xleft + rect_size$x
  nodes$ytop <- nodes$ybottom + rect_size$y
  nodes
}
if(bw) {
  col <- rep(NA, 2)
  bcol <- rep(1, 2)
  dens <- 20
} else {
  bcol <- col <- color
  dens <- NULL
}
nodes <- chr_rect(300, 385, rect_size, offset)
rect(nodes$xleft, nodes$ybottom, nodes$xright, nodes$ytop,
     lend=1, ljoin=1, lwd = 2, col=col[1], border=bcol[1])
nodes <- chr_rect(526, 385, rect_size, offset)
rect(nodes$xleft, nodes$ybottom, nodes$xright, nodes$ytop,
     lend=1, ljoin=1, lwd = 2, col=col[2], border=bcol[2],
     angle = 45, density = dens)

points(432,440,pch=4,cex=2.5, lwd=2)
arrows(432, 400, 432, 300+18, len=0.1, lwd=2)

text(300-25,(480+385)/2,expression(P[1]),cex=2,adj=c(1,0.5))
text(564+25,(480+385)/2,expression(P[2]),cex=2,adj=c(0,0.5))

nodes <- chr_rect(413, 18+192, rect_size, offset, FALSE)
rect(nodes$xleft, nodes$ybottom, nodes$xright, nodes$ytop,
     lend=1, ljoin=1, lwd = 2, col=col[1], border=bcol[1])
nodes <- chr_rect(441, 18+192, rect_size, offset, FALSE)
rect(nodes$xleft, nodes$ybottom, nodes$xright, nodes$ytop,
     lend=1, ljoin=1, lwd = 2, col=col[2], border=bcol[2],
     angle = 45, density = dens)
segments(432,207,432,147, lwd=2)
segments(57,147,849,147, lwd=2)
arrows(seq(57,849,by=88),rep(147,10),seq(57,849,by=88),rep(107,10),len=0.1, lwd=2)
text(451+25,18+(287+192)/2,expression(F[1]),cex=2,adj=c(0,0.5))

f1 <- create.par(100,c(1,2))
set.seed(99019)
f2 <- vector("list",10)
for(i in 1:10) f2[[i]] <- cross(f1,f1,m=10,obl=TRUE)

xloc <- 38
mult <- 95/f2[[1]]$mat[1,ncol(f2[[1]]$mat)]
recomb_size <- rect_size
for(i in 1:10) {
  nodes <- chr_rect(xloc, 0, rect_size, offset)
  rect(nodes$xleft, nodes$ybottom, nodes$xright, nodes$ytop,
       lend=1, ljoin=1, lwd = 2, col=col[1], border=bcol[1])
  f2p <- f2[[i]]$pat
  for(j in 2:ncol(f2p)) {
    if(f2p[2,j]==2) {
      recomb_size$y <- -diff(f2p[1,j-c(0,1)])*mult
      nodes <- chr_rect(xloc, f2p[1,j-1]*mult, recomb_size, offset)
      rect(nodes$xleft, nodes$ybottom, nodes$xright, nodes$ytop,
           lend=1, ljoin=1, lwd = 2, col=col[2], border=bcol[2],
           angle = 45, density = dens)
    }
  }
  xloc <- xloc+38+50
}
text(38-25,95/2,"DH",cex=2,adj=c(1,0.5))
dev.off()




# end of backcross_fig.R
