

plotnsum <- function (nsum, ny = 5, by.y = FALSE, xlim = c(0.75,120), ylim = xlim, 
					  label = "", labelx = 0.15, type = c("Fletcher", "Wedderburn"), ...) {
	type <- match.arg(type)
	varnmat <- matrix(nsum[,'varration'], nrow = ny)
	if (type=="Fletcher")
		cmat <- matrix(nsum[,"chatF"], nrow = ny)
	else 
		cmat <- matrix(nsum[,"chatW"], nrow = ny)
	nx <- nrow(nsum)/ny
	ticks <- c(1,2,5,10,20,50,100)
	plot(1,1,type='n', log = 'xy', xlim=xlim, ylim=ylim, 
		 xlab = expression(paste("var(", italic(n), ") / E[var(", italic(n), ")]")), 
		 ylab = "", axes = F)
	axis(1, at=ticks, labels = ticks)
	axis(2, at=ticks, labels = ticks, las=1)
	box()
	if (by.y) 
		for (i in 1:nx) points(varnmat[,i], cmat[,i], type = 'o', ...)
	else 
		for (i in 1:ny) points(varnmat[i,], cmat[i,], type = 'o', ...)
	mtext(side = 2, line = 2.8, las=1, expression(paste(italic(hat(c)))), cex = 1)
	abline(0,1,lty=2)
	text (40, 60, 'y = x', srt=45, cex = 1)
	text(labelx, 240, label, xpd = TRUE, cex = 1.3, adj = 0)
}

