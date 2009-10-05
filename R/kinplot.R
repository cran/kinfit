kinplot <- function(kinobject, 
	xlab = "Time [days]", ylab = "Parent [% of applied radioactivity]",
        ylim = c("auto", "auto"),
	lpos = "topright")
{
	kindata <- na.omit(kinobject$data)
	kinfits <- kinobject$fits
        if (ylim[1] == "auto") ylim[1] <- 0
        if (ylim[2] == "auto") ylim[2] <- max(kindata$parent)
        ylim <- as.numeric(ylim)
        
	plot(kindata$t, kindata$parent,
	  xlab = xlab,
	  ylab = ylab,
	  ylim = ylim
        )
	n.m <- length(kinfits)
	colors <- ltys <- 1:n.m
	names(colors) <- names(ltys) <- names(kinfits)
        ltext <- paste(kinobject$parent, "measured")
	for (kinmodel in names(kinfits))
	{
		m = kinfits[[kinmodel]]
		if(class(m) == "nls") {
			switch(kinmodel,
				SFO = curve(SFO(x, 
					coef(m)[["parent.0"]], 
					coef(m)[["k"]]),
					from = min(kindata$t), to = max(kindata$t), add=TRUE,
					col = colors[[kinmodel]],
					lty = ltys[[kinmodel]]),
				FOMC = curve(FOMC(x, 
					coef(m)[["parent.0"]],
					coef(m)[["alpha"]],
					coef(m)[["beta"]]),
					from = min(kindata$t), to = max(kindata$t), add=TRUE,
					col = colors[[kinmodel]],
					lty = ltys[[kinmodel]]),
				HS = curve(HS(x, 
					coef(m)[["parent.0"]], 
					coef(m)[["k1"]],
					coef(m)[["k2"]],
					coef(m)[["tb"]]),
					from = min(kindata$t), to = max(kindata$t), add=TRUE,
					col = colors[[kinmodel]],
					lty = ltys[[kinmodel]]),
				DFOP = curve(DFOP(x, 
					coef(m)[["parent.0"]], 
					coef(m)[["k1"]],
					coef(m)[["k2"]],
					coef(m)[["g"]]),
					from = min(kindata$t), to = max(kindata$t), add=TRUE,
					col = colors[[kinmodel]],
					lty = ltys[[kinmodel]]))
                        ltext <- c(ltext, paste("Fitted", kinmodel, "model"))
		} else {
                        ltext <- c(ltext, paste(kinmodel, "model failed"))
                        ltys[[kinmodel]] <- NA
		} 
	}
	legend(lpos, bty="n", inset = 0.05, 
		legend = ltext,
		pch = c(1, rep(NA, n.m)),
		lty = c(NA, ltys),
		col = c(1, colors))
}
