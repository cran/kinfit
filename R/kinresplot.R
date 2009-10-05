kinresplot <- function(kinobject, kinmodel,
	xlab = "Time [days]", ylab = "Residual [% of applied radioactivity]",
	maxabs = "auto")
{
	m <- kinobject$fits[[kinmodel]]
	t <- m$model$t
	residuals <- residuals(m)
	if (maxabs == "auto") maxabs = max(abs(residuals))
	plot(t, residuals,
		xlab = xlab,
		ylab = ylab,
		ylim = c( -1.2 * maxabs, 1.2 * maxabs))
	title(paste("Residuals of", kinmodel, "fit"), font.main = 1)
}
