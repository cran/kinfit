kinerrmin <- function(kinfits, kinmodel = "SFO", alpha = 0.05)
{
	m = kinfits[[kinmodel]]

	kindata <- data.frame(t = kinfits[[kinmodel]]$model$t, 
		parent = kinfits[[kinmodel]]$model$parent)
        kindata.means <- aggregate(kindata, list(kindata$t), mean)
	kindata.means.mean <- mean(kindata.means$parent, na.rm=TRUE)

	n.parms = length(coef(m))
	df = length(kindata.means$parent) - n.parms
	kindata.means$est <- predict(m, kindata.means)

	f <- function(err)
	{
		(sum((kindata.means$parent - kindata.means$est)^2/((err*kindata.means.mean)^2)) - 
		 qchisq(1 - alpha,df))^2
	}
	err.min <- optimize(f, c(0.01,0.9))$minimum
	return(err.min)
}
