kinresults <- function(kinfits, alpha = 0.05, SFORB=TRUE)
{
	kindata <- data.frame(t = kinfits[[1]]$model$t, parent = kinfits[[1]]$model$parent)
        kindata.means <- aggregate(kindata, list(kindata$t), mean)
	kindata.means.mean <- mean(kindata.means$parent, na.rm=TRUE)
	n.times <- length(kindata.means$parent)
	parms <- list()
	df <- err.min <- RSS <- vector()
	DT50 <- DT90 <- vector()
	f <- list()
	for (kinmodel in names(kinfits))
	{
		m = kinfits[[kinmodel]]
		if(class(m) == "nls") {
			kindata.means$est <- predict(m, kindata.means)
			parms[[kinmodel]] <- switch(kinmodel,
				SFO = list(parent.0 = coef(m)[["parent.0"]], 
                                    k = coef(m)[["k"]]),
				FOMC = list(parent.0 = coef(m)[["parent.0"]],
                                    alpha = coef(m)[["alpha"]],
				    beta = coef(m)[["beta"]]),
				HS = list(parent.0 = coef(m)[["parent.0"]], 
                                    k1 = coef(m)[["k1"]],
				    k2 = coef(m)[["k2"]], 
                                    tb = coef(m)[["tb"]]),
				DFOP = list(parent.0 = coef(m)[["parent.0"]],
                                    k1 = coef(m)[["k1"]],
				    k2 = coef(m)[["k2"]], 
                                    g = coef(m)[["g"]]))
			if(kinmodel == "DFOP" & SFORB) {
				k1 = coef(m)[["k1"]]
				k2 = coef(m)[["k2"]]
				g = coef(m)[["g"]]
				parms[["SFORB"]] = 
                                    list(parent.0 = coef(m)[["parent.0"]],
					k1out = g * k1 + (1 - g) * k2,
					k21 = k1 * k2 / (g * k1 + (1 - g) * k2),
					k12 = (g * (1 - g) * (k1 - k2)^2) / (g * k1 + (1 - g) * k2))
			}
			n.parms = length(coef(m))
			f[[kinmodel]] = switch(kinmodel,
				HS = function(t, x) {
					(HS(t, coef(m)[["parent.0"]], 
						coef(m)[["k1"]], coef(m)[["k2"]], coef(m)[["tb"]]) - 
					(1 - x/100) * coef(m)[["parent.0"]])^2
				},
				DFOP = function(t, x) {
					(DFOP(t, coef(m)[["parent.0"]], 
						coef(m)[["k1"]], coef(m)[["k2"]], coef(m)[["g"]]) - 
					(1 - x/100) * coef(m)[["parent.0"]])^2
				}
			)
			coef(m)

			df[[kinmodel]] = n.times - n.parms
			RSS[[kinmodel]] = sum(summary(m)$residuals^2)
			DT50[[kinmodel]] = switch(kinmodel,
					SFO = log(2)/coef(m)[["k"]],
				FOMC = coef(m)[["beta"]] * (2^(1/coef(m)[["alpha"]]) - 1),
				HS = optimize(f[[kinmodel]], c(0, max(kindata$t)), x=50)$minimum,
				DFOP = optimize(f[[kinmodel]], c(0, max(kindata$t)), x=50)$minimum)
			DT90[[kinmodel]] = switch(kinmodel,
				SFO = log(10)/coef(m)[["k"]],
				FOMC = coef(m)[["beta"]] * (10^(1/coef(m)[["alpha"]]) - 1),
				HS = optimize(f[[kinmodel]], c(0, max(kindata$t)), x=90)$minimum,
				DFOP = optimize(f[[kinmodel]], c(0, max(kindata$t)), x=90)$minimum)
			err.min[[kinmodel]] <- kinerrmin(kinfits, kinmodel)
		}
	}
	stats <- data.frame(n.times = n.times, df = df, mean.means = kindata.means.mean, 
		RSS = RSS, err.min = err.min)
	results <- data.frame(DT50 = DT50, DT90 = DT90)
	list(parms = parms, stats = stats, results = results)
}
