kinfit <- function(kindata, kinmodels = c("SFO"), 
	parent.0.user = NA, 
	start.SFO = list(parent.0 = NA, k = NA), 
	start.FOMC = list(parent.0 = NA, alpha = NA, beta = NA), 
	start.DFOP = list(parent.0 = NA, k1 = NA, k2 = NA, g = NA),
	start.HS = list(parent.0 = NA, k1 = NA, k2 = NA, tb = NA),
        algorithm = "port")
{
	kindata <- subset(kindata, !is.na(kindata$parent))
	kinfits <- list()

	if (!is.na(parent.0.user)) {
		start.SFO$parent.0 = parent.0.user
		start.FOMC$parent.0 = parent.0.user
	}

	lmlogged = lm(log(parent) ~ t, data = kindata)

	for (kinmodel in kinmodels)
	{

		if (kinmodel == "SFO") {
			if (is.na(start.SFO$parent.0)) {
                                start.SFO$parent.0 = max(kindata$parent)
			}
			if (is.na(start.SFO$k)) {
				start.SFO$k = - coef(lmlogged)[["t"]]
			}
			kinfits[[kinmodel]] = try(
				nls(parent ~ SFO(t, parent.0, k),
					data = kindata, model = TRUE,
					start = start.SFO,
                                        algorithm = algorithm), silent=TRUE)
		}	
		k.est = ifelse(is.na(coef(kinfits$SFO)[["k"]]),
			-coef(lmlogged)[["t"]],
			coef(kinfits$SFO)[["k"]])
		if (kinmodel == "FOMC") {
			if (is.na(start.FOMC$parent.0)) {
                                start.FOMC$parent.0 = max(kindata$parent)
			}
			if (is.na(start.FOMC$alpha)) {
				start.FOMC$alpha = 1
			}
			if (is.na(start.FOMC$beta)) {
				start.FOMC$beta = start.FOMC$alpha / k.est 
			}
			kinfits[[kinmodel]] = try(
				nls(parent ~ FOMC(t, parent.0, alpha, beta),
					data = kindata, model = TRUE,
					start = start.FOMC,
                                        algorithm = algorithm), silent=TRUE)
		}	
		if (kinmodel == "DFOP") {
			if (is.na(start.DFOP$parent.0)) {
                                start.DFOP$parent.0 = max(kindata$parent)
			}
			if (is.na(start.DFOP$k1)) {
				start.DFOP$k1 = k.est * 2
			}
			if (is.na(start.DFOP$k2)) {
				start.DFOP$k2 = k.est / 2
			}
			if (is.na(start.DFOP$g)) {
				start.DFOP$g = 0.5
			}
			kinfits[[kinmodel]] = try(
				nls(parent ~ DFOP(t, parent.0, k1, k2, g),
					data = kindata, model = TRUE,
					start = start.DFOP,
                                        algorithm = algorithm), silent=TRUE)
		}	
		if (kinmodel == "HS") {
			if (is.na(start.HS$parent.0)) {
                                start.HS$parent.0 = max(kindata$parent)
			}
			if (is.na(start.HS$k1)) {
				start.HS$k1 = k.est
			}
			if (is.na(start.HS$k2)) {
				start.HS$k2 = k.est / 10
			}
			if (is.na(start.HS$tb)) {
				start.HS$tb = 0.05 * max(kindata$t)
			}
			kinfits[[kinmodel]] = try(
				nls(parent ~ HS(t, parent.0, k1, k2, tb),
					data = kindata, model = TRUE,
					start = start.HS,
                                        algorithm = algorithm), silent=TRUE)
		}	
	}
	return(kinfits)		
}
