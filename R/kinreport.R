kinreport <- function(kinobject, file = NA, vcov = FALSE, endpoint.digits = 1)
{
	if (!is.na(file)) {
		sink(file, split=TRUE)
	}

	cat("Parent compound: ", kinobject$parent, "\n")
        if (!is.null(kinobject$label)) cat("Label position:\t\t", kinobject$label, "\n")
	cat("Study type:      ", kinobject$type, "\n")
	cat("System:          ", kinobject$system, "\n")
	if (!is.null(kinobject$source)) {
          cat("Source:          ", kinobject$source, "\n")
        }
	cat("\n")
	fit.names <- names(kinobject$fits)
	for (kinmodel in fit.names)
	{
                m <- kinobject$fits[[kinmodel]]
                if (!(class(m) == "try-error")) {
                    cat("\n\n---\n")
                    cat("Nonlinear least squares fit of the", kinmodel, "model\n\n")
                    cat("Parameter estimation:\t")
                    s <- summary(m)
                    df <- s$df[2]
                    p <- 1 - pt(s$parameters[,3], df = df)
                    parms <- cbind(s$parameters[,c(1,2,3)], "Pr(>t)" = p)
                    cat("\n")
                    print(parms, digits=3)
                    cat("\n")
                    if(vcov)
                    {
                        cat("Variance-covariance matrix:\n")
                        print(vcov(m))
                        cat("\n")
                    }
                    cat("Chi2 error estimation:\t", 
                            round(100 * kinobject$results$stats[kinmodel, "err.min"], digits=2), 
                            " %\n", sep="")
                    cat("\n")
                }
	}
	cat("\n\n---\n")
	cat("Endpoint estimates\n\n")
	print(round(kinobject$results$results, digits=endpoint.digits))

	if (!is.na(file)) sink()
}
