kinwrite.KinGUI <- function(kinobject, file, comment=NA)
{
	sink(file)
	cat("Version:\t1.1\n")
	cat("Project:\t", kinobject$parent, "\n", sep = "")
	cat("Testsystem:\t", kinobject$type, "\n", sep = "")
	cat("Comment:\t", comment, "\n", sep = "")
	write.table(kinobject$data, sep = "\t", na = "NaN", 
                quote = FALSE, row.names = FALSE)
	sink()
}
