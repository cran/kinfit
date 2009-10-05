kinobject <- function(parent, type, system, 
        layers = NA, sampling_times = NA)
{
        kinobject <- list(parent = parent, 
                type = type, system = system)
        if (!is.na(layers[1])) kinobject$layers = layers
        if (!is.na(sampling_times[1])) {
                kinobject$sampling_times = layers
        }
        return(kinobject)
}
