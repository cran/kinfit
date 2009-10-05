kinobjects<- function(parent, type, systems,
        layers = NA, sampling_times = NA)
{
        kinobjects <- list()
        for (system in systems) {
            kinobjects[[system]] <- kinobject(parent = parent, 
                    type = type, system = system)
            if (!is.na(layers[1])) kinobjects[[system]]$layers = layers
            if (!is.na(sampling_times[1])) {
                    kinobjects[[system]]$sampling_times = layers
            }
        }
        return(kinobjects)
}
