misclassification.function <- function(times, status, screening_time, screening_crit, use_surrogate_times=TRUE) {
    event_times <- times
    if(use_surrogate_times) {
        if(screening_crit=="low risk") {
            ## Pessimistically assume individual died when censored
            ## Noting to do
        } else {
            ## Pessimistically assume individual never died
            event_times[status==0] <- Inf
        }
    }
    if(screening_crit=="low risk") {
        errors <- (event_times <= screening_time)
    } else {
        errors <- (event_times > screening_time)
    }
    return(errors)
}
