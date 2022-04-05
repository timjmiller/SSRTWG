make_basic_info <- function(base_years = 1982:2021, ages = 1:10, Fhist = "updown", n_feedback_years = 0) { #changed years
    info <- list()
    info$ages <- ages
    info$years <- as.integer(base_years[1] - 1 + 1:(length(base_years) + n_feedback_years))
    info$n_fleets <- 1 
    info$n_indices <- 1
    na <- length(info$ages)
    ny <- length(info$years)
    
    nby <- length(base_years)
    mid <- floor(nby/2)
    #up then down
    if(Fhist == "updown") info$F <- matrix(0.2 + c(seq(0,0.4,length.out = mid),seq(0.4,0,length.out=nby-mid)),nby, info$n_fleets)

    #down then up
    if(Fhist == "downup") info$F <- matrix(0.2 + c(seq(0.4,0,length.out = mid),seq(0,0.4,length.out=nby-mid)),nby, info$n_fleets)

    if(n_feedback_years>0) info$F <- rbind(info$F, info$F[rep(nby, n_feedback_years),, drop = F]) #same F as terminal year for feedback period

    info$catch_cv <- matrix(0.1, ny, info$n_fleets)
    info$catch_Neff <- matrix(200, ny, info$n_fleets)
    
    info$index_cv <- matrix(0.3, ny, info$n_indices)
    info$index_Neff <- matrix(100, ny, info$n_indices)
    info$fracyr_indices <- matrix(0.5, ny, info$n_indices)
    info$index_units <- rep(1, length(info$n_indices)) #biomass
    info$index_paa_units <- rep(2, length(info$n_indices)) #abundance
    
    info$maturity <- t(matrix(1/(1 + exp(-1*(1:na - na/2))), na, ny))

    L <- 100*(1-exp(-0.3*(1:na - 0)))
    W <- exp(-11)*L^3
    nwaa <- info$n_indices + info$n_fleets + 2
    info$waa <- array(NA, dim = c(nwaa, ny, na))
    for(i in 1:nwaa) info$waa[i,,] <- t(matrix(W, na, ny))

    info$fracyr_SSB <- rep(0.25,ny)
    info$q <- rep(0.3, info$n_indices)

    info$selblock_pointer_fleets <- t(matrix(1:info$n_fleets, info$n_fleets, ny))
    info$selblock_pointer_indices <- t(matrix(info$n_fleets + 1:info$n_indices, info$n_indices, ny))

    #Don't bias correct anything
    info$bias_correct_process = FALSE
    info$bias_correct_observation = FALSE
    return(info)
}
