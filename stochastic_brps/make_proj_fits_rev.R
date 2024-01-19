devtools::load_all("c:/work/wham_master/wham")
source("results_plotting_functions.R")
#remotes::install_github("timjmiller/wham", dependencies=TRUE, ref = "devel", INSTALL_opts=c("--no-multiarch"))
#library(wham) #from devel branch

# Load dat
asap3 <- read_asap3_dat("WGOM_COD_ASAP_2023_SEL3_2021_c.dat")
sel <- list(model = c("logistic", "logistic", "age-specific", "logistic", "logistic", "age-specific",
    rep("age-specific", 10)),
    re = c(rep("none", 6), 
    rep("none", 10)))

sel$initial_pars = list(c(rep(9/2, 2)), 
    c(rep(9/2, 2)), 
    c(0.5, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0, 0.5, 0.5),
    c(rep(9/2, 2)), 
    c(rep(9/2, 2)), 
    c(0.5, 0.5, 1.0, 1.0, 1.0, 1.0, 1.0, 0.5, 0.5),
    c(rep(0.5, 4), rep(1, 5)), #<< start surveys
    c(rep(0.5, 4), rep(1, 5)),
    c(1, 0.5, 0.5, rep(0.5, 3),
     rep(0, 3)),
    c(1, 0.5, 0.5, 1, 1, 1,
     0, 0, 0),
    c(0.5, 0.5, 0.5, rep(1, 6)),
    c(rep(0.5, 4), rep(1, 5)),
    c(rep(0.5, 3), rep(1, 6)),
    c(rep(0.5, 2), rep(1, 7)),
    c(rep(0.5, 7), rep(1, 2)),
    c(rep(0.5, 7), rep(1, 2)))
sel$fix_pars = list(NULL, 
   NULL,  
   c(5:7), 
   NULL,
   NULL,  
   c(3:7), 
   5:9,  #<< start surveys
   5:9, 
   c(1, 6:9),
   c(1, 4:9),
   c(4:9),
   5:9,
   4:9,
   3:9,
   8:9,
   8:9)


# Setup input
input <- prepare_wham_input(asap3, selectivity = sel,
    recruit_model = 2, 
    model_name = "GOM cod 2023 RT sel3ahl Run",
    NAA_re = list(sigma="rec+1", cor = "2dar1"),
    age_comp = list(fleets = rep("dirichlet-miss0", 2), 
    indices = rep("dirichlet-miss0", 10)),
    basic_info = list(simulate_process_error = rep(FALSE, 5)#,
    #XSPR_R_opt = 1 #1: annual R estimates
    ))
input$data$XSRP_R_opt <- 1
input$data$XSRP_R_opt #2 means average R over some time period (which would only be up to the last year of the unprojected model)
input$data$simulate_state[] <- 1
input$data$simulate_period[1] <- 0
# Fit model
mod_bc   <- fit_wham(input, do.osa = F, do.retro = F)
saveRDS(mod, "mod_bc.RDS")
mod_bc$input$data$simulate_state[] <- 1
mod_bc$input$data$simulate_period[1] <- 0
mod_bc_proj  <- project_wham(mod_bc, proj.opts  = list(
    n.yrs = 100, 
    use.FXSPR = T),
    save.sdrep = F,
    do.sdrep = F,
)
saveRDS(mod_bc_proj, "mod_bc_proj.RDS")

mod_bc_proj$input$data$simulate_state
mod_bc_proj$input$data$simulate_period

mod_bc_proj$rep$NAA[40:45,1]
mod_bc_proj$simulate()$NAA[40:45,1]

mod_bc_proj$rep$NAA[140,1]
mod_bc_proj$simulate()$NAA[140,1]

mod_bc_F0_proj  <- project_wham(mod_bc, proj.opts  = list(
    n.yrs = 100, 
    proj.F = rep(0.00001,100),
    use.FXSPR = F),
    save.sdrep = F,
    do.sdrep = F,
)
saveRDS(mod_bc_F0_proj, "mod_bc_F0_proj.RDS")

mod_bc_proj_sdrep_bc <- TMB::sdreport(mod_bc_proj, bias.correct = TRUE)#, bias.correct.control = list(sd = TRUE, split = NULL, nsplit = NULL))
mod_bc_F0_proj_sdrep_bc <- TMB::sdreport(mod_bc_F0_proj, bias.correct = TRUE)

saveRDS(mod_bc_proj_sdrep_bc, "mod_bc_proj_sdrep_bc.RDS")
saveRDS(mod_bc_F0_proj_sdrep_bc, "mod_bc_F0_proj_sdrep_bc.RDS")

input_nbc <- input
input_nbc$data$bias_correct_pe[] <- 0
input_nbc$data$bias_correct_oe[] <- 0

mod_nbc   <- fit_wham(input_nbc, do.osa = F, do.retro = F)
saveRDS(mod, "mod_nbc.RDS")

mod_nbc_proj  <- project_wham(mod_nbc, proj.opts  = list(
    n.yrs = 100, 
    use.FXSPR = T),
    save.sdrep = F,
    do.sdrep = F,
)
saveRDS(mod_nbc_proj, "mod_nbc_proj.RDS")

mod_nbc_F0_proj  <- project_wham(mod_nbc, proj.opts  = list(
    n.yrs = 100, 
    proj.F = rep(0.00001,100),
    use.FXSPR = F),
    save.sdrep = F,
    do.sdrep = F,
)
saveRDS(mod_nbc_F0_proj, "mod_nbc_F0_proj.RDS")

mod_nbc_proj_sdrep_bc <- TMB::sdreport(mod_nbc_proj, bias.correct = TRUE)#, bias.correct.control = list(sd = TRUE, split = NULL, nsplit = NULL))
mod_nbc_F0_proj_sdrep_bc <- TMB::sdreport(mod_nbc_F0_proj, bias.correct = TRUE)

saveRDS(mod_nbc_proj_sdrep_bc, "mod_nbc_proj_sdrep_bc.RDS")
saveRDS(mod_nbc_F0_proj_sdrep_bc, "mod_nbc_F0_proj_sdrep_bc.RDS")


types <- c("bc")
mods <- paste0("mod_", types, "_proj")
mods <- c(mods, paste0("mod_", types, "_F0_proj"))

set.seed(123)
seeds <- sample(-1e6:1e6, 1000)
proj_sims <- lapply(mods, function(z){
    sims <- lapply(1:1000, function(x) {
        print(x)
        set.seed(seeds[x])
        temp <- get(z)$simulate(complete = TRUE)
        temp$pred_SSB <- get_SSB(temp, pred=TRUE)
        return(temp[c("NAA","SSB","pred_NAA","pred_SSB")])
    })
    return(sims)
})
names(proj_sims) <- mods
saveRDS(proj_sims, "mod_bc_proj_sims.RDS")
get(mods[1])$rep$NAA[1:5,1]
get(mods[1])$simulate()$NAA[1:5,1]
get(mods[1])$simulate()$sims[1:5,1]
get(mods[1])$rep$pred_NAA[1:5,1]
get(mods[1])$simulate()$pred_NAA[1:5,1]
get(mods[1])$simulate()$sims[1:5,10]
get(mods[1])$simulate()$NAA_devs[1:5,1]


x <- cbind(
  exp(as.list(mod_bc_proj_sdrep_bc, report = TRUE, "Est. (bias.correct)")$log_NAA[,1]),
  as.list(mod_bc_proj_sdrep_bc, report = TRUE, "Est. (bias.correct)")$NAA[,1]
)

n_ages <- input$data$n_ages
proj_median_NAA <- lapply(proj_sims, function(x) sapply(1:n_ages, function(z) apply(sapply(x, function(y) y$NAA[,z]),1,median)))
proj_median_R <- lapply(proj_sims, function(x) apply(sapply(x, function(y) y$NAA[,1]),1,median))
proj_mean_R <- lapply(proj_sims, function(x) apply(sapply(x, function(y) y$NAA[,1]),1,mean))
proj_mean_NAA <- lapply(proj_sims, function(x) sapply(1:n_ages, function(z) apply(sapply(x, function(y) y$NAA[,z]),1,mean)))

yrs <- 1:41
plot(mod_bc_proj$years_full[yrs], x[yrs,1], type = 'n', ylim = c(0,max(x[yrs,])))
lines(mod_bc_proj$years_full[yrs], x[yrs,1], lwd = 2)
lines(mod_bc_proj$years_full[yrs], x[yrs,2], lwd = 2, col = 'red')
lines(mod_bc_proj$years_full[yrs], mod_bc_proj$rep$pred_NAA[yrs,1], lwd = 2, col = 'green')
lines(mod_bc_proj$years_full[yrs], mod_bc_proj$rep$NAA[yrs,1], lwd = 2, col = 'brown')

x <- cbind(
  exp(as.list(mod_bc_proj_sdrep_bc, report = TRUE, "Est. (bias.correct)")$log_SSB),
  as.list(mod_bc_proj_sdrep_bc, report = TRUE, "Est. (bias.correct)")$SSB
)

plot(mod_bc_proj$years_full, x[,1], type = 'n', ylim = c(0,max(x)))
lines(mod_bc_proj$years_full, x[,1], lwd = 2)
lines(mod_bc_proj$years_full, x[,2], lwd = 2, col = 'red')
lines(mod_bc_proj$years_full, proj_median_NAA[[1]][,1], lwd = 2, col = 'blue')
lines(mod_bc_proj$years_full, proj_mean_NAA[[1]][,1], lwd = 2, col = 'orange')

proj_sims[[1]][[1]]$NAA[,1]
proj_sims[[1]][[2]]$NAA[,1]
yrs <- 1:41
plot(mod_bc_proj$years, x[yrs,2] - x[yrs,1], type = 'n')
lines(mod_bc_proj$years, x[yrs,2] - x[yrs,1], lwd = 2)
lines(mod_bc_proj$years, x[yrs,2], lwd = 2, col = 'red')

length(unique(sapply(proj_sims[[1]], function(x) x$NAA[3,1])))
length(mod_bc_proj$years)


