# misspecification_example.R
# demonstrate how to use bias_data() function to misspecify a state-space model

library(wham)
library(tidyr)
library(dplyr)
library(here)
library(ggplot2)
source(file.path(here(), "common_code", "make_basic_info.R"))
source(file.path(here(), "common_code", "bias_data.R"))

# create directory for analysis, e.g.
write.dir <- file.path(here(),"misspecification_study", "results")

if(!exists("write.dir")) write.dir = getwd()
if(!dir.exists(write.dir)) dir.create(write.dir, recursive = T)
setwd(write.dir)

# how many simulations
nsim <- 10

# create input
groundfish_info <- make_basic_info()

selectivity_om = list(
  model = c(rep("logistic", groundfish_info$n_fleets),rep("logistic", groundfish_info$n_indices)),
  initial_pars = rep(list(c(5,1)), groundfish_info$n_fleets + groundfish_info$n_indices)) #fleet, index

M_om = list(initial_means = rep(0.2, length(groundfish_info$ages)))

NAA_re_om = list(
  N1_pars = exp(10)*exp(-(0:(length(groundfish_info$ages)-1))*M_om$initial_means[1]),
  sigma = "rec", #random about mean
  cor="iid", #random effects are independent
  use_steepness = 0,
  recruit_model = 2, #random effects with a constant mean
  recruit_pars = exp(10)
)

input <- prepare_wham_input(basic_info = groundfish_info, NAA_re = NAA_re_om, selectivity = selectivity_om, M = M_om)

# run starter input
om <- fit_wham(input, do.fit = FALSE, MakeADFun.silent = TRUE)

#fit a SCAA model instead of Recruits as random effects?
NAA_re_em <- list(
  sigma = "rec", #random about mean
  cor="iid", #random effects are independent
  use_steepness = 0,
  recruit_model = 2, #random effects with a constant mean
  recruit_pars = exp(10),
  N1_pars = exp(10)*exp(-(0:(length(groundfish_info$ages)-1))*M_om$initial_means[1])
)
em_input <- prepare_wham_input(basic_info = groundfish_info, NAA_re = NAA_re_em, selectivity = selectivity_om, M = M_om)
agg_catch_multiplier <- create_agg_catch_multiplier(em_input)

#simulate data from operating model
#all RE and data are simulated
sim_input <- list()
# sim_input[[1]] has no data modification
set.seed(8675309) #use same seed for all operating models?
sim_input[[1]] = lapply(1:nsim, function(x) {
  input_i = em_input #changed to not estimate recruits as RE
  sim = om$simulate(complete=TRUE)
  input_i$data = sim
  return(input_i)
})


# sim_input[[2]] has catch data modification
set.seed(8675309) #use same seed for all operating models?
sim_input[[2]] = lapply(1:nsim, function(x) {
  input_i = em_input #changed to not estimate recruits as RE
  sim = om$simulate(complete=TRUE)
  sim <- bias_data(sim, multiply_agg_catch_flag=TRUE, agg_catch_multiplier=agg_catch_multiplier)
  input_i$data = sim
  return(input_i)
})

# see if it worked - it did!
df <- tibble(Source = character(),
             Year = integer(),
             Sim = integer(),
             agg_catch = double())
for(i in 1:nsim){
  thisdf1 <- tibble(Source = "orig",
                    Year = input$years,
                    Sim = i,
                    agg_catch = sim_input[[1]][[i]]$data$agg_catch)
  thisdf2 <- tibble(Source = "modified",
                    Year = input$years,
                    Sim = i,
                    agg_catch = sim_input[[2]][[i]]$data$agg_catch)
  df <- rbind(df, thisdf1, thisdf2)
}
ggplot(df, aes(x=Year, y=agg_catch, color=Source)) +
  geom_point() +
  geom_line() +
  facet_wrap(~Sim, ncol=2) +
  ylim(0, NA) +
  theme_bw()


# run the models
em_fits = list()
for(m in 1:2){
  em_fits[[m]] = lapply(1:nsim, function(x){
    cat(paste("model:",m, "fit:", x, "start \n"))
    out = fit_wham(sim_input[[m]][[x]], do.osa = FALSE, MakeADFun.silent = TRUE, retro.silent = TRUE, save.sdrep = FALSE)
    cat(paste("model:",m, "fit:", x, "done \n"))
    return(out)
  })
}
saveRDS(em_fits, file.path(write.dir, "em_fits.RDS"))
sim_input[[1]][[1]]$data$obsvec[sim_input[[1]][[1]]$data$keep_C+1]/sim_input[[2]][[1]]$data$obsvec[sim_input[[2]][[1]]$data$keep_C+1]
temp = sim_input[[1]][[1]]
temp$par$F_devs[] = 0
out = fit_wham(temp, do.osa = FALSE, do.retro = FALSE)
mohns_rho(out)
temp = sim_input[[2]][[1]]
temp$par$F_devs[] = 0
out2 = fit_wham(temp, do.osa = FALSE, do.retro = FALSE)

out = fit_wham(sim_input[[1]][[1]], do.osa = FALSE, do.retro = FALSE)
out2 = fit_wham(sim_input[[2]][[1]], do.osa = FALSE, do.retro = FALSE)
# hmmm, something did not work, getting same results from original and modified catch
mohns_rho(em_fits[[1]][[4]])
mohns_rho(em_fits[[2]][[4]])

em_fits[[1]][[5]]$report()$SSB
em_fits[[2]][[5]]$report()$SSB

