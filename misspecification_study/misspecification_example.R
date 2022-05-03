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
nsim <- 4

# create input
groundfish_info <- make_basic_info()

gf_selectivity = list(
  model = c(rep("logistic", groundfish_info$n_fleets),rep("logistic", groundfish_info$n_indices)),
  initial_pars = rep(list(c(5,1)), groundfish_info$n_fleets + groundfish_info$n_indices)) #fleet, index

gf_M = list(initial_means = rep(0.2, length(groundfish_info$ages)))

gf_NAA_re = list(
  N1_pars = exp(10)*exp(-(0:(length(groundfish_info$ages)-1))*gf_M$initial_means[1]),
  sigma = "rec", #random about mean
  cor="iid", #random effects are independent
  use_steepness = 0,
  recruit_model = 2, #random effects with a constant mean
  recruit_pars = exp(10)
)

input <- prepare_wham_input(basic_info = groundfish_info, selectivity = gf_selectivity, NAA_re = gf_NAA_re, M= gf_M)

# run starter input
om <- fit_wham(input, do.fit = FALSE, MakeADFun.silent = TRUE)

# set up estimation model
em_input <- input # self-test

#simulate data from operating model
#all RE and data are simulated
sim_input <- list()
# sim_input[[1]] has no data modification
set.seed(8675309) #use same seed for all operating models?
sim_input[[1]] = lapply(1:nsim, function(x) {
  input_i = em_input #maybe changed to not estimate recruits as RE
  sim = om$simulate(complete=TRUE)
  input_i$data = sim
  return(input_i)
})


# sim_input[[2]] has catch data modification
set.seed(8675309) #use same seed for all operating models?
agg_catch_multiplier <- create_agg_catch_multiplier(input, multiplier=0.25)
sim_input[[2]] = lapply(1:nsim, function(x) {
  input_i = em_input #maybe changed to not estimate recruits as RE
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

####################################
# confirm input data and exp(obsvec) the same
data.check <- matrix(NA, nrow=nsim, ncol=2)
for (i in 1:nsim){
  data.check[i,1] <- all.equal(as.vector(sim_input[[1]][[i]]$data$agg_catch), exp(sim_input[[1]][[i]]$data$obsvec[sim_input[[1]][[i]]$data$keep_C+1]))
  data.check[i,2] <- all.equal(as.vector(sim_input[[2]][[i]]$data$agg_catch), exp(sim_input[[2]][[i]]$data$obsvec[sim_input[[2]][[i]]$data$keep_C+1]))
}
data.check
#test fitting of one data set
#look at ratio of osbvec agg_catches
exp(sim_input[[2]][[1]]$data$obsvec[sim_input[[2]][[1]]$data$keep_C+1]-sim_input[[1]][[1]]$data$obsvec[sim_input[[1]][[1]]$data$keep_C+1])  
#out = fit_wham(sim_input[[1]][[1]], do.osa = FALSE, do.retro = FALSE)
#out2 = fit_wham(sim_input[[2]][[1]], do.osa = FALSE, do.retro = FALSE)
#cbind(out$years, out$rep$SSB, out2$rep$SSB) # getting different results, so working
####################################

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


# compare SSB time series
df2 <- tibble(Source = character(),
              Year = integer(),
              Sim = integer(),
              SSB = double())

for(i in 1:nsim){
  thisdf1 <- tibble(Source = "orig",
                    Year = input$years,
                    Sim = i,
                    SSB = em_fits[[1]][[i]]$report()$SSB)
  thisdf2 <- tibble(Source = "modified",
                    Year = input$years,
                    Sim = i,
                    SSB = em_fits[[2]][[i]]$report()$SSB)
  df2 <- rbind(df2, thisdf1, thisdf2)
}
ggplot(df2, aes(x=Year, y=SSB, color=Source)) +
  geom_point() +
  geom_line() +
  facet_wrap(~Sim, ncol=2) +
  ylim(0, NA) +
  theme_bw()

rhoSSB <- list()
for (m in 1:2){
  rhoSSB[[m]] <- lapply(1:nsim, function(x){
    out <- mohns_rho(em_fits[[m]][[x]])[1]
  return(out)
  })
}
rhoSSB
