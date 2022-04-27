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
input <- prepare_wham_input(basic_info = groundfish_info)

# run starter input
om <- fit_wham(input, do.fit = FALSE, MakeADFun.silent = TRUE)

#simulate data from operating model
#all RE and data are simulated
sim_input <- list()
# sim_input[[1]] has no data modification
set.seed(8675309) #use same seed for all operating models?
sim_input[[1]] = lapply(1:nsim, function(x) {
  input_i = input
  sim = om$simulate(complete=TRUE)
  input_i$data = sim
  return(input_i)
})
# sim_input[[2]] has catch data modification
set.seed(8675309) #use same seed for all operating models?
agg_catch_multiplier <- create_agg_catch_multiplier(input)
sim_input[[2]] = lapply(1:nsim, function(x) {
  input_i = input
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
