

## if(file.exists("c:/Users/timothy.j.miller")) {
##   library(wham, lib.loc = "c:/work/wham/old_packages/77bbd94")
## } else library(wham) #make sure to use the right version of wham

## requires special growth branch
## devtools::install_github("GiancarloMCorrea/wham", ref='growth')
## devtools::install_github("timjmiller/wham", dependencies=TRUE, ref="77bbd94")
## devtools::install_github("timjmiller/wham", dependencies=TRUE, ref="growth")
library(tidyr)
library(dplyr)
library(here)
library(ggplot2)
theme_set(theme_bw())
library(wham)
source(file.path(here(), "common_code", "make_basic_info.R"))
source(file.path(here(), "common_code", "set_NAA.R"))
source(file.path(here(), "common_code", "set_M.R"))
source(file.path(here(), "common_code", "set_q.R"))
source(file.path(here(), "common_code", "set_ecov.R"))
source(file.path(here(), "common_code", "set_selectivity.R"))
source(file.path(here(), "common_code", "set_simulation_options.R"))
source(file.path(here(), "common_code", "get_FMSY.R"))
source(file.path(here(), "Ecov_study", "growth", "code", "make_om.R"))
source(file.path(here(), "Ecov_study", "growth", "code", "sim_management.R"))
source("make_om.R")

## verify_version()

write.dir <- file.path(here(),"Ecov_study", "growth", "inputs")

#### This breaks b/c the input are missing some growth updates,
#### so just manually specify SRab which apparently all these are
#### used for?
## naa_om_inputs = readRDS(file.path(here::here(),"Project_0","inputs", "NAA_om_inputs.RDS"))
## ##SR parameters are the same for all naa_om models
## temp = fit_wham(naa_om_inputs[[1]], do.fit = FALSE, MakeADFun.silent = TRUE)
## SRab = exp(c(temp$rep$log_SR_a[1], temp$rep$log_SR_b[1]))
SRab <- c(5.955694e-01, 2.404283e-05)


#if(!exists("write.dir")) write.dir = getwd()
if(!dir.exists(write.dir)) dir.create(write.dir, recursive = T)
#setwd(write.dir)


#number of simulations for each scenario
#nsim = 1000
#nsim = 100

#Operating model factors
#NAA sigmas for each scenario
R_sig <- c(0.5)
NAA_sig <- c(0.3)
#F time series
M_sig <- 0.3
M_cor <- 0
Ecov_obs_sig <- c(0.1, 0.5)[1]
Ecov_re_sig <- c(0.1,0.5)[2]
# For example. 0 == (sigma = 1 internally)
Ecov_re_cor <- c(0, 0.9)[2] # This parameter is not being specified correctly in the OM. Make sure you use the right scale
Ecov_effect <- c(0, 0.2, 0.4) # This parameter is not being specified correctly in the OM
Fhist = c("H-MSY","MSY")[2]
#how much observation error
obs_error = c("L", "H")
NAA_M_re <- c("rec","rec+1", "rec+M")[1]
df.oms <- expand.grid(NAA_M_re = NAA_M_re,
  Ecov_obs_sig=Ecov_obs_sig, Ecov_re_sig=Ecov_re_sig, Ecov_re_cor=Ecov_re_cor, Ecov_effect = Ecov_effect,
  Fhist = Fhist, obs_error = obs_error, stringsAsFactors = FALSE)
#logistic-normal age comp SDs for L/H observation error (both indices AND CATCH!!!????)
L_N_sigma = c(L = 0.3, H = 1.5)
#(log) index SDs for L/H observation error
index_sigma = c(L = 0.1, H = 0.4)
index_NeffL <- c(H=50, L=200)

n.mods <- dim(df.oms)[1] #288 operating model scenarios
df.oms$Model <- paste0("om_",1:n.mods)
df.oms <- df.oms %>% select(Model, everything()) # moves Model to first col
# look at model table
df.oms
saveRDS(df.oms, file.path(here(),"Ecov_study", "growth", "inputs", "df.oms.RDS"))


gf_info = make_basic_info()
#selectivity is not changing
gf_selectivity = list(
  ## model = c(rep("logistic", gf_info$n_fleets),rep("logistic", gf_info$n_indices)),
  model = c("logistic", "logistic", "len-logistic"),
  initial_pars = list(c(5,1), c(5,1), c(30,4))) # original value: 65
#M set is not changing
gf_M = list(model = "age-specific",
  initial_means = rep(0.2, length(gf_info$ages))#,
  #re = "ar1_y"#, # This is needed to set up operating models with a single annual re, iid or ar1_y
  #sigma_vals = 0.1,
  #cor_vals = 0
  )
#NAA_re set up that can be changed for each OM scenario
gf_NAA_re = list(
  N1_pars = exp(10)*exp(-(0:(length(gf_info$ages)-1))*gf_M$initial_means[1]),
  #sigma = "rec", #random about mean
  cor="iid", #random effects are independent
  use_steepness = 0,
  #recruit_model = 2, #random effects with a constant mean
  recruit_model = 3, #B-H
  recruit_pars = SRab) #defined above from naa_om_inputs
gf_ecov <- list(
  label = "Ecov",
  process_model = "ar1",
  lag = 0,
  mean = cbind(rep(0, length(gf_info$years))),
  year = gf_info$years,
  ages = list(1:10),
  use_obs = cbind(rep(1, length(gf_info$years))),
  # process_mean_vals = 0, # what is this?
  where = list('growth'),
  where_subindex = 3 # on L1
)

### ------------------------------------------------------------
## Setup the growth structure of the model. Assume survey two
## collects marginal length compositions and uses logistic selex
## at length. Growth (k, Linf, CV1, CV10) is estimated
## internally. WAA is then derived from weight-length
## relationship (not estimated).

## These are are the growth values used to create WAA in the
## other models, so use them here
Linf <- 85
k <- 0.3
t0 <- 0
a_LW <- exp(-12.1)
b_LW <- 3.2
L <- Linf*(1-exp(-k*(1:10 - t0)))
W <- a_LW*L^b_LW
L[1] ## length at reference age (1)
CV <- .1
## growth CVs at age 1 and 10
CV*L[1]
CV*L[10]
gf_growth <- list(model='vB_classic', init_vals=c(k, Linf, L[1]),
                  est_pars=1:2, SD_vals=c(CV*L[1], CV*L[10]),
                  SD_est=1:2)
gf_LW <- list(init_vals=c(a_LW, b_LW))
## end of growth changes
### --------------------------------------------------

beta_vals <- list(rep(list(matrix(0,1,length(gf_info$ages))), 8))
# base_om = make_om(Fhist = "Fmsy", N1_state = "overfished", selectivity = gf_selectivity,
#     M = gf_M, NAA_re = gf_NAA_re, age_comp = "logistic-normal-miss0", brp_year = 1, eq_F_init = 0.3,
#     om_input = TRUE, max_mult_Fmsy = 1, min_mult_Fmsy = 1)

#make inputs for operating model (smaller objects to save, can recreate simulated data sets)
om_inputs = list()
for(i in 1:NROW(df.oms)){
  print(paste0("row ", i))
  NAA_re = gf_NAA_re
  #M
  M_i = gf_M
  if(df.oms$NAA_M_re[i] == "rec"){
    NAA_re$sigma = "rec"
    NAA_re$sigma_vals = R_sig
  }
  if(df.oms$NAA_M_re[i] == "rec+1"){
    NAA_re$sigma = "rec+1"
    NAA_re$sigma_vals = c(R_sig, NAA_sig)
  }
  if(df.oms$NAA_M_re[i] == "rec+M"){
    NAA_re$sigma = "rec"
    NAA_re$sigma_vals = R_sig
    ## M_i$re = "ar1_y"
    ## M_i$sigma_vals = M_sig * sqrt(1-M_cor^2) #defining marginal variance, but wham estimates conditional var.
    ## M_i$cor_vals = M_cor
  }
  ecov_i = gf_ecov
  ecov_i$logsigma = cbind(rep(log(df.oms$Ecov_obs_sig[i]), length(ecov_i$year)))
  ecov_i$process_sig_vals = df.oms$Ecov_re_sig[i]
  ecov_i$process_cor_vals = df.oms$Ecov_re_cor[i]*100
  if(df.oms$Ecov_effect[i] < 1e-7){
    ecov_i$how = 0
    ecov_i$where = "none"
  } else {
    ecov_i$how = 1
    ecov_i$where = "growth"
    # beta_vals_i = beta_vals
    # beta_vals_i[[1]][[5]][] <- df.oms$Ecov_effect[i]
    # ecov_i$beta_vals = beta_vals_i # wrong way to do it
  }
  Fhist. = "Fmsy"
  if(df.oms$Fhist[i] == "H-MSY") Fhist. = "H-L"
  max_mult = 2.5 # fishing at 2.5 x Fmsy
  if(Fhist. == "Fmsy") max_mult = 1
  min_mult = 1 # fishing at Fmsy
  om_inputs[[i]] <-
    make_om(Fhist = Fhist., N1_state = "overfished", selectivity = gf_selectivity,
            M = M_i, NAA_re = NAA_re, ecov = ecov_i,
            growth=gf_growth, LW=gf_LW,
            age_comp = "logistic-normal-miss0", brp_year = 1, eq_F_init = 0.3,
            om_input = TRUE, max_mult_Fmsy = max_mult, min_mult_Fmsy = min_mult,
            df.oms = df.oms[i,]) # I added the dfOM here, there are some things that you are specifying but cannot be passed to the
            # parameter section using prepare_wham_input, so they will need to be modified in the make_om function.
  #turn off bias correction
  om_inputs[[i]] = set_simulation_options(om_inputs[[i]], simulate_data = TRUE, simulate_process = TRUE, simulate_projection = TRUE,
    bias_correct_pe = FALSE, bias_correct_oe = FALSE)
  #set L-N SD parameters for catch and index age comp
  om_inputs[[i]]$par$catch_paa_pars[,1] = log(L_N_sigma[df.oms$obs_error[i]])
  om_inputs[[i]]$par$index_paa_pars[,1] = log(L_N_sigma[df.oms$obs_error[i]])
  om_inputs[[i]]$data$agg_catch_sigma[] = 0.1
  om_inputs[[i]]$data$agg_index_sigma[] = index_sigma[df.oms$obs_error[i]]
  #change Neff so scalar doesn't affect L-N SD
  om_inputs[[i]]$data$catch_Neff[] = 1
  om_inputs[[i]]$data$index_Neff[] = 1
  ## length comps can only use multinomial so set depending on
  ## the data scenario
  om_inputs[[i]]$data$index_NeffL[,2] <- index_NeffL[df.oms$obs_error[i]]
   ## this didn't get updated so do it manually here. THis
  ## shouldn't be necessary b/c internal to make_om it should get
  ## updated via set_ecov but that doesn't work right now
  ## om_inputs[[i]]$par$Ecov_beta[5,1,1,] <-
  ## df.oms$Ecov_effect[i]
  ## manually set the Ecov process pars based on AR1 model
   om_inputs[[i]]$par$Ecov_process_pars[,1] <- c(0,log(df.oms$Ecov_re_sig[i]), -log(2/(1+df.oms$Ecov_re_cor[i])-1))
}

## ## Not sure the ecov process pars are right??
## x <- -log(2/(1+df.oms$Ecov_re_cor[i])-1)
## -1+2/(1+exp(-x))
## lapply(om_inputs, function(x) x$par$Ecov_process_pars)
## om_inputs[[i]]$par$Ecov_process_pars[,1]

### test the new growth stuff is working in the OM
png('../plots/OM_basics.png', width=7, height=5, units='in', res=500)
x <- om_inputs[[1]]
log(gf_NAA_re$recruit_pars)
x$par$mean_rec_pars
test <- fit_wham(input=x, do.fit=FALSE)
par(mfrow=c(2,3), mar=c(2.5,2.5,.5,.5), mgp=c(1.5,.25,0), tck=-.02)
plot(x$data$lengths, test$rep$selAL[[3]][1,], type='b',
     xlab='length', ylab='selex')
## derived selex at age
selL <- test$rep$selAL[[3]][1,]
phi <- test$rep$phi_mat[1,,]
selA <- as.numeric(selL %*% phi)
plot(1:10, selA/max(selA), ylim=c(0,1), xlab='age', ylab='selex', type='b')
legend('bottomright', legend=c('survey 1', 'survey 2 derived'),
       col=2:1, lty=1)
lines(1:10, test$rep$selAL[[2]][1,], col=2, type='b')
plot(1:10, L, type='l', xlab='age', ylab='length')
lines(1:10, test$rep$LAA[1,], type='l', col=2)
legend('bottomright', legend=c('default', 'growth'), lty=1, col=1:2)
## matplot(cbind(test$rep$pred_waa[5,1,],W))
sim <- test$simulate(complete=TRUE)
plot(1:50, sim$index_pal[2,1,], type='l', xlab='length',
     ylab='proportion')
mtext(line=-1.5, text='Data L')
lines(test$rep$pred_index_pal[1,2,], type='b', col=2)
##legend('topleft', legend=c('expected', 'simulated'), lty=1, col=2:1)
## with more data
test2 <- fit_wham(input=om_inputs[[4]], do.fit=FALSE)
sim2 <- test2$simulate(complete=TRUE)
plot(1:50, sim2$index_pal[2,1,], type='l', xlab='length',
     ylab='proportion')
mtext(line=-1.5, text='Data H')
lines(test2$rep$pred_index_pal[1,2,], type='b', col=2)
##legend('topleft', legend=c('expected', 'simulated'), lty=1, col=2:1)
plot(1:10, test$rep$pred_waa[5,1,], type='b', col=2, xlab='age', ylab='weight')
lines(1:10, W, type='b')
legend('topleft', legend=c('default', 'growth'), lty=1, col=1:2)
dev.off()


## test that time-varying growth works
fits <- lapply(1:6, function(i){
  fit <- fit_wham(input=om_inputs[[i]], do.fit=FALSE)
  fit$simulate(complete=TRUE)
})
get_laa <- function(fit, om){
  x <- fit$LAA %>% reshape2::melt() %>% setNames(c('year', 'age', 'length')) %>%
    mutate(OM=om, year=year+1981, cohort=year-age)
  merge(x, df.oms, by.x='OM', by.y='Model') %>%
    mutate(Ecov_effectf=paste0('Ecov_effect=', Ecov_effect),
           Ecov_re_sigf=paste0('Ecov_re_sig=', Ecov_re_sig))
}
g <- lapply(1:6, function(i) get_laa(fits[[i]], df.oms$Model[i])) %>%
  bind_rows %>% filter(cohort>1981) %>%
  ggplot(aes(age, length, group=cohort)) + geom_line() +
  facet_grid(Ecov_effectf~obs_error)
ggsave('../plots/OM_laa.png', g, width=6, height=4)
get_waa <- function(fit, om){
  x <- fit$pred_waa[1,,] %>% reshape2::melt() %>% setNames(c('year', 'age', 'weight')) %>%
    mutate(OM=om, year=year+1981, cohort=year-age)
  merge(x, df.oms, by.x='OM', by.y='Model') %>%
    mutate(Ecov_effectf=paste0('Ecov_effect=', Ecov_effect),
           Ecov_re_sigf=paste0('Ecov_re_sig=', Ecov_re_sig))
}
g <- lapply(1:6, function(i) get_waa(fits[[i]], df.oms$Model[i])) %>%
  bind_rows %>% filter(cohort>1981) %>%
  ggplot(aes(age, weight, group=cohort)) + geom_line() +
  facet_grid(Ecov_effectf~obs_error)
ggsave('../plots/OM_waa.png', g, width=6, height=4)
get_ecov <- function(fit, om){
  x <- data.frame(Ecov=fit$Ecov_x[,1], OM=om, year=gf_info$years) %>%
    merge(df.oms, by.x='OM', by.y='Model')  %>%
        mutate(Ecov_effectf=paste0('Ecov_effect=', Ecov_effect),
               Ecov_re_sigf=paste0('Ecov_re_sig=', Ecov_re_sig))
}
g <- lapply(1:6, function(i) get_ecov(fits[[i]], df.oms$Model[i])) %>%
  bind_rows %>%
  ggplot(aes(year, Ecov, group=Ecov_effect)) + geom_line() +
  facet_grid(.~obs_error)
ggsave('../plots/OM_ecov.png', g, width=6, height=4)
print(g)


## ##start out at MSY and continue
## om_msy = make_om(Fhist = "Fmsy", N1_state = "overfished", selectivity = gf_selectivity,
##   M = gf_M, NAA_re = NAA_re, brp_year = 1, eq_F_init = 0.3, max_mult_Fmsy = 1)
## input = om_msy#$input
## input$par$log_NAA_sigma[] = -100 #no process error
## temp <- fit_wham(input, do.fit = FALSE, MakeADFun.silent = TRUE)
## sim = temp$simulate(complete=TRUE)
## sim$NAA #NAA should stay the same throughout time series
## #start out overfished and continue for 20 years
## om_msy = make_om(Fhist = "H-L", N1_state = "overfished", selectivity = gf_selectivity,
##   M = gf_M, NAA_re = NAA_re, brp_year = 1, eq_F_init = 0.3, max_mult_Fmsy = 2.5)
## input = om_msy#$input
## input$par$log_NAA_sigma[] = -100 #no process error
## temp <- fit_wham(input, do.fit = FALSE, MakeADFun.silent = TRUE)
## sim = temp$simulate(complete=TRUE)
## sim$NAA #NAA should stay the same throughout time series

saveRDS(om_inputs, file.path(here(), "Ecov_study", "growth", "inputs", "om_inputs.RDS"))

#I don't think we want to use the same (e.g. 1000) seeds for everything.
set.seed(8675309)
seeds = sample(x = (-1e9):(1e9), size = NROW(df.oms)*1000, replace = FALSE)
seeds <- lapply(1:NROW(df.oms), function(x) seeds[(1:1000) + 1000*(x-1)])
saveRDS(seeds, file.path(here(), "Ecov_study", "growth", "inputs","seeds.RDS"))
seeds = readRDS(file.path(here::here(), "Ecov_study", "growth", "inputs","seeds.RDS"))
