# devtools::install_github("timjmiller/wham", dependencies=TRUE, ref="devel")
library(wham)
library(here)
source(file.path(here(),"common_code", "set_selectivity.R")) # overwrite wham's set_selectivity() with ours
write.dir <- file.path(here(),"selectivity", "results")

if(!exists("write.dir")) write.dir = getwd()
if(!dir.exists(write.dir)) dir.create(write.dir, recursive = T)
setwd(write.dir)

#number of simulations for each scenario
#nsim = 1000
nsim = 10

#NAA sigmas for each scenario
R_sig <- c(0.4,0.8,1.2)
NAA_sig <- c(0.1,0.3,0.5)
ar1_y <- c(NA,0.4,0.8)
ar1_a <- c(NA,0.4,0.8)

df.mods <- expand.grid(Recruitment = 2, R_sig = R_sig, NAA_sig = NA, ar1_y = ar1_y, ar1_a = NA)
df.mods <- rbind(df.mods, expand.grid(Recruitment = 2, R_sig = R_sig, NAA_sig = NAA_sig, ar1_y = ar1_y, ar1_a = ar1_a))

n.mods <- dim(df.mods)[1] #90 scenarios!
df.mods$Model <- paste0("m_",1:n.mods)
df.mods <- df.mods %>% select(Model, everything()) # moves Model to first col

# look at model table
df.mods

#make a list of input components that prepare_wham_input can use to generate an input for fit_wham
#this will be the generic flatfish/groundfish life histor information
source(file.path(here(),"common_code", "make_basic_info.R"))
groundfish_info <- make_basic_info()

selectivity = list(model = c(rep("logistic", groundfish_info$n_fleets),rep("logistic", groundfish_info$n_indices)),
    initial_pars = rep(list(c(5,1)), groundfish_info$n_fleets + groundfish_info$n_indices)) #fleet, index

M = list(initial_means = rep(0.2, length(groundfish_info$ages)))

#initial numbers at age

NAA_re = list(N1_pars = exp(10)*exp(-(0:(length(groundfish_info$ages)-1))*M$initial_means[1]))

input = prepare_wham_input(basic_info = groundfish_info, selectivity = selectivity, M = M, NAA_re = NAA_re)


sim_input = list()
# sim the 90 models!
for(m in 1:n.mods){

    #initial numbers at age
    NAA_re = list(N1_pars = exp(10)*exp(-(0:(length(groundfish_info$ages)-1))*M$initial_means[1]))
    input = prepare_wham_input(basic_info = groundfish_info, selectivity = selectivity, M = M, NAA_re = NAA_re)
    
    NAA_re$recruit_model = df.mods$Recruitment[m] #random effects with a constant mean
    NAA_re$recruit_pars = exp(10) #mean recruitment, if recruit model changes, need to change this
    
    #whether RE on recruitment or on all NAA
    if(is.na(df.mods$NAA_sig[m])) {
      NAA_re$sigma = "rec" #random about mean
      NAA_re$sigma_vals = df.mods$R_sig[m]
    } else {
      NAA_re$sigma = "rec+1"
      NAA_re$sigma_vals = c(df.mods$R_sig[m], df.mods$NAA_sig[m])
    }
    
    #what type of correlation structure
    if(!is.na(df.mods$ar1_a[m]) & !is.na(df.mods$ar1_y[m])) {
        NAA_re$cor = "2dar1"
        NAA_re$cor_vals = c(df.mods$ar1_a[m],df.mods$ar1_y[m])
    }
    if(is.na(df.mods$ar1_a[m]) & !is.na(df.mods$ar1_y[m])) {
        NAA_re$cor = "ar1_y"
        NAA_re$cor_vals = df.mods$ar1_y[m]
    }
    if(!is.na(df.mods$ar1_a[m]) & is.na(df.mods$ar1_y[m])) {
        NAA_re$cor = "ar1_a"
        NAA_re$cor_vals = df.mods$ar1_a[m]
    }
    if(is.na(df.mods$ar1_a[m]) & is.na(df.mods$ar1_y[m])) {
        NAA_re$cor = "iid"
    }
    input = set_NAA(input, NAA_re)

    print(paste0("m: ", m))
    om = fit_wham(input, do.fit = FALSE, MakeADFun.silent = TRUE)


    #simulate data from operating model
    set.seed(8675309) #use same seed for all operating models?
    #all RE and data are simulated
    sim_input[[m]] = lapply(1:nsim, function(x) {
        input_i = input
        sim = om$simulate(complete=TRUE)
        input_i$data = sim
        return(input_i)
    })

}

#save simulated data sets including true population, catches, etc. to google drive?
saveRDS(sim_input, file.path(write.dir, "om_sim_data.RDS"))
#sim_input <- readRDS(file.path(write.dir, "om_sim_data.RDS"))

#This will now generate a list of inputs to estimate models that match the operating model.
#initial values are ommitted to start at "generic" starting values for estimation.
#this could be modified to estimate models that do not match the operating model.
em_input = list()
# sim the 90 models!
for(m in 1:n.mods){

    #initial numbers at age
    NAA_re = list()
    input = prepare_wham_input(basic_info = groundfish_info, selectivity = selectivity, M = M, NAA_re = NAA_re)
    
    NAA_re$recruit_model = df.mods$Recruitment[m] #random effects with a constant mean
    #whether RE on recruitment or on all NAA
    if(is.na(df.mods$NAA_sig[m])) {
      NAA_re$sigma = "rec" #random about mean
    } else {
      NAA_re$sigma = "rec+1"
    }
    
    #what type of correlation structure
    if(!is.na(df.mods$ar1_a[m]) & !is.na(df.mods$ar1_y[m])) {
        NAA_re$cor = "2dar1"
    }
    if(is.na(df.mods$ar1_a[m]) & !is.na(df.mods$ar1_y[m])) {
        NAA_re$cor = "ar1_y"
    }
    if(!is.na(df.mods$ar1_a[m]) & is.na(df.mods$ar1_y[m])) {
        NAA_re$cor = "ar1_a"
    }
    if(is.na(df.mods$ar1_a[m]) & is.na(df.mods$ar1_y[m])) {
        NAA_re$cor = "iid"
    }
    input = set_NAA(input, NAA_re)
    
    #put in data simulated from operating model
    em_input[[m]] = lapply(1:nsim, function(x) {
        input_i = input
        input_i$data = sim_input[[m]][[x]]$data #put in simulated operating model data
        return(input_i)
    })
}
#save 
saveRDS(em_input, file.path(write.dir, "em_input.RDS"))
#em_input <- readRDS(file.path(write.dir, "em_input.RDS"))

#NOT RUN YET
#em_input <- readRDS(file.path(write.dir, "em_input.RDS"))
em_fits = list()
#for(m in 1:n.mods){ #just do the first scenario until we are ready to do all of them.
for(m in 1){
    em_fits[[m]] = lapply(1:nsim, function(x){
        cat(paste("model:",m, "fit:", x, "start \n"))
        out = fit_wham(em_input[[m]][[x]], do.osa = FALSE, MakeADFun.silent = TRUE, retro.silent = TRUE, save.sdrep = FALSE)
        cat(paste("model:",m, "fit:", x, "done \n"))
        return(out)
        })
}

sapply(sim_fits[[1]], function(x) x$parList$log_NAA_sigma[1])


