
# btime <- Sys.time()
# wham.dir <- find.package("wham")
# etime <- Sys.time()
# runtime = etime - btime

# devtools::load_all("~/work/wham/wham")
# devtools::install_github("timjmiller/wham", dependencies=TRUE, ref="devel")
library(wham)
library(ggplot2)
library(tidyr)
library(dplyr)
library(here)
# create directory for analysis, e.g.

write.dir <- file.path(here(),"results","flatfish_NAA")

if(!exists("write.dir")) write.dir = getwd()
if(!dir.exists(write.dir)) dir.create(write.dir, recursive = T)
setwd(write.dir)

#number of simulations for each scenario
nsim = 100

#NAA sigmas for each scenario
sigs = cbind(R_sig = c(0.3,0.5,0.7), NAA_sig = c(0.1,0.3,0.5))
sigs_use = expand.grid(sig_ind = 1:NROW(sigs), ar1_y = c(NA,0.25,0.5), ar1_a = c(NA,0.25,0.5))

df.mods = expand.grid(Recruitment = 2, R_sig = c(0.3,0.5,0.7), NAA_sig = c(NA), ar1_y = c(NA,0.25,0.5), ar1_a = NA)
#df.mods = rbind(df.mods, cbind(Recruitment = 2, sigs[sigs_use$sig_ind,], sigs_use[,2:3]))
df.mods = rbind(df.mods, expand.grid(Recruitment = 2, R_sig = c(0.3,0.5,0.7), NAA_sig = c(0.1,0.3,0.5), ar1_y = c(NA,0.25,0.5), ar1_a = c(NA,0.25,0.5)))

n.mods <- dim(df.mods)[1] #90 scenarios!
df.mods$Model <- paste0("m_",1:n.mods)
df.mods <- df.mods %>% select(Model, everything()) # moves Model to first col

# look at model table
df.mods

#taken from example 10 script in WHAM
#make a list of input components that prepare_wham_input can use to generate an input for fit_wham
make_digifish <- function(years = 1975:2014) {
    digifish = list()
    digifish$ages = 1:10
    digifish$years = years
    na = length(digifish$ages)
    ny = length(digifish$years)

    digifish$n_fleets = 1
    digifish$catch_cv = matrix(0.1, ny, digifish$n_fleets)
    digifish$catch_Neff = matrix(200, ny, digifish$n_fleets)
    digifish$n_indices = 1
    digifish$index_cv = matrix(0.3, ny, digifish$n_indices)
    digifish$index_Neff = matrix(100, ny, digifish$n_indices)
    digifish$fracyr_indices = matrix(0.5, ny, digifish$n_indices)
    digifish$index_units = rep(1, length(digifish$n_indices)) #biomass
    digifish$index_paa_units = rep(2, length(digifish$n_indices)) #abundance
    digifish$maturity = t(matrix(1/(1 + exp(-1*(1:na - na/2))), na, ny))

    L = 100*(1-exp(-0.3*(1:na - 0)))
    W = exp(-11)*L^3
    nwaa = digifish$n_indices + digifish$n_fleets + 2
    digifish$waa = array(NA, dim = c(nwaa, ny, na))
    for(i in 1:nwaa) digifish$waa[i,,] = t(matrix(W, na, ny))

    digifish$fracyr_SSB = rep(0.25,ny)
    digifish$q = rep(0.3, digifish$n_indices)
    mid = floor(ny/2)
    #up then down
    digifish$F = matrix(0.2 + c(seq(0,0.4,length.out = ny/2),seq(0.4,0,length.out=ny-mid)),ny, digifish$n_fleets)

    digifish$selblock_pointer_fleets = t(matrix(1:digifish$n_fleets, digifish$n_fleets, ny))
    digifish$selblock_pointer_indices = t(matrix(digifish$n_fleets + 1:digifish$n_indices, digifish$n_indices, ny))
    return(digifish)
}
digifish = make_digifish()

selectivity = list(model = c(rep("logistic", digifish$n_fleets),rep("logistic", digifish$n_indices)),
    initial_pars = rep(list(c(5,1)), digifish$n_fleets + digifish$n_indices)) #fleet, index

M = list(initial_means = rep(0.2, length(digifish$ages)))

#initial numbers at age
NAA_re = list(N1_pars = exp(10)*exp(-(0:(length(digifish$ages)-1))*M$initial_means[1]))

sim_input = list()
# sim the 90 models!
for(m in 1:n.mods){

    if(is.na(df.mods$NAA_sig[m])) {
        NAA_re$sigma = "rec" #random about mean
    } else NAA_re$sigma = "rec+1"
    

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
    
    NAA_re$use_steepness = 0
    NAA_re$recruit_model = df.mods$Recruitment[m] #random effects with a constant mean
    NAA_re$recruit_pars = exp(10) #mean recruitment, if recruit model changes, need to change this

    input = prepare_wham_input(basic_info = digifish, selectivity = selectivity, M = M, NAA_re = NAA_re)
    input$par$log_NAA_sigma[1] = log(df.mods$R_sig[m])
    #input$par$trans_NAA_rho[2] = log(df.mods$ar1_y[m]+1) - log(1-df.mods$ar1_y[m])
    if(!is.na(df.mods$NAA_sig[m])) {
        input$par$log_NAA_sigma[2] = log(df.mods$NAA_sig[m])
    #    input$par$trans_NAA_rho[1] = log(df.mods$ar1_a[m]+1) - log(1-df.mods$ar1_a[m])
    }

    if(!is.na(df.mods$ar1_a[m])) {
        input$par$trans_NAA_rho[1] = log(df.mods$ar1_a[m]+1) - log(1-df.mods$ar1_a[m])
    }
    if(!is.na(df.mods$ar1_y[m])) {
        input$par$trans_NAA_rho[2] = log(df.mods$ar1_y[m]+1) - log(1-df.mods$ar1_y[m])
    }

    om = fit_wham(input, do.fit = FALSE)


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

#save simulated data sets to google drive?

#for(m in 1:n.mods){
m = 1
    sim_fits = lapply(1:nsim, function(x){
        input = sim_input[[1]][[x]]
        cat(paste("sim:", x))
        fit_wham(sim_input[[1]][[1]], do.osa = FALSE, MakeADFun.silent = TRUE, retro.silent = TRUE, save.sdrep = FALSE)
        cat(paste("model:",m, "fit:", x, "done \n"))
        })
#}

temp = fit_wham(sim_input[[1]][[1]], do.osa = FALSE, MakeADFun.silent = TRUE, retro.silent = TRUE, save.sdrep = FALSE)
plot_wham_output(temp, out.type = "html")

