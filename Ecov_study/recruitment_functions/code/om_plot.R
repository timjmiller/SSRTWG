library(here)
library(wham) #make sure to use the right version of wham
source(file.path(here(), "Ecov_study", "mortality", "code", "sim_management.R"))
verify_version()

file_loc  <- 'Ecov_study/recruitment_functions/om_inputs_06_13_2023'
om_inputs <- readRDS(file.path(here(),file_loc,'om_inputs.rds'))
df.oms    <- readRDS(file.path(here(),file_loc,'df.oms.rds'))

#####################################################################
##--Recruitment parameter set in simulation--#########
#GLB: check this
a <- 5.955694e-01 
b <- 2.404283e-05

ncols <- 7
nrows <- 7

ylims <- c(0,1E5)
xlims <- c(0,5E5)

ssb0 <- seq(0,xlims[2],length.out=1000)
r0   <- (a*ssb0)/(1+b*ssb0)

#Subset stocks
lls <- which(df.oms$Ecov_how%in%c(1,2,4) & 
               df.oms$Fhist%in%c('H-MSY','MSY') & 
               df.oms$Ecov_obs_sig==min(df.oms$Ecov_obs_sig) &
               df.oms$Ecov_effect==max(df.oms$Ecov_effect) & 
               df.oms$obs_error=='L' & 
               df.oms$R_sig==max(df.oms$R_sig) &
               df.oms$NAA_cor==min(df.oms$NAA_cor) & 
               df.oms$Ecov_re_cor==min(df.oms$Ecov_re_cor))

#PLOT
set.seed(1234)

make_plots(lls=lls,om_inputs=om_inputs,
           ncols=ncols,nrows=nrows,ylims=ylims,xlims=xlims)


