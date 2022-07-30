make_OSA_residuals = function(mod,osa.opts = list(method="cdf", parallel=TRUE)){
  # one-step-ahead residuals
  if(mod$is_sdrep){ # only do OSA residuals if sdrep ran
    cat("Doing OSA residuals...\n");
    input = mod$input
    full_set = 1:length(input$data$obsvec)
    input$data$obs$residual = NA
    if(!is.null(input$data$condition_no_osa)) cat("OSA not available for some age comp likelihoods...\n")
    #first do continuous obs, condition on obs without osa (probably none)
    subset. = setdiff(full_set, c(input$data$subset_discrete_osa, input$data$conditional_no_osa))
    print("subset.")
    print(subset.)
    print(osa.opts$method)
    print(osa.opts$parallel)
    mod$env$data$do_osa = 1
    mod$env$retape()
    x = TMB::oneStepPredict(mod, method = "cdf", subset = subset.,
    observation.name = "obsvec", data.term.indicator="keep", discrete=F, parallel= T)
    print(head(x))
    stop()
#
    OSA.continuous <- suppressWarnings(TMB::oneStepPredict(obj=mod, observation.name="obsvec",
                                data.term.indicator="keep",
                                method=osa.opts$method,
                                discrete=FALSE, parallel=osa.opts$parallel,
                                subset = subset., conditional = input$data$conditional_no_osa))
    input$data$obs$residual[subset.] <- OSA.continuous$residual;
    print(OSA.continuous$residual)
    stop()
    mod$OSA.continuous = OSA.continuous
    if(!is.null(input$data$subset_discrete_osa)) {
      cat("Doing OSA for discrete age comp likelihoods...\n")
      conditional = union(input$data$condition_no_osa, subset.) #all with continuous and without osa 
      subset. = input$data$subset_discrete_osa
      #first do continuous
      OSA.discrete <- suppressWarnings(TMB::oneStepPredict(obj=mod, observation.name="obsvec",
                                  data.term.indicator="keep",
                                  method= osa.opts$method,
                                  discrete=TRUE, parallel=osa.opts$parallel,
                                  conditional = conditional))
      input$data$obs$residual[subset.] <- OSA.discrete$residual
      mod$OSA.discrete = OSA.discrete
    }
    mod$osa <- input$data$obs
    mod$env$data$do_osa = 0 #set this back to not using OSA likelihoods
    mod$env$retape()
    return(mod)
  }
}
full_set = 1:length(tfit$input$data$obsvec)
if(!is.null(tfit$input$data$condition_no_osa)) cat("OSA not available for some age comp likelihoods...\n")
#first do continuous obs, condition on obs without osa (probably none)
subset. = setdiff(full_set, c(tfit$input$data$subset_discrete_osa, tfit$input$data$conditional_no_osa))
sum(subset. != first)
length(tfit$input$data$obsvec)
tfit <- fit_wham(sim_input[[1]][[1]], do.osa=F, do.retro=F)
tfit_alt <- fit_wham(sim_input[[1]][[1]], do.osa=F, do.retro=F)
y = make_OSA_residuals(tfit_alt)
x = TMB::oneStepPredict(tfit, method = "cdf", subset = first,
  observation.name = "obsvec", data.term.indicator="keep", discrete=F, parallel= T)
#x = make_OSA_residuals(em_fits[[1]][[1]])

xfirst = TMB::oneStepPredict(obj = em_fits[[1]][[1]], observation.name="obsvec",
                                data.term.indicator="keep",
                                method="cdf",
                                discrete=FALSE, parallel=TRUE,
                                subset = first, conditional = input$data$conditional_no_osa)
xfirst = TMB::oneStepPredict(obj = em_fits[[1]][[1]], observation.name="obsvec",
                                data.term.indicator="keep",
                                method="cdf",
                                discrete=FALSE,
                                subset = first)
temp = TMB::oneStepPredict(em_fits[[1]][[1]], method = "cdf", subset = first,
  observation.name = "obsvec", data.term.indicator="keep", discrete=F, parallel= T)
temp2 = TMB::oneStepPredict(em_fits[[1]][[1]], method = "cdf", conditional = first, 
  observation.name = "obsvec", data.term.indicator="keep", discrete=T, parallel= T)
temp = TMB::oneStepPredict(em_fits[[1]][[1]], method = "oneStepGaussianOffMode", subset = first,
  observation.name = "obsvec", data.term.indicator="keep", discrete=F, parallel= T)
