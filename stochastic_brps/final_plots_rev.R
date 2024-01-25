source("results_plotting_functions.R")
###########################################################
#rec+1, 2dar1


# mod_nbc_proj <- readRDS("msel3ahl_nbc_proj.RDS")
# mod_nbcpe_proj <- readRDS("msel3ahl_nbcpe_proj.RDS")
# mod_nbcoe_proj <- readRDS("msel3ahl_nbcoe_proj.RDS")
#mod_bc_proj <- readRDS("msel3ahl_bc_proj.RDS")
# mod_nbc_proj_sdrep_bc <- readRDS("msel3ahl_nbc_proj_sdrep_bc.RDS")
#mod_bc_proj_sdrep_bc <- readRDS("msel3ahl_bc_proj_sdrep_bc.RDS")
# mod_nbcpe_proj_sdrep_bc <- readRDS("msel3ahl_nbcpe_proj_sdrep_bc.RDS")
# mod_nbcoe_proj_sdrep_bc <- readRDS("msel3ahl_nbcoe_proj_sdrep_bc.RDS")

mod_bc_proj <- readRDS("mod_bc_proj.RDS")
mod_bc_proj_sdrep_bc <- readRDS("mod_bc_proj_sdrep_bc.RDS")
mod_bc_F0_proj <- readRDS("mod_bc_F0_proj.RDS")
mod_bc_F0_proj_sdrep_bc <- readRDS("mod_bc_F0_proj_sdrep_bc.RDS")

# mods <- c("mod_bc_proj","mod_nbcpe_proj","mod_nbcoe_proj","mod_nbc_proj")
# sdrep_bcs <- c("mod_bc_proj_sdrep_bc","mod_nbcpe_proj_sdrep_bc","mod_nbcoe_proj_sdrep_bc","mod_nbc_proj_sdrep_bc")
# titles <- c("bias-correct all", "bias-correct observations", "bias-correct process errors", "no bias-correction")

mods <- c("mod_bc_proj", 'mod_bc_F0_proj')
sdrep_bcs <- c("mod_bc_proj_sdrep_bc","mod_bc_F0_proj_sdrep_bc")
titles <- c("F = F40", "F = 0")

SSBAA_eq <- list()
for(i in 1:length(mods)){
    rep <- get(mods[i])$rep
    dat <- get(mods[i])$input$data
    par <- get(mods[i])$parList
    temp <- c(rep,dat,par)
    SSBAA_eq[[i]] <- get_SSBAA_eq(temp)
}

SSBAA <- list()
for(i in 1:length(mods)){
    rep <- get(mods[i])$rep
    dat <- get(mods[i])$input$data
    temp <- c(rep,dat)
    SSBAA[[i]] <- get_SSB(temp, total = FALSE)
}
n_ages <- 9
proj_sims <- readRDS("mod_bc_proj_sims.RDS")
proj_median_pred_NAA <- lapply(proj_sims, function(x) sapply(1:n_ages, function(z) apply(sapply(x, function(y) y$pred_NAA[,z]),1,median)))
proj_median_NAA <- lapply(proj_sims, function(x) sapply(1:n_ages, function(z) apply(sapply(x, function(y) y$NAA[,z]),1,median)))
proj_mean_NAA <- lapply(proj_sims, function(x) sapply(1:n_ages, function(z) apply(sapply(x, function(y) y$NAA[,z]),1,mean)))
names(proj_median_pred_NAA) <- mods

proj_median_SSB <- lapply(proj_sims, function(x) apply(sapply(x, function(y) y$SSB),1,median))
names(proj_median_SSB) <- mods
proj_mean_SSB <- lapply(proj_sims, function(x) apply(sapply(x, function(y) y$SSB),1,mean))
names(proj_mean_SSB) <- mods

median_SSBAA <- list()
for(i in 1:length(mods)){
    rep <- get(mods[i])$rep
    dat <- get(mods[i])$input$data
    temp <- c(rep,dat)
    median_SSBAA[[i]] <- get_SSB(temp, median_NAA = proj_median_NAA[[i]], total = FALSE)
}

median_cumsum_SSBAA <- lapply(1:length(mods), function(x) sapply(1:n_ages, function(a) {
    apply(sapply(proj_sims[[x]], function(y) {
        res <- get(mods[x])
        zaa <- res$rep$FAA_tot + res$rep$MAA
        S <- exp(-zaa * res$input$data$fracyr_SSB)
        if(is.null(res$input$data$n_ages_pop)) res$input$data$n_ages_pop <- res$input$data$n_ages
        ind <- c(1:res$input$data$n_ages, rep(res$input$data$n_ages, res$input$data$n_ages_pop-res$input$data$n_ages))
        ssbaa <- y$NAA * S *res$input$data$waa[res$input$data$waa_pointer_ssb,,ind] * res$input$data$mature[,ind]
        out <- apply(cbind(ssbaa[,1:a]),1,sum)
        return(out)
    }),1,median)
}))

pred_SSB <- list()
for(i in 1:length(mods)){
    rep <- get(mods[i])$rep
    dat <- get(mods[i])$input$data
    temp <- c(rep,dat)
    pred_SSB[[i]] <- get_SSB(temp, pred = TRUE, total = TRUE)
}

png("compare_recruitment.png", width=8, height=8, units="in", res=300)
yrs <- 40:141
ymax <- max(sapply(proj_median_NAA, function(x) return(x[yrs,1])))
ymax <- max(sapply(proj_mean_NAA, function(x) return(x[yrs,1])))
ymax <- max(c(ymax, max(sapply(mods, function(x) get(x)$rep$pred_NAA[yrs,1]))))
ymax <- max(c(ymax, max(sapply(mods, function(x) get(x)$rep$NAA[yrs,1]))))
par(mfrow = c(1,1), mar = c(2,2,3,1), oma = c(3,3,0,0))
for(i in 1:1){
    plot(1980 + yrs, proj_median_NAA[[1]][yrs,1], type = "n", ylim = c(0,ymax), ylab = "", xlab = "")
    lines(1980 + yrs, proj_median_NAA[[i]][yrs,1], col = 1, lwd = 2)
    lines(1980 + yrs, proj_mean_NAA[[i]][yrs,1], col = 2, lwd = 2)
    lines(1980 + yrs, get(mods[i])$rep$pred_NAA[yrs,1], col = 2, lty = 2, lwd = 2)
    lines(1980 + yrs, get(mods[i])$rep$NAA[yrs,1], col = 1, lty = 2, lwd = 2)
    lines(1980 + yrs, exp(as.list(get(sdrep_bcs[i]), "Est. (bias.correct)", report=T)$log_NAA)[yrs,1], col = 3, lty = 3, lwd = 2)    
    lines(1980 + yrs, as.list(get(sdrep_bcs[i]), "Est. (bias.correct)", report=T)$NAA[yrs,1], col = 4, lty = 3, lwd = 2)    
    # title(titles[i], line = 1)
    # lines(1980 + yrs, proj_median_R[[mods[i+4]]][yrs], col = 1, lty = 2)
    # lines(1980 + yrs, jitter(get(mods[i+4])$rep$pred_NAA[yrs,1], factor = 0.5), col = 2, lty = 2)
    # lines(1980 + yrs, jitter(get(mods[i+4])$rep$NAA[yrs,1], factor = 0.5), col = 3, lty = 2)
    if(i == 1) legend("bottomright", legend = paste(c("simulation median", "simulation mean", "hat E(R)", "hat R", "exp(hat log R (TMB b-c))","hat R (TMB b-c)")),
        lty = c(rep(1:2, each = 2),3,3), col = c(1,2,2,1,3,4), lwd = 2)
}
mtext(side = 2, "Recruitment (1000s)", line = 1, outer = TRUE)
mtext(side = 1, "Year", line = 1, outer = TRUE)
dev.off()

png("compare_ssb.png", width=12, height=8, units="in", res=300)
yrs <- 40:141
ymax <- max(sapply(proj_median_SSB, function(x) return(x[yrs])))
ymax <- max(sapply(proj_mean_SSB, function(x) return(x[yrs])))
ymax <- max(c(ymax, max(sapply(mods, function(x) get(x)$rep$SSB[yrs]))))
ymax <- max(c(ymax, max(sapply(1:length(mods), function(x) pred_SSB[[x]][yrs]))))
ymax <- 200000
par(mfrow = c(1,2), mar = c(2,2,3,1), oma = c(3,3,0,0))
for(i in 1:length(mods)){
    plot(1980 + yrs, proj_median_SSB[[1]][yrs], type = "n", ylim = c(0,ymax), ylab = "", xlab = "")
    lines(1980 + yrs, proj_median_SSB[[mods[i]]][yrs], col = 1, lwd = 2)
    lines(1980 + yrs, proj_mean_SSB[[mods[i]]][yrs], col = 2, lwd = 2)
    lines(1980 + yrs, pred_SSB[[i]][yrs], col = 2, lty = 2, lwd = 2)
    lines(1980 + yrs, get(mods[i])$rep$SSB[yrs], col = 1, lty = 2, lwd = 2)
    # if(i==1) lines(1980 + yrs, get(mods[i])$rep$NAA[yrs,1] * exp(get(mods[i])$rep$log_SPR_FXSPR[yrs]), col = 1, lty = 4, lwd = 2)
    # if(i==2) lines(1980 + yrs, get(mods[i])$rep$NAA[yrs,1] * exp(get(mods[i])$rep$log_SPR0[yrs]), col = 1, lty = 4, lwd = 2)
    lines(1980 + yrs, exp(as.list(get(sdrep_bcs[i]), "Est. (bias.correct)", report=T)$log_SSB)[yrs], col = 3, lty = 3, lwd = 2)    
    lines(1980 + yrs, as.list(get(sdrep_bcs[i]), "Est. (bias.correct)", report=T)$SSB[yrs], col = 4, lty = 3, lwd = 2)    
    #lines(1980 + yrs, pred_SSB_eq[[i]][yrs], col = 2, lty = 3)
    title(titles[i], line = 1)
    if(i == 2) legend("topright", legend = paste(c("simulation median", "simulation mean", "f(hat E(NAA))", "f(hat NAA)", "exp(hat log SSB (TMB b-c))","hat SSB (TMB b-c)")),
        lty = c(rep(1:2, each = 2),3,3), col = c(1,2,2,1,3,4), lwd = 2)
}
#mtext(side = 2, "SSB(F=F40) (mt)", line = 1, outer = TRUE)
mtext(side = 2, "SSB (mt)", line = 1, outer = TRUE)
mtext(side = 1, "Year", line = 1, outer = TRUE)
dev.off()

png(paste0("compare_SSBAA_F40_at_age.png"), width=12, height=8, units="in", res=300)
plot.cols <- rep(c("grey","red", "blue"),2)
plot.cols <- sapply(1:length(plot.cols), function(x) adjustcolor(plot.cols[x], alpha = rep(c(0.3,0.8),each = 3)[x]))
par(mfrow = c(3,3), mar = c(2,2,3,1), oma = c(3,3,0,0))
for(a in 1:9){
  yrs <- 42:141
  ymax <- max(sapply(median_SSBAA[1:2], function(x) return(x[yrs,a])))
  #if(a != 20) ymax <- max(c(ymax, sapply(proj_mean_NAA[2], function(x) return(x[yrs,a]))))
  ymax <- max(ymax, sapply(SSBAA[1:2], function(x) return(x[yrs,a])))
  ymax <- max(ymax, sapply(SSBAA_eq[1:2], function(x) return(x[yrs,a])))
  #ymax <- 300000
  plot(1980 + yrs, median_SSBAA[[1]][yrs,a], type = "n", ylim = c(0,ymax), ylab = "", xlab = "")
  for(i in 1:2){
    lines(1980 + yrs, median_SSBAA[[i]][yrs,a], col = i, lwd = 2)
    lines(1980 + yrs, SSBAA[[i]][yrs,a], col = i, lty = 2, lwd = 2)
    lines(1980 + yrs, SSBAA_eq[[i]][yrs,a], col = i, lty = 3, lwd = 2)
  }
  #title(titles[i-4], line = 1)
  legend("bottomright", legend = paste0("Age ", a))
  if(a == 1) legend("bottom", legend = paste(c("simulation median", "estimated", "Rhat x SSB/R"), rep(c("F=F40", "F = 0"), each = 3)),
      lty = c(1,2,3), col = rep(1:2, each = 3), lwd = 2)
}
mtext(side = 2, paste0("SSB at age"), line = 1, outer = TRUE)
mtext(side = 1, "Year", line = 1, outer = TRUE)
dev.off()


png(paste0("compare_cumulative_SSBAA_F40_at_age.png"), width=12, height=8, units="in", res=300)
par(mfrow = c(3,3), mar = c(2,2,3,1), oma = c(3,3,0,0))
for(a in 1:9){
  print(a)
  yrs <- 42:141
  ymax <- max(sapply(median_SSBAA[1:2], function(x) return(apply(x[yrs,1:a, drop = F],1,sum))))
  #if(a != 20) ymax <- max(c(ymax, sapply(proj_mean_NAA[2], function(x) return(x[yrs,a]))))
  ymax <- max(ymax, sapply(SSBAA[1:2], function(x) return(apply(x[yrs,1:a, drop = F],1,sum))))
  ymax <- max(ymax, sapply(SSBAA_eq[1:2], function(x) return(apply(x[yrs,1:a, drop = F],1,sum))))
  ymax <- max(ymax, sapply(median_cumsum_SSBAA[1:2], function(x) return(x[yrs,a])))
  #ymax <- 300000
  plot(1980 + yrs, apply(median_SSBAA[[1]][yrs,1:a, drop = F],1,sum), type = "n", ylim = c(0,ymax), ylab = "", xlab = "")
  for(i in 1:2){
    lines(1980 + yrs, apply(median_SSBAA[[i]][yrs,1:a, drop = F],1,sum), col = i, lwd = 2)
    lines(1980 + yrs, apply(SSBAA[[i]][yrs,1:a, drop = F],1,sum), col = i, lty = 2, lwd = 2)
    lines(1980 + yrs, apply(SSBAA_eq[[i]][yrs,1:a, drop = F],1,sum), col = i, lty = 3, lwd = 2)
    lines(1980 + yrs, median_cumsum_SSBAA[[i]][yrs,a, drop = F], col = i, lty = 4, lwd = 2)
  }
  #title(titles[i-4], line = 1)
  legend("bottomright", legend = paste0("Age ", a))
  if(a == 1) legend("bottom", legend = paste(c("sum(median)", "sum(estimated)", "sum(Rhat x SSB/R)", "median(sum)"), rep(c("F=F40", "F=0"), each = 4)),
      lty = c(1,2,3,4), col = rep(1:2, each = 4), lwd = 2)
}
mtext(side = 2, paste0("Cumulative SSB at age"), line = 1, outer = TRUE)
mtext(side = 1, "Year", line = 1, outer = TRUE)
dev.off()


#############################################
#rec, ar1

mod_nbcpe_proj <- readRDS("rec_nbcpe_proj.RDS")
mod_bc_proj <- readRDS("rec_bc_proj.RDS")
mod_nbcpe_F0_proj <- readRDS("rec_nbcpe_F0_proj.RDS")
mod_bc_F0_proj <- readRDS("rec_bc_F0_proj.RDS")

types <- c("bc","nbcpe")
mods <- paste0("mod_", types, "_proj")
mods <- c(mods, paste0("mod_", types, "_F0_proj"))
sdrep_bcs <- c("mod_bc_proj_sdrep_bc","mod_nbcpe_proj_sdrep_bc","mod_bc_F0_proj_sdrep_bc","mod_nbcpe_F0_proj_sdrep_bc")

proj_sims <- readRDS("rec_proj_sims.RDS")
proj_median_pred_NAA <- lapply(proj_sims, function(x) sapply(1:9, function(z) apply(sapply(x, function(y) y$pred_NAA[,z]),1,median)))
names(proj_median_pred_NAA) <- mods
proj_median_NAA <- lapply(proj_sims, function(x) sapply(1:9, function(z) apply(sapply(x, function(y) y$NAA[,z]),1,median)))
names(proj_median_NAA) <- mods
#mods <- c("mod_bc_proj","mod_nbcpe_proj","mod_nbcoe_proj","mod_nbc_proj")

median_SSBAA <- list()
for(i in 1:length(mods)){
    rep <- get(mods[i])$rep
    dat <- get(mods[i])$input$data
    temp <- c(rep,dat)
    median_SSBAA[[i]] <- get_SSB(temp, median_NAA = proj_median_NAA[[i]], total = FALSE)
}

SSBAA_eq <- list()
for(i in 1:length(mods)){
    rep <- get(mods[i])$rep
    dat <- get(mods[i])$input$data
    par <- get(mods[i])$parList
    temp <- c(rep,dat,par)
    SSBAA_eq[[i]] <- get_SSBAA_eq(temp)
}

SSBAA <- list()
for(i in 1:length(mods)){
    rep <- get(mods[i])$rep
    dat <- get(mods[i])$input$data
    temp <- c(rep,dat)
    SSBAA[[i]] <- get_SSB(temp, total = FALSE)
}

png(paste0("compare_SSBAA_F40_at_age_rec.png"), width=12, height=8, units="in", res=300)
plot.cols <- rep(c("grey","red", "blue"),2)
plot.cols <- sapply(1:length(plot.cols), function(x) adjustcolor(plot.cols[x], alpha = rep(c(0.3,0.8),each = 3)[x]))
par(mfrow = c(3,3), mar = c(2,2,3,1), oma = c(3,3,0,0))
for(a in 1:9){
  yrs <- 42:141
  ymax <- max(sapply(median_SSBAA[1:2], function(x) return(x[yrs,a])))
  #if(a != 20) ymax <- max(c(ymax, sapply(proj_mean_NAA[2], function(x) return(x[yrs,a]))))
  ymax <- max(ymax, sapply(SSBAA[1:2], function(x) return(x[yrs,a])))
  ymax <- max(ymax, sapply(SSBAA_eq[1:2], function(x) return(x[yrs,a])))
  #ymax <- 300000
  plot(1980 + yrs, median_SSBAA[[1]][yrs,a], type = "n", ylim = c(0,ymax), ylab = "", xlab = "")
  for(i in 1:2){
    lines(1980 + yrs, median_SSBAA[[i]][yrs,a], col = i, lwd = 2)
    lines(1980 + yrs, SSBAA[[i]][yrs,a], col = i, lty = 2, lwd = 2)
    lines(1980 + yrs, SSBAA_eq[[i]][yrs,a], col = i, lty = 3, lwd = 2)
  }
  #title(titles[i-4], line = 1)
  legend("bottomright", legend = paste0("Age ", a))
  if(a == 1) legend("bottom", legend = paste(c("simulation median", "estimated", "Rhat x SSB/R"), rep(c("bias-correct", "no bias-correct"), each = 3)),
      lty = c(1,2,3), col = rep(1:2, each = 3), lwd = 2)
}
mtext(side = 2, paste0("SSB at age (F=F40)"), line = 1, outer = TRUE)
mtext(side = 1, "Year", line = 1, outer = TRUE)
dev.off()

na <- 9
median_cumsum_SSBAA <- lapply(1:4, function(x) sapply(1:na, function(a) {
  print(x)
    apply(sapply(proj_sims[[x]], function(y) {
        res <- get(mods[x])
        zaa <- res$rep$FAA_tot + res$rep$MAA
        S <- exp(-zaa * res$input$data$fracyr_SSB)
        if(is.null(res$input$data$n_ages_pop)) res$input$data$n_ages_pop <- res$input$data$n_ages
        ind <- c(1:res$input$data$n_ages, rep(res$input$data$n_ages, res$input$data$n_ages_pop-res$input$data$n_ages))
        ssbaa <- y$NAA * S *res$input$data$waa[res$input$data$waa_pointer_ssb,,ind] * res$input$data$mature[,ind]
        out <- apply(cbind(ssbaa[,1:a]),1,sum)
        return(out)
    }),1,median)
}))

png(paste0("compare_cumulative_SSBAA_F40_at_age_rec.png"), width=12, height=8, units="in", res=300)
par(mfrow = c(3,3), mar = c(2,2,3,1), oma = c(3,3,0,0))
for(a in 1:9){
  print(a)
  yrs <- 42:141
  ymax <- max(sapply(median_SSBAA[1:2], function(x) return(apply(x[yrs,1:a, drop = F],1,sum))))
  #if(a != 20) ymax <- max(c(ymax, sapply(proj_mean_NAA[2], function(x) return(x[yrs,a]))))
  ymax <- max(ymax, sapply(SSBAA[1:2], function(x) return(apply(x[yrs,1:a, drop = F],1,sum))))
  ymax <- max(ymax, sapply(SSBAA_eq[1:2], function(x) return(apply(x[yrs,1:a, drop = F],1,sum))))
  ymax <- max(ymax, sapply(median_cumsum_SSBAA[1:2], function(x) return(x[yrs,a])))
  #ymax <- 300000
  plot(1980 + yrs, apply(median_SSBAA[[1]][yrs,1:a, drop = F],1,sum), type = "n", ylim = c(0,ymax), ylab = "", xlab = "")
  for(i in 1:2){
    lines(1980 + yrs, apply(median_SSBAA[[i]][yrs,1:a, drop = F],1,sum), col = i, lwd = 2)
    lines(1980 + yrs, apply(SSBAA[[i]][yrs,1:a, drop = F],1,sum), col = i, lty = 2, lwd = 2)
    lines(1980 + yrs, apply(SSBAA_eq[[i]][yrs,1:a, drop = F],1,sum), col = i, lty = 3, lwd = 2)
    lines(1980 + yrs, median_cumsum_SSBAA[[i]][yrs,a, drop = F], col = i, lty = 4, lwd = 2)
  }
  #title(titles[i-4], line = 1)
  legend("bottomright", legend = paste0("Age ", a))
  if(a == 1) legend("bottom", legend = paste(c("sum(median)", "sum(estimated)", "sum(Rhat x SSB/R)", "median(sum)"), rep(c("bias-correct", "no bias-correct"), each = 4)),
      lty = c(1,2,3,4), col = rep(1:2, each = 4), lwd = 2)
}
mtext(side = 2, paste0("cumulative sum of SSB at age (F=F40)"), line = 1, outer = TRUE)
mtext(side = 1, "Year", line = 1, outer = TRUE)
dev.off()


#############################################
#rec, iid

mod_nbcpe_proj <- readRDS("rec_iid_nbcpe_proj.RDS")
mod_bc_proj <- readRDS("rec_iid_bc_proj.RDS")
mod_nbcpe_F0_proj <- readRDS("rec_iid_nbcpe_F0_proj.RDS")
mod_bc_F0_proj <- readRDS("rec_iid_bc_F0_proj.RDS")

types <- c("bc","nbcpe")
mods <- paste0("mod_", types, "_proj")
mods <- c(mods, paste0("mod_", types, "_F0_proj"))
sdrep_bcs <- c("mod_bc_proj_sdrep_bc","mod_nbcpe_proj_sdrep_bc","mod_bc_F0_proj_sdrep_bc","mod_nbcpe_F0_proj_sdrep_bc")

proj_sims <- readRDS("rec_iid_proj_sims.RDS")
proj_median_pred_NAA <- lapply(proj_sims, function(x) sapply(1:9, function(z) apply(sapply(x, function(y) y$pred_NAA[,z]),1,median)))
names(proj_median_pred_NAA) <- mods
proj_median_NAA <- lapply(proj_sims, function(x) sapply(1:9, function(z) apply(sapply(x, function(y) y$NAA[,z]),1,median)))
names(proj_median_NAA) <- mods
#mods <- c("mod_bc_proj","mod_nbcpe_proj","mod_nbcoe_proj","mod_nbc_proj")

median_SSBAA <- list()
for(i in 1:length(mods)){
    rep <- get(mods[i])$rep
    dat <- get(mods[i])$input$data
    temp <- c(rep,dat)
    median_SSBAA[[i]] <- get_SSB(temp, median_NAA = proj_median_NAA[[i]], total = FALSE)
}

SSBAA_eq <- list()
for(i in 1:length(mods)){
    rep <- get(mods[i])$rep
    dat <- get(mods[i])$input$data
    par <- get(mods[i])$parList
    temp <- c(rep,dat,par)
    SSBAA_eq[[i]] <- get_SSBAA_eq(temp)
}

SSBAA <- list()
for(i in 1:length(mods)){
    rep <- get(mods[i])$rep
    dat <- get(mods[i])$input$data
    temp <- c(rep,dat)
    SSBAA[[i]] <- get_SSB(temp, total = FALSE)
}

SSBAA_sdrep_bc <- list()
for(i in 1:length(mods)){
    rep <- get(mods[i])$rep
    dat <- get(mods[i])$input$data
    temp <- c(rep,dat)
    NAA <- exp(as.list(get(sdrep_bcs[i]), "Est. (bias.correct)", report=T)$log_NAA)
    SSBAA_sdrep_bc[[i]] <- get_SSB(temp, median_NAA = NAA, total = FALSE)
}

png(paste0("compare_SSBAA_F40_at_age_rec_iid.png"), width=12, height=8, units="in", res=300)
#plot.cols <- viridisLite::viridis(n = 3, begin = 0, end = 1, alpha = 0.4, option = "turbo")
plot.cols <- rep(c("grey","red", "blue", "darkgreen"),2)
plot.cols <- sapply(1:length(plot.cols), function(x) adjustcolor(plot.cols[x], alpha = rep(c(0.3,0.8),each = 4)[x]))
par(mfrow = c(3,3), mar = c(2,2,3,1), oma = c(3,3,0,0))
for(a in 1:9){
  yrs <- 42:141
  ymax <- max(sapply(median_SSBAA[1:2], function(x) return(x[yrs,a])))
  #if(a != 20) ymax <- max(c(ymax, sapply(proj_mean_NAA[2], function(x) return(x[yrs,a]))))
  ymax <- max(ymax, sapply(SSBAA[1:2], function(x) return(x[yrs,a])))
  ymax <- max(ymax, sapply(SSBAA_eq[1:2], function(x) return(x[yrs,a])))
  ymax <- max(ymax, sapply(SSBAA_sdrep_bc[1:2], function(x) return(x[yrs,a])))
  plot(1980 + yrs, median_SSBAA[[1]][yrs,a], type = "n", ylim = c(0,ymax), ylab = "", xlab = "")
  k <- 0
  for(i in 1:2){
    lines(1980 + yrs, median_SSBAA[[i]][yrs,a], col = plot.cols[k+1], lty = i, lwd = 2)
    lines(1980 + yrs, SSBAA[[i]][yrs,a], col = plot.cols[k+2], lty = i, lwd = 2)
    lines(1980 + yrs, SSBAA_eq[[i]][yrs,a], col = plot.cols[k+3], lty = i, lwd = 2)
    lines(1980 + yrs, SSBAA_sdrep_bc[[i]][yrs,a], col = plot.cols[k+4], lty = i, lwd = 2)    
    k <- k + 4
  }
  legend("bottomright", legend = paste0("Age ", a))
  if(a == 1) legend("bottom", legend = paste(c("simulation median", "estimated", "Rhat x SSB/R", "TMB bias.correct"), rep(c("bias-correct", "no bias-correct"), each = 4)),
      lty = rep(1:2, each = 4), col = plot.cols, lwd = 2)
}
mtext(side = 2, paste0("SSB at age (F=F40)"), line = 1, outer = TRUE)
mtext(side = 1, "Year", line = 1, outer = TRUE)
dev.off()

na <- 9
median_cumsum_SSBAA <- lapply(1:4, function(x) sapply(1:na, function(a) {
  print(x)
    apply(sapply(proj_sims[[x]], function(y) {
        res <- get(mods[x])
        zaa <- res$rep$FAA_tot + res$rep$MAA
        S <- exp(-zaa * res$input$data$fracyr_SSB)
        if(is.null(res$input$data$n_ages_pop)) res$input$data$n_ages_pop <- res$input$data$n_ages
        ind <- c(1:res$input$data$n_ages, rep(res$input$data$n_ages, res$input$data$n_ages_pop-res$input$data$n_ages))
        ssbaa <- y$NAA * S *res$input$data$waa[res$input$data$waa_pointer_ssb,,ind] * res$input$data$mature[,ind]
        out <- apply(cbind(ssbaa[,1:a]),1,sum)
        return(out)
    }),1,median)
}))

png(paste0("compare_cumulative_SSBAA_F40_at_age_rec_iid.png"), width=12, height=8, units="in", res=300)
plot.cols <- rep(c("grey","red", "blue", "darkgreen"),2)
plot.cols <- sapply(1:length(plot.cols), function(x) adjustcolor(plot.cols[x], alpha = rep(c(0.3,0.8),each = 4)[x]))
par(mfrow = c(3,3), mar = c(2,2,3,1), oma = c(3,3,0,0))
for(a in 1:9){
  print(a)
  yrs <- 42:141
  ymax <- max(sapply(median_SSBAA[1:2], function(x) return(apply(x[yrs,1:a, drop = F],1,sum))))
  #if(a != 20) ymax <- max(c(ymax, sapply(proj_mean_NAA[2], function(x) return(x[yrs,a]))))
  ymax <- max(ymax, sapply(SSBAA[1:2], function(x) return(apply(x[yrs,1:a, drop = F],1,sum))))
  ymax <- max(ymax, sapply(SSBAA_eq[1:2], function(x) return(apply(x[yrs,1:a, drop = F],1,sum))))
  ymax <- max(ymax, sapply(median_cumsum_SSBAA[1:2], function(x) return(x[yrs,a])))
  #ymax <- 300000
  plot(1980 + yrs, apply(median_SSBAA[[1]][yrs,1:a, drop = F],1,sum), type = "n", ylim = c(0,ymax), ylab = "", xlab = "")
  k <- 0  
  for(i in 1:2){
    lines(1980 + yrs, apply(median_SSBAA[[i]][yrs,1:a, drop = F],1,sum), col = plot.cols[k+1], lty = i, lwd = 2)
    lines(1980 + yrs, apply(SSBAA[[i]][yrs,1:a, drop = F],1,sum), col = plot.cols[k+2], lty = i, lwd = 2)
    lines(1980 + yrs, apply(SSBAA_eq[[i]][yrs,1:a, drop = F],1,sum), col = plot.cols[k+3], lty = i, lwd = 2)
    lines(1980 + yrs, median_cumsum_SSBAA[[i]][yrs,a, drop = F], col = plot.cols[k+4], lty = i, lwd = 2)
    k <- k + 4
  }
  #title(titles[i-4], line = 1)
  legend("bottomright", legend = paste0("Age ", a))
  if(a == 1) legend("bottom", legend = paste(c("sum(median)", "sum(estimated)", "sum(Rhat x SSB/R)", "median(sum)"), rep(c("bias-correct", "no bias-correct"), each = 4)),
      lty = rep(1:2, each = 4), col = plot.cols, lwd = 2)
}
mtext(side = 2, paste0("Cumulative sum of SSB at age (F=F40)"), line = 1, outer = TRUE)
mtext(side = 1, "Year", line = 1, outer = TRUE)
dev.off()

