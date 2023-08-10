# plot NAA (one multipanel)

# # testing
# library(wham)
# library(here)
# library(tidyverse)
# library(ggplotFL)
# library(ggsci)
# library(cowplot)
# library(simhelpers)
# # stock.id = "SNEMAYT"
# # stock.id = "butterfish"
# stock.id = "ICEherring"
# re="NAA"
# bc.type = 2
# ty = sim.types = 2
# n.mods=4
# n.sim=100
# inv.rho.trans <- function(x) return(2/(1 + exp(-2*x)) - 1)
# res_dir=file.path("/home/bstock/Documents/ms/wham-sim/results")
# simdata_dir=file.path("/home/bstock/Documents/ms/wham-sim/data/simdata")
# plots_dir=file.path("/home/bstock/Documents/ms/wham-sim/plots")
# source("/home/bstock/Documents/ms/wham-sim/code/get_results.R")
# results <- get_results(stock.id, re, bc.type)

# plot_NAA(stock.id=stock.id, bc.type=bc.type, sim.types=sim.types, n.mods=4, n.sim=100, multipanel=FALSE, plot.eps=FALSE)

# bc.type     bias corrected obs only (= 1) or obs + process (= 2)
# sim.types   simulated obs only (= 1) or obs + process (= 2)
# n.mods      4 if all NAA models converged
# multipanel  TRUE makes a 5-panel plot (B, F, relB, relF, recruit), FALSE makes individual plots
# plot.eps    FALSE does not plot the TMB epsilon results, TRUE does
plot_rel_err <- function(results, om.id="R+S", multipanel=TRUE, plot.eps=FALSE,
                      res_dir=file.path(getwd(),"results"), 
                      simdata_dir=file.path(getwd(),"data","simdata"),
                      plots_dir=file.path(getwd(),"plots")){ 
  id <- paste0(stock.id)#,"_",re)
  if(bc.type == 1){
    bc <- "bias_correct_oe"
    plots_dir <- file.path(plots_dir, bc, id)
  }
  if(bc.type == 2){
    bc <- "bias_correct_oepe"
    plots_dir <- file.path(plots_dir, bc, id)
  }
  dir.create(plots_dir, showWarnings = FALSE)
  types <- c("OE","OEPE") # simulation type, not bias correction type

  # 95% CI bounds for median (sims in a given year)
  n.mods <- length(table(results$om))
  n.sim <- length(table(results$sim))
  bnds <- qbinom(c(0.025,0.975), n.sim, 0.5)/n.sim 
  # https://www-users.york.ac.uk/~mb55/intro/cicent.htm
  # https://stats.stackexchange.com/questions/122001/confidence-intervals-for-median
  # http://www.jstor.com/stable/2957563
  
  # boxplot whiskers = 95% CI for median to match time-series plots
  # boxplot boxes = 80% CI for median
  # n.yrs <- length(table(results$year))
  #n.yrs <- 1
  custom_boxplot_stat <- function(x, n.yrs = 1, n.sim = 100) {
    bnds95 <- qbinom(c(0.025,0.975), n.sim*n.yrs, 0.5)/(n.sim*n.yrs) # 95% CI bounds for median (sims x years)
    bnds80 <- qbinom(c(0.1,0.9), n.sim*n.yrs, 0.5)/(n.sim*n.yrs) # 80% CI bounds for median (sims x years)
    r <- quantile(x, probs = c(bnds95[1], bnds80[1], 0.5, bnds80[2], bnds95[2]))
    names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
    r
  }
  
  # ---------------------------------------------------------------------------------
  # Multipanel plots
  if(multipanel){
    for(ty in sim.types){
    # collapse across years, group by om/em
    	df.plot <- filter(results, type==levels(results$type)[ty]) %>%
    	  group_by(om2, em.x, em2, sim) %>%
    	  # group_by(om2, em.x, em2, year) %>%
    	  summarize(SSB.rel=median(SSB.rel), F.rel=median(F.rel), relB.rel=median(relB.rel, na.rm=T),
    	            relF.rel=median(relF.rel,na.rm=T), R.rel=median(R.rel), .groups = 'drop_last') %>%
    	  select(om2, em2, em.x, SSB.rel, F.rel, relB.rel, relF.rel, R.rel) %>% as.data.frame
    	df.plot <- do.call(data.frame, lapply(df.plot, function(x) replace(x, is.infinite(x), NA)))
    	# df.plot[df.plot == 0] = NA
    	df.plot <- df.plot %>%
    	  pivot_longer(-c(om2,em2,em.x), names_to = "variable", values_to = "val") %>%
    	  group_by(om2, em.x, em2, variable) %>%
    	  summarize(rel_err_var = var(val, na.rm=T), rel_err_mean = mean(val, na.rm=T)-1, K=sum(!is.na(val)), .groups = 'keep') %>%
    	  mutate(rel_err_se = sqrt(rel_err_var/K)) %>%
    	  mutate(rel_err_lo = rel_err_mean - 1.96*rel_err_se, rel_err_hi = rel_err_mean + 1.96*rel_err_se)       
    	df.plot$variable <- factor(df.plot$variable, levels=c("SSB.rel", "F.rel", "relB.rel", "relF.rel", "R.rel"), 
    	                       labels=c("SSB", "F", expression(B/B[40]["%"]), expression(F/F[40]["%"]), "Recruitment"))

    	 p <- ggplot(df.plot, aes(x=em.x, y=rel_err_mean)) +
    	  xlab("Estimation model") +
    	  ylab("Relative error") +
    	  geom_hline(yintercept = 0, linetype=2, color='black') +
    	  geom_linerange(aes(ymin=rel_err_lo, ymax=rel_err_hi), size=.5) +
    	  geom_point(aes(fill=em2), shape=21, size=3) +
    	   # scale_fill_jco(name="", labels=lapply(levels(df.plot$em2), function(x) parse(text=x))) +
    	   scale_fill_jco() +
    	   coord_cartesian(ylim=c(-.5,.5)) +
    	  facet_grid(rows=vars(variable), cols=vars(om2), labeller = label_parsed) +
    	  theme_bw() +
    	  # theme(legend.position="bottom",
    	  theme(legend.position="none",
    	        strip.text.y = element_text(size = 10), strip.text.x = element_text(size = 10, margin = margin(3,1,1,1, "pt")),
    	        axis.title.x = element_text(size = 12), axis.text.x = element_text(size = 8),
    	        legend.text = element_text(margin = margin(r = 6, l=1,unit = "pt"), hjust = 0, size=10), 
    	        legend.box.margin = margin(0,0,0,0), legend.margin = margin(0,0,0,0))
    	title <- ggdraw() + draw_label("Operating model", hjust = 0.3, vjust=1) + theme(plot.margin = margin(0, 0, 0, 0))
      p1 <- plot_grid(title, p, ncol = 1, rel_heights = c(0.045, 1))

      if(n.mods == 5) png(file.path(plots_dir, paste0("0_",id,"_medianCI_",types[ty],".png")), height=8, width=7, units='in', res=300)
      if(n.mods == 4) png(file.path(plots_dir, paste0("0_",id,"_medianCI_",types[ty],".png")), height=8, width=6.5, units='in', res=300)
      if(n.mods == 3) png(file.path(plots_dir, paste0("0_",id,"_medianCI_",types[ty],".png")), height=8, width=5, units='in', res=300)
      if(n.mods == 2) png(file.path(plots_dir, paste0("0_",id,"_medianCI_",types[ty],".png")), height=8, width=4, units='in', res=300)
      if(n.mods == 1) png(file.path(plots_dir, paste0("0_",id,"_medianCI_",types[ty],".png")), height=8, width=3, units='in', res=300)
      print(p1)
      dev.off()
    }
  } else { # individual plots -----------------------------------------------------------------------------
    for(ty in sim.types){
      df.plot <- filter(results, type==levels(results$type)[ty])
      
      # 1. SSB
      p <- ggplot(df.plot, aes(x=year, y=SSB.rel-1)) +
        stat_flquantiles(probs=bnds, alpha=0.7, fill="grey", geom="ribbon") +
        stat_summary(fun = "median", geom = "line", color = "red") +
        coord_cartesian(ylim=c(-.5,.5)) +
        xlab("Year") +
        ylab("Relative error in SSB") +
        geom_hline(yintercept = 0, linetype=2, color='black') +
        facet_grid(rows=vars(em2), cols=vars(om2), labeller = label_parsed) +
        theme_bw() +
        theme(axis.text.x = element_text(size=8), plot.margin = unit(c(0.3,0.3,0.1,0.1), "in"))
      png(file.path(plots_dir,paste0("1_ssb_",types[ty],".png")), width=7, height=7, units='in',res=100)
      print(p)
      grid::grid.text(unit(0.98,"npc"),0.5, label = 'Estimation model', rot = 270) # right
      grid::grid.text(unit(0.5,"npc"),unit(.98,'npc'), label = 'Operating model', rot = 0)   # top)
      dev.off()
      
      # # v2 = mean and 95% CI(mean)
      # # https://meghapsimatrix.github.io/simhelpers/articles/MCSE.html#relative-criteria
      # p <- df.plot %>% 
      #   group_by(om2, em2, year) %>%
      #   summarize(Tbar = mean(SSB_fit, na.rm=T), theta=mean(SSB_sim), Tvar=var(SSB_fit, na.rm=T), K=sum(!is.na(SSB_fit)), .groups = 'keep') %>%
      #   mutate(rel_err_mean = Tbar/theta-1, rel_err_se = sqrt(Tvar/(K*theta^2))) %>%
      #   mutate(rel_err_lo = rel_err_mean - 1.96*rel_err_se, rel_err_hi = rel_err_mean + 1.96*rel_err_se) %>%
      #   mutate(rel_err_lo2 = rel_err_mean - 0.674*rel_err_se, rel_err_hi2 = rel_err_mean + 0.674*rel_err_se) %>%
      #   # group_modify(~ calc_relative(.x, estimates = SSB_fit, true_param = SSB_sim, perfm_criteria="relative bias")) %>%
      #   # mutate(rel_err = rel_bias-1) %>%
      #   # mutate(rel_err_lo = rel_err - 1.96*rel_bias_mcse, rel_err_hi = rel_err + 1.96*rel_bias_mcse) %>%
      #   ggplot(aes(x=year, y=rel_err_mean)) +
      #   geom_ribbon(aes(ymin=rel_err_lo2, ymax=rel_err_hi2), alpha=.5, fill='grey') +
      #   geom_ribbon(aes(ymin=rel_err_lo, ymax=rel_err_hi), alpha=.35, fill='grey') +
      #   geom_line(color='red') +
      #   coord_cartesian(ylim=c(-.5,.5)) +
      #   xlab("Year") +
      #   ylab("Relative error in SSB") +
      #   geom_hline(yintercept = 0, linetype=2, color='black') +
      #   facet_grid(rows=vars(em2), cols=vars(om2), labeller = label_parsed) +
      #   theme_bw() +
      #   theme(axis.text.x = element_text(size=8), plot.margin = unit(c(0.3,0.3,0.1,0.1), "in"))
      # png(file.path(plots_dir,paste0("1_ssb_",types[ty],"_medianCI.png")), width=7, height=7, units='in',res=100)
      # print(p)
      # grid::grid.text(unit(0.98,"npc"),0.5, label = 'Estimation model', rot = 270) # right
      # grid::grid.text(unit(0.5,"npc"),unit(.98,'npc'), label = 'Operating model', rot = 0)   # top)
      # dev.off()      

      # png(file.path(plots_dir,paste0("1_ssb_boxplots",types[ty],".png")), width=8, height=3, units='in',res=100)
      # print(ggplot(df.plot, aes(x=em.x, y=SSB.rel-1)) +
      #   geom_boxplot(aes(fill=em), outlier.shape = NA) +
      #   scale_fill_jco(name="Estimation model") +
      #   coord_cartesian(ylim=c(-.5,.5)) +
      #   xlab("Estimation model") +
      #   ylab("Relative error in SSB") +
      #   labs(title="Operating model") +
      #   geom_hline(yintercept = 0, linetype=2, color='black') +
      #   facet_wrap(vars(om2), nrow=1, labeller = label_parsed) +
      #   theme_bw() +
      #   theme(plot.title = element_text(hjust = 0.5)))
      # dev.off()
      
      png(file.path(plots_dir,paste0("1_ssb_boxplots",types[ty],".png")), width=8, height=3, units='in',res=100)
      print(ggplot(df.plot, aes(x=em.x, y=SSB.rel-1)) +
              # geom_boxplot(aes(fill=em), outlier.shape = NA) +
              stat_summary(aes(fill=em), fun.data=custom_boxplot_stat, geom="boxplot") + 
              scale_fill_jco(name="Estimation model") +
              coord_cartesian(ylim=c(-.5,.5)) +
              xlab("Estimation model") +
              ylab("Relative error in SSB") +
              labs(title="Operating model") +
              geom_hline(yintercept = 0, linetype=2, color='black') +
              facet_wrap(vars(om2), nrow=1, labeller = label_parsed) +
              theme_bw() +
              theme(plot.title = element_text(hjust = 0.5)))
      dev.off()
      
      p <- df.plot %>% 
        group_by(om2, em.x, em, sim) %>%
        summarize(rel_err_median = median(SSB.rel, na.rm=T), .groups = 'keep') %>%
        group_by(om2, em.x, em) %>%
        summarize(rel_err_var = var(rel_err_median), rel_err_mean = mean(rel_err_median)-1, K=sum(!is.na(rel_err_median)), .groups = 'keep') %>%
        mutate(rel_err_se = sqrt(rel_err_var/K)) %>%
        mutate(rel_err_lo = rel_err_mean - 1.96*rel_err_se, rel_err_hi = rel_err_mean + 1.96*rel_err_se) %>%         
        # mutate(rel_err_lo2 = rel_err_mean - 0.674*rel_err_se, rel_err_hi2 = rel_err_mean + 0.674*rel_err_se) %>%
        ggplot(aes(x=em.x, y=rel_err_mean)) +
        geom_hline(yintercept = 0, linetype=2, color='black') +
        geom_linerange(aes(ymin=rel_err_lo, ymax=rel_err_hi), size=.5) +
        # geom_linerange(aes(ymin=rel_err_lo2, ymax=rel_err_hi2), size=1.3) +
        geom_point(aes(fill=em), shape=21, size=3) +
        scale_fill_jco(name="Estimation model") +
        coord_cartesian(ylim=c(-.5,.5)) +
        xlab("Estimation model") +
        ylab("Relative error in SSB") +
        labs(title="Operating model") +
        facet_wrap(vars(om2), nrow=1, labeller = label_parsed) +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5))
      png(file.path(plots_dir,paste0("1_ssb_boxplots",types[ty],"_medianCI.png")), width=8, height=3, units='in',res=100)
      print(p)
      dev.off()
    
      # 2. Fishing mortality
      p <- ggplot(df.plot, aes(x=year, y=F.rel-1)) +
        stat_flquantiles(probs=bnds, alpha=0.7, fill="grey", geom="ribbon") +
        stat_summary(fun = "median", geom = "line", color = "red") +
        coord_cartesian(ylim=c(-.5,.5)) +
        xlab("Year") +
        ylab("Relative error in F") +
        geom_hline(yintercept = 0, linetype=2, color='black') +
        facet_grid(rows=vars(em2), cols=vars(om2), labeller = label_parsed) +
        theme_bw() +
        theme(axis.text.x = element_text(size=8), plot.margin = unit(c(0.3,0.3,0.1,0.1), "in"))
      png(file.path(plots_dir,paste0("2_F_",types[ty],".png")), width=7, height=7, units='in',res=100)
      print(p)
      grid::grid.text(unit(0.98,"npc"),0.5, label = 'Estimation model', rot = 270) # right
      grid::grid.text(unit(0.5,"npc"),unit(.98,'npc'), label = 'Operating model', rot = 0)   # top)
      dev.off()
      
      
      png(file.path(plots_dir,paste0("2_F_boxplots",types[ty],".png")), width=8, height=3, units='in',res=100)
      print(ggplot(df.plot, aes(x=em.x, y=F.rel-1)) +
              # geom_boxplot(aes(fill=em), outlier.shape = NA) +
              stat_summary(aes(fill=em), fun.data=custom_boxplot_stat, geom="boxplot") + 
        scale_fill_jco(name="Estimation model") +
        coord_cartesian(ylim=c(-.5,.5)) +
        xlab("Estimation model") +
        ylab("Relative error in F") +
        labs(title="Operating model") +
        geom_hline(yintercept = 0, linetype=2, color='black') +
        facet_wrap(vars(om2), nrow=1, labeller = label_parsed) +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5)))
      dev.off()   
      
      p <- df.plot %>% 
        filter(!is.na(F_sim) & is.finite(F_sim) & !is.nan(F_sim)) %>% 
        group_by(om2, em.x, em, sim) %>%
        summarize(rel_err_median = median(F.rel, na.rm=T), .groups = 'keep') %>%
        group_by(om2, em.x, em) %>%
        summarize(rel_err_var = var(rel_err_median), rel_err_mean = mean(rel_err_median)-1, K=sum(!is.na(rel_err_median)), .groups = 'keep') %>%
        mutate(rel_err_se = sqrt(rel_err_var/K)) %>%
        mutate(rel_err_lo = rel_err_mean - 1.96*rel_err_se, rel_err_hi = rel_err_mean + 1.96*rel_err_se) %>%         
        # mutate(rel_err_lo2 = rel_err_mean - 0.674*rel_err_se, rel_err_hi2 = rel_err_mean + 0.674*rel_err_se) %>%
        ggplot(aes(x=em.x, y=rel_err_mean)) +
        geom_hline(yintercept = 0, linetype=2, color='black') +
        geom_linerange(aes(ymin=rel_err_lo, ymax=rel_err_hi), size=.5) +
        # geom_linerange(aes(ymin=rel_err_lo2, ymax=rel_err_hi2), size=1.3) +
        geom_point(aes(fill=em), shape=21, size=3) +
        scale_fill_jco(name="Estimation model") +
        coord_cartesian(ylim=c(-.5,.5)) +
        xlab("Estimation model") +
        ylab("Relative error in F") +
        labs(title="Operating model") +
        facet_wrap(vars(om2), nrow=1, labeller = label_parsed) +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5))
      png(file.path(plots_dir,paste0("2_F_boxplots",types[ty],"_medianCI.png")), width=8, height=3, units='in',res=100)
      print(p)
      dev.off()

      # 3. SSB / SSB_40
      p <- df.plot %>%
        filter(!is.na(relB_sim) & is.finite(relB_sim) & !is.nan(relB_sim) & !is.na(relF_sim) & is.finite(relF_sim) & !is.nan(relF_sim)) %>% 
        ggplot(aes(x=year, y=relB.rel-1)) +
        stat_flquantiles(probs=bnds, alpha=0.7, fill="grey", geom="ribbon") +
        stat_summary(fun = "median", geom = "line", color = "red") +
        coord_cartesian(ylim=c(-.5,.5)) +
        xlab("Year") +
        ylab(expression("Relative error in"~SSB/SSB[40]["%"])) +
        geom_hline(yintercept = 0, linetype=2, color='black') +
        facet_grid(rows=vars(em2), cols=vars(om2), labeller = label_parsed) +
        theme_bw() +
        theme(axis.text.x = element_text(size=8), plot.margin = unit(c(0.3,0.3,0.1,0.1), "in"))
      png(file.path(plots_dir,paste0("3_relB_",types[ty],".png")), width=7, height=7, units='in',res=100)
      print(p)
      grid::grid.text(unit(0.98,"npc"),0.5, label = 'Estimation model', rot = 270) # right
      grid::grid.text(unit(0.5,"npc"),unit(.98,'npc'), label = 'Operating model', rot = 0)   # top)
      dev.off()
      
      png(file.path(plots_dir,paste0("3_relB_boxplots",types[ty],".png")), width=8, height=3, units='in',res=100)
      print(ggplot(df.plot, aes(x=em.x, y=relB.rel-1)) +
              # geom_boxplot(aes(fill=em), outlier.shape = NA) +
              stat_summary(aes(fill=em), fun.data=custom_boxplot_stat, geom="boxplot") + 
        scale_fill_jco(name="Estimation model") +
        coord_cartesian(ylim=c(-.5,.5)) +
        xlab("Estimation model") +
        ylab(expression("Relative error in"~SSB/SSB[40]["%"])) +
        labs(title="Operating model") +
        geom_hline(yintercept = 0, linetype=2, color='black') +
        facet_wrap(vars(om2), nrow=1, labeller = label_parsed) +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5)))
      dev.off() 
      
      p <- df.plot %>% 
        filter(!is.na(relB_sim) & is.finite(relB_sim) & !is.nan(relB_sim) & !is.na(relF_sim) & is.finite(relF_sim) & !is.nan(relF_sim)) %>% 
        group_by(om2, em.x, em, sim) %>%
        summarize(rel_err_median = median(relB.rel, na.rm=T), .groups = 'keep') %>%
        group_by(om2, em.x, em) %>%
        summarize(rel_err_var = var(rel_err_median), rel_err_mean = mean(rel_err_median)-1, K=sum(!is.na(rel_err_median)), .groups = 'keep') %>%
        mutate(rel_err_se = sqrt(rel_err_var/K)) %>%
        mutate(rel_err_lo = rel_err_mean - 1.96*rel_err_se, rel_err_hi = rel_err_mean + 1.96*rel_err_se) %>%         
        # mutate(rel_err_lo2 = rel_err_mean - 0.674*rel_err_se, rel_err_hi2 = rel_err_mean + 0.674*rel_err_se) %>%
        ggplot(aes(x=em.x, y=rel_err_mean)) +
        geom_hline(yintercept = 0, linetype=2, color='black') +
        geom_linerange(aes(ymin=rel_err_lo, ymax=rel_err_hi), size=.5) +
        # geom_linerange(aes(ymin=rel_err_lo2, ymax=rel_err_hi2), size=1.3) +
        geom_point(aes(fill=em), shape=21, size=3) +
        scale_fill_jco(name="Estimation model") +
        coord_cartesian(ylim=c(-.5,.5)) +
        xlab("Estimation model") +
        ylab(expression("Relative error in"~SSB/SSB[40]["%"])) +
        labs(title="Operating model") +
        facet_wrap(vars(om2), nrow=1, labeller = label_parsed) +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5))
      png(file.path(plots_dir,paste0("3_relB_boxplots",types[ty],"_medianCI.png")), width=8, height=3, units='in',res=100)
      print(p)
      dev.off()

      # 4. F / F40
      p <- df.plot %>%
        filter(!is.na(relB_sim) & is.finite(relB_sim) & !is.nan(relB_sim) & !is.na(relF_sim) & is.finite(relF_sim) & !is.nan(relF_sim)) %>% 
        ggplot(aes(x=year, y=relF.rel-1)) +
        stat_flquantiles(probs=bnds, alpha=0.7, fill="grey", geom="ribbon") +
        stat_summary(fun = "median", geom = "line", color = "red") +
        coord_cartesian(ylim=c(-.5,.5)) +
        xlab("Year") +
        ylab(expression("Relative error in"~F/F[40]["%"])) +
        geom_hline(yintercept = 0, linetype=2, color='black') +
        facet_grid(rows=vars(em2), cols=vars(om2), labeller = label_parsed) +
        theme_bw() +
        theme(axis.text.x = element_text(size=8), plot.margin = unit(c(0.3,0.3,0.1,0.1), "in"))
      png(file.path(plots_dir,paste0("4_relF_",types[ty],".png")), width=7, height=7, units='in',res=100)
      print(p)
      grid::grid.text(unit(0.98,"npc"),0.5, label = 'Estimation model', rot = 270) # right
      grid::grid.text(unit(0.5,"npc"),unit(.98,'npc'), label = 'Operating model', rot = 0)   # top)
      dev.off()
      
      png(file.path(plots_dir,paste0("4_relF_boxplots",types[ty],".png")), width=8, height=3, units='in',res=100)
      print(ggplot(df.plot, aes(x=em.x, y=relF.rel-1)) +
              # geom_boxplot(aes(fill=em), outlier.shape = NA) +
              stat_summary(aes(fill=em), fun.data=custom_boxplot_stat, geom="boxplot") + 
        scale_fill_jco(name="Estimation model") +
        coord_cartesian(ylim=c(-.5,.5)) +
        xlab("Estimation model") +
        ylab(expression("Relative error in"~F/F[40]["%"])) +
        labs(title="Operating model") +
        geom_hline(yintercept = 0, linetype=2, color='black') +
        facet_wrap(vars(om2), nrow=1, labeller = label_parsed) +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5)))
      dev.off()  
      
      p <- df.plot %>% 
        filter(!is.na(relB_sim) & is.finite(relB_sim) & !is.nan(relB_sim) & !is.na(relF_sim) & is.finite(relF_sim) & !is.nan(relF_sim)) %>% 
        group_by(om2, em.x, em, sim) %>%
        summarize(rel_err_median = median(relF.rel, na.rm=T), .groups = 'keep') %>%
        group_by(om2, em.x, em) %>%
        summarize(rel_err_var = var(rel_err_median), rel_err_mean = mean(rel_err_median)-1, K=sum(!is.na(rel_err_median)), .groups = 'keep') %>%
        mutate(rel_err_se = sqrt(rel_err_var/K)) %>%
        mutate(rel_err_lo = rel_err_mean - 1.96*rel_err_se, rel_err_hi = rel_err_mean + 1.96*rel_err_se) %>%         
        # mutate(rel_err_lo2 = rel_err_mean - 0.674*rel_err_se, rel_err_hi2 = rel_err_mean + 0.674*rel_err_se) %>%
        ggplot(aes(x=em.x, y=rel_err_mean)) +
        geom_hline(yintercept = 0, linetype=2, color='black') +
        geom_linerange(aes(ymin=rel_err_lo, ymax=rel_err_hi), size=.5) +
        # geom_linerange(aes(ymin=rel_err_lo2, ymax=rel_err_hi2), size=1.3) +
        geom_point(aes(fill=em), shape=21, size=3) +
        scale_fill_jco(name="Estimation model") +
        coord_cartesian(ylim=c(-.5,.5)) +
        xlab("Estimation model") +
        ylab(expression("Relative error in"~F/F[40]["%"])) +
        labs(title="Operating model") +
        facet_wrap(vars(om2), nrow=1, labeller = label_parsed) +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5))
      png(file.path(plots_dir,paste0("4_relF_boxplots",types[ty],"_medianCI.png")), width=8, height=3, units='in',res=100)
      print(p)
      dev.off()

      # 5. Catch
      p <- ggplot(df.plot, aes(x=year, y=catch.rel-1)) +
        stat_flquantiles(probs=bnds, alpha=0.7, fill="grey", geom="ribbon") +
        stat_summary(fun = "median", geom = "line", color = "red") +
        coord_cartesian(ylim=c(-.5,.5)) +
        xlab("Year") +
        ylab("Relative error in Catch") +
        geom_hline(yintercept = 0, linetype=2, color='black') +
        facet_grid(rows=vars(em2), cols=vars(om2), labeller = label_parsed) +
        theme_bw() +
        theme(axis.text.x = element_text(size=8), plot.margin = unit(c(0.3,0.3,0.1,0.1), "in"))
      png(file.path(plots_dir,paste0("5_catch_",types[ty],".png")), width=7, height=7, units='in',res=100)
      print(p)
      grid::grid.text(unit(0.98,"npc"),0.5, label = 'Estimation model', rot = 270) # right
      grid::grid.text(unit(0.5,"npc"),unit(.98,'npc'), label = 'Operating model', rot = 0)   # top)
      dev.off()
      
      # boxplots (collapse time series)
      png(file.path(plots_dir,paste0("5_catch_boxplots",types[ty],".png")), width=8, height=3, units='in',res=100)
      print(ggplot(df.plot, aes(x=em.x, y=catch.rel-1)) +
              # geom_boxplot(aes(fill=em), outlier.shape = NA) +
              stat_summary(aes(fill=em), fun.data=custom_boxplot_stat, geom="boxplot") + 
        scale_fill_jco(name="Estimation model") +
        coord_cartesian(ylim=c(-.5,.5)) +
        xlab("Estimation model") +
        ylab("Relative error in Catch") +
        labs(title="Operating model") +
        geom_hline(yintercept = 0, linetype=2, color='black') +
        facet_wrap(vars(om2), nrow=1, labeller = label_parsed) +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5)))
      dev.off()  
      
      p <- df.plot %>% 
        filter(!is.na(relB_sim) & is.finite(relB_sim) & !is.nan(relB_sim) & !is.na(relF_sim) & is.finite(relF_sim) & !is.nan(relF_sim)) %>% 
        group_by(om2, em.x, em, sim) %>%
        summarize(rel_err_median = median(catch.rel, na.rm=T), .groups = 'keep') %>%
        group_by(om2, em.x, em) %>%
        summarize(rel_err_var = var(rel_err_median), rel_err_mean = mean(rel_err_median)-1, K=sum(!is.na(rel_err_median)), .groups = 'keep') %>%
        mutate(rel_err_se = sqrt(rel_err_var/K)) %>%
        mutate(rel_err_lo = rel_err_mean - 1.96*rel_err_se, rel_err_hi = rel_err_mean + 1.96*rel_err_se) %>%         
        # mutate(rel_err_lo2 = rel_err_mean - 0.674*rel_err_se, rel_err_hi2 = rel_err_mean + 0.674*rel_err_se) %>%
        ggplot(aes(x=em.x, y=rel_err_mean)) +
        geom_hline(yintercept = 0, linetype=2, color='black') +
        geom_linerange(aes(ymin=rel_err_lo, ymax=rel_err_hi), size=.5) +
        # geom_linerange(aes(ymin=rel_err_lo2, ymax=rel_err_hi2), size=1.3) +
        geom_point(aes(fill=em), shape=21, size=3) +
        scale_fill_jco(name="Estimation model") +
        coord_cartesian(ylim=c(-.5,.5)) +
        xlab("Estimation model") +
        ylab("Relative error in Catch") +
        labs(title="Operating model") +
        facet_wrap(vars(om2), nrow=1, labeller = label_parsed) +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5))
      png(file.path(plots_dir,paste0("5_catch_boxplots",types[ty],"_medianCI.png")), width=8, height=3, units='in',res=100)
      print(p)
      dev.off()

      # 6. Recruitment
      p <- ggplot(df.plot, aes(x=year, y=R.rel-1)) +
        stat_flquantiles(probs=bnds, alpha=0.7, fill="grey", geom="ribbon") +
        stat_summary(fun = "median", geom = "line", color = "red") +
        coord_cartesian(ylim=c(-.5,.5)) +
        xlab("Year") +
        ylab("Relative error in Recruitment") +
        geom_hline(yintercept = 0, linetype=2, color='black') +
        facet_grid(rows=vars(em2), cols=vars(om2), labeller = label_parsed) +
        theme_bw() +
        theme(axis.text.x = element_text(size=8), plot.margin = unit(c(0.3,0.3,0.1,0.1), "in"))
      png(file.path(plots_dir,paste0("6_R_",types[ty],".png")), width=7, height=7, units='in',res=100)
      print(p)
      grid::grid.text(unit(0.98,"npc"),0.5, label = 'Estimation model', rot = 270) # right
      grid::grid.text(unit(0.5,"npc"),unit(.98,'npc'), label = 'Operating model', rot = 0)   # top)
      dev.off()

      png(file.path(plots_dir,paste0("6_R_boxplots",types[ty],".png")), width=8, height=3, units='in',res=100)
      print(ggplot(df.plot, aes(x=em.x, y=R.rel-1)) +
              # geom_boxplot(aes(fill=em), outlier.shape = NA) +
              stat_summary(aes(fill=em), fun.data=custom_boxplot_stat, geom="boxplot") + 
        scale_fill_jco(name="Estimation model") +
        coord_cartesian(ylim=c(-.5,.5)) +
        xlab("Estimation model") +
        ylab("Relative error in Recruitment") +
        labs(title="Operating model") +
        geom_hline(yintercept = 0, linetype=2, color='black') +
        facet_wrap(vars(om2), nrow=1, labeller = label_parsed) +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5)))
      dev.off()           
      
      p <- df.plot %>% 
        # filter(!is.na(relB_sim) & is.finite(relB_sim) & !is.nan(relB_sim) & !is.na(relF_sim) & is.finite(relF_sim) & !is.nan(relF_sim)) %>% 
        group_by(om2, em.x, em, sim) %>%
        summarize(rel_err_median = median(R.rel, na.rm=T), .groups = 'keep') %>%
        group_by(om2, em.x, em) %>%
        summarize(rel_err_var = var(rel_err_median), rel_err_mean = mean(rel_err_median)-1, K=sum(!is.na(rel_err_median)), .groups = 'keep') %>%
        mutate(rel_err_se = sqrt(rel_err_var/K)) %>%
        mutate(rel_err_lo = rel_err_mean - 1.96*rel_err_se, rel_err_hi = rel_err_mean + 1.96*rel_err_se) %>%         
        # mutate(rel_err_lo2 = rel_err_mean - 0.674*rel_err_se, rel_err_hi2 = rel_err_mean + 0.674*rel_err_se) %>%
        ggplot(aes(x=em.x, y=rel_err_mean)) +
        geom_hline(yintercept = 0, linetype=2, color='black') +
        geom_linerange(aes(ymin=rel_err_lo, ymax=rel_err_hi), size=.5) +
        # geom_linerange(aes(ymin=rel_err_lo2, ymax=rel_err_hi2), size=1.3) +
        geom_point(aes(fill=em), shape=21, size=3) +
        scale_fill_jco(name="Estimation model") +
        coord_cartesian(ylim=c(-.5,.5)) +
        xlab("Estimation model") +
        ylab("Relative error in Recruitment") +
        labs(title="Operating model") +
        facet_wrap(vars(om2), nrow=1, labeller = label_parsed) +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5))
      png(file.path(plots_dir,paste0("6_R_boxplots",types[ty],"_medianCI.png")), width=8, height=3, units='in',res=100)
      print(p)
      dev.off()

# ---------------------------------------------------------------------------
      # now with TMB bias correction
      if(plot.eps){
        # 1. SSB
        p <- ggplot(df.plot, aes(x=year, y=SSB.rel.bc-1)) +
          stat_flquantiles(probs=bnds, alpha=0.7, fill="grey", geom="ribbon") +
          stat_summary(fun = "median", geom = "line", color = "red") +
          coord_cartesian(ylim=c(-.5,.5)) +
          xlab("Year") +
          ylab("Relative error in SSB") +
          geom_hline(yintercept = 0, linetype=2, color='black') +
          facet_grid(rows=vars(em2), cols=vars(om2), labeller = label_parsed) +
          theme_bw() +
          theme(axis.text.x = element_text(size=8), plot.margin = unit(c(0.3,0.3,0.1,0.1), "in"))

        png(file.path(plots_dir,paste0("1_ssb_",types[ty],"_bc.png")), width=7, height=7, units='in',res=100)
        print(p)
        grid::grid.text(unit(0.98,"npc"),0.5, label = 'Estimation model', rot = 270) # right
        grid::grid.text(unit(0.5,"npc"),unit(.98,'npc'), label = 'Operating model', rot = 0)   # top)
        dev.off()

        png(file.path(plots_dir,paste0("1_ssb_boxplots",types[ty],"_bc.png")), width=8, height=3, units='in',res=100)
        print(ggplot(df.plot, aes(x=em.x, y=SSB.rel.bc-1)) +
          geom_boxplot(aes(fill=em), outlier.shape = NA) +
          scale_fill_jco(name="Estimation model") +
          coord_cartesian(ylim=c(-.5,.5)) +
          xlab("Estimation model") +
          ylab("Relative error in SSB") +
          labs(title="Operating model") +
          geom_hline(yintercept = 0, linetype=2, color='black') +
          facet_wrap(vars(om2), nrow=1, labeller = label_parsed) +
          theme_bw() +
          theme(plot.title = element_text(hjust = 0.5)))
        dev.off() 

        # 2. Fishing mortality
        p <- ggplot(df.plot, aes(x=year, y=F.rel.bc-1)) +
          stat_flquantiles(probs=bnds, alpha=0.7, fill="grey", geom="ribbon") +
          stat_summary(fun = "median", geom = "line", color = "red") +
          coord_cartesian(ylim=c(-.5,.5)) +
          xlab("Year") +
          ylab("Relative error in F") +
          geom_hline(yintercept = 0, linetype=2, color='black') +
          facet_grid(rows=vars(em2), cols=vars(om2), labeller = label_parsed) +
          theme_bw() +
          theme(axis.text.x = element_text(size=8), plot.margin = unit(c(0.3,0.3,0.1,0.1), "in"))

        png(file.path(plots_dir,paste0("2_F_",types[ty],"_bc.png")), width=7, height=7, units='in',res=100)
        print(p)
        grid::grid.text(unit(0.98,"npc"),0.5, label = 'Estimation model', rot = 270) # right
        grid::grid.text(unit(0.5,"npc"),unit(.98,'npc'), label = 'Operating model', rot = 0)   # top)
        dev.off()

        png(file.path(plots_dir,paste0("2_F_boxplots",types[ty],"_bc.png")), width=8, height=3, units='in',res=100)
        print(ggplot(df.plot, aes(x=em.x, y=F.rel.bc-1)) +
          geom_boxplot(aes(fill=em), outlier.shape = NA) +
          scale_fill_jco(name="Estimation model") +
          coord_cartesian(ylim=c(-.5,.5)) +
          xlab("Estimation model") +
          ylab("Relative error in F") +
          labs(title="Operating model") +
          geom_hline(yintercept = 0, linetype=2, color='black') +
          facet_wrap(vars(om2), nrow=1, labeller = label_parsed) +
          theme_bw() +
          theme(plot.title = element_text(hjust = 0.5)))
        dev.off()     

        # 3. SSB / SSB40
        p <- ggplot(df.plot, aes(x=year, y=relB.rel.bc-1)) +
          stat_flquantiles(probs=bnds, alpha=0.7, fill="grey", geom="ribbon") +
          stat_summary(fun = "median", geom = "line", color = "red") +
          coord_cartesian(ylim=c(-.5,.5)) +
          xlab("Year") +
          ylab(expression("Relative error in"~SSB/SSB[40]["%"])) +
          geom_hline(yintercept = 0, linetype=2, color='black') +
          facet_grid(rows=vars(em2), cols=vars(om2), labeller = label_parsed) +
          theme_bw() +
          theme(axis.text.x = element_text(size=8), plot.margin = unit(c(0.3,0.3,0.1,0.1), "in"))

        png(file.path(plots_dir,paste0("3_relB_",types[ty],"_bc.png")), width=7, height=7, units='in',res=100)
        print(p)
        grid::grid.text(unit(0.98,"npc"),0.5, label = 'Estimation model', rot = 270) # right
        grid::grid.text(unit(0.5,"npc"),unit(.98,'npc'), label = 'Operating model', rot = 0)   # top)
        dev.off()

        png(file.path(plots_dir,paste0("3_relB_boxplots",types[ty],"_bc.png")), width=8, height=3, units='in',res=100)
        print(ggplot(df.plot, aes(x=em.x, y=relB.rel.bc-1)) +
          geom_boxplot(aes(fill=em), outlier.shape = NA) +
          scale_fill_jco(name="Estimation model") +
          coord_cartesian(ylim=c(-.5,.5)) +
          xlab("Estimation model") +
          ylab(expression("Relative error in"~SSB/SSB[40]["%"])) +
          labs(title="Operating model") +
          geom_hline(yintercept = 0, linetype=2, color='black') +
          facet_wrap(vars(om2), nrow=1, labeller = label_parsed) +
          theme_bw() +
          theme(plot.title = element_text(hjust = 0.5)))
        dev.off()   

        # 4. F / F40
        p <- ggplot(df.plot, aes(x=year, y=relF.rel.bc-1)) +
          stat_flquantiles(probs=bnds, alpha=0.7, fill="grey", geom="ribbon") +
          stat_summary(fun = "median", geom = "line", color = "red") +
          coord_cartesian(ylim=c(-.5,.5)) +
          xlab("Year") +
          ylab(expression("Relative error in"~F/F[40]["%"])) +
          geom_hline(yintercept = 0, linetype=2, color='black') +
          facet_grid(rows=vars(em2), cols=vars(om2), labeller = label_parsed) +
          theme_bw() +
          theme(axis.text.x = element_text(size=8), plot.margin = unit(c(0.3,0.3,0.1,0.1), "in"))
        png(file.path(plots_dir,paste0("4_relF_",types[ty],"_bc.png")), width=7, height=7, units='in',res=100)
        print(p)
        grid::grid.text(unit(0.98,"npc"),0.5, label = 'Estimation model', rot = 270) # right
        grid::grid.text(unit(0.5,"npc"),unit(.98,'npc'), label = 'Operating model', rot = 0)   # top)
        dev.off()

        png(file.path(plots_dir,paste0("4_relF_boxplots",types[ty],"_bc.png")), width=8, height=3, units='in',res=100)
        print(ggplot(df.plot, aes(x=em.x, y=relF.rel.bc-1)) +
          geom_boxplot(aes(fill=em), outlier.shape = NA) +
          scale_fill_jco(name="Estimation model") +
          coord_cartesian(ylim=c(-.5,.5)) +
          xlab("Estimation model") +
          ylab(expression("Relative error in"~F/F[40]["%"])) +
          labs(title="Operating model") +
          geom_hline(yintercept = 0, linetype=2, color='black') +
          facet_wrap(vars(om2), nrow=1, labeller = label_parsed) +
          theme_bw() +
          theme(plot.title = element_text(hjust = 0.5)))
        dev.off()   

        # 5. Catch
        p <- ggplot(df.plot, aes(x=year, y=catch.rel.bc-1)) +
          stat_flquantiles(probs=bnds, alpha=0.7, fill="grey", geom="ribbon") +
          stat_summary(fun = "median", geom = "line", color = "red") +
          coord_cartesian(ylim=c(-.5,.5)) +
          xlab("Year") +
          ylab("Relative error in Catch") +
          geom_hline(yintercept = 0, linetype=2, color='black') +
          facet_grid(rows=vars(em2), cols=vars(om2), labeller = label_parsed) +
          theme_bw() +
          theme(axis.text.x = element_text(size=8), plot.margin = unit(c(0.3,0.3,0.1,0.1), "in"))
        png(file.path(plots_dir,paste0("5_catch_",types[ty],"_bc.png")), width=7, height=7, units='in',res=100)
        print(p)
        grid::grid.text(unit(0.98,"npc"),0.5, label = 'Estimation model', rot = 270) # right
        grid::grid.text(unit(0.5,"npc"),unit(.98,'npc'), label = 'Operating model', rot = 0)   # top)
        dev.off()

        png(file.path(plots_dir,paste0("5_catch_boxplots",types[ty],"_bc.png")), width=8, height=3, units='in',res=100)
        print(ggplot(df.plot, aes(x=em.x, y=catch.rel.bc-1)) +
          geom_boxplot(aes(fill=em), outlier.shape = NA) +
          scale_fill_jco(name="Estimation model") +
          coord_cartesian(ylim=c(-.5,.5)) +
          xlab("Estimation model") +
          ylab("Relative error in Catch") +
          labs(title="Operating model") +
          geom_hline(yintercept = 0, linetype=2, color='black') +
          facet_wrap(vars(om2), nrow=1, labeller = label_parsed) +
          theme_bw() +
          theme(plot.title = element_text(hjust = 0.5)))
        dev.off() 

        # 6. Recruitment                 
        p <- ggplot(df.plot, aes(x=year, y=R.rel.bc-1)) +
          stat_flquantiles(probs=bnds, alpha=0.7, fill="grey", geom="ribbon") +
          stat_summary(fun = "median", geom = "line", color = "red") +
          coord_cartesian(ylim=c(-.5,.5)) +
          xlab("Year") +
          ylab("Relative error in Recruitment") +
          geom_hline(yintercept = 0, linetype=2, color='black') +
          facet_grid(rows=vars(em2), cols=vars(om2), labeller = label_parsed) +
          theme_bw() +
          theme(axis.text.x = element_text(size=8), plot.margin = unit(c(0.3,0.3,0.1,0.1), "in"))
        png(file.path(plots_dir,paste0("6_R_",types[ty],"_bc.png")), width=7, height=7, units='in',res=100)
        print(p)
        grid::grid.text(unit(0.98,"npc"),0.5, label = 'Estimation model', rot = 270) # right
        grid::grid.text(unit(0.5,"npc"),unit(.98,'npc'), label = 'Operating model', rot = 0)   # top)
        dev.off()

        png(file.path(plots_dir,paste0("6_R_boxplots",types[ty],"_bc.png")), width=8, height=3, units='in',res=100)
        print(ggplot(df.plot, aes(x=em.x, y=R.rel.bc-1)) +
          geom_boxplot(aes(fill=em), outlier.shape = NA) +
          scale_fill_jco(name="Estimation model") +
          coord_cartesian(ylim=c(-.5,.5)) +
          xlab("Estimation model") +
          ylab("Relative error in Recruitment") +
          labs(title="Operating model") +
          geom_hline(yintercept = 0, linetype=2, color='black') +
          facet_wrap(vars(om2), nrow=1, labeller = label_parsed) +
          theme_bw() +
          theme(plot.title = element_text(hjust = 0.5)))
        dev.off()         
      }
    }
  }
}
