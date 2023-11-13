
#' @param convergedONLY A boolean, when TRUE plots only show converged runs, default = TRUE.
#' @param outdir A string for the directory where a plot folder will be generated


# results <- perfMet
# convergedONLY = TRUE

plotResults <- function(results = NULL, convergedONLY = TRUE, outdir = here::here()){
  
  # Simulation summary
  results %>% count(OMshortName, EMshortName) %>% mutate(nsim = n/40) # /40 years so nsim = number of full simulations for each OM/EM
  # results %>% count(OMshortName) %>% mutate(nsim = n/40/4) # Full sims of all 4 correctly specified EMs
  # results %>% filter(OMshortName == 15) %>% select(EMshortName) %>% unique()
  
  # setdiff(subsetOM$OMname, results$OMshortName) # Check what simulations still not run from subsetOM (in Catchability_ecov_sims.Rmd)
  
  if(convergedONLY == TRUE){
    results <- results %>% filter(EM_converged == TRUE)
  }
  
  # Pull out results that don't have time series
  #!!! problem with sim == 351 (has both H-L and Fmsy )
  # results <- results %>% filter(sim %in% c(351, 352, 227, 228, 231, 234, 235, 236, 237, 238, 239, 240, 241, 243, 249, 250, 101, 102) == FALSE)
  singleResults <- results %>% group_by(sim) %>% dplyr::summarize(seed = unique(seed), Fhist = unique(F_hist), ageComp_sig = unique(ageComp_sig), log_catch_sig = unique(log_catch_sig), log_index_sig = unique(log_index_sig),
                                                                                 OM_ecov_effect = unique(OM_ecov_effect), OM_ecov_process_cor = unique(OM_ecov_process_cor), OM_ecov_process_obs_sig = unique(OM_ecov_process_obs_sig), OM_ecov_process_sig = unique(OM_ecov_process_sig), OMshortName = unique(OMshortName),
                                                                                                         EM_miss_q = unique(EM_miss_q), EM_miss_season = unique(EM_miss_season), EMshortName = unique(EMshortName),
                                                                                 EM_MohnsRho_SSB = unique(EM_MohnsRho_SSB), EM_MohnsRho_F = unique(EM_MohnsRho_F), EM_MohnsRho_R = unique(EM_MohnsRho_R), EM_converged = unique(EM_converged), 
                                                                                 relSSBMSY = unique(relSSBMSY), relFMSY = unique(relFMSY), relMSY = unique(relMSY), 
                                                                                 
                                                                                 relselInd1_1 = unique(relselInd1_1),
                                                                                 relselInd1_2 = unique(relselInd1_2),
                                                                                 relselInd1_3 = unique(relselInd1_3),
                                                                                 relselInd1_4 = unique(relselInd1_4),
                                                                                 relselInd1_5 = unique(relselInd1_5),
                                                                                 relselInd1_6 = unique(relselInd1_6),
                                                                                 relselInd1_7 = unique(relselInd1_7),
                                                                                 relselInd1_8 = unique(relselInd1_8),
                                                                                 relselInd1_9 = unique(relselInd1_9),
                                                                                 relselInd1_10 = unique(relselInd1_10),
                                                                                 
                                                                                 relselInd2_1 = unique(relselInd2_1),
                                                                                 relselInd2_2 = unique(relselInd2_2),
                                                                                 relselInd2_3 = unique(relselInd2_3),
                                                                                 relselInd2_4 = unique(relselInd2_4),
                                                                                 relselInd2_5 = unique(relselInd2_5),
                                                                                 relselInd2_6 = unique(relselInd2_6),
                                                                                 relselInd2_7 = unique(relselInd2_7),
                                                                                 relselInd2_8 = unique(relselInd2_8),
                                                                                 relselInd2_9 = unique(relselInd2_9),
                                                                                 relselInd2_10 = unique(relselInd2_10))
    
                                                                                 
  # Time series include status SSB/F/Y, relSSB, relF, relR, relq both indices, relFAA, relNAA, relCAA
  
  # MohnsRho_SSB
  mohnsRho_SSB_series <- singleResults %>% 
    ggplot() + 
    geom_point(aes(x=sim, y=EM_MohnsRho_SSB, color = OM_ecov_effect)) # Plot by sim
  ggsave(mohnsRho_SSB_series, filename = here::here(outfile, "mohnsRho_SSB_series.png"), width = 10)
  
  mohnsRho_SSB <- singleResults %>%
    ggplot()+
    geom_boxplot(aes(x=EMshortName, y=EM_MohnsRho_SSB)) + facet_grid(. ~ OM_ecov_process_sig + OM_ecov_process_cor + OM_ecov_process_obs_sig + ageComp_sig + log_index_sig) + # Facet by OM settings
    geom_hline(aes(yintercept = 0), color="red")
    # ggtitle("OM_ecov_process_sig")
  ggsave(mohnsRho_SSB, filename = here::here(outfile, "mohnsRho_SSB.png"), width = 20)
  
  # MohnsRho_F
  mohnsRho_F <- singleResults %>%
    ggplot()+
    geom_boxplot(aes(x=EMshortName, y=EM_MohnsRho_F)) + facet_grid(. ~ OM_ecov_process_sig + OM_ecov_process_cor + OM_ecov_process_obs_sig + ageComp_sig + log_index_sig) + # Facet by OM settings
    geom_hline(aes(yintercept = 0), color="red")
  ggsave(mohnsRho_F, filename = here::here(outfile, "mohnsRho_F.png"), width = 20)
  
  # MohnsRho_R
  mohnsRho_R <- singleResults %>%
    ggplot()+
    geom_boxplot(aes(x=EMshortName, y=EM_MohnsRho_R)) + facet_grid(. ~ OM_ecov_process_sig + OM_ecov_process_cor + OM_ecov_process_obs_sig + ageComp_sig + log_index_sig) + # Facet by OM settings
    geom_hline(aes(yintercept = 0), color="red")
  ggsave(mohnsRho_R, filename = here::here(outfile, "mohnsRho_R.png"), width = 20)
  
  # Relative SSBmsy
  relSSBMSY <- singleResults %>%
    ggplot()+
    geom_boxplot(aes(x=EMshortName, y=relSSBMSY)) + facet_grid(. ~ OM_ecov_process_sig + OM_ecov_process_cor + OM_ecov_process_obs_sig + ageComp_sig + log_index_sig) + # Facet by OM settings
    geom_hline(aes(yintercept = 1), color="red")
  ggsave(relSSBMSY, filename = here::here(outfile, "relSSBMSY.png"), width = 20)
  
  # Relative Fmsy
  relFMSY <- singleResults %>%
    ggplot()+
    geom_boxplot(aes(x=EMshortName, y=relFMSY)) + facet_grid(. ~ OM_ecov_process_sig + OM_ecov_process_cor + OM_ecov_process_obs_sig + ageComp_sig + log_index_sig) + # Facet by OM settings
    geom_hline(aes(yintercept = 1), color="red")
  ggsave(relFMSY, filename = here::here(outfile, "relFMSY.png"), width = 20)
  
  # Relative MSY
  relMSY <- singleResults %>%
    ggplot()+
    geom_boxplot(aes(x=EMshortName, y=relMSY)) + facet_grid(. ~ OM_ecov_process_sig + OM_ecov_process_cor + OM_ecov_process_obs_sig + ageComp_sig + log_index_sig) + # Facet by OM settings
    geom_hline(aes(yintercept = 1), color="red")
  ggsave(relMSY, filename = here::here(outfile, "relMSY.png"), width = 20)
  
  # Relative Index 1 
  relSelInd1 <- singleResults %>% 
    pivot_longer(., cols = c(relselInd1_1, relselInd1_2, relselInd1_3, relselInd1_4, relselInd1_5, relselInd1_6, relselInd1_7, relselInd1_8, relselInd1_9, relselInd1_10), names_to = "age", names_prefix = "relselInd1_", values_to = "sel")  
  relSelInd1$age <- factor(relSelInd1$age, levels = c(1:10))
  relSelInd1plot <- relSelInd1 %>% 
    ggplot()+
    geom_boxplot(aes(x=EMshortName, y=sel, fill = age)) + facet_grid(. ~ OM_ecov_process_sig + OM_ecov_process_cor + OM_ecov_process_obs_sig + ageComp_sig + log_index_sig) + # Facet by OM settings
    geom_hline(aes(yintercept =1), color="red")
  ggsave(relSelInd1plot, filename = here::here(outfile, "relSelInd1.png"), width = 20)
  
  # Relative Index 2 
  relSelInd2 <- singleResults %>% 
    pivot_longer(., cols = c(relselInd2_1, relselInd2_2, relselInd2_3, relselInd2_4, relselInd2_5, relselInd2_6, relselInd2_7, relselInd2_8, relselInd2_9, relselInd2_10), names_to = "age", names_prefix = "relselInd2_", values_to = "sel") 
  relSelInd2$age <- factor(relSelInd2$age, levels = c(1:10))
  relSelInd2plot <- relSelInd2 %>%
    ggplot()+
    geom_boxplot(aes(x=EMshortName, y=sel, fill = age)) + facet_grid(. ~ OM_ecov_process_sig + OM_ecov_process_cor + OM_ecov_process_obs_sig + ageComp_sig + log_index_sig) + # Facet by OM settings
    geom_hline(aes(yintercept =1), color="red")
     # ageComp_sig impacts spread of results across all ages
  ggsave(relSelInd2plot, filename = here::here(outfile, "relSelInd2.png"), width = 20)
    
  ### Time series plots ###
  #include status SSB/F/Y, relSSB, relF, relR, relq both indices, relFAA, relNAA, relCAA
  
  # Relative SSB
  #   # Time series
  # results %>% 
  #   ggplot()+
  #   geom_line(aes(x=Year, y = relSSB, group = seed), alpha = 0.5) + # Could color lines by F history
  #   geom_hline(aes(yintercept = 1.0), color="red")
    # Terminal year
  relSSB <- results %>% filter(Year == 40) %>%
    ggplot()+
    geom_boxplot(aes(x=EMshortName, y=relSSB)) + facet_grid(. ~ OM_ecov_process_sig + OM_ecov_process_cor + OM_ecov_process_obs_sig + ageComp_sig + log_index_sig) + # Facet by OM settings
    geom_hline(aes(yintercept = 1), color="red")
  ggsave(relSSB, filename = here::here(outfile, "relSSB.png"), width = 10)
  
  # # EM SSB
  # results %>%
  #   ggplot()+
  #   geom_line(aes(x=Year, y = EM_SSB, group = seed, color = OM_ecov_effect), alpha = 0.5)  # Could color lines by F history
  # # EM F
  # results %>%
  #   ggplot()+
  #   geom_line(aes(x=Year, y = EM_F, group = seed, color = OM_ecov_effect), alpha = 0.5)  # Could color lines by F history
  # # EM R
  # results %>%
  #   ggplot()+
  #   geom_line(aes(x=Year, y = EM_R, group = seed, color = OM_ecov_effect), alpha = 0.5)  # Could color lines by F history
  # 
  # # Look at subset of results where F was estimated low for first half of time series
  # res_index <- results[which(results$Year == 1 & results$EM_F < 0.6),]$seed
  # results %>% filter(seed %in% res_index) %>%
  #   summary()
  #   ggplot() + geom_line(aes(x=Year, y = EM_F, group = seed), alpha = 0.5)
  
  # Relative F
  #   # Time series
  # results %>% 
  #   ggplot()+
  #   geom_line(aes(x=Year, y = relF, group = seed), alpha = 0.5) + # Could color lines by F history
  #   geom_hline(aes(yintercept = 1.0), color="red")
  # Terminal year
  relativeF_terminal <- results %>% filter(Year == 40) %>%
    ggplot()+
    geom_boxplot(aes(x=EMshortName, y=relF)) + facet_grid(. ~ OM_ecov_process_sig + OM_ecov_process_cor + OM_ecov_process_obs_sig + ageComp_sig + log_index_sig) + # Facet by OM settings
    geom_hline(aes(yintercept = 1), color="red")
  ggsave(relativeF_terminal, filename = here::here(outfile, "relativeF_terminal.png"), width = 20)
  
  # Relative R
  #   # Time series
  # results %>% 
  #   ggplot()+
  #   geom_line(aes(x=Year, y = relR, group = seed, color = F_hist), alpha = 0.5) + # Could color lines by F history
  #   geom_hline(aes(yintercept = 1.0), color="red")
    # Terminal year
  relativeR_terminal <- results %>% filter(Year == 40) %>%
    ggplot()+
    geom_boxplot(aes(x=EMshortName, y=relR)) + facet_grid(. ~ OM_ecov_process_sig + OM_ecov_process_cor + OM_ecov_process_obs_sig + ageComp_sig + log_index_sig) + # Facet by OM settings
    geom_hline(aes(yintercept = 1), color="red")
  ggsave(relativeR_terminal, filename = here::here(outfile, "relativeR_terminal.png"), width = 20)
  
  # Stock status SSB/SSBmsy
  statusSSB <- results %>%
    ggplot()+
    geom_line(aes(x=Year, y=status_SSB, group = seed, color = F_hist)) +
    geom_hline(aes(yintercept = 1), color = "red") +
    geom_hline(aes(yintercept = 0.5), color = "red", linetype = "dashed") # Check this threshold definition
  ggsave(statusSSB, filename = here::here(outfile, "status_SSB.png"), width = 20)
  
  # Stock status F/Fmsy
  statusF <- results %>%
    ggplot()+
    geom_line(aes(x=Year, y=status_F, group = seed, color = F_hist))
  ggsave(statusF, filename = here::here(outfile, "status_F.png"), width = 20)
  
  # Stock status Y/MSY
  statusY <- results %>%
    ggplot()+
    geom_line(aes(x=Year, y=status_Y, group = seed, color = F_hist))
  ggsave(statusY, filename = here::here(outfile, "status_Y.png"), width = 20)
  
  # Relative q index 1
  relQInd1 <- results %>% 
    ggplot()+
    geom_line(aes(x=Year, y = relq_index1, group = seed), alpha = 0.5) + # Could color lines by F history
    geom_hline(aes(yintercept = 1.0), color="red")
  ggsave(relQInd1, filename = here::here(outfile, "relq_index1.png"), width = 20)
  
  # Relative q index 2
  relQInd2 <- results %>% 
    ggplot()+
    geom_line(aes(x=Year, y = relq_index2, group = seed), alpha = 0.5) + # Could color lines by F history
    geom_hline(aes(yintercept = 1.0), color="red")
  ggsave(relQInd2, filename = here::here(outfile, "relq_index2.png"), width = 20)
  
  # Plot terminal year relative SSB vs. relative F
  SSB_F_terminal <- results %>% filter(Year == 40) %>%
    ggplot()+
    geom_point(aes(x=relSSB, y =relF))
  ggsave(SSB_F_terminal, filename = here::here(outfile, "SSB_F_terminal.png"), width = 20)
  
  # Terminal year kobe plot
  kobe_terminal <- results %>% filter(Year == 40) %>%
    ggplot() + 
    geom_hline(aes(yintercept = 1)) +
    geom_vline(aes(xintercept = 1)) +
    #geom_point(aes(x = status_SSB, y = status_F, color = EMshortName)) # By EM
    #geom_point(aes(x = status_SSB, y = status_F, color = OM_ecov_process_sig)) # OM ecov process sig
    #geom_point(aes(x = status_SSB, y = status_F, color = OM_ecov_process_cor)) # OM ecov process correlation
    #geom_point(aes(x = status_SSB, y = status_F, color = OM_ecov_process_obs_sig)) # OM ecov process obs sig
    #geom_point(aes(x = status_SSB, y = status_F, color = ageComp_sig)) # age comp sigma
    geom_point(aes(x = status_SSB, y = status_F, color = exp(log_index_sig))) # log index sigma
  ggsave(kobe_terminal, filename = here::here(outfile, "kobe_terminal.png"), width = 10)
    
  
}








