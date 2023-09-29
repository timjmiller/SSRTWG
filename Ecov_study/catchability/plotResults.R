
#' @param convergedONLY A boolean, when TRUE plots only show converged runs, default = TRUE.
#' @param outdir A string for the directory where a plot folder will be generated


results <- perfMet
convergedONLY = TRUE

plotResults <- function(results = NULL, convergedONLY = TRUE, outdir = here::here()){
  
  # Simulation summary
  results %>% count(OMshortName, EMshortName) %>% mutate(nsim = n/40) # /40 years so nsim = number of full simulations for each OM/EM
  # results %>% count(OMshortName) %>% mutate(nsim = n/40/4) # Full sims of all 4 correctly specified EMs
  # results %>% filter(OMshortName == 15) %>% select(EMshortName) %>% unique()
  
  setdiff(subsetOM$OMname, results$OMshortName) # Check what simulations still not run from subsetOM (in Catchability_ecov_sims.Rmd)
  
  if(convergedONLY == TRUE){
    results <- results %>% filter(EM_converged == TRUE)
  }
  
  # Pull out results that don't have time series
  singleResults <- results %>% group_by(sim) %>% dplyr::summarize(seed = unique(seed), F_hist = unique(F_hist), ageComp_sig = unique(ageComp_sig), log_catch_sig = unique(log_catch_sig), log_index_sig = unique(log_index_sig),
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
  singleResults %>% 
    ggplot() + 
    geom_point(aes(x=sim, y=EM_MohnsRho_SSB, color = OM_ecov_effect)) # Plot by sim
  
  singleResults %>%
    ggplot()+
    geom_boxplot(aes(x=EMshortName, y=EM_MohnsRho_SSB)) + facet_grid(. ~ OM_ecov_process_sig + OM_ecov_process_cor + OM_ecov_process_obs_sig + ageComp_sig + log_index_sig) + # Facet by OM settings
    geom_hline(aes(yintercept = 0), color="red")
    # ggtitle("OM_ecov_process_sig")
  
  # MohnsRho_F
  singleResults %>%
    ggplot()+
    geom_boxplot(aes(x=EMshortName, y=EM_MohnsRho_F)) + facet_grid(. ~ OM_ecov_process_sig + OM_ecov_process_cor + OM_ecov_process_obs_sig + ageComp_sig + log_index_sig) + # Facet by OM settings
    geom_hline(aes(yintercept = 0), color="red")
  
  # MohnsRho_R
  singleResults %>%
    ggplot()+
    geom_boxplot(aes(x=EMshortName, y=EM_MohnsRho_R)) + facet_grid(. ~ OM_ecov_process_sig + OM_ecov_process_cor + OM_ecov_process_obs_sig + ageComp_sig + log_index_sig) + # Facet by OM settings
    geom_hline(aes(yintercept = 0), color="red")
  
  # Relative SSBmsy
  singleResults %>%
    ggplot()+
    geom_boxplot(aes(x=EMshortName, y=relSSBMSY)) + facet_grid(. ~ OM_ecov_process_sig + OM_ecov_process_cor + OM_ecov_process_obs_sig + ageComp_sig + log_index_sig) + # Facet by OM settings
    geom_hline(aes(yintercept = 1), color="red")
  
  # Relative Fmsy
  singleResults %>%
    ggplot()+
    geom_boxplot(aes(x=EMshortName, y=relFMSY)) + facet_grid(. ~ OM_ecov_process_sig + OM_ecov_process_cor + OM_ecov_process_obs_sig + ageComp_sig + log_index_sig) + # Facet by OM settings
    geom_hline(aes(yintercept = 1), color="red")
  
  # Relative MSY
  singleResults %>%
    ggplot()+
    geom_boxplot(aes(x=EMshortName, y=relMSY)) + facet_grid(. ~ OM_ecov_process_sig + OM_ecov_process_cor + OM_ecov_process_obs_sig + ageComp_sig + log_index_sig) + # Facet by OM settings
    geom_hline(aes(yintercept = 1), color="red")
  
  # Relative Index 1 !!! may go back and fix so ages are 1-10 not 1, 10, 2-9
  singleResults %>% 
    pivot_longer(., cols = c(relselInd1_1, relselInd1_2, relselInd1_3, relselInd1_4, relselInd1_5, relselInd1_6, relselInd1_7, relselInd1_8, relselInd1_9, relselInd1_10), names_to = "age", values_to = "sel") %>%
    ggplot()+
    geom_boxplot(aes(x=EMshortName, y=sel, fill = age)) + facet_grid(. ~ OM_ecov_process_sig + OM_ecov_process_cor + OM_ecov_process_obs_sig + ageComp_sig + log_index_sig) + # Facet by OM settings
    geom_hline(aes(yintercept =1), color="red")
  
  # Relative Index 2 !!! may go back and fix so ages are 1-10 not 1, 10, 2-9
  singleResults %>% 
    pivot_longer(., cols = c(relselInd2_1, relselInd2_2, relselInd2_3, relselInd2_4, relselInd2_5, relselInd2_6, relselInd2_7, relselInd2_8, relselInd2_9, relselInd2_10), names_to = "age", values_to = "sel") %>%
    ggplot()+
    geom_boxplot(aes(x=EMshortName, y=sel, fill = age)) + facet_grid(. ~ OM_ecov_process_sig + OM_ecov_process_cor + OM_ecov_process_obs_sig + ageComp_sig + log_index_sig) + # Facet by OM settings
    geom_hline(aes(yintercept =1), color="red")
     # ageComp_sig impacts spread of results across all ages
    
  ### Time series plots ###
  #include status SSB/F/Y, relSSB, relF, relR, relq both indices, relFAA, relNAA, relCAA
  
  # Relative SSB
  #   # Time series
  # results %>% 
  #   ggplot()+
  #   geom_line(aes(x=Year, y = relSSB, group = seed), alpha = 0.5) + # Could color lines by F history
  #   geom_hline(aes(yintercept = 1.0), color="red")
    # Terminal year
  results %>% filter(Year == 40) %>%
    ggplot()+
    geom_boxplot(aes(x=EMshortName, y=relSSB)) + facet_grid(. ~ OM_ecov_process_sig + OM_ecov_process_cor + OM_ecov_process_obs_sig + ageComp_sig + log_index_sig) + # Facet by OM settings
    geom_hline(aes(yintercept = 1), color="red")
  
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
  results %>% filter(Year == 40) %>%
    ggplot()+
    geom_boxplot(aes(x=EMshortName, y=relF)) + facet_grid(. ~ OM_ecov_process_sig + OM_ecov_process_cor + OM_ecov_process_obs_sig + ageComp_sig + log_index_sig) + # Facet by OM settings
    geom_hline(aes(yintercept = 1), color="red")
  
  
  # Relative R
  #   # Time series
  # results %>% 
  #   ggplot()+
  #   geom_line(aes(x=Year, y = relR, group = seed, color = F_hist), alpha = 0.5) + # Could color lines by F history
  #   geom_hline(aes(yintercept = 1.0), color="red")
    # Terminal year
  results %>% filter(Year == 40) %>%
    ggplot()+
    geom_boxplot(aes(x=EMshortName, y=relR)) + facet_grid(. ~ OM_ecov_process_sig + OM_ecov_process_cor + OM_ecov_process_obs_sig + ageComp_sig + log_index_sig) + # Facet by OM settings
    geom_hline(aes(yintercept = 1), color="red")
  
  # Stock status SSB/SSBmsy
  results %>%
    ggplot()+
    geom_line(aes(x=Year, y=status_SSB, group = seed, color = F_hist)) +
    geom_hline(aes(yintercept = 1), color = "red") +
    geom_hline(aes(yintercept = 0.5), color = "red", linetype = "dashed") # Check this threshold definition
  
  # Stock status F/Fmsy
  results %>%
    ggplot()+
    geom_line(aes(x=Year, y=status_F, group = seed, color = F_hist))
  
  # Stock status Y/MSY
  results %>%
    ggplot()+
    geom_line(aes(x=Year, y=status_Y, group = seed, color = F_hist))
  
  # Relative q index 1
  results %>% 
    ggplot()+
    geom_line(aes(x=Year, y = relq_index1, group = seed), alpha = 0.5) + # Could color lines by F history
    geom_hline(aes(yintercept = 1.0), color="red")
  
  # Relative q index 2
  results %>% 
    ggplot()+
    geom_line(aes(x=Year, y = relq_index2, group = seed), alpha = 0.5) + # Could color lines by F history
    geom_hline(aes(yintercept = 1.0), color="red")
  
  # Plot terminal year relative SSB vs. relative F
  results %>% filter(Year == 40) %>%
    ggplot()+
    geom_point(aes(x=relSSB, y =relF))
  
  # Terminal year kobe plot
  results %>% filter(Year == 40) %>%
    ggplot() + 
    geom_hline(aes(yintercept = 1)) +
    geom_vline(aes(xintercept = 1)) +
    #geom_point(aes(x = status_SSB, y = status_F, color = EMshortName)) # By EM
    #geom_point(aes(x = status_SSB, y = status_F, color = OM_ecov_process_sig)) # OM ecov process sig
    #geom_point(aes(x = status_SSB, y = status_F, color = OM_ecov_process_cor)) # OM ecov process correlation
    #geom_point(aes(x = status_SSB, y = status_F, color = OM_ecov_process_obs_sig)) # OM ecov process obs sig
    #geom_point(aes(x = status_SSB, y = status_F, color = ageComp_sig)) # age comp sigma
    geom_point(aes(x = status_SSB, y = status_F, color = exp(log_index_sig))) # log index sigma
    
    
  
}








