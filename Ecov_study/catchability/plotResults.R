
#' @param convergedONLY A boolean, when TRUE plots only show converged runs, default = TRUE.
#' @param outfile A string for the directory where a plot folder will be generated


# results <- perfMet
# convergedONLY = TRUE

plotResults <- function(results = NULL, convergedONLY = TRUE, outfile = here::here()){
  
  # Simulation summary
  results %>% count(OMshortName, EMshortName) %>% mutate(nsim = n/40) %>% write.csv(., paste0(outfile,"/simSummary2.csv")) # /40 years so nsim = number of full simulations for each OM/EM
  
  # setdiff(subsetOM$OMname, results$OMshortName) # Check what simulations still not run from subsetOM (in Catchability_ecov_sims.Rmd)
  
  if(convergedONLY == TRUE){
    # Filter out converged runs
    results <- results %>% filter(EM_converged == TRUE)
    # Calculate number of converged runs
    convergeCount <- results %>% count(OMshortName, EMshortName) %>% mutate(nsim = n/40) # /40 years so nsim = number of full simulations for each OM/EM
    #convergeCount$OMshortName <- convergeCount$OMshortName %>% as.numeric()
    totalCount <- read.csv(paste0(outfile,"/simSummary2.csv")) # Created before results subset by convergence status
    totalCount$OMshortName <- totalCount$OMshortName %>% as.character()
    
    convergeRate <- full_join(convergeCount, totalCount, by = c("OMshortName", "EMshortName")) %>% mutate(convergeRate = nsim.x/nsim.y) %>% drop_columns(c("n.x", "nsim.x", "X", "n.y", "nsim.y"))
  }
  
  numericIndex <- which(colnames(results) %in% c("OMshortName", 'EMshortName', "EM_miss_season", "EM_miss_q", "F_hist", "Converged") == FALSE)
  results[,numericIndex] <- sapply(results[,numericIndex], as.numeric) # This introduces NAs by coercion since using NAs as placeholders
  
  # Pull out results that don't have time series
  # singleResults <- results %>% group_by(sim) %>% dplyr::summarize(seed = unique(seed), Fhist = unique(F_hist), ageComp_sig = unique(ageComp_sig), log_catch_sig = unique(log_catch_sig), log_index_sig = unique(log_index_sig),
  #                                                                                OM_ecov_effect = unique(OM_ecov_effect), OM_ecov_process_cor = unique(OM_ecov_process_cor), OM_ecov_process_obs_sig = unique(OM_ecov_process_obs_sig), OM_ecov_process_sig = unique(OM_ecov_process_sig), OMshortName = unique(OMshortName),
  #                                                                                                        EM_miss_q = unique(EM_miss_q), EM_miss_season = unique(EM_miss_season), EMshortName = unique(EMshortName),
  #                                                                                EM_MohnsRho_SSB = unique(EM_MohnsRho_SSB), EM_MohnsRho_F = unique(EM_MohnsRho_F), EM_MohnsRho_R = unique(EM_MohnsRho_R), EM_converged = unique(EM_converged), 
  #                                                                                relSSBMSY = unique(relSSBMSY), relFMSY = unique(relFMSY), relMSY = unique(relMSY), 
  #                                                                                
  #                                                                                relselInd1_1 = unique(relselInd1_1),
  #                                                                                relselInd1_2 = unique(relselInd1_2),
  #                                                                                relselInd1_3 = unique(relselInd1_3),
  #                                                                                relselInd1_4 = unique(relselInd1_4),
  #                                                                                relselInd1_5 = unique(relselInd1_5),
  #                                                                                relselInd1_6 = unique(relselInd1_6),
  #                                                                                relselInd1_7 = unique(relselInd1_7),
  #                                                                                relselInd1_8 = unique(relselInd1_8),
  #                                                                                relselInd1_9 = unique(relselInd1_9),
  #                                                                                relselInd1_10 = unique(relselInd1_10),
  #                                                                                
  #                                                                                relselInd2_1 = unique(relselInd2_1),
  #                                                                                relselInd2_2 = unique(relselInd2_2),
  #                                                                                relselInd2_3 = unique(relselInd2_3),
  #                                                                                relselInd2_4 = unique(relselInd2_4),
  #                                                                                relselInd2_5 = unique(relselInd2_5),
  #                                                                                relselInd2_6 = unique(relselInd2_6),
  #                                                                                relselInd2_7 = unique(relselInd2_7),
  #                                                                                relselInd2_8 = unique(relselInd2_8),
  #                                                                                relselInd2_9 = unique(relselInd2_9),
  #                                                                                relselInd2_10 = unique(relselInd2_10),
  #                                                                                
  #                                                                 relselCat_1 = unique(relselCat_1),
  #                                                                 relselCat_2 = unique(relselCat_2),
  #                                                                 relselCat_3 = unique(relselCat_3),
  #                                                                 relselCat_4 = unique(relselCat_4),
  #                                                                 relselCat_5 = unique(relselCat_5),
  #                                                                 relselCat_6 = unique(relselCat_6),
  #                                                                 relselCat_7 = unique(relselCat_7),
  #                                                                 relselCat_8 = unique(relselCat_8),
  #                                                                 relselCat_9 = unique(relselCat_9),
  #                                                                 relselCat_10 = unique(relselCat_10),
  #                                                                 
  #                                                                 med_relNAA_1 = median(relNAA_1), # Median EM:OM ratio for NAA in ages 1-10+
  #                                                                 med_relNAA_2 = median(relNAA_2),
  #                                                                 med_relNAA_3 = median(relNAA_3),
  #                                                                 med_relNAA_4 = median(relNAA_4),
  #                                                                 med_relNAA_5 = median(relNAA_5),
  #                                                                 med_relNAA_6 = median(relNAA_6),
  #                                                                 med_relNAA_7 = median(relNAA_7),
  #                                                                 med_relNAA_8 = median(relNAA_8),
  #                                                                 med_relNAA_9 = median(relNAA_9),
  #                                                                 med_relNAA_10 = median(relNAA_10),
  #                                                                 
  #                                                                 med_relFAA_1 = median(relFAA_1), # Median EM:OM ratio for FAA in ages 1-10+
  #                                                                 med_relFAA_2 = median(relFAA_2),
  #                                                                 med_relFAA_3 = median(relFAA_3),
  #                                                                 med_relFAA_4 = median(relFAA_4),
  #                                                                 med_relFAA_5 = median(relFAA_5),
  #                                                                 med_relFAA_6 = median(relFAA_6),
  #                                                                 med_relFAA_7 = median(relFAA_7),
  #                                                                 med_relFAA_8 = median(relFAA_8),
  #                                                                 med_relFAA_9 = median(relFAA_9),
  #                                                                 med_relFAA_10 = median(relFAA_10),
  #                                                                 
  #                                                                 med_relCAA_1 = median(relCAA_1), # Median EM:OM ratio for CAA in ages 1-10+
  #                                                                 med_relCAA_2 = median(relCAA_2),
  #                                                                 med_relCAA_3 = median(relCAA_3),
  #                                                                 med_relCAA_4 = median(relCAA_4),
  #                                                                 med_relCAA_5 = median(relCAA_5),
  #                                                                 med_relCAA_6 = median(relCAA_6),
  #                                                                 med_relCAA_7 = median(relCAA_7),
  #                                                                 med_relCAA_8 = median(relCAA_8),
  #                                                                 med_relCAA_9 = median(relCAA_9),
  #                                                                 med_relCAA_10 = median(relCAA_10),
  #                                                                 EM_ecovBeta_ind1 = unique(EM_ecovBeta_ind1), # EM beta parameter for index 1
  #                                                                 EM_ecovBeta_ind2 = unique(EM_ecovBeta_ind2), # EM beta parameter for index 2
  #                                                                 med_relCatch = median(EM_Catch/OM_Catch)) # Median EM:OM ratio for catch
  #   
  if(convergedONLY == TRUE){
    singleResults <- results %>% group_by(sim) %>%
      dplyr::summarize(seed = unique(seed), Fhist = unique(F_hist), ageComp_sig = unique(ageComp_sig), log_catch_sig = unique(log_catch_sig), log_index_sig = unique(log_index_sig),
                       OM_ecov_effect = unique(OM_ecov_effect), OM_ecov_process_cor = unique(OM_ecov_process_cor), OM_ecov_process_obs_sig = unique(OM_ecov_process_obs_sig), OM_ecov_process_sig = unique(OM_ecov_process_sig), OMshortName = unique(OMshortName),
                       EM_miss_q = unique(EM_miss_q), EM_miss_season = unique(EM_miss_season), EMshortName = unique(EMshortName),
                       EM_MohnsRho_SSB = unique(EM_MohnsRho_SSB), EM_MohnsRho_F = unique(EM_MohnsRho_F), EM_MohnsRho_R = unique(EM_MohnsRho_R), EM_converged = unique(EM_converged),
                       relFMSY = unique(relFMSY),
                       relMSY = unique(relMSY),
                       relSSBMSY = unique(relSSBMSY))
      
    singleResults_converge <- full_join(convergeRate, singleResults, by = c("OMshortName", "EMshortName", "ageComp_sig", "log_index_sig", "log_catch_sig")) %>% # remember that some runs only had NoEcov and qRand since the Ecov not used so Ecov and qRandEcov have fewer samples in these plots
    group_by(OMshortName, EMshortName) %>% #group_by(EMshortName) %>%
      dplyr::summarize(OMshortName = unique(OMshortName), 
                       EMshortName = unique(EMshortName),
                       convergeRate = unique(convergeRate),
                       OM_ecov_process_sig = unique(OM_ecov_process_sig),
                       OM_ecov_process_cor = unique(OM_ecov_process_cor),
                       OM_ecov_process_obs_sig = unique(OM_ecov_process_obs_sig),
                       OM_ecov_effect = unique(OM_ecov_effect),
                       ageComp_sig = unique(ageComp_sig),
                       log_index_sig = unique(log_index_sig),
                       Fhist = unique(Fhist),
                       relFMSY = unique(relFMSY),
                       relMSY = unique(relMSY),
                       relSSBMSY = unique(relSSBMSY))
    # Plot convergence rate
    convergePlot <- singleResults_converge %>% 
      ggplot()+
      geom_boxplot(aes(x=EMshortName, y=convergeRate)) + facet_grid(. ~ OM_ecov_process_sig + OM_ecov_process_cor + OM_ecov_process_obs_sig + ageComp_sig + log_index_sig + OM_ecov_effect)  # Facet by OM settings
    ggsave(convergePlot, filename = here::here(outfile, "plots/convergeRate_boxplot.png"), width = 20) # Formerly saved as convergeRate.png
    # Label order: Ecov, NoEcov, qRand, qRandEcov
    
    ## Extra plots looking at results across F history
    Fhist <- singleResults_converge %>%
      ggplot() +
      geom_boxplot(aes(x=EMshortName, y = convergeRate)) +
      facet_grid(. ~ Fhist)
    ggsave(Fhist, filename = here::here(outfile, paste0("plots/Fhist_summary.png")), width = 10)
    
    Fhist_relFMSY <- singleResults_converge %>%
      ggplot() +
      geom_boxplot(aes(x=EMshortName, y = relFMSY)) +
      facet_grid(. ~ Fhist)
    ggsave(Fhist_relFMSY, filename = here::here(outfile, paste0("plots/Fhist_relFMSY_summary.png")), width = 10)
    
    Fhist_relMSY <- singleResults_converge %>%
      ggplot() +
      geom_boxplot(aes(x=EMshortName, y = relMSY)) +
      facet_grid(. ~ Fhist)
    ggsave(Fhist_relMSY, filename = here::here(outfile, paste0("plots/Fhist_relMSY_summary.png")), width = 10)
    
    Fhist_relSSBMSY <- singleResults_converge %>%
      ggplot() +
      geom_boxplot(aes(x=EMshortName, y = relSSBMSY)) +
      facet_grid(. ~ Fhist)
    ggsave(Fhist_relSSBMSY, filename = here::here(outfile, paste0("plots/Fhist_relSSBMSY_summary.png")), width = 10)
  }
                                                                                 
  # Time series include status SSB/F/Y, relSSB, relF, relR, relq both indices, relFAA, relNAA, relCAA
  
  # MohnsRho_SSB
  mohnsRho_results <- results %>% group_by(sim) %>%
    dplyr::summarize(seed = unique(seed), Fhist = unique(F_hist), ageComp_sig = unique(ageComp_sig), log_catch_sig = unique(log_catch_sig), log_index_sig = unique(log_index_sig),
                     OM_ecov_effect = unique(OM_ecov_effect), OM_ecov_process_cor = unique(OM_ecov_process_cor), OM_ecov_process_obs_sig = unique(OM_ecov_process_obs_sig), OM_ecov_process_sig = unique(OM_ecov_process_sig), OMshortName = unique(OMshortName),
                     EM_miss_q = unique(EM_miss_q), EM_miss_season = unique(EM_miss_season), EMshortName = unique(EMshortName),
                     EM_MohnsRho_SSB = unique(EM_MohnsRho_SSB), EM_MohnsRho_F = unique(EM_MohnsRho_F), EM_MohnsRho_R = unique(EM_MohnsRho_R), EM_converged = unique(EM_converged)) 
  mohnsRho_SSB_series <- mohnsRho_results %>%
    ggplot() + 
    geom_point(aes(x=sim, y=EM_MohnsRho_SSB, color = OM_ecov_effect)) +  # Plot by sim
    ggtitle("EM Mohn's Rho Time Series: SSB")
  ggsave(mohnsRho_SSB_series, filename = here::here(outfile, "plots/mohnsRho_SSB_series.png"), width = 10)
  
  mohnsRho_SSB <- mohnsRho_results %>%
    ggplot()+
    geom_boxplot(aes(x=EMshortName, y=EM_MohnsRho_SSB, fill = EMshortName)) +
    facet_grid(cols = vars(OM_ecov_process_sig, OM_ecov_process_cor, OM_ecov_process_obs_sig, ageComp_sig, log_index_sig), rows = vars(OM_ecov_effect, Fhist)) + # Facet by OM settings
    geom_hline(aes(yintercept = 0), color="red") +
    scale_fill_grey(start = 0.45, end = 1.0) +
    theme(axis.text.x = element_blank()) +
    labs(fill = "EM",
         x = "EM",
         y = "Mohn's Rho: SSB") +
    ggtitle("EM Mohn's Rho: SSB")
  ggsave(mohnsRho_SSB, filename = here::here(outfile, "plots/mohnsRho_SSB_boxplot.png"), width = 20) # Formerly saved as mohnsRho_SSB.png
  
  # MohnsRho_F
  mohnsRho_F <-  mohnsRho_results %>%
    ggplot()+
    geom_boxplot(aes(x=EMshortName, y=EM_MohnsRho_F, fill = EMshortName)) + 
    facet_grid(cols = vars(OM_ecov_process_sig, OM_ecov_process_cor, OM_ecov_process_obs_sig, ageComp_sig, log_index_sig), rows = vars(OM_ecov_effect, Fhist)) + # Facet by OM settings
    geom_hline(aes(yintercept = 0), color="red") +
    scale_fill_grey(start = 0.45, end = 1.0) +
    theme(axis.text.x = element_blank()) +
    labs(fill = "EM",
         x = "EM",
         y = "Mohn's Rho: F") +
    ggtitle("EM Mohn's Rho: F")
  ggsave(mohnsRho_F, filename = here::here(outfile, "plots/mohnsRho_F_boxplot.png"), width = 20) # Formerly saved as mohnsRho_F.png
  
  # MohnsRho_R
  mohnsRho_R <- mohnsRho_results %>%
    ggplot()+
    geom_boxplot(aes(x=EMshortName, y=EM_MohnsRho_R, fill = EMshortName)) +
    facet_grid(cols = vars(OM_ecov_process_sig, OM_ecov_process_cor, OM_ecov_process_obs_sig, ageComp_sig, log_index_sig), rows = vars(OM_ecov_effect, Fhist)) + # Facet by OM settings
    geom_hline(aes(yintercept = 0), color="red") +
    scale_fill_grey(start = 0.45, end = 1.0) +
    theme(axis.text.x = element_blank()) +
    labs(fill = "EM", 
         x = "EM",
         y = "Mohn's Rho: R") +
    ggtitle("EM Mohn's Rho: R")
  ggsave(mohnsRho_R, filename = here::here(outfile, "plots/mohnsRho_R_boxplot.png"), width = 20) # Formerly saved as mohnsRho_R.png
  
  # Relative SSBmsy
  refpts_results <- results %>% group_by(sim) %>%
    dplyr::summarize(seed = unique(seed), Fhist = unique(F_hist), ageComp_sig = unique(ageComp_sig), log_catch_sig = unique(log_catch_sig), log_index_sig = unique(log_index_sig),
                     OM_ecov_effect = unique(OM_ecov_effect), OM_ecov_process_cor = unique(OM_ecov_process_cor), OM_ecov_process_obs_sig = unique(OM_ecov_process_obs_sig), OM_ecov_process_sig = unique(OM_ecov_process_sig), OMshortName = unique(OMshortName),
                     EM_miss_q = unique(EM_miss_q), EM_miss_season = unique(EM_miss_season), EMshortName = unique(EMshortName),
                     relSSBMSY = unique(relSSBMSY), relFMSY = unique(relFMSY), relMSY = unique(relMSY))
  relSSBMSY <- refpts_results %>%
    ggplot()+
    geom_boxplot(aes(x=EMshortName, y=relSSBMSY, fill = EMshortName)) + 
    facet_grid(cols = vars(OM_ecov_process_sig, OM_ecov_process_cor, OM_ecov_process_obs_sig, ageComp_sig, log_index_sig), rows = vars(OM_ecov_effect, Fhist)) + # Facet by OM settings
    geom_hline(aes(yintercept = 1), color="red") +
    scale_fill_grey(start = 0.45, end = 1.0) +
    theme(axis.text.x = element_blank()) +
    labs(fill = "EM", 
         x = "EM",
         y = "EM:OM SSBmsy ratio") +
    ggtitle("EM:OM ratio of SSBmsy")
  ggsave(relSSBMSY, filename = here::here(outfile, "plots/relSSBMSY_boxplot.png"), width = 20) # Formerly saved as relSSBMSY.png
  
  # Relative Fmsy
  relFMSY <- refpts_results %>%
    ggplot()+
    geom_boxplot(aes(x=EMshortName, y=relFMSY, fill = EMshortName)) + 
    facet_grid(cols = vars(OM_ecov_process_sig, OM_ecov_process_cor, OM_ecov_process_obs_sig, ageComp_sig, log_index_sig), rows = vars(OM_ecov_effect, Fhist)) + # Facet by OM settings
    geom_hline(aes(yintercept = 1), color="red") +
    scale_fill_grey(start = 0.45, end = 1.0) +
    theme(axis.text.x = element_blank()) +
    labs(fill = "EM", 
         x = "EM",
         y = "EM:OM FMSY ratio") +
    ggtitle("EM:OM ratio of FMSY")
  ggsave(relFMSY, filename = here::here(outfile, "plots/relFMSY_boxplot.png"), width = 20) # Formerly saved as relFMSY.png
  
  # Relative MSY
  relMSY <- refpts_results %>%
    ggplot()+
    geom_boxplot(aes(x=EMshortName, y=relMSY, fill = EMshortName)) + 
    facet_grid(cols = vars(OM_ecov_process_sig, OM_ecov_process_cor, OM_ecov_process_obs_sig, ageComp_sig, log_index_sig), rows = vars(OM_ecov_effect, Fhist)) + # Facet by OM settings
    geom_hline(aes(yintercept = 1), color="red") +
    scale_fill_grey(start = 0.45, end = 1.0) +
    theme(axis.text.x = element_blank()) +
    labs(fill = "EM", 
         x = "EM",
         y = "EM:OM ratio of MSY") + 
    ggtitle("EM:OM ratio of MSY")
  ggsave(relMSY, filename = here::here(outfile, "plots/relMSY_boxplot.png"), width = 20) # Formerly saved as relMSY.png
  
  # Relative Index 1 
  selectivity_results <- results %>% group_by(sim) %>%
    dplyr::summarize(seed = unique(seed), Fhist = unique(F_hist), ageComp_sig = unique(ageComp_sig), log_catch_sig = unique(log_catch_sig), log_index_sig = unique(log_index_sig),
                     OM_ecov_effect = unique(OM_ecov_effect), OM_ecov_process_cor = unique(OM_ecov_process_cor), OM_ecov_process_obs_sig = unique(OM_ecov_process_obs_sig), OM_ecov_process_sig = unique(OM_ecov_process_sig), OMshortName = unique(OMshortName),
                     EM_miss_q = unique(EM_miss_q), EM_miss_season = unique(EM_miss_season), EMshortName = unique(EMshortName),
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
                     relselInd2_10 = unique(relselInd2_10),
                     
                     relselCat_1 = unique(relselCat_1),
                     relselCat_2 = unique(relselCat_2),
                     relselCat_3 = unique(relselCat_3),
                     relselCat_4 = unique(relselCat_4),
                     relselCat_5 = unique(relselCat_5),
                     relselCat_6 = unique(relselCat_6),
                     relselCat_7 = unique(relselCat_7),
                     relselCat_8 = unique(relselCat_8),
                     relselCat_9 = unique(relselCat_9),
                     relselCat_10 = unique(relselCat_10))
                     
  relSelInd1 <- selectivity_results %>% 
    drop_columns(c("relselInd2_1", "relselInd2_2", "relselInd2_3", "relselInd2_4", "relselInd2_5", "relselInd2_6", "relselInd2_7", "relselInd2_8", "relselInd2_9", "relselInd2_10")) %>%
    drop_columns(c("relselCat_1", "relselCat_2", "relselCat_3", "relselCat_4", "relselCat_5", "relselCat_6", "relselCat_7", "relselCat_8", "relselCat_9", "relselCat_10")) %>%
    pivot_longer(., cols = c(relselInd1_1, relselInd1_2, relselInd1_3, relselInd1_4, relselInd1_5, relselInd1_6, relselInd1_7, relselInd1_8, relselInd1_9, relselInd1_10), names_to = "age", names_prefix = "relselInd1_", values_to = "sel")  
  relSelInd1$age <- factor(relSelInd1$age, levels = c(1:10))
  relSelInd1plot <- relSelInd1 %>% 
    ggplot()+
    geom_boxplot(aes(x=EMshortName, y=sel, fill = EMshortName)) + 
    facet_grid(cols = vars(OM_ecov_process_sig, OM_ecov_process_cor, OM_ecov_process_obs_sig, ageComp_sig, log_index_sig), rows = vars(OM_ecov_effect, Fhist)) + # Facet by OM settings
    geom_hline(aes(yintercept = 1), color="red") +
    scale_fill_grey(start = 0.45, end = 1.0) +
    theme(axis.text.x = element_blank()) +
    labs(fill = "EM", 
         x = "EM",
         y = "Selectivity") +
    ggtitle("EM:OM ratio of spring index selectivity")
  ggsave(relSelInd1plot, filename = here::here(outfile, "plots/relSelInd1_boxplot.png"), width = 20) # Formerly saved as relSelInd1.png
  
  relSelInd1 <- selectivity_results %>% 
    drop_columns(c("relselInd2_1", "relselInd2_2", "relselInd2_3", "relselInd2_4", "relselInd2_5", "relselInd2_6", "relselInd2_7", "relselInd2_8", "relselInd2_9", "relselInd2_10")) %>%
    drop_columns(c("relselCat_1", "relselCat_2", "relselCat_3", "relselCat_4", "relselCat_5", "relselCat_6", "relselCat_7", "relselCat_8", "relselCat_9", "relselCat_10")) %>%
    pivot_longer(., cols = c(relselInd1_1, relselInd1_2, relselInd1_3, relselInd1_4, relselInd1_5, relselInd1_6, relselInd1_7, relselInd1_8, relselInd1_9, relselInd1_10), names_to = "age", names_prefix = "relselInd1_", values_to = "sel")  
  relSelInd1$age <- factor(relSelInd1$age, levels = c(1:10))
  for(iage in 1:10){
    relSelInd1plot <- relSelInd1 %>% filter(age == iage) %>%
      ggplot()+
      geom_boxplot(aes(x=EMshortName, y=sel, fill = EMshortName)) + 
      facet_grid(cols = vars(OM_ecov_process_sig, OM_ecov_process_cor, OM_ecov_process_obs_sig, ageComp_sig, log_index_sig), rows = vars(OM_ecov_effect, Fhist)) + # Facet by OM settings
      geom_hline(aes(yintercept = 1), color="red") +
      scale_fill_grey(start = 0.45, end = 1.0) +
      theme(axis.text.x = element_blank()) +
      labs(fill = "EM", 
           x = "EM",
           y = "EM:OM ratio of selectivity") +
      ggtitle(paste0("EM:OM ratio of spring index selectivity ", iage)) +
      ylim(0, 2.7)
    ggsave(relSelInd1plot, filename = here::here(outfile, paste0("plots/relSelInd1_", iage, "_boxplot.png")), width = 20) # Formerly saved as relSelInd1_'iage'.png
  }
  
  # Relative Index 2 
  relSelInd2 <- selectivity_results %>% 
    drop_columns(c("relselInd1_1", "relselInd1_2", "relselInd1_3", "relselInd1_4", "relselInd1_5", "relselInd1_6", "relselInd1_7", "relselInd1_8", "relselInd1_9", "relselInd1_10")) %>%
    drop_columns(c("relselCat_1", "relselCat_2", "relselCat_3", "relselCat_4", "relselCat_5", "relselCat_6", "relselCat_7", "relselCat_8", "relselCat_9", "relselCat_10")) %>%
    pivot_longer(., cols = c(relselInd2_1, relselInd2_2, relselInd2_3, relselInd2_4, relselInd2_5, relselInd2_6, relselInd2_7, relselInd2_8, relselInd2_9, relselInd2_10), names_to = "age", names_prefix = "relselInd2_", values_to = "sel") 
  relSelInd2$age <- factor(relSelInd2$age, levels = c(1:10))
  relSelInd2plot <- relSelInd2 %>%
    ggplot()+
    geom_boxplot(aes(x=EMshortName, y=sel, fill = EMshortName)) + 
    facet_grid(cols = vars(OM_ecov_process_sig, OM_ecov_process_cor, OM_ecov_process_obs_sig, ageComp_sig, log_index_sig), rows = vars(OM_ecov_effect, Fhist)) + # Facet by OM settings
    geom_hline(aes(yintercept = 1), color="red") +
    scale_fill_grey(start = 0.45, end = 1.0) +
    theme(axis.text.x = element_blank()) +
    labs(fill = "EM", 
         x = "EM",
         y = "EM:OM ratio of selectivity") +
    ggtitle("EM:OM ratio of fall index selectivity")
     # ageComp_sig impacts spread of results across all ages
  ggsave(relSelInd2plot, filename = here::here(outfile, "plots/relSelInd2_boxplot.png"), width = 20) # Formerly saved as relSelInd2.png
  
  relSelInd2 <- selectivity_results %>% 
    drop_columns(c("relselInd1_1", "relselInd1_2", "relselInd1_3", "relselInd1_4", "relselInd1_5", "relselInd1_6", "relselInd1_7", "relselInd1_8", "relselInd1_9", "relselInd1_10")) %>%
    drop_columns(c("relselCat_1", "relselCat_2", "relselCat_3", "relselCat_4", "relselCat_5", "relselCat_6", "relselCat_7", "relselCat_8", "relselCat_9", "relselCat_10")) %>%
    pivot_longer(., cols = c(relselInd2_1, relselInd2_2, relselInd2_3, relselInd2_4, relselInd2_5, relselInd2_6, relselInd2_7, relselInd2_8, relselInd2_9, relselInd2_10), names_to = "age", names_prefix = "relselInd2_", values_to = "sel") 
  relSelInd2$age <- factor(relSelInd2$age, levels = c(1:10))
  for(iage in 1:10){
    relSelInd2plot <- relSelInd2 %>% filter(age == iage) %>%
      ggplot()+
      geom_boxplot(aes(x=EMshortName, y=sel, fill = EMshortName)) + 
      facet_grid(cols = vars(OM_ecov_process_sig, OM_ecov_process_cor, OM_ecov_process_obs_sig, ageComp_sig, log_index_sig), rows = vars(OM_ecov_effect, Fhist)) + # Facet by OM settings
      geom_hline(aes(yintercept = 1), color="red") +
      scale_fill_grey(start = 0.45, end = 1.0) +
      theme(axis.text.x = element_blank()) +
      labs(fill = "EM", 
           x = "EM",
           y = "EM:OM ratio of selectivity") +
      ggtitle(paste0("EM:OM ratio of fall index selectivity ", iage)) +
      ylim(0, 3)
    # ageComp_sig impacts spread of results across all ages
    ggsave(relSelInd2plot, filename = here::here(outfile, paste0("plots/relSelInd2_", iage, "_boxplot.png")), width = 20) # Formerly saved as relSelInd2_'iage'.png
  }
  
  # Relative selectivity fleet
  relSelFleet <- selectivity_results %>% 
    drop_columns(c("relselInd1_1", "relselInd1_2", "relselInd1_3", "relselInd1_4", "relselInd1_5", "relselInd1_6", "relselInd1_7", "relselInd1_8", "relselInd1_9", "relselInd1_10")) %>%
    drop_columns(c("relselInd2_1", "relselInd2_2", "relselInd2_3", "relselInd2_4", "relselInd2_5", "relselInd2_6", "relselInd2_7", "relselInd2_8", "relselInd2_9", "relselInd2_10")) %>%
    pivot_longer(., cols = c(relselCat_1, relselCat_2, relselCat_3, relselCat_4, relselCat_5, relselCat_6, relselCat_7, relselCat_8, relselCat_9, relselCat_10), names_to = "age", names_prefix = "relselCat_", values_to = "sel")
  relSelFleet$age <- factor(relSelFleet$age, levels = c(1:10))
  for(iage in 1:10){
    relSelFleetplot <- relSelFleet %>% filter(age == iage) %>%
      ggplot() +
      geom_boxplot(aes(x = EMshortName, y = sel, fill = EMshortName)) + 
      facet_grid(cols = vars(OM_ecov_process_sig, OM_ecov_process_cor, OM_ecov_process_obs_sig, ageComp_sig, log_index_sig), rows = vars(OM_ecov_effect, Fhist)) + # Facet by OM settings
      geom_hline(aes(yintercept = 1), color="red") +
      scale_fill_grey(start = 0.45, end = 1.0) +
      theme(axis.text.x = element_blank()) +
      labs(fill = "EM", 
           x = "EM",
           y = "EM:OM ratio of selectivity") +
      ggtitle(paste0("Fleet relative selectivity ", iage)) +
      ylim(0,2.7) +
      ggtitle(paste0("EM:OM ratio of fleet selectivity at age ", iage))
    ggsave(relSelFleetplot, filename = here::here(outfile, paste0("plots/relSelFleet_", iage, "_boxplot.png")), width = 20) # Formerly saved as relSelFleet_'iage'.png
  }
  
  # Median (across years) relative FAA
  ageEcov_results <- results %>% group_by(sim) %>% dplyr::summarize(seed = unique(seed), Fhist = unique(F_hist), ageComp_sig = unique(ageComp_sig), log_catch_sig = unique(log_catch_sig), log_index_sig = unique(log_index_sig),
                                                                    OM_ecov_effect = unique(OM_ecov_effect), OM_ecov_process_cor = unique(OM_ecov_process_cor), OM_ecov_process_obs_sig = unique(OM_ecov_process_obs_sig), OM_ecov_process_sig = unique(OM_ecov_process_sig), OMshortName = unique(OMshortName),
                                                                    EM_miss_q = unique(EM_miss_q), EM_miss_season = unique(EM_miss_season), EMshortName = unique(EMshortName),
                                                                    med_relNAA_1 = median(relNAA_1), # Median EM:OM ratio for NAA in ages 1-10+
                                                                    med_relNAA_2 = median(relNAA_2),
                                                                    med_relNAA_3 = median(relNAA_3),
                                                                    med_relNAA_4 = median(relNAA_4),
                                                                    med_relNAA_5 = median(relNAA_5),
                                                                    med_relNAA_6 = median(relNAA_6),
                                                                    med_relNAA_7 = median(relNAA_7),
                                                                    med_relNAA_8 = median(relNAA_8),
                                                                    med_relNAA_9 = median(relNAA_9),
                                                                    med_relNAA_10 = median(relNAA_10),
                                                                    
                                                                    med_relFAA_1 = median(relFAA_1), # Median EM:OM ratio for FAA in ages 1-10+
                                                                    med_relFAA_2 = median(relFAA_2),
                                                                    med_relFAA_3 = median(relFAA_3),
                                                                    med_relFAA_4 = median(relFAA_4),
                                                                    med_relFAA_5 = median(relFAA_5),
                                                                    med_relFAA_6 = median(relFAA_6),
                                                                    med_relFAA_7 = median(relFAA_7),
                                                                    med_relFAA_8 = median(relFAA_8),
                                                                    med_relFAA_9 = median(relFAA_9),
                                                                    med_relFAA_10 = median(relFAA_10),
                                                                    
                                                                    med_relCAA_1 = median(relCAA_1), # Median EM:OM ratio for CAA in ages 1-10+
                                                                    med_relCAA_2 = median(relCAA_2),
                                                                    med_relCAA_3 = median(relCAA_3),
                                                                    med_relCAA_4 = median(relCAA_4),
                                                                    med_relCAA_5 = median(relCAA_5),
                                                                    med_relCAA_6 = median(relCAA_6),
                                                                    med_relCAA_7 = median(relCAA_7),
                                                                    med_relCAA_8 = median(relCAA_8),
                                                                    med_relCAA_9 = median(relCAA_9),
                                                                    med_relCAA_10 = median(relCAA_10),
                                                                    EM_ecovBeta_ind1 = unique(EM_ecovBeta_ind1), # EM beta parameter for index 1
                                                                    EM_ecovBeta_ind2 = unique(EM_ecovBeta_ind2), # EM beta parameter for index 2
                                                                    med_relCatch = median(EM_Catch/OM_Catch)) # Median EM:OM ratio for catch
  
  med_relFAA <- ageEcov_results %>%
    drop_columns(c("med_relNAA_1", "med_relNAA_2", "med_relNAA_3", "med_relNAA_4", "med_relNAA_5", "med_relNAA_6", "med_relNAA_7", "med_relNAA_8", "med_relNAA_9", "med_relNAA_10")) %>%
    drop_columns(c("med_relCAA_1", "med_relCAA_2", "med_relCAA_3", "med_relCAA_4", "med_relCAA_5", "med_relCAA_6", "med_relCAA_7", "med_relCAA_8", "med_relCAA_9", "med_relCAA_10")) %>%
    pivot_longer(., cols = c(med_relFAA_1, med_relFAA_2, med_relFAA_3, med_relFAA_4, med_relFAA_5, med_relFAA_6, med_relFAA_7, med_relFAA_8, med_relFAA_9, med_relFAA_10), names_to = "age", names_prefix = "med_relFAA_", values_to = "relFAA")
  med_relFAA$age <- factor(med_relFAA$age, levels = c(1:10))
  for(iage in 1:10){
    med_relFAAplot <- med_relFAA %>% filter(age == iage) %>%
      ggplot() +
      geom_boxplot(aes(x = EMshortName, y = relFAA, fill = EMshortName)) + 
      facet_grid(cols = vars(OM_ecov_process_sig, OM_ecov_process_cor, OM_ecov_process_obs_sig, ageComp_sig, log_index_sig), rows = vars(OM_ecov_effect, Fhist)) + # Facet by OM settings
      geom_hline(aes(yintercept = 1), color="red") +
      scale_fill_grey(start = 0.45, end = 1.0) +
      theme(axis.text.x = element_blank()) +
      labs(fill = "EM", 
           x = "EM",
           y = "EM:OM ratio of FAA") +
      ggtitle(paste0("Median EM:OM ratio of F at age ", iage)) +
      ylim(-5,5)
    ggsave(med_relFAAplot, filename = here::here(outfile, paste0("plots/relFAA_", iage, "_boxplot.png")), width = 20) # Formerly saved as relFAA_'iage'.png
  }
  
  # Median (across years) relative CAA
  med_relCAA <- ageEcov_results %>%
    drop_columns(c("med_relNAA_1", "med_relNAA_2", "med_relNAA_3", "med_relNAA_4", "med_relNAA_5", "med_relNAA_6", "med_relNAA_7", "med_relNAA_8", "med_relNAA_9", "med_relNAA_10")) %>%
    drop_columns(c("med_relFAA_1", "med_relFAA_2", "med_relFAA_3", "med_relFAA_4", "med_relFAA_5", "med_relFAA_6", "med_relFAA_7", "med_relFAA_8", "med_relFAA_9", "med_relFAA_10")) %>%
    pivot_longer(., cols = c(med_relCAA_1, med_relCAA_2, med_relCAA_3, med_relCAA_4, med_relCAA_5, med_relCAA_6, med_relCAA_7, med_relCAA_8, med_relCAA_9, med_relCAA_10), names_to = "age", names_prefix = "med_relCAA_", values_to = "relCAA")
  med_relCAA$age <- factor(med_relCAA$age, levels = c(1:10))
  for(iage in 1:10){
    med_relCAAplot <- med_relCAA %>% filter(age == iage) %>%
      ggplot() +
      geom_boxplot(aes(x = EMshortName, y = relCAA, fill = EMshortName)) + 
      facet_grid(cols = vars(OM_ecov_process_sig, OM_ecov_process_cor, OM_ecov_process_obs_sig, ageComp_sig, log_index_sig), rows = vars(OM_ecov_effect, Fhist)) + # Facet by OM settings
      geom_hline(aes(yintercept = 1), color="red") +
      scale_fill_grey(start = 0.45, end = 1.0) +
      theme(axis.text.x = element_blank()) +
      labs(fill = "EM", 
           x = "EM",
           y = "EM:OM ratio of CAA") +
      ggtitle(paste0("Median EM:OM ratio of Catch at age ", iage)) +
      ylim(-5,5)
    ggsave(med_relCAAplot, filename = here::here(outfile, paste0("plots/relCAA_", iage, "_boxplot.png")), width = 20) # Formerly saved as relCAA_'iage'.png
  }
  
  # Median (across years) relative NAA
  med_relNAA <- ageEcov_results %>%
    drop_columns(c("med_relFAA_1", "med_relFAA_2", "med_relFAA_3", "med_relFAA_4", "med_relFAA_5", "med_relFAA_6", "med_relFAA_7", "med_relFAA_8", "med_relFAA_9", "med_relFAA_10")) %>%
    drop_columns(c("med_relCAA_1", "med_relCAA_2", "med_relCAA_3", "med_relCAA_4", "med_relCAA_5", "med_relCAA_6", "med_relCAA_7", "med_relCAA_8", "med_relCAA_9", "med_relCAA_10")) %>%
    pivot_longer(., cols = c(med_relNAA_1, med_relNAA_2, med_relNAA_3, med_relNAA_4, med_relNAA_5, med_relNAA_6, med_relNAA_7, med_relNAA_8, med_relNAA_9, med_relNAA_10), names_to = "age", names_prefix = "med_relNAA_", values_to = "relNAA")
  med_relNAA$age <- factor(med_relNAA$age, levels = c(1:10))
  for(iage in 1:10){
    med_relNAAplot <- med_relNAA %>% filter(age == iage) %>%
      ggplot() +
      geom_boxplot(aes(x = EMshortName, y = relNAA, fill = EMshortName)) + 
      facet_grid(cols = vars(OM_ecov_process_sig, OM_ecov_process_cor, OM_ecov_process_obs_sig, ageComp_sig, log_index_sig), rows = vars(OM_ecov_effect, Fhist)) + # Facet by OM settings
      geom_hline(aes(yintercept = 1), color="red") +
      scale_fill_grey(start = 0.45, end = 1.0) +
      theme(axis.text.x = element_blank()) +
      labs(fill = "EM", 
           x = "EM",
           y = "EM:OM ratio of NAA") +
      ggtitle(paste0("Median EM:OM ratio of N at age ", iage)) +
      ylim(-5,5)
    ggsave(med_relNAAplot, filename = here::here(outfile, paste0("plots/relNAA_", iage, "_boxplot.png")), width = 20) # Formerly saved as relNAA_'iage'.png
  }
  
  ## Extra plots small ecov process sig & high ecov process obs sig
  ecov_process_relSSBMSY <- refpts_results %>%
    filter(OM_ecov_process_sig == 0.1) %>%
    ggplot()+
    geom_boxplot(aes(x=EMshortName, y = relSSBMSY)) +
    facet_grid(. ~ OM_ecov_process_obs_sig) + 
    ggtitle("EM:OM ratio of SSBmsy when ecov process sigma = 0.1")
  ggsave(ecov_process_relSSBMSY, filename = here::here(outfile, paste0("plots/ecov_process_relSSBMSY.png")), width = 20)
  
  ecov_process_relFMSY <- refpts_results %>%
    filter(OM_ecov_process_sig == 0.1) %>%
    ggplot()+
    geom_boxplot(aes(x=EMshortName, y = relFMSY)) +
    facet_grid(. ~ OM_ecov_process_obs_sig) +
    ggtitle("EM:OM ratio of Fmsy when ecov process sigma = 0.1")
  ggsave(ecov_process_relFMSY, filename = here::here(outfile, paste0("plots/ecov_process_relFMSY.png")), width = 20)
  
  ecov_process_relMSY <- refpts_results %>%
    filter(OM_ecov_process_sig == 0.1) %>%
    ggplot()+
    geom_boxplot(aes(x=EMshortName, y = relMSY)) +
    facet_grid(. ~ OM_ecov_process_obs_sig) +
    ggtitle("EM:OM ratio of MSY when ecov process sigma = 0.1")
  ggsave(ecov_process_relMSY, filename = here::here(outfile, paste0("plots/ecov_process_relMSY.png")), width = 20)
  
  ecov_process_MohnsRho_SSB <- mohnsRho_results %>%
    filter(OM_ecov_process_sig == 0.1) %>%
    ggplot()+
    geom_boxplot(aes(x=EMshortName, y = EM_MohnsRho_SSB)) +
    facet_grid(. ~ OM_ecov_process_obs_sig) +
    ggtitle("EM Mohn's rho for SSB when ecov process sigma = 0.1")
  ggsave(ecov_process_MohnsRho_SSB, filename = here::here(outfile, paste0("plots/ecov_process_MohnsRho_SSB.png")), width = 10)
  
  ecov_process_MohnsRho_R <- mohnsRho_results %>%
    filter(OM_ecov_process_sig == 0.1) %>%
    ggplot()+
    geom_boxplot(aes(x=EMshortName, y = EM_MohnsRho_R)) +
    facet_grid(. ~ OM_ecov_process_obs_sig + OM_ecov_effect) +
    ggtitle("EM Mohn's rho for R when ecov process sigma = 0.1")
  ggsave(ecov_process_MohnsRho_R, filename = here::here(outfile, paste0("plots/ecov_process_MohnsRho_R.png")), width = 10)
    
  ecov_process_MohnsRho_F <- mohnsRho_results %>%
    filter(OM_ecov_process_sig == 0.1) %>%
    ggplot()+
    geom_boxplot(aes(x=EMshortName, y = EM_MohnsRho_F)) +
    facet_grid(. ~ OM_ecov_process_obs_sig) +
    ggtitle("EM Mohn's rho for F when ecov process sigma = 0.1")
  ggsave(ecov_process_MohnsRho_F, filename = here::here(outfile, paste0("plots/ecov_process_MohnsRho_F.png")), width = 10)
  
  ## Extra plots focused on qRand and qRandEcov models
  # qrand <- singleResults_converge %>%
  #   filter(EMshortName  =="EM_NONE_qRand" | EMshortName == "EM_NONE_qRandEcov") %>%
  #   ggplot() +
  #   geom_boxplot(aes(x=EMshortName, y = convergeRate)) +
  #   facet_grid(. ~ OM_ecov_process_sig + OM_ecov_process_cor + OM_ecov_process_obs_sig + ageComp_sig + log_index_sig)
  # ggsave(qrand, filename = here::here(outfile, paste0("qrand.png")), width = 20)
  
  # Extra plots breaking down differences in results between Ecov beta parameters
  beta_relSSBMSY <- refpts_results %>% 
    ggplot() +
    geom_boxplot(aes(x=EMshortName, y = relSSBMSY, fill = EMshortName)) +
    facet_grid(. ~ OM_ecov_effect) +
    scale_fill_grey(start = 0.45, end = 1.0) +
    theme(axis.text.x = element_blank()) +
    labs(fill = "EM", 
         x = "EM",
         y = "EM:OM ratio of SSBmsy") +
    ggtitle("EM:OM ratio of SSBMSY across environmental effect sizes")
  ggsave(beta_relSSBMSY, filename = here::here(outfile, "plots/beta_relSSBMSY_boxplot.png"))
  
  # plot OM ecov beta vs. EM beta by EM
  beta_Ind2 <- ageEcov_results %>%
    ggplot() +
    geom_point(aes(x=OM_ecov_effect, y = EM_ecovBeta_ind2)) + #!!! compares OM_ecov_effect to the correctly implemented index
    geom_abline() +
    facet_grid(.~EMshortName) +
    ggtitle("Environmental effect on fall index q")
  ggsave(beta_Ind2, filename = here::here(outfile, "plots/beta_Ind2.png"), width = 20)
  
  beta_Ind2_box <- ageEcov_results %>%
    ggplot() +
    geom_boxplot(aes(x=EMshortName, y = EM_ecovBeta_ind2/OM_ecov_effect)) +
    ggtitle("EM:OM ratio of environmental effect size")
  ggsave(beta_Ind2_box, filename = here::here(outfile, "plots/beta_Ind2_boxplot.png"), width = 5) # Formerly saved as beta_Ind2_box.png
  
  beta_Ind1 <- ageEcov_results %>%
    ggplot() +
    geom_point(aes(x=OM_ecov_effect, y = EM_ecovBeta_ind1)) +
    facet_grid(.~EMshortName) +
    ggtitle("Environmental effect on spring index q")
  ggsave(beta_Ind1, filename = here::here(outfile, "plots/beta_Ind1.png"), width = 20) 
  
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
    geom_boxplot(aes(x=EMshortName, y=relSSB, fill = EMshortName)) + 
    facet_grid(cols = vars(OM_ecov_process_sig, OM_ecov_process_cor, OM_ecov_process_obs_sig, ageComp_sig, log_index_sig), rows = vars(OM_ecov_effect, F_hist)) + # Facet by OM settings
    geom_hline(aes(yintercept = 1), color="red") +
    scale_fill_grey(start = 0.45, end = 1.0) +
    theme(axis.text.x = element_blank()) +
    labs(fill = "EM", 
         x = "EM",
         y = "EM:OM ratio of SSB") +
    ggtitle("EM:OM ratio of terminal SSB time series")
  ggsave(relSSB, filename = here::here(outfile, "plots/relativeSSB_terminal.png"), width = 10)
  
  # Extra plot of terminal year SSB facet by ecov_beta
  relSSB_ecovBeta <- results %>% filter(Year == 40) %>%
    ggplot()+
    geom_boxplot(aes(x=EMshortName, y=relSSB, fill = EMshortName)) + 
    facet_grid(. ~ OM_ecov_effect) + # Facet by OM ecov beta (effect size)
    geom_hline(aes(yintercept = 1), color="red") +
    scale_fill_grey(start = 0.45, end = 1.0) +
    theme(axis.text.x = element_blank()) +
    labs(fill = "EM", 
         x = "EM",
         y = "EM:OM ratio of SSB") +
    ggtitle("EM:OM ratio of terminal SSB time series by ecov effect size")
  ggsave(relSSB_ecovBeta, filename = here::here(outfile, "plots/relativeSSB_ecovBeta_terminal.png"), width = 10)
  
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
    geom_boxplot(aes(x=EMshortName, y=relF, fill = EMshortName)) + 
    facet_grid(cols = vars(OM_ecov_process_sig, OM_ecov_process_cor, OM_ecov_process_obs_sig, ageComp_sig, log_index_sig), rows = vars(OM_ecov_effect, F_hist)) + # Facet by OM settings
    geom_hline(aes(yintercept = 1), color="red") +
    scale_fill_grey(start = 0.45, end = 1.0) +
    theme(axis.text.x = element_blank()) +
    labs(fill = "EM", 
         x = "EM",
         y = "EM:OM ratio of F") +
    ggtitle("EM:OM ratio of terminal F")
  ggsave(relativeF_terminal, filename = here::here(outfile, "plots/relativeF_terminal.png"), width = 20)
  
  # Relative R
  #   # Time series
  # results %>% 
  #   ggplot()+
  #   geom_line(aes(x=Year, y = relR, group = seed, color = F_hist), alpha = 0.5) + # Could color lines by F history
  #   geom_hline(aes(yintercept = 1.0), color="red")
    # Terminal year
  relativeR_terminal <- results %>% filter(Year == 40) %>%
    ggplot()+
    geom_boxplot(aes(x=EMshortName, y=relR, fill = EMshortName)) + 
    facet_grid(cols = vars(OM_ecov_process_sig, OM_ecov_process_cor, OM_ecov_process_obs_sig, ageComp_sig, log_index_sig), rows = vars(OM_ecov_effect, F_hist)) + # Facet by OM settings
    geom_hline(aes(yintercept = 1), color="red") +
    scale_fill_grey(start = 0.45, end = 1.0) +
    theme(axis.text.x = element_blank()) +
    labs(fill = "EM", 
         x = "EM",
         y = "EM:OM ratio of R") +
    ggtitle("EM:OM ratio of terminal R")
  ggsave(relativeR_terminal, filename = here::here(outfile, "plots/relativeR_terminal.png"), width = 20)
  
  # Stock status SSB/SSBmsy
  statusSSB <- results %>%
    ggplot()+
    geom_line(aes(x=Year, y=status_SSB, group = seed, color = F_hist)) +
    geom_hline(aes(yintercept = 1), color = "red") +
    geom_hline(aes(yintercept = 0.5), color = "red", linetype = "dashed") + # Check this threshold definition
    ggtitle("Stock status time series: SSB/SSBmsy")
  ggsave(statusSSB, filename = here::here(outfile, "plots/status_SSB.png"), width = 20)
  
  # Stock status F/Fmsy
  statusF <- results %>%
    ggplot()+
    geom_line(aes(x=Year, y=status_F, group = seed, color = F_hist)) +
    ggtitle("Stock status time series: F/Fmsy")
  ggsave(statusF, filename = here::here(outfile, "plots/status_F.png"), width = 20)
  
  # Stock status Y/MSY
  statusY <- results %>%
    ggplot()+
    geom_line(aes(x=Year, y=status_Y, group = seed, color = F_hist)) +
    ggtitle("Stock status time series: Y/MSY")
  ggsave(statusY, filename = here::here(outfile, "plots/status_Y.png"), width = 20)
  
  # Relative q index 1
  relQInd1 <- results %>% 
    ggplot()+
    geom_line(aes(x=Year, y = relq_index1, group = seed), alpha = 0.5) + # Could color lines by F history
    geom_hline(aes(yintercept = 1.0), color="red") +
    facet_grid(.~ OM_ecov_effect + EMshortName) +
    ggtitle("Time series of EM:OM spring index catchability")
  ggsave(relQInd1, filename = here::here(outfile, "plots/relq_index1.png"), width = 20)
  
  # Relative q index 2
  relQInd2 <- results %>% 
    ggplot()+
    geom_line(aes(x=Year, y = relq_index2, group = seed), alpha = 0.5) + # Could color lines by F history
    geom_hline(aes(yintercept = 1.0), color="red") + 
    facet_grid(.~ OM_ecov_effect + EMshortName) +
    ggtitle("Time series of EM:OM fall index catchability")
  ggsave(relQInd2, filename = here::here(outfile, "plots/relq_index2.png"), width = 20)
  
  relQInd2_box <- results %>% 
    ggplot()+
    geom_boxplot(aes(x=EMshortName, y = relq_index2), alpha = 0.5) + # Could color lines by F history
    facet_grid(.~ OM_ecov_effect) +
    geom_hline(aes(yintercept = 1), color="red") +
    scale_fill_grey(start = 0.45, end = 1.0) +
    theme(axis.text.x = element_blank()) +
    labs(fill = "EM", 
         x = "EM",
         y = "EM:OM ratio of fall index q") +
    ggtitle("EM:OM fall index catchability")
  ggsave(relQInd2_box, filename = here::here(outfile, "plots/relq_index2_boxplot.png"), width = 20) # Formerly saved as relq_index2_box.png
  
  # Plot index 2 q for OM only - pick handful of examples
  q_ind2_OM <- results %>% filter(OM_ecov_effect > 0) %>% #filter(sim %in% c(201, 18718, 18905, 9126)) %>%
    ggplot()+
    geom_line(aes(x=Year, y = OM_q_index2, group = seed), alpha = 0.5) +
    #geom_line(aes(x=Year, y=seq(0:5,39)), colour = "red")+
    facet_grid(. ~ OM_ecov_effect + EMshortName) +
    ggtitle("OM fall catchability when environmental effect > 0")
  ggsave(q_ind2_OM, filename = here::here(outfile, "plots/q_ind2_OM.png"), width = 20)
  
  # Plot index 2 q for EM only - pick handful of examples
  q_ind2_EM <- results %>% filter(OM_ecov_effect > 0) %>% filter(EM_miss_season== "NONE") %>% #filter(sim %in% c(201, 18718, 18905, 9126)) %>%
    ggplot()+
    geom_line(aes(x=Year, y = EM_q_index2, group = seed), alpha = 0.5) +
    #geom_line(aes(x=Year, y=seq(0:5,39)), colour = "red")+
    facet_grid(. ~ OM_ecov_effect + EMshortName) +
    ggtitle("EM fall catchability when environmental effect > 0")
  ggsave(q_ind2_EM, filename = here::here(outfile, "plots/q_ind2_EM.png"), width = 20)
  
  # Plot environmental covariate that drives index 2 for q in OM
  ecov_OM <- results %>% filter(OM_ecov_effect > 0) %>%
    ggplot() +
    geom_line(aes(x=Year, y= OM_Ecov_obs)) + # If broken try V85 (85th column in post-processed results)
    facet_grid(.~OM_ecov_effect) +
    ggtitle("OM environmental covariate time series")
  ggsave(ecov_OM, filename = here::here(outfile, "plots/ecov_OM.png"))
  
  # Plot predicted environmental covariate from EM
  ecov_EM <- results %>% filter(OM_ecov_effect > 0) %>%
    ggplot() +
    geom_line(aes(x=Year, y=EM_Ecov_pred)) +
    facet_grid(.~OM_ecov_effect) +
    ggtitle("EM environmental covariate time series")
  ggsave(ecov_EM, filename = here::here(outfile, "plots/ecov_EM.png"))

  # Plot relative recovery of ecov observations
  ecov_relObs <- results %>% filter(OM_ecov_effect > 0) %>%
    ggplot() +
    geom_boxplot(aes(x=EMshortName, y=EM_Ecov_pred/OM_ecov_effect, fill = EMshortName)) +
    facet_grid(.~OM_ecov_effect) +
    scale_fill_grey(start = 0.45, end = 1.0) +
    theme(axis.text.x = element_blank()) +
    labs(fill = "EM", 
         x = "EM",
         y = "EM:OM ratio of environmental covariate") +
    ggtitle("EM:OM ratio of the environmental covariate time series")
  ggsave(ecov_relObs, filename = here::here(outfile, "plots/ecov_relObs.png"))
  
  # Plot terminal year relative SSB vs. relative F
  SSB_F_terminal <- results %>% filter(Year == 40) %>%
    ggplot()+
    geom_point(aes(x=relSSB, y =relF)) + 
    facet_grid(cols = vars(OM_ecov_process_sig, OM_ecov_process_cor, OM_ecov_process_obs_sig, ageComp_sig, log_index_sig), rows = vars(OM_ecov_effect, Fhist)) + # Facet by OM settings
    geom_hline(yintercept = 1, color = "red") + geom_vline(xintercept = 1, color = "red") +
    ggtitle("EM:OM terminal SSB vs. F")
  ggsave(SSB_F_terminal, filename = here::here(outfile, "plots/SSB_F_terminal.png"), width = 20)
  
  # Plot terminal year EM:OM relative catch
  relCatch  <- results %>% filter(Year == 40) %>%
    ggplot() +
    geom_line(aes(x=EMshortName, y = EM_Catch/OM_Catch, fill = EMshortName)) +
    facet_grid(cols = vars(OM_ecov_process_sig, OM_ecov_process_cor, OM_ecov_process_obs_sig, ageComp_sig, log_index_sig), rows = vars(OM_ecov_effect, Fhist)) + # Facet by OM settings
    scale_fill_grey(start = 0.45, end = 1.0) +
    theme(axis.text.x = element_blank()) +
    labs(fill = "EM", 
         x = "EM",
         y = "EM:OM ratio of catch") +
    ggtitle("Terminal EM:OM catch")
  ggsave(relCatch, filename = here::here(outfile, "plots/relCatch_terminal.png"))
  
  # Plot median year EM:OM relative catch
  relCatch_med <- singleResults %>% 
    ggplot() +
    geom_boxplot(aes(x=EMshortName, y = med_relCatch, fill = EMshortName)) +
    facet_grid(cols = vars(OM_ecov_process_sig, OM_ecov_process_cor, OM_ecov_process_obs_sig, ageComp_sig, log_index_sig), rows = vars(OM_ecov_effect, Fhist)) + # Facet by OM settings
    geom_hline(aes(yintercept = 1), color="red") +
    scale_fill_grey(start = 0.45, end = 1.0) +
    theme(axis.text.x = element_blank()) +
    labs(fill = "EM", 
         x = "EM",
         y = "EM:OM ratio of catch") +
    ggtitle("Median EM:OM catch") # Median over years in each simulation
  ggsave(relCatch_med, filename = here::here(outfile, "plots/relCatch_boxplot.png"))
  
  # Plot time series of EM:OM relative catch
  relCatch_series  <- results %>% 
    ggplot() +
    geom_line(aes(x=EMshortName, y = EM_Catch/OM_Catch, fill = EMshortName)) +
    facet_grid(cols = vars(OM_ecov_process_sig, OM_ecov_process_cor, OM_ecov_process_obs_sig, ageComp_sig, log_index_sig), rows = vars(OM_ecov_effect, Fhist)) + # Facet by OM settings
    scale_fill_grey(start = 0.45, end = 1.0) +
    theme(axis.text.x = element_blank()) +
    labs(fill = "EM", 
         x = "EM",
         y = "EM:OM ratio of Catch") +
    ggtitle("EM:OM catch time series")
  ggsave(relCatch_series, filename = here::here(outfile, "plots/relCatch_series.png"))
  
  # Plot indices: EM_agg_ind1, EM_agg_ind2, EM_ind1_paa, EM_ind2_paa, OM_agg_ind1, OM_agg_ind2, OM_ind1_paa, OM_ind2_paa
  
  
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
    geom_point(aes(x = status_SSB, y = status_F, color = exp(log_index_sig))) + # log index sigma
    ggtitle("Terminal year stock status")
  ggsave(kobe_terminal, filename = here::here(outfile, "plots/kobe_terminal.png"), width = 10)
    

  
}








