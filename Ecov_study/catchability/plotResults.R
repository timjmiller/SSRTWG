
#' @param convergedONLY A boolean, when TRUE plots only show converged runs, default = TRUE.
#' @param outfile A string for the directory where a plot folder will be generated


# results <- perfMet
# convergedONLY = TRUE

plotResults <- function(results = NULL, convergedONLY = TRUE, outfile = here::here()){
  
  # Simulation summary
  results %>% count(OMshortName, EMshortName) %>% mutate(nsim = n/40) %>% write.csv(., paste0(outfile,"/simSummary.csv")) # /40 years so nsim = number of full simulations for each OM/EM
  # results %>% count(OMshortName) %>% mutate(nsim = n/40/4) # Full sims of all 4 correctly specified EMs
  # results %>% filter(OMshortName == 15) %>% select(EMshortName) %>% unique()
  
  # setdiff(subsetOM$OMname, results$OMshortName) # Check what simulations still not run from subsetOM (in Catchability_ecov_sims.Rmd)
  
  if(convergedONLY == TRUE){
    # Filter out converged runs
    results <- results %>% filter(EM_converged == TRUE)
    # Calculate number of converged runs
    convergeCount <- results %>% count(OMshortName, EMshortName) %>% mutate(nsim = n/40) # /40 years so nsim = number of full simulations for each OM/EM
    #convergeCount$OMshortName <- convergeCount$OMshortName %>% as.numeric()
    totalCount <- read.csv(paste0(outfile,"/simSummary.csv"))
    totalCount$OMshortName <- totalCount$OMshortName %>% as.character()
    
    convergeRate <- full_join(convergeCount, totalCount, by = c("OMshortName", "EMshortName")) %>% mutate(convergeRate = nsim.x/nsim.y) %>% drop_columns(c("n.x", "nsim.x", "X", "n.y", "nsim.y"))
  }
  
  numericIndex <- which(colnames(results) %in% c("OMshortName", 'EMshortName', "EM_miss_season", "EM_miss_q", "F_hist", "Converged") == FALSE)
  results[,numericIndex] <- sapply(results[,numericIndex], as.numeric) # This introduces NAs by coercion since using NAs as placeholders
  
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
                                                                  relselCat_10 = unique(relselCat_10),
                                                                  
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
                                                                  EM_ecovBeta_ind2 = unique(EM_ecovBeta_ind2)) # EM beta parameter for index 2
    
  if(convergedONLY == TRUE){
    singleResults_converge <- full_join(convergeRate, singleResults, by = c("OMshortName", "EMshortName")) %>% # remember that some runs only had NoEcov and qRand since the Ecov not used so Ecov and qRandEcov have fewer samples in these plots
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
    ggsave(convergePlot, filename = here::here(outfile, "convergeRate.png"), width = 20)
    # Label order: Ecov, NoEcov, qRand, qRandEcov
    
    ## Extra plots looking at results across F history
    Fhist <- singleResults_converge %>%
      ggplot() +
      geom_boxplot(aes(x=EMshortName, y = convergeRate)) +
      facet_grid(. ~ Fhist)
    ggsave(Fhist, filename = here::here(outfile, paste0("Fhist_summary.png")), width = 10)
    
    Fhist_relFMSY <- singleResults_converge %>%
      ggplot() +
      geom_boxplot(aes(x=EMshortName, y = relFMSY)) +
      facet_grid(. ~ Fhist)
    ggsave(Fhist_relFMSY, filename = here::here(outfile, paste0("Fhist_relFMSY_summary.png")), width = 10)
    
    Fhist_relMSY <- singleResults_converge %>%
      ggplot() +
      geom_boxplot(aes(x=EMshortName, y = relMSY)) +
      facet_grid(. ~ Fhist)
    ggsave(Fhist_relMSY, filename = here::here(outfile, paste0("Fhist_relMSY_summary.png")), width = 10)
    
    Fhist_relSSBMSY <- singleResults_converge %>%
      ggplot() +
      geom_boxplot(aes(x=EMshortName, y = relSSBMSY)) +
      facet_grid(. ~ Fhist)
    ggsave(Fhist_relSSBMSY, filename = here::here(outfile, paste0("Fhist_relSSBMSY_summary.png")), width = 10)
  }
                                                                                 
  # Time series include status SSB/F/Y, relSSB, relF, relR, relq both indices, relFAA, relNAA, relCAA
  
  # MohnsRho_SSB
  mohnsRho_SSB_series <- singleResults %>% 
    ggplot() + 
    geom_point(aes(x=sim, y=EM_MohnsRho_SSB, color = OM_ecov_effect))  # Plot by sim
  ggsave(mohnsRho_SSB_series, filename = here::here(outfile, "mohnsRho_SSB_series.png"), width = 10)
  
  mohnsRho_SSB <- singleResults %>%
    ggplot()+
    geom_boxplot(aes(x=EMshortName, y=EM_MohnsRho_SSB)) + facet_grid(. ~ OM_ecov_process_sig + OM_ecov_process_cor + OM_ecov_process_obs_sig + ageComp_sig + log_index_sig + OM_ecov_effect) + # Facet by OM settings
    geom_hline(aes(yintercept = 0), color="red")
    # ggtitle("OM_ecov_process_sig")
  ggsave(mohnsRho_SSB, filename = here::here(outfile, "mohnsRho_SSB.png"), width = 20)
  
  # MohnsRho_F
  mohnsRho_F <- singleResults %>%
    ggplot()+
    geom_boxplot(aes(x=EMshortName, y=EM_MohnsRho_F)) + facet_grid(. ~ OM_ecov_process_sig + OM_ecov_process_cor + OM_ecov_process_obs_sig + ageComp_sig + log_index_sig + OM_ecov_effect) + # Facet by OM settings
    geom_hline(aes(yintercept = 0), color="red")
  ggsave(mohnsRho_F, filename = here::here(outfile, "mohnsRho_F.png"), width = 20)
  
  # MohnsRho_R
  mohnsRho_R <- singleResults %>%
    ggplot()+
    geom_boxplot(aes(x=EMshortName, y=EM_MohnsRho_R)) + facet_grid(. ~ OM_ecov_process_sig + OM_ecov_process_cor + OM_ecov_process_obs_sig + ageComp_sig + log_index_sig + OM_ecov_effect) + # Facet by OM settings
    geom_hline(aes(yintercept = 0), color="red")
  ggsave(mohnsRho_R, filename = here::here(outfile, "mohnsRho_R.png"), width = 20)
  
  # Relative SSBmsy
  relSSBMSY <- singleResults %>%
    ggplot()+
    geom_boxplot(aes(x=EMshortName, y=relSSBMSY)) + facet_grid(. ~ OM_ecov_process_sig + OM_ecov_process_cor + OM_ecov_process_obs_sig + ageComp_sig + log_index_sig + OM_ecov_effect) + # Facet by OM settings
    geom_hline(aes(yintercept = 1), color="red")
  ggsave(relSSBMSY, filename = here::here(outfile, "relSSBMSY.png"), width = 20)
  
  # Relative Fmsy
  relFMSY <- singleResults %>%
    ggplot()+
    geom_boxplot(aes(x=EMshortName, y=relFMSY)) + facet_grid(. ~ OM_ecov_process_sig + OM_ecov_process_cor + OM_ecov_process_obs_sig + ageComp_sig + log_index_sig + OM_ecov_effect) + # Facet by OM settings
    geom_hline(aes(yintercept = 1), color="red")
  ggsave(relFMSY, filename = here::here(outfile, "relFMSY.png"), width = 20)
  
  # Relative MSY
  relMSY <- singleResults %>%
    ggplot()+
    geom_boxplot(aes(x=EMshortName, y=relMSY)) + facet_grid(. ~ OM_ecov_process_sig + OM_ecov_process_cor + OM_ecov_process_obs_sig + ageComp_sig + log_index_sig + OM_ecov_effect) + # Facet by OM settings
    geom_hline(aes(yintercept = 1), color="red")
  ggsave(relMSY, filename = here::here(outfile, "relMSY.png"), width = 20)
  
  # Relative Index 1 
  relSelInd1 <- singleResults %>% 
    pivot_longer(., cols = c(relselInd1_1, relselInd1_2, relselInd1_3, relselInd1_4, relselInd1_5, relselInd1_6, relselInd1_7, relselInd1_8, relselInd1_9, relselInd1_10), names_to = "age", names_prefix = "relselInd1_", values_to = "sel")  
  relSelInd1$age <- factor(relSelInd1$age, levels = c(1:10))
  relSelInd1plot <- relSelInd1 %>% 
    ggplot()+
    geom_boxplot(aes(x=EMshortName, y=sel, fill = age)) + facet_grid(. ~ OM_ecov_process_sig + OM_ecov_process_cor + OM_ecov_process_obs_sig + ageComp_sig + log_index_sig + OM_ecov_effect) + # Facet by OM settings
    geom_hline(aes(yintercept =1), color="red") +
    ggtitle("Index 1 relative selectivity")
  ggsave(relSelInd1plot, filename = here::here(outfile, "relSelInd1.png"), width = 20)
  
  for(iage in 1:10){
    relSelInd1 <- singleResults %>% 
      pivot_longer(., cols = c(relselInd1_1, relselInd1_2, relselInd1_3, relselInd1_4, relselInd1_5, relselInd1_6, relselInd1_7, relselInd1_8, relselInd1_9, relselInd1_10), names_to = "age", names_prefix = "relselInd1_", values_to = "sel")  
    relSelInd1$age <- factor(relSelInd1$age, levels = c(1:10))
    relSelInd1plot <- relSelInd1 %>% filter(age == iage) %>%
      ggplot()+
      geom_boxplot(aes(x=EMshortName, y=sel)) + facet_grid(. ~ OM_ecov_process_sig + OM_ecov_process_cor + OM_ecov_process_obs_sig + ageComp_sig + log_index_sig + OM_ecov_effect) + # Facet by OM settings
      geom_hline(aes(yintercept =1), color="red") +
      ggtitle(paste0("Index 1 relative selectivity ", iage)) +
      ylim(0, 2.7)
    ggsave(relSelInd1plot, filename = here::here(outfile, paste0("relSelInd1_", iage, ".png")), width = 20)
  }
  
  # Relative Index 2 
  relSelInd2 <- singleResults %>% 
    pivot_longer(., cols = c(relselInd2_1, relselInd2_2, relselInd2_3, relselInd2_4, relselInd2_5, relselInd2_6, relselInd2_7, relselInd2_8, relselInd2_9, relselInd2_10), names_to = "age", names_prefix = "relselInd2_", values_to = "sel") 
  relSelInd2$age <- factor(relSelInd2$age, levels = c(1:10))
  relSelInd2plot <- relSelInd2 %>%
    ggplot()+
    geom_boxplot(aes(x=EMshortName, y=sel, fill = age)) + facet_grid(. ~ OM_ecov_process_sig + OM_ecov_process_cor + OM_ecov_process_obs_sig + ageComp_sig + log_index_sig + OM_ecov_effect) + # Facet by OM settings
    geom_hline(aes(yintercept =1), color="red") +
    ggtitle("Index 2 relative selectivity")
     # ageComp_sig impacts spread of results across all ages
  ggsave(relSelInd2plot, filename = here::here(outfile, "relSelInd2.png"), width = 20)
  
  for(iage in 1:10){
    relSelInd2 <- singleResults %>% 
      pivot_longer(., cols = c(relselInd2_1, relselInd2_2, relselInd2_3, relselInd2_4, relselInd2_5, relselInd2_6, relselInd2_7, relselInd2_8, relselInd2_9, relselInd2_10), names_to = "age", names_prefix = "relselInd2_", values_to = "sel") 
    relSelInd2$age <- factor(relSelInd2$age, levels = c(1:10))
    relSelInd2plot <- relSelInd2 %>% filter(age == iage) %>%
      ggplot()+
      geom_boxplot(aes(x=EMshortName, y=sel)) + facet_grid(. ~ OM_ecov_process_sig + OM_ecov_process_cor + OM_ecov_process_obs_sig + ageComp_sig + log_index_sig + OM_ecov_effect) + # Facet by OM settings
      geom_hline(aes(yintercept =1), color="red") +
      ggtitle(paste0("Index 2 relative selectivity ", iage)) +
      ylim(0, 3)
    # ageComp_sig impacts spread of results across all ages
    ggsave(relSelInd2plot, filename = here::here(outfile, paste0("relSelInd2_", iage, ".png")), width = 20)
  }
  
  # Relative selectivity fleet
  relSelFleet <- singleResults %>%
    pivot_longer(., cols = c(relselCat_1, relselCat_2, relselCat_3, relselCat_4, relselCat_5, relselCat_6, relselCat_7, relselCat_8, relselCat_9, relselCat_10), names_to = "age", names_prefix = "relselCat_", values_to = "sel")
  relSelFleet$age <- factor(relSelFleet$age, levels = c(1:10))
  for(iage in 1:10){
    relSelFleetplot <- relSelFleet %>% filter(age == iage) %>%
      ggplot() +
      geom_boxplot(aes(x = EMshortName, y = sel)) + facet_grid(. ~ OM_ecov_process_sig + OM_ecov_process_cor + OM_ecov_process_obs_sig + ageComp_sig + log_index_sig + OM_ecov_effect) + # Facet by OM settings
      geom_hline(aes(yintercept = 1), color="red")+
      ggtitle(paste0("Fleet relative selectivity ", iage)) +
      ylim(0,2.7)
    ggsave(relSelFleetplot, filename = here::here(outfile, paste0("relSelFleet_", iage, ".png")), width = 20)
  }
  
  # Median (across years) relative FAA
  med_relFAA <- singleResults %>%
    pivot_longer(., cols = c(med_relFAA_1, med_relFAA_2, med_relFAA_3, med_relFAA_4, med_relFAA_5, med_relFAA_6, med_relFAA_7, med_relFAA_8, med_relFAA_9, med_relFAA_10), names_to = "age", names_prefix = "med_relFAA_", values_to = "relFAA")
  med_relFAA$age <- factor(med_relFAA$age, levels = c(1:10))
  for(iage in 1:10){
    med_relFAAplot <- med_relFAA %>% filter(age == iage) %>%
      ggplot() +
      geom_boxplot(aes(x = EMshortName, y = relFAA)) + 
      facet_grid(. ~ OM_ecov_process_sig + OM_ecov_process_cor + OM_ecov_process_obs_sig + ageComp_sig + log_index_sig + OM_ecov_effect) + # Facet by OM settings
      geom_hline(aes(yintercept = 1), color="red") +
      ggtitle(paste0("Median relative F at age ", iage)) +
      ylim(-5,5)
    ggsave(med_relFAAplot, filename = here::here(outfile, paste0("relFAA_", iage, ".png")), width = 20)
  }
  
  # Median (across years) relative CAA
  med_relCAA <- singleResults %>%
    pivot_longer(., cols = c(med_relCAA_1, med_relCAA_2, med_relCAA_3, med_relCAA_4, med_relCAA_5, med_relCAA_6, med_relCAA_7, med_relCAA_8, med_relCAA_9, med_relCAA_10), names_to = "age", names_prefix = "med_relCAA_", values_to = "relCAA")
  med_relCAA$age <- factor(med_relCAA$age, levels = c(1:10))
  for(iage in 1:10){
    med_relCAAplot <- med_relCAA %>% filter(age == iage) %>%
      ggplot() +
      geom_boxplot(aes(x = EMshortName, y = relCAA)) + 
      facet_grid(. ~ OM_ecov_process_sig + OM_ecov_process_cor + OM_ecov_process_obs_sig + ageComp_sig + log_index_sig + OM_ecov_effect) + # Facet by OM settings
      geom_hline(aes(yintercept = 1), color="red") +
      ggtitle(paste0("Median relative Catch at age ", iage)) +
      ylim(-5,5)
    ggsave(med_relCAAplot, filename = here::here(outfile, paste0("relCAA_", iage, ".png")), width = 20)
  }
  
  # Median (across years) relative NAA
  med_relNAA <- singleResults %>%
    pivot_longer(., cols = c(med_relNAA_1, med_relNAA_2, med_relNAA_3, med_relNAA_4, med_relNAA_5, med_relNAA_6, med_relNAA_7, med_relNAA_8, med_relNAA_9, med_relNAA_10), names_to = "age", names_prefix = "med_relNAA_", values_to = "relNAA")
  med_relNAA$age <- factor(med_relNAA$age, levels = c(1:10))
  for(iage in 1:10){
    med_relNAAplot <- med_relNAA %>% filter(age == iage) %>%
      ggplot() +
      geom_boxplot(aes(x = EMshortName, y = relNAA)) + 
      facet_grid(. ~ OM_ecov_process_sig + OM_ecov_process_cor + OM_ecov_process_obs_sig + ageComp_sig + log_index_sig + OM_ecov_effect) + # Facet by OM settings
      geom_hline(aes(yintercept = 1), color="red") +
      ggtitle(paste0("Median relative N at age ", iage)) +
      ylim(-5,5)
    ggsave(med_relNAAplot, filename = here::here(outfile, paste0("relNAA_", iage, ".png")), width = 20)
  }
  
  ## Extra plots small ecov process sig & high ecov process obs sig
  ecov_process_relSSBMSY <- singleResults %>%
    filter(OM_ecov_process_sig == 0.1) %>%
    ggplot()+
    geom_boxplot(aes(x=EMshortName, y = relSSBMSY)) +
    facet_grid(. ~ OM_ecov_process_obs_sig)
  ggsave(ecov_process_relSSBMSY, filename = here::here(outfile, paste0("ecov_process_relSSBMSY.png")), width = 20)
  
  ecov_process_relFMSY <- singleResults %>%
    filter(OM_ecov_process_sig == 0.1) %>%
    ggplot()+
    geom_boxplot(aes(x=EMshortName, y = relFMSY)) +
    facet_grid(. ~ OM_ecov_process_obs_sig)
  ggsave(ecov_process_relFMSY, filename = here::here(outfile, paste0("ecov_process_relFMSY.png")), width = 20)
  
  ecov_process_relMSY <- singleResults %>%
    filter(OM_ecov_process_sig == 0.1) %>%
    ggplot()+
    geom_boxplot(aes(x=EMshortName, y = relMSY)) +
    facet_grid(. ~ OM_ecov_process_obs_sig)
  ggsave(ecov_process_relMSY, filename = here::here(outfile, paste0("ecov_process_relMSY.png")), width = 20)
  
  ecov_process_MohnsRho_SSB <- singleResults %>%
    filter(OM_ecov_process_sig == 0.1) %>%
    ggplot()+
    geom_boxplot(aes(x=EMshortName, y = EM_MohnsRho_SSB)) +
    facet_grid(. ~ OM_ecov_process_obs_sig)
  ggsave(ecov_process_MohnsRho_SSB, filename = here::here(outfile, paste0("ecov_process_MohnsRho_SSB.png")), width = 10)
  
  ecov_process_MohnsRho_R <- singleResults %>%
    filter(OM_ecov_process_sig == 0.1) %>%
    ggplot()+
    geom_boxplot(aes(x=EMshortName, y = EM_MohnsRho_R)) +
    facet_grid(. ~ OM_ecov_process_obs_sig + OM_ecov_effect)
  ggsave(ecov_process_MohnsRho_R, filename = here::here(outfile, paste0("ecov_process_MohnsRho_R.png")), width = 10)
    
  ecov_process_MohnsRho_F <- singleResults %>%
    filter(OM_ecov_process_sig == 0.1) %>%
    ggplot()+
    geom_boxplot(aes(x=EMshortName, y = EM_MohnsRho_F)) +
    facet_grid(. ~ OM_ecov_process_obs_sig)
  ggsave(ecov_process_MohnsRho_F, filename = here::here(outfile, paste0("ecov_process_MohnsRho_F.png")), width = 10)
  
  ## Extra plots focused on qRand and qRandEcov models
  # qrand <- singleResults_converge %>%
  #   filter(EMshortName  =="EM_NONE_qRand" | EMshortName == "EM_NONE_qRandEcov") %>%
  #   ggplot() +
  #   geom_boxplot(aes(x=EMshortName, y = convergeRate)) +
  #   facet_grid(. ~ OM_ecov_process_sig + OM_ecov_process_cor + OM_ecov_process_obs_sig + ageComp_sig + log_index_sig)
  # ggsave(qrand, filename = here::here(outfile, paste0("qrand.png")), width = 20)
  
  # Extra plots breaking down differences in results between Ecov beta parameters
  singleResults %>% 
    ggplot() +
    geom_boxplot(aes(x=EMshortName, y = relSSBMSY)) +
    facet_grid(. ~ OM_ecov_effect)
  
  # plot OM ecov beta vs. EM beta by EM
  beta_Ind2 <- singleResults %>%
    ggplot() +
    geom_point(aes(x=OM_ecov_effect, y = EM_ecovBeta_ind2)) +
    geom_abline() +
    facet_grid(.~EMshortName)
  ggsave(beta_Ind2, filename = here::here(outfile, "beta_Ind2.png"), width = 20)
  
  beta_Ind2_box <- singleResults %>%
    ggplot() +
    geom_boxplot(aes(x=EMshortName, y = EM_ecovBeta_ind2/OM_ecov_effect)) 
  ggsave(beta_Ind2_box, filename = here::here(outfile, "beta_Ind2_box.png"), width = 5)
  
  beta_Ind1 <- singleResults %>%
    ggplot() +
    geom_point(aes(x=OM_ecov_effect, y = EM_ecovBeta_ind1)) +
    facet_grid(.~EMshortName)
  ggsave(beta_Ind1, filename = here::here(outfile, "beta_Ind1.png"), width = 20)
  
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
  ggsave(relSSB, filename = here::here(outfile, "relativeSSB_terminal.png"), width = 10)
  
  # Extra plot of terminal year SSB facet by ecov_beta
  relSSB_ecovBeta <- results %>% filter(Year == 40) %>%
    ggplot()+
    geom_boxplot(aes(x=EMshortName, y=relSSB)) + facet_grid(. ~ OM_ecov_effect) + # Facet by OM ecov beta (effect size)
    geom_hline(aes(yintercept = 1), color="red")
  ggsave(relSSB_ecovBeta, filename = here::here(outfile, "relativeSSB_ecovBeta_terminal.png"), width = 10)
  
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
    geom_hline(aes(yintercept = 1.0), color="red") + 
    facet_grid(.~ OM_ecov_effect + EMshortName)
  ggsave(relQInd2, filename = here::here(outfile, "relq_index2.png"), width = 20)
  
  relQInd2_box <- results %>% 
    ggplot()+
    geom_boxplot(aes(x=EMshortName, y = relq_index2), alpha = 0.5) + # Could color lines by F history
    geom_hline(aes(yintercept = 1.0), color="red") + 
    facet_grid(.~ OM_ecov_effect)
  ggsave(relQInd2_box, filename = here::here(outfile, "relq_index2_box.png"), width = 20)
  
  # Plot index 2 q for OM only - pick handful of examples
  q_ind2_OM <- results %>% filter(OM_ecov_effect > 0) %>% #filter(sim %in% c(201, 18718, 18905, 9126)) %>%
    ggplot()+
    geom_line(aes(x=Year, y = OM_q_index2, group = seed), alpha = 0.5) +
    #geom_line(aes(x=Year, y=seq(0:5,39)), colour = "red")+
    facet_grid(. ~ OM_ecov_effect + EMshortName)
  ggsave(q_ind2_OM, filename = here::here(outfile, "q_ind2_OM.png"), width = 20)
  
  # Plot index 2 q for EM only - pick handful of examples
  q_ind2_EM <- results %>% filter(OM_ecov_effect > 0) %>% filter(EM_miss_season== "NONE") %>% #filter(sim %in% c(201, 18718, 18905, 9126)) %>%
    ggplot()+
    geom_line(aes(x=Year, y = EM_q_index2, group = seed), alpha = 0.5) +
    #geom_line(aes(x=Year, y=seq(0:5,39)), colour = "red")+
    facet_grid(. ~ OM_ecov_effect + EMshortName)
  ggsave(q_ind2_EM, filename = here::here(outfile, "q_ind2_EM.png"), width = 20)
  
  # Plot environmental covariate that drives index 2 for q in OM
  ecov_OM <- results %>% filter(OM_ecov_effect > 0) %>%
    ggplot() +
    geom_line(aes(x=Year, y= OM_Ecov_obs)) + # If broken try V85 (85th column in post-processed results)
    facet_grid(.~OM_ecov_effect)
  ggsave(ecov_OM, filename = here::here(outfile, "ecov_OM.png"))
  
  # # Plot predicted environmental covariate from EM
  # ecov_EM <- results %>% filter(OM_ecov_effect > 0) %>%
  #   ggplot() +
  #   geom_line(aes(x=Year, y=EM_Ecov_pred)) +
  #   facet_grid(.~OM_ecov_effect)
  # ggsave(ecov_EM, filename = here::here(outfile, "ecov_EM.png"))
  # 
  # # Plot relative recovery of ecov observations
  # ecov_relObs <- results %>% filter(OM_ecov_effect > 0) %>%
  #   ggplot() +
  #   geom_boxplot(aes(x=EMshortName, y=EM_Ecov_pred/OM_ecov_effect)) %>%
  #   facet_grid(.~OM_ecov_effect)
  # ggsave(ecov_relObs, filename = here::here(outfile, "ecov_relObs.png"))
  
  # Plot terminal year relative SSB vs. relative F
  SSB_F_terminal <- results %>% filter(Year == 40) %>%
    ggplot()+
    geom_point(aes(x=relSSB, y =relF)) + 
    facet_grid(. ~ OM_ecov_process_sig + OM_ecov_process_cor + OM_ecov_process_obs_sig + ageComp_sig + log_index_sig) + # Facet by OM settings
    geom_hline(yintercept = 1, color = "red") + geom_vline(xintercept = 1, color = "red")
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








