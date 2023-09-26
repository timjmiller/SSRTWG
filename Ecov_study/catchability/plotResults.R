



results <- perfMet

plotResults <- function(results = NULL){
  
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
  
  # Plot MohnsRho_SSB
  singleResults %>% 
    ggplot() + 
    geom_point(aes(x=sim, y=EM_MohnsRho_SSB, color = OM_ecov_effect)) # Plot by sim
  
  singleResults %>%
    ggplot()+
    geom_boxplot(aes(x=EMshortName, y=EM_MohnsRho_SSB)) + facet_grid(. ~ OM_ecov_process_sig + OM_ecov_process_cor) + # Facet by OM settings
    ggtitle("OM_ecov_process_sig")
    
  # Plot SSB vs. F
  
  
}








