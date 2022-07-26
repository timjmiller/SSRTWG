# functions to calculate mohn's rho for wham parameter estimates and make plots
# liz brooks


# function calc.rho.diffs ====

calc.rho.diffs <- function(mod, param, peel0, od, save.name, ...)  {
  # mod is a fitted model
  # param is the parameter to calculate rho (q, index, MAA, NAA.devs, or Fleet.Sel
  # peel0 is the object with values for original model fit ("peel 0")
  # od is output directory (passed from plot.retro.pars function)
  # save.name is text prepended to each saved filename (passed from plot.retro.pars function)
  
  
  npeels <- length(mod$peels) 
  nyears <- length(mod$years)
  nindices <-  mod$input$data$n_indices
  nages <- mod$input$data$n_ages
  ind.sel.num <- mod$input$data$selblock_pointer_indices[nyears,]
  last.sel.row <- sapply(mod$rep$selAA[ind.sel.num], function(x)  x[nyears,])
  fleet.sel.num <- mod$input$data$selblock_pointer_fleet[nyears,]        #fleet sel blk in final year
  fleet.pointer.yrs <- mod$input$data$selblock_pointer_fleets #vector of fleet sel blk pointers 
  n.sigma <- mod$input$data$n_NAA_sigma 
  
  
  # CATCHABILITY ====
  if(param =="q")  {
    
    rho.mat <- matrix(NA, nrow=1, ncol=nindices)
    colnames(rho.mat) <- paste0(rep(param, nindices), seq(1,nindices))
    
    for (i in 1:nindices) {
      tmp.mat <- matrix(NA, nrow=npeels, ncol=4 )
      colnames(tmp.mat) <- c("Peel", param, "Diff", "Rel.diff")
      
      for (p in 1:npeels)  {
        
        tmp.mat[p,1] <- p
        tmp.mat[p,2] <- mod$peels[[p]]$rep$q[(nyears-p),i]
        tmp.mat[p,3] <- (tmp.mat[p,2]-peel0[(nyears-p),(1+i)]) 
        tmp.mat[p,4] <- (tmp.mat[p,2]-peel0[(nyears-p),(1+i)])/peel0[(nyears-p),(1+i)]
        
      } # end loop over p
      
      write.csv(tmp.mat, file=file.path(od, paste0(save.name,"_Rel.Diffs_" ,param, "_", i, ".csv") ), row.names=F)
      
      rho.mat[1,i] <- mean(tmp.mat[,4])
      
      
    }# end loop over i (number of indices)
    
    write.csv(rho.mat, file=file.path(od, paste0(save.name,"_Rho_", param,   ".csv") ), row.names=F)
    
    
  } # end test for q
  
  
  # INDEX ====
  if(param =="index")  {
    
    rho.mat <- matrix(NA, nrow=nages, ncol=nindices)
    colnames(rho.mat) <- paste0(rep(param, nindices), seq(1,nindices))
    
    for (i in 1:nindices) {
      
      for (a in 1:nages)  {
        tmp.mat <- matrix(NA, nrow=npeels, ncol=4 )
        colnames(tmp.mat) <- c("Peel", param, "Diff", "Rel.diff")
        
        
        for (p in 1:npeels)  {
          
          tmp.mat[p,1] <- p
          tmp.mat[p,2] <- sapply(mod$peels[[p]]$rep$selAA[ind.sel.num[i]], function(x)  x[(nyears-p),]) [a,1]
          tmp.mat[p,3] <- (tmp.mat[p,2]-peel0[a,(1+i)]) 
          tmp.mat[p,4] <- (tmp.mat[p,2]-peel0[a,(1+i)])/peel0[a,(1+i)]
          
          
        } # end loop over p
        
        write.csv(tmp.mat, file=file.path(od, paste0(save.name,"_Rel.Diffs_" ,param, "_", i, "_age_", a, ".csv") ), row.names=F)
        
        rho.mat[a,i] <- mean(tmp.mat[,4])
        
      } # end loop over a
      
    }# end loop over i (number of indices)
    
    write.csv(rho.mat, file=file.path(od, paste0(save.name,"_Rho_", param,  ".csv") ), row.names=F)
    
  } # end test for index
  
  
  
  # FLEET SELECTIVITY ====
  if(param =="Fleet.Sel")  {
    
    rho.mat <- matrix(NA, nrow=nages, ncol=1)
    colnames(rho.mat) <- param
    
    for (a in 1:nages)  {
      tmp.mat <- matrix(NA, nrow=npeels, ncol=4 )
      colnames(tmp.mat) <- c("Peel", param, "Diff", "Rel.diff")
      
      
      for (p in 1:npeels)  {
        
        tmp.mat[p,1] <- p
        tmp.mat[p,2] <- mod$peels[[p]]$rep$selAA[[fleet.pointer.yrs[(nyears-p)] ]] [(nyears-p),a]
        
        tmp.mat[p,3] <- (tmp.mat[p,2]-peel0[(nyears-p),(a+1)]) 
        tmp.mat[p,4] <- (tmp.mat[p,2]-peel0[(nyears-p),(a+1)])/peel0[(nyears-p),(a+1)]
        
        
      } # end loop over p
      
      write.csv(tmp.mat, file=file.path(od, paste0(save.name,"_Rel.Diffs_" ,param,  "_age_", a, ".csv") ), row.names=F)
      
      rho.mat[a,1] <- mean(tmp.mat[,4])
      
      
    } # end loop over a
    
    write.csv(rho.mat, file=file.path(od, paste0(save.name,"_Rho_", param,  ".csv") ), row.names=F)
    
    
  } # end test for Fleet.Sel
  
  
  
  
  # NAA DEVS ====
  if(param =="NAA.devs")  {
    
    rho.mat <- matrix(NA, nrow=nages, ncol=1)
    colnames(rho.mat) <- param
    nages.devs <- ifelse(n.sigma==1, 1, nages)
    
    for (a in 1:nages.devs)  {
      tmp.mat <- matrix(NA, nrow=npeels, ncol=4 )
      colnames(tmp.mat) <- c("Peel", param, "Diff", "Rel.diff")
      
      
      for (p in 1:npeels)  {
        
        tmp.mat[p,1] <- p
        tmp.mat[p,2] <- mod$peels[[p]]$rep$NAA_devs[(nyears-p-1),a]
        
        tmp.mat[p,3] <- (tmp.mat[p,2]-peel0[(nyears-p-1),(a+1)]) 
        tmp.mat[p,4] <- (tmp.mat[p,2]-peel0[(nyears-p-1),(a+1)])/peel0[(nyears-p-1),(a+1)]
        
        
      } # end loop over p
      
      write.csv(tmp.mat, file=file.path(od, paste0(save.name,"_Rel.Diffs_" ,param,  "_age_", a, ".csv") ), row.names=F)
      
      rho.mat[a,1] <- mean(tmp.mat[,4])
      
      
      
    } # end loop over a
    
    write.csv(rho.mat, file=file.path(od, paste0(save.name,"_Rho_", param, ".csv") ), row.names=F)
    
    
  } # end test for NAA DEVS
  
  
  
  
  # MAA ====
  if(param =="MAA")  {
    
    rho.mat <- matrix(NA, nrow=nages, ncol=1)
    colnames(rho.mat) <- param
    
    for (a in 1:nages)  {
      tmp.mat <- matrix(NA, nrow=npeels, ncol=4 )
      colnames(tmp.mat) <- c("Peel", param, "Diff", "Rel.diff")
      
      
      for (p in 1:npeels)  {
        
        tmp.mat[p,1] <- p
        tmp.mat[p,2] <- mod$peels[[p]]$rep$MAA[(nyears-p),a]
        
        tmp.mat[p,3] <- (tmp.mat[p,2]-peel0[(nyears-p),(a+1)]) 
        tmp.mat[p,4] <- (tmp.mat[p,2]-peel0[(nyears-p),(a+1)])/peel0[(nyears-p),(a+1)]
        
        
      } # end loop over p
      
      write.csv(tmp.mat, file=file.path(od, paste0(save.name,"_Rel.Diffs_" ,param,  "_age_", a, ".csv") ), row.names=F)
      
      rho.mat[a,1] <- mean(tmp.mat[,4])
      
      
    } # end loop over a
    
    write.csv(rho.mat, file=file.path(od, paste0(save.name,"_Rho_", param, ".csv") ), row.names=F)
    
    
  } # end test for MAA 
  
  
  
  
  
  return(rho.mat)
  
}  #end function calc.rho.diffs


# function plot.retro.pars ====

# plot generated (2 each, one is full time series, second is shorter time series starting at yr.plot.start):
# M by peel  (could also make this plot by age?)
# q by peel
# index sel by peel
# fleet sel by peel (needs to work w or w/o RE)
# RE in NAA by peel (will have to do by age)




plot.retro.pars <- function (m, yr.plot.start, od, save.name, plot.f='png', ...){
  # m is a wham fitted model
  # yr.plot.start is the first year plotted for the "short" plots ("full" plots use all years)
  # od is the output directory
  # save.name is prepended to the pdf (comprising all plots) and individual plot names
  # plot.f is the graphics format (png, jpg, etc.)
  
  # Note: not tested with >1 fleet
  
if (dir.exists(od)==F)  dir.create(od)
  
npeels <- length(m$peels)  
nyears <- m$input$data$n_years_model
nindices <-  m$input$data$n_indices
nages <- m$input$data$n_ages

y1 <- m$input$data$year1_model 
years <- seq(y1, (y1+nyears-1))

q.peel <- cbind(years, m$rep$q, rep(0, nyears))
colnames(q.peel) <- c("Year", paste0(rep("q", nindices), seq(1, nindices) ) , "Peel")
q.peel0 <- q.peel

n.sel.blks <- m$input$data$n_selblocks

ind.sel.num <- m$input$data$selblock_pointer_indices[nyears,]
last.sel.row <- sapply(m$rep$selAA[ind.sel.num], function(x)  x[nyears,])
ind.sel.peel <- cbind(seq(1,nages), last.sel.row,  rep(0, nages))
colnames(ind.sel.peel) <- c("Age", paste0(rep("Ind_", nindices), seq(1, nindices) ) , "Peel")
ind.sel.peel0 <- ind.sel.peel

fleet.sel.num <- m$input$data$selblock_pointer_fleet[nyears,]
fleet.block.yrs <- m$input$data$selblock_years[, c(1:fleet.sel.num)]

if (fleet.sel.num==1) tmp.sel <- m$rep$selAA[[1]]
  
if (fleet.sel.num>1) {
for (f in 1:fleet.sel.num) {
  tmp.sel.mat <- m$rep$selAA[[f]]
  if(f==1) tmp.sel <- tmp.sel.mat[(fleet.block.yrs[,f]==1), ]
  if (f>1)  tmp.sel <- rbind(tmp.sel, tmp.sel.mat[(fleet.block.yrs[,f]==1), ])
  }
}

fleet.sel <- cbind(years, tmp.sel, rep(0, nyears))
colnames(fleet.sel) <- c("Year", paste0(rep("Age_", nages), seq(1, nages) ) , "Peel")
fleet.sel.peel0 <- fleet.sel
                   

naa.devs <- cbind(years[2:nyears], m$rep$NAA_devs, rep(0, (nyears-1)))
colnames(naa.devs) <- c("Year", paste0(rep("Age_", nages), seq(1, nages) ) , "Peel")
naa.devs.peel0 <- naa.devs

maa <- cbind(years, m$rep$MAA, rep(0, nyears))
colnames(maa) <- c("Year", paste0(rep("Age_", nages), seq(1, nages) ) , "Peel")
maa.peel0 <- maa



for (p in 1:npeels)  {

  tmp.nyrs <- nyears-p 
  tmp.yrs <-  seq(y1, (y1+tmp.nyrs-1))
 q.peel <- rbind(q.peel, cbind(tmp.yrs, m$peels[[p]]$rep$q, rep(p, tmp.nyrs))  )
 
 tmp.ind.sel.peel <- sapply(m$peels[[p]]$rep$selAA[ind.sel.num], function(x)  x[tmp.nyrs,])
 ind.sel.peel <- rbind(ind.sel.peel, cbind(seq(1,nages), tmp.ind.sel.peel,  rep(p, nages)))
 
 naa.devs <- rbind(naa.devs, cbind(years[2:tmp.nyrs], m$peels[[p]]$rep$NAA_devs, rep(p, (tmp.nyrs-1) ) ))
 
 maa <- rbind(maa , cbind(tmp.yrs, m$peels[[p]]$rep$MAA, rep(p, tmp.nyrs)))
 
 if (fleet.sel.num==1)  tmp.sel <- m$peels[[p]]$rep$selAA[[1]]
   
 
   
 if (fleet.sel.num>1) { 
 for (f in 1:fleet.sel.num) {
   tmp.sel.mat <- m$peels[[p]]$rep$selAA[[f]]
   tmp.fleet.block.yrs <- fleet.block.yrs[(1:(nyears-p)),]
   if(f==1) tmp.sel <- tmp.sel.mat[(tmp.fleet.block.yrs[,f]==1), ]
   if (f>1)  tmp.sel <- rbind(tmp.sel, tmp.sel.mat[(tmp.fleet.block.yrs[,f]==1), ])
 }
}
 
 fleet.sel <- rbind(fleet.sel, cbind(tmp.yrs, tmp.sel, rep(p, tmp.nyrs))  )
 
  
}# end p-loop over npeels
  

# check if NAA devs were estimated and if so, whether input$data$n_NAA_sigma is <2  ====
if(length(which(m$input$random %in% "log_NAA"))==0)  naa.devs[, (2:(2+nages-1))] <- NA  # no NAA devs
if(length(which(m$input$random %in% "log_NAA"))>0) {
  if(m$input$data$n_NAA_sigma <2)  naa.devs[ , (3:(3+nages-2))] <- NA   #only recr Devs
}

# set up tibbles ====

maa.tib <- as_tibble(maa) %>%
  pivot_longer(cols=starts_with("Age_"), names_prefix="Age_", names_to="Age", values_to="M") %>%
  mutate(Peel=as.factor(Peel))


naa.devs.tib <- as_tibble(naa.devs) %>%
  pivot_longer(cols=starts_with("Age_"), names_prefix="Age_", names_to="Age", values_to="Dev") %>%
  mutate(Peel=as.factor(Peel))


ind.sel.tib <- as_tibble(ind.sel.peel) %>%
  pivot_longer(cols=starts_with("Ind_"), names_prefix="Ind_", names_to="Index", values_to="Index_Selectivity") %>%
  mutate(Peel=as.factor(Peel))

q.peel.tib <- as_tibble(q.peel) %>% 
  pivot_longer(cols=starts_with("q"), names_prefix="q", names_to="Index", values_to="Catchability") %>%
  mutate(Peel=as.factor(Peel))
 
fleet.sel.tib <- as_tibble(fleet.sel) %>%
  pivot_longer(cols=starts_with("Age_"), names_prefix="Age_", names_to="Age", values_to="Fleet_Selectivity") %>%
  mutate(Peel=as.factor(Peel))



# make plots ====  

## get facet labels (mohn's rho values) ====
q.rho <- calc.rho.diffs(mod=m, param="q", peel0=q.peel0, od=od, save.name=save.name)

q.labels <- data.frame(Index= as.character(seq(1, nindices) ), 
                       rho=as.character(round(q.rho,3))  ) 
 
q.peel.tib <- q.peel.tib %>%
  left_join(q.labels) %>%
  mutate(facet_label = paste0(Index, " (rho = ", rho, ")"))
# can't seem to get rho to show up as greek letter -- failed attempts in 2 lines below
# mutate(facet_label = paste0(Index, " (", expression(~rho) , " = ", rho, ")"))
# mutate(facet_label = expression(paste0(Index, " (", rho , " = ", rho, ")" )))


index.rho <- calc.rho.diffs(mod=m, param="index", peel0=ind.sel.peel0, od=od, save.name=save.name)
index.labels <- data.frame(Index=as.character( rep(seq(1, nindices), each=nages) ), 
                           rho=as.character(round(index.rho,3)), 
                           Age= as.double(rep(seq(1, nages), nindices) ) )



fleet.sel.rho <- calc.rho.diffs(mod=m, param="Fleet.Sel", peel0=fleet.sel.peel0, od=od, save.name=save.name)
fleet.sel.labels <- data.frame(rho = as.character(round(fleet.sel.rho,3)),
                               Age= as.character( seq(1, nages) ))
fleet.sel.tib <- fleet.sel.tib %>%
  left_join(fleet.sel.labels)  %>%
  mutate(facet_label =  paste0(Age, " ( rho" , " = ", rho, ")"  ) )


naa.dev.rho <- calc.rho.diffs(mod=m, param="NAA.devs", peel0=naa.devs.peel0, od=od, save.name=save.name)
naa.dev.labels <- data.frame(rho = as.character(round(naa.dev.rho,3)),
                             Age= as.character( seq(1, nages) ))
                               
naa.devs.tib <- naa.devs.tib %>%
  left_join(naa.dev.labels)  %>%
  mutate(facet_label =  paste0(Age, " ( rho" , " = ", rho, ")"  ) )


maa.rho <- calc.rho.diffs(mod=m, param="MAA", peel0=maa.peel0, od=od, save.name=save.name)
maa.labels <- data.frame(rho = as.character(round(maa.rho,3)),
                             Age= as.character( seq(1, nages) ))

maa.tib <- maa.tib %>%
  left_join(maa.labels)  %>%
  mutate(facet_label =  paste0(Age, " ( rho" , " = ", rho, ")"  ) )

## plots ====


fleet.sel.full.plot <- ggplot(fleet.sel.tib, aes(x=Year, y=Fleet_Selectivity, col=Peel)) +
  facet_wrap(~facet_label, nrow=3) +
  geom_line() +
  geom_point()
ggsave(fleet.sel.full.plot, filename=here(od,paste(save.name,"_fleet.sel.full.plot", plot.f, sep=".")), 
       device=plot.f,  height=8, width=8) 


fleet.sel.short.plot <- ggplot(fleet.sel.tib, aes(x=Year, y=Fleet_Selectivity, col=Peel)) +
  facet_wrap(~facet_label, nrow=3) +
  geom_line() +
  geom_point() +
  xlim(yr.plot.start, years[nyears])
ggsave(fleet.sel.short.plot, filename=here(od,paste(save.name,"_fleet.sel.short.plot", plot.f, sep=".")), 
       device=plot.f,  height=8, width=8) 


maa.full.plot <- ggplot(maa.tib, aes(x=Year, y=M, col=Peel)) +
  facet_wrap(~facet_label, nrow=3) +
  geom_line() +
  geom_point()
ggsave(maa.full.plot, filename=here(od,paste(save.name,"_maa.full.plot", plot.f, sep=".")), 
       device=plot.f,  height=8, width=8) 

maa.short.plot <- ggplot(maa.tib, aes(x=Year, y=M, col=Peel)) +
  facet_wrap(~facet_label, nrow=3) +
  geom_line() +
  geom_point() +
  xlim(yr.plot.start, years[nyears])
ggsave(maa.short.plot, filename=here(od,paste(save.name,"_maa.short.plot", plot.f, sep=".")), 
       device=plot.f,  height=8, width=8) 


if (length(which(m$input$random %in% "log_NAA"))>0 ) {
naa.devs.full.plot <- ggplot(naa.devs.tib, aes(x=Year, y=Dev, col=Peel)) +
  facet_wrap(~facet_label, nrow=3, scales="free_y") +
  geom_hline(yintercept=0, col='grey25', linetype=2) +
  geom_line() +
  geom_point() +
  ylab("NAA log_Deviation")
ggsave(naa.devs.full.plot, filename=here(od,paste(save.name,"_naa.devs.full.plot", plot.f, sep=".")), 
       device=plot.f,  height=8, width=8) 


naa.devs.short.plot <- ggplot(naa.devs.tib, aes(x=Year, y=Dev, col=Peel)) +
  facet_wrap(~facet_label, nrow=3, scales="free_y") +
  geom_hline(yintercept=0, col='grey25', linetype=2) +
  geom_line() +
  geom_point() +
  xlim(yr.plot.start, years[nyears]) +
  ylab("NAA log_Deviation")
ggsave(naa.devs.short.plot, filename=here(od,paste(save.name,"_naa.devs.short.plot", plot.f, sep=".")), 
       device=plot.f,  height=8, width=8) 

}


ind.sel.full.plot <- ggplot(ind.sel.tib, aes(x=Age, y=Index_Selectivity, col=Peel)) +
  facet_wrap(~Index, nrow=3) +
  geom_hline(yintercept=1, col='grey25', linetype=2) +
  geom_line() +
  geom_point() +
  ylim(0, 1.25) +
  geom_text(data=index.labels, x=index.labels$Age, y=1.07, aes(label=rho),  col='black')  +
  geom_text(data=index.labels, x=as.double(m$ages.lab[1]), y=1.17, aes(label="rho="),  col='black')
ggsave(ind.sel.full.plot, filename=here(od,paste(save.name,"_ind.sel.full.plot", plot.f, sep=".")), 
       device=plot.f,  height=8, width=8) 


q.short.plot <- ggplot(q.peel.tib, aes(x=Year, y=Catchability, col=Peel)) +
  facet_wrap(~facet_label, nrow=3) +
  geom_line() +
  geom_point() +
  xlim(yr.plot.start, years[nyears])  
ggsave(q.short.plot, filename=here(od,paste(save.name,"_q.short.plot", plot.f, sep=".")), 
       device=plot.f,  height=8, width=8) 


q.full.plot <- ggplot(q.peel.tib, aes(x=Year, y=Catchability, col=Peel)) +
  facet_wrap(~facet_label, nrow=3) +
  geom_line() +
  geom_point() 
ggsave(q.full.plot, filename=here(od,paste( save.name,"_q.full.plot", plot.f, sep=".")), 
       device=plot.f,  height=8, width=8) 



# plot pdf ====
grDevices::cairo_pdf(filename = here(od,paste0(save.name,"_Retro_Est_Peels.pdf") ), 
                     family = "sans", height = 10, width = 10, onefile = TRUE)

print(fleet.sel.full.plot)
print(fleet.sel.short.plot)
print(ind.sel.full.plot)
print(q.full.plot)
print(q.short.plot)
if (length(which(m$input$random %in% "log_NAA"))>0 ) {
print(naa.devs.full.plot)
print(naa.devs.short.plot)
}
print(maa.full.plot)
print(maa.short.plot)

dev.off()
} # end plot.retro.pars

