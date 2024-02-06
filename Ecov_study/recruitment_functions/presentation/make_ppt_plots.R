library(patchwork)
library(tidyverse)
library(wham)
naa_om_inputs <- readRDS(file.path(here::here(),"Project_0","inputs", "NAA_om_inputs.RDS"))
om.data <- naa_om_inputs[[1]]$data
mat.waa <- as_tibble(cbind(Age=seq(1, om.data$n_ages), WAA=om.data$waa[1,1,], Mat = om.data$mature[1,] ))
temp <- fit_wham(naa_om_inputs[[1]], do.fit = FALSE, MakeADFun.silent = TRUE)
SRab <- exp(c(temp$rep$log_SR_a[1], temp$rep$log_SR_b[1]))
srr <- as_tibble(cbind(SSB=seq(0,3.5e5, length=1000), Recruits= (SRab[1]*seq(0,3.5e5, length=1000)/(1+SRab[2]*seq(0,3.5e5, length=1000)) )  ))

mat.plot <- ggplot(mat.waa, aes(x=Age, y=Mat)) +
  geom_line(linewidth=1) +
  geom_point(size=2.5) +
  ylab('Maturity at age') +
  scale_x_continuous(breaks=mat.waa$Age, labels=as.character(mat.waa$Age))+
  theme(axis.text.x = element_text(size = 16))   + 
  theme(axis.text.y = element_text(size = 16)) +
  theme(axis.title.x = element_text(size = 16))   + 
  theme(axis.title.y = element_text(size = 16))  
ggsave(mat.plot, filename=file.path(here::here(), 'Ecov_study','recruitment_functions', 'presentation','plots', 
                                    'mat.plot.png'), height=3, width=4)
waa.plot <- ggplot(mat.waa, aes(x=Age, y=WAA)) +
  geom_line(linewidth=1, color='blue') +
  geom_point(size=2.5, color='blue') +
  ylab('WAA (SSB and Catch)') +
  scale_x_continuous(breaks=mat.waa$Age, labels=as.character(mat.waa$Age))  +
  theme(axis.text.x = element_text(size = 16))   + 
  theme(axis.text.y = element_text(size = 16))  + 
  theme(axis.title.y = element_text(size = 16))  + 
  theme(axis.title.x = element_text(size = 16))  
ggsave(waa.plot, filename=file.path(here::here(), 'Ecov_study','recruitment_functions', 'presentation','plots', 
                                    'waa.plot.png'), height=3, width=4)

srr.plot <- ggplot(srr, aes(x=SSB, y=Recruits)) +
  geom_line(linewidth=1, color='red')  +
  theme(axis.text.x = element_text(size = 16))   + 
  theme(axis.text.y = element_text(size = 16))  + 
  theme(axis.title.y = element_text(size = 16))  + 
  theme(axis.title.x = element_text(size = 16))  
# geom_point(size=2.5, color='blue') +
# ylab('Weight at age (SSB and Catch)')
ggsave(srr.plot, filename=file.path(here::here(), 'Ecov_study','recruitment_functions', 'presentation','plots', 
                                    'srr.plot.png'), height=5, width=8)


bio.plot <- (mat.plot|waa.plot)/srr.plot
ggsave(bio.plot, filename=file.path(here::here(), 'Ecov_study','recruitment_functions', 'presentation','plots', 'bio.plot.png'),  height=6, width=9)
