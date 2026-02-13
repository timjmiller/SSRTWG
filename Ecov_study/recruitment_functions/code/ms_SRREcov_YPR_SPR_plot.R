
# plot selectivity to see how it looks  ====
library(dplyr)
library(ggplot2)
library(reshape2)
library(tibble) 
library(tidyr)
library(patchwork)  

scriptdir <- dirname(rstudioapi::getSourceEditorContext()$path) # get path where this script is
source(file.path(scriptdir, 'fn_msy_ypr_ssbpr_calcs.R'))
setwd(scriptdir)
setwd('../..') #back up to SSRTWG main directory
getwd()

plot.dir <- file.path(here::here(), 'Ecov_study', 'recruitment_functions', 'manuscript_plots')
if(!dir.exists(plot.dir))  dir.create(plot.dir, recursive=TRUE)

library(here) # load this after establishing root directory


# bio params ====
ages <- seq(1,10)
nages <- length(ages)
naa_om_inputs <- readRDS(file.path(here::here(),"Project_0","inputs", "NAA_om_inputs.RDS"))
waa.fleet <- naa_om_inputs[[24]]$data$waa[4,1,]
waa.ssb <- naa_om_inputs[[24]]$data$waa[5,1,]  #looks like these are the same

# stock recruit params ====
R0 = exp(10)
steep = 0.69
ahat = 4*steep/(1-steep)
spr0 <- s.per.recr(nages=nages, fec.age=waa.ssb, mat.age=maturity, mm= 0.2, 
                                  F.mult=0, sel.age=rep(0,nages), spawn.time=0.25) 


sel.orig.list <-  list('Fleet-1' = c(5,1), 'Index-1' = c(5, 1), 'Index-2' = c(5, 1)  ) #fleet, index1, index2
sel.orig <- lapply(sel.orig.list, function(x) {
  tmp=logistic.fn(a50=x[1], k=x[2], ages=ages)
  tmp2 = tmp/max(tmp)
  tmp2
  }
  )
 
sel.new.sel.list <-list('Fleet-1' = c(2.89, 0.88), 'Index-1' = c(0.8, 3), 'Index-2' = c(3.5, 0.7)  ) #fleet, index1, index2
sel.new <- lapply(sel.new.sel.list, function(x) {
  tmp=logistic.fn(a50=x[1], k=x[2], ages=ages)
  tmp2 = tmp/max(tmp)
  tmp2
}
)

maturity.list <- list('Maturity' = c(2.89, 0.88))
maturity1 <- logistic.fn(a50=maturity.list$Maturity[1], k=maturity.list$Maturity[2], ages=ages)
maturity <- maturity1/max(maturity1)
mat.tib <- as_tibble(cbind(Maturity=maturity, Age=ages) )

sel.orig.tib <- as_tibble(melt(sel.orig)) %>%
  rename(Selectivity=value, Source=L1) %>%
  mutate(Case='Original')
sel.new.tib <- as_tibble(melt(sel.new))   %>%
          rename(Selectivity=value, Source=L1) %>%
  mutate(Case='Sensitivity')

sel.tib <- sel.orig.tib %>%
  add_row(sel.new.tib) %>%
  mutate(Age=rep(ages,6))

# get msy (run 3 times to check solution) ====  none of these gives 0.348 for Fmsy, they all give 0.379)
msy.orig1 <- as.list( get.MSY(nages=nages, steep=0.69, R0=exp(10), fec.age=waa.ssb, mat.age=maturity, M.age=rep(0.2,10), 
                    fleet.waa=waa.fleet, sel.age=sel.tib$Selectivity[sel.tib$Source=='Fleet-1' & sel.tib$Case=='Original'], 
                    spawn.time=0.25,
                    F.start=0.1) )

msy.orig2 <- get.MSY(nages=nages, steep=0.69, R0=exp(10), fec.age=waa.ssb, mat.age=maturity, M.age=rep(0.2,10), 
                     fleet.waa=waa.fleet, sel.age=sel.tib$Selectivity[sel.tib$Source=='Fleet-1' & sel.tib$Case=='Original'], 
                     spawn.time=0.25,
                     F.start=0.5) 
msy.orig3 <- get.MSY(nages=nages, steep=0.69, R0=exp(10), fec.age=waa.ssb, mat.age=maturity, M.age=rep(0.2,10), 
                     fleet.waa=waa.fleet, sel.age=sel.tib$Selectivity[sel.tib$Source=='Fleet-1' & sel.tib$Case=='Original'], 
                     spawn.time=0.25,
                     F.start=0.348) 
# ----> all 3 solutions agree

ypr.half.msy.orig <- ypr(nages=nages, wgt.age=waa.fleet, mm=0.2, F.mult=0.5*msy.orig1$Fmsy, 
                    sel.age=sel.tib$Selectivity[sel.tib$Source=='Fleet-1' & sel.tib$Case=='Original'])
ypr.2msy.orig <- ypr(nages=nages, wgt.age=waa.fleet, mm=0.2, F.mult=2.0*msy.orig1$Fmsy, 
                sel.age=sel.tib$Selectivity[sel.tib$Source=='Fleet-1' & sel.tib$Case=='Original'])

ssbpr.half.msy.orig <- s.per.recr(nages=nages, fec.age=waa.ssb, mat.age=maturity, mm= 0.2, 
                                  F.mult=0.5*msy.orig1$Fmsy, sel.age=sel.tib$Selectivity[sel.tib$Source=='Fleet-1' & sel.tib$Case=='Original'], 
                                  spawn.time=0.25) 
ssbpr.2msy.orig <- s.per.recr(nages=nages, fec.age=waa.ssb, mat.age=maturity, mm= 0.2, 
                                  F.mult=2*msy.orig1$Fmsy, sel.age=sel.tib$Selectivity[sel.tib$Source=='Fleet-1' & sel.tib$Case=='Original'], 
                              spawn.time=0.25) 
ssbeq.half.msy.orig <- ssb.eq( alpha.BH=ahat, R0.BH=R0, spr=ssbpr.half.msy.orig, spr0=spr0 )
ssbeq.2msy.orig <- ssb.eq( alpha.BH=ahat, R0.BH=R0, spr=ssbpr.2msy.orig, spr0=spr0 )



# get msy for sensitivity (run 3 times to check solution) ====  none of these gives 0.348 for Fmsy, they all give 0.379)
msy.new1 <- as.list( get.MSY(nages=nages, steep=0.69, R0=exp(10), fec.age=waa.ssb, mat.age=maturity, M.age=rep(0.2,10), 
                     fleet.waa=waa.fleet, sel.age=sel.tib$Selectivity[sel.tib$Source=='Fleet-1' & sel.tib$Case=='Sensitivity'], 
                     spawn.time=0.25,
                     F.start=0.1)  )

msy.new2 <- get.MSY(nages=nages, steep=0.69, R0=exp(10), fec.age=waa.ssb, mat.age=maturity, M.age=rep(0.2,10), 
                     fleet.waa=waa.fleet, sel.age=sel.tib$Selectivity[sel.tib$Source=='Fleet-1' & sel.tib$Case=='Sensitivity'], 
                     spawn.time=0.25,
                     F.start=0.5) 
msy.new3 <- get.MSY(nages=nages, steep=0.69, R0=exp(10), fec.age=waa.ssb, mat.age=maturity, M.age=rep(0.2,10), 
                     fleet.waa=waa.fleet, sel.age=sel.tib$Selectivity[sel.tib$Source=='Fleet-1' & sel.tib$Case=='Sensitivity'], 
                     spawn.time=0.25,
                     F.start=0.348) 
# ----> all 3 solutions agree

ypr.half.msy.new <- ypr(nages=nages, wgt.age=waa.fleet, mm=0.2, F.mult=0.5*msy.new1$Fmsy, 
                        sel.age=sel.tib$Selectivity[sel.tib$Source=='Fleet-1' & sel.tib$Case=='Sensitivity']) 
ypr.2msy.new <- ypr(nages=nages, wgt.age=waa.fleet, mm=0.2, F.mult=2*msy.new1$Fmsy, 
                        sel.age=sel.tib$Selectivity[sel.tib$Source=='Fleet-1' & sel.tib$Case=='Sensitivity']) 
ssbpr.half.msy.new <- s.per.recr(nages=nages, fec.age=waa.ssb, mat.age=maturity, mm= 0.2, 
                                  F.mult=0.5*msy.new1$Fmsy, sel.age=sel.tib$Selectivity[sel.tib$Source=='Fleet-1' & sel.tib$Case=='Sensitivity'], 
                                  spawn.time=0.25) 
ssbpr.2msy.new <- s.per.recr(nages=nages, fec.age=waa.ssb, mat.age=maturity, mm= 0.2, 
                              F.mult=2*msy.new1$Fmsy, sel.age=sel.tib$Selectivity[sel.tib$Source=='Fleet-1' & sel.tib$Case=='Sensitivity'], 
                              spawn.time=0.25) 
ssbeq.half.msy.new <- ssb.eq( alpha.BH=ahat, R0.BH=R0, spr=ssbpr.half.msy.new, spr0=spr0 )
ssbeq.2msy.new <- ssb.eq( alpha.BH=ahat, R0.BH=R0, spr=ssbpr.2msy.new, spr0=spr0 )



# calculate ypr ====
F.seq <- seq(0, 1.5, by=0.01)
ypr.tab.orig <-as_tibble(cbind(F.seq=F.seq, YPR.original=unlist(lapply( F.seq, function(x)  ypr(nages=nages, wgt.age=waa.fleet, mm=0.2, F.mult=x, 
                                              sel.age=sel.tib$Selectivity[sel.tib$Source=='Fleet-1' & sel.tib$Case=='Original'])  ) ) )  )

ypr.tab.new <- as_tibble(cbind(F.seq=F.seq, YPR.maturity=unlist(lapply( F.seq, function(x)  ypr(nages=nages, wgt.age=waa.fleet, mm=0.2, F.mult=x, 
                                            sel.age=sel.tib$Selectivity[sel.tib$Source=='Fleet-1' & sel.tib$Case=='Sensitivity'])  )  ) )  )
fmax.new.row <- which(ypr.tab.new$YPR.maturity==max(ypr.tab.new$YPR.maturity)) 
fmax.orig.row <- which(ypr.tab.orig$YPR.original==max(ypr.tab.orig$YPR.original)) 
f.max.new <- F.seq[fmax.new.row]  #0.4
f.max.orig <- F.seq[fmax.orig.row] #0.85


ypr.tibble <- (ypr.tab.orig) %>%
  left_join(ypr.tab.new ) %>%
  pivot_longer(cols=starts_with("YPR"), names_prefix="YPR.", names_to="Selectivity" , values_to="YPR")


# ypr plot ====
ypr.plot <- ggplot(ypr.tibble, aes(x=F.seq, y=YPR, col=Selectivity, fill=Selectivity)) +
    geom_line( linewidth=1) +
    #Fmax
    annotate("point", x = F.seq[fmax.orig.row], y = ypr.tab.orig$YPR.original[fmax.orig.row], colour = "#2233ff", size = 4, shape=2) +
    annotate("point", x = F.seq[fmax.new.row], y = ypr.tab.new$YPR.maturity[fmax.new.row], colour = "#2233ff", size = 4, shape=2) +
    #Fmsy
    annotate("point", x = msy.orig1$Fmsy, y = msy.orig1$YPRmsy, colour = "#787171", size = 4, shape='$') +
    annotate("point", x = msy.new1$Fmsy, y = msy.new1$YPRmsy, colour = "#BF13BF", size = 4, shape='$') +
    #0.5*Fmsy
    annotate("point", x = 0.5*msy.orig1$Fmsy, y = ypr.half.msy.orig, colour = "#787171", size = 2.5, shape=16) +
    annotate("point", x = 0.5*msy.new1$Fmsy, y = ypr.half.msy.new, colour = "#BF13BF", size = 2.5, shape=16) +
    #2*Fmsy
    annotate("point", x = 2.0*msy.orig1$Fmsy, y = ypr.2msy.orig, colour = "#787171", size = 2.5, shape=16) +
    annotate("point", x = 2.0*msy.new1$Fmsy, y = ypr.2msy.new, colour = "#BF13BF", size = 2.5, shape=16) +
    
    xlab('F') +
    theme_light()  +
    theme(strip.background =element_rect(fill="white", color="grey65"))+
    theme(strip.text = element_text(colour = 'black', size=12)) +
    theme(axis.text.x = element_text(size = 12))   +
    theme(axis.text.y = element_text(size = 13))  +
    theme(axis.title.y = element_text(size = 14))  +
    theme(axis.title.x = element_text(size = 14)) +
    theme(legend.position = "bottom") +
    theme(legend.title = element_blank())  +
    scale_fill_manual(values=c("#BF13BF66", "#78717166")) +
    scale_color_manual(values=c("#BF13BF66", "#78717166"))
  ggsave(ypr.plot, filename = file.path(plot.dir, "ypr.plot.png"),
         height=6, width=9) 
  
  
# SRR plot with equilibrium (deterministic) SSB range ====
SSB0 <- spr0*R0  
ssb.seq <- seq(0,SSB0, length=3000 )
r.seq <- bev.holt.alpha(S=ssb.seq,R0=R0,recr.par=ahat,spr0=spr0, is.steepness=FALSE)
srr.tib <- as_tibble(cbind(SSB=ssb.seq, R=r.seq))

srr.plot <- ggplot(srr.tib, aes(x=SSB, y=R)) +
    geom_line( linewidth=1, color='black') +

    #Fmsy
    annotate("point", x = msy.orig1$SSBmsy, y = msy.orig1$Rmsy, colour = "#787171", size = 4, shape='$') +
    annotate("point", x = msy.new1$SSBmsy, y = msy.new1$Rmsy, colour = "#BF13BF", size = 4, shape='$') +
    #0.5*Fmsy
    annotate("point", x = ssbeq.half.msy.orig, y = bev.holt.alpha(S=ssbeq.half.msy.orig,R0=R0,recr.par=ahat,spr0=spr0, is.steepness=FALSE), 
             colour = "#787171", size = 3.5, shape=16) +
    annotate("point", x = ssbeq.half.msy.new, y = bev.holt.alpha(S=ssbeq.half.msy.new,R0=R0,recr.par=ahat,spr0=spr0, is.steepness=FALSE), 
             colour = "#BF13BF", size = 3.5, shape=16) +
    #2*Fmsy
    annotate("point", x = ssbeq.2msy.orig, y = bev.holt.alpha(S=ssbeq.2msy.orig,R0=R0,recr.par=ahat,spr0=spr0, is.steepness=FALSE), 
             colour = "#787171", size = 3.5, shape=16) +
    annotate("point", x = ssbeq.2msy.new, y = bev.holt.alpha(S=ssbeq.2msy.new,R0=R0,recr.par=ahat,spr0=spr0, is.steepness=FALSE), 
             colour = "#BF13BF", size = 3.5, shape=16) +
    
    theme_light()  +
    theme(strip.background =element_rect(fill="white", color="grey65"))+
    theme(strip.text = element_text(colour = 'black', size=12)) +
    theme(axis.text.x = element_text(size = 12))   +
    theme(axis.text.y = element_text(size = 13))  +
    theme(axis.title.y = element_text(size = 14))  +
    theme(axis.title.x = element_text(size = 14)) +
    theme(legend.position = "bottom") +
    theme(legend.title = element_blank())  +
    scale_fill_manual(values=c("#BF13BF66", "#78717166")) +
    scale_color_manual(values=c("#BF13BF66", "#78717166"))
  ggsave(srr.plot, filename = file.path(plot.dir, "srr.plot.png"),
         height=6, width=9) 
  
  

  fig_ypr_srr <- (ypr.plot / srr.plot)  +
    plot_annotation(tag_levels = list(c('(a)', '(b)' )) )
  ggsave(fig_ypr_srr, filename=file.path(plot.dir, "fig_ypr_srr.pdf"),  
         units="in", height=9, width=7, device = "pdf", dpi=300)
  
