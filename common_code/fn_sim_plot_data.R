# functions to simulate data from om_inputs 
#     and plot om specs, time series, and SR 
# liz brooks
# feb 21, 2023




########### simulate_data fn ====

simulate_data <- function(seed.rds=NULL, same.seeds=TRUE,
                          om.inputs.rds=NULL, df.oms.rds=NULL,
                          om.mod.nums=NULL, sim.nums=NULL,
                          simdata_dir=NULL,
                          write.sim.data=FALSE)   {  
 
  if (!dir.exists(simdata_dir))  dir.create(simdata_dir)
  
  df.oms <- readRDS(df.oms.rds)
  om_inputs <- readRDS(om.inputs.rds)
  
  # same.seeds==TRUE will use the same 1000 seeds for all oms
  seeds = readRDS(seed.rds)
  if (same.seeds==TRUE) {
      seeds <- rep(seeds[[1]], times= NROW(df.oms) )
      seeds <- lapply(1:NROW(df.oms), function(x) seeds[(1:1000) + 1000*(x-1)])
                   
  }
  
#specify number of data sets to generate and number of fits to make ====
om_spec <- om.mod.nums   #which oms to generate simulations (select from 1-288)
n_om <- length(om_spec)
sim_spec <- sim.nums     #which simulations to execute (select from 1-1000)
n_sim <- length(sim_spec)

# create objects to store data
om_fit <- list()
om_sim <- list()
sim_data <-list(list())
tmp <- list()

# generate data sets   ====
for (iom in 1:n_om)  {
  
  tmp.om <- om_spec[iom]
  om_fit[[iom]] <- fit_wham(om_inputs[[tmp.om]], do.fit = FALSE, MakeADFun.silent = TRUE)
  om_fit[[iom]]$Model <- df.oms$Model[tmp.om]
  
  sim_data[[iom]] <- df.oms$Model[om_spec[iom]]
  
  for (isim in 1:n_sim)  {
    
    tmp.sim <- sim_spec[isim]
    set.seed(seeds[[tmp.om]][tmp.sim])
    sim_data[[iom]][isim] <- list(om_fit[[iom]]$simulate(complete=TRUE))
    

  } # end isim loop

  names(sim_data[[iom]]) <- paste0(rep("sim", n_sim), sim_spec)
  
} # end iom loop  

  names(sim_data) <- df.oms$Model[om_spec]
  
  if (write.sim.data==T) saveRDS(sim_data, 
                            file=file.path(simdata_dir, paste0("sim_data_oms_", 
                                                               paste0(as.character(om.mod.nums), sep=".", collapse=""), "RDS")  ) )
  return(sim_data)
  
}  #end function simulate_data   ====





############# make plots for the simulated data ====

plot_sim_data <- function(ts.obj.names=NULL,   # vector with any of the above obj names (haven't tested any other output)
                          sim_data=NULL, simdata_dir=NULL,  #list of simulated data; directory for saving output
                          df.oms.rds=NULL,              #matrix with om_specs
                          pheight=6, pwidth=6, text.size=4)   {  
  
  
  df.oms <- readRDS(df.oms.rds)
  
  
  # plot the om specs (text plot) ====
  
  plot_om_specs(df.oms=df.oms, 
                sim_data=sim_data, 
                simdata_dir=simdata_dir,
                pheight=pheight, pwidth=pwidth, text.size=text.size)
  
  
  
  # plot the simulated time series objects ====
  vec.ts.obj.names <- ts.obj.names
  if (length(vec.ts.obj.names)>1)  sapply(vec.ts.obj.names, function(x) plot_ts_obj(ts.obj.names=x, 
                                                                  sim_data=sim_data, 
                                                                  simdata_dir=simdata_dir,
                                                                  pheight=pheight, pwidth=pwidth) )
  
  if (length(vec.ts.obj.names)==1)  plot_ts_obj(ts.obj.names=ts.obj.names, 
                                                sim_data=sim_data, 
                                                simdata_dir=simdata_dir,
                                                pheight=pheight, pwidth=pwidth)  # call to plot_ts_obj
  
    
    
    
    
  # plot the SR objects ====
  
  plot_SR(sim_data=sim_data, 
          simdata_dir=simdata_dir,
          pheight=pheight, pwidth=pwidth)
  
  
  
} #end function plot_sim_data   ====





############# plot the om specs (text plot) fn ====

plot_om_specs <- function(df.oms=NULL, sim_data=NULL, simdata_dir=NULL, 
                        pheight=6, pwidth=6, text.size=text.size) {   #ignoring these for now; width needs to be >height (6x10 look ok?)


  if (!dir.exists(simdata_dir))  dir.create(simdata_dir)
  
  n_om <- length(sim_data)
  om.nums <-  as.integer(sapply( names(sim_data), function(x)  substr(x, 4, nchar(x)) ) )
  prefix.plot <- paste0(om.nums, sep=".", collapse="")
  
  
  
  head.chars <- nchar(paste(names(df.oms)))
  head.char.place <- cumsum( (c(0, head.chars)+max(head.chars+7))   ) [1:length(head.chars)]
  temp.head <- paste(names(df.oms), collapse="    ")
  ncols.df.oms <- NCOL(df.oms)
  
  df.oms.nums <- as_tibble(df.oms[om.nums,])  %>%
    mutate_all(as.character )  %>%
    pivot_longer(cols=everything(), names_to="Factor", values_to="Level")  %>%
    mutate(xval=rep(head.char.place, n_om),
           yval=rev( rep(1*seq(1,n_om), each=length(head.chars) )  ) ,
           yval.factor = (rep(1*rep(n_om,n_om), each=length(head.chars)) +1 )  )
  
  
  
  om.text.size <- text.size
  om.width <- NCOL(df.oms) +3
  om.height <- round(0.5*NCOL(df.oms),0)
  
  ytop <- max(df.oms.nums$yval)+2
  xmax <- head.char.place[length(head.char.place)] + head.chars[length(head.chars)]
  
  om.spec.plot <- ggplot(df.oms.nums, aes(x=as.integer(xval), y=as.integer(yval), label=Level)) +
    geom_point(size=0.05, color="grey90")  +
    xlim(0, xmax) +
    ylim(0, ytop)  + 
    geom_text(hjust="left", size=om.text.size)  +  
    geom_point(data=df.oms.nums[df.oms.nums$yval==2, ] , aes(x=xval, y=yval.factor), size=0.05, color="grey90" )  +
    geom_text(data=df.oms.nums[df.oms.nums$yval==2, ] , 
              aes(x=(xval-2), y=yval.factor, label=Factor, hjust=0, fontface=2 ), size=om.text.size)  +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
          axis.text.y = element_blank(), axis.ticks.y = element_blank())  +
    ylab("") +
    xlab("")  +
    theme(panel.background = element_rect(fill="white", color="white" ))
  #theme_void()
  ggsave(om.spec.plot, filename = file.path(simdata_dir, paste0("om_",  prefix.plot, "_A.om.spec.plot", ".png")),
         height=6, width=10)
  
  
}  #end function plot_om_specs   ====   





############ plot_ts_obj fn ====

plot_ts_obj <- function(ts.obj.names=NULL, sim_data=NULL, simdata_dir=NULL ,
                        pheight=6, pwidth=6) {  
  
 # objects that can be plotted (that i have tested so far):
  #"agg_catch"  # "agg_indices"  "Ecov_x"  "Ecov_obs"  "Ecov_obs_sigma"  "SSB"  "F"
  
  if (!dir.exists(simdata_dir))  dir.create(simdata_dir)
  
  n_om <- length(sim_data)
  n_sim <- length(sim_data[[1]])
  
  om.nums <-  as.integer(sapply( names(sim_data), function(x)  substr(x, 4, nchar(x)) ) )
  sim.nums <-  as.integer(sapply( names(sim_data[[1]]), function(x)  substr(x, 4, nchar(x)) ) )
  
  
  obj.slot <- which((names(sim_data[[1]][[1]])==ts.obj.names))
  

  if( ts.obj.names %in% c("MAA", "NAA")  ) {
    if (ts.obj.names=="MAA" ) {
      ts.obj <-  lapply(sim_data, function(x) lapply(x, function(y)  get_recr(y, obj.slot) )   )  #this grabs first column
      ts.obj.names <- "M[1]"
    }
    if (ts.obj.names=="NAA") {
        ts.obj <-  lapply(sim_data, function(x) lapply(x, function(y)  get_total(y, obj.slot) )   )  #this grabs first column
        ts.obj.names <- "Total N"
      } # if NAA
  }   else {
ts.obj <-  lapply(sim_data, function(x) lapply(x, function(y)  y[obj.slot] )   ) 
  } #end if-else to check for MAA or NAA

ts.tib <-  as_tibble(lapply(ts.obj, function(x) as_tibble(unlist(x) )) )

n.obj <- dim(as.data.frame(ts.obj[[1]][1]) ) [2]
len.obj <- dim(as.data.frame(ts.obj[[1]][1]) ) [1]
years <- seq(1, len.obj)

ts.tib1 <- ts.tib %>%
  mutate(obj= rep(paste0(rep(ts.obj.names, n.obj*len.obj),  rep(seq(1, n.obj), each=len.obj)  ), n_sim) ) %>%
  mutate(sim = rep(rep(sim.nums, each=len.obj*n.obj)))  %>%
  mutate( Year = rep(years, n.obj*n_sim)) %>%
  pivot_longer(cols=starts_with("om"), names_to="OM", values_to="Value")  %>%
  mutate(obj=as.factor(obj), OM=as.factor(OM))  %>%
  mutate(obj.value = unlist(Value))

prefix.plot <- paste0(om.nums, sep=".", collapse="")

obj.plot <- ggplot(ts.tib1, 
                   aes(x=Year, y=obj.value, col=as.factor(sim) ))  +
  facet_grid( OM ~obj) +
  geom_line(alpha=0.45)  + 
  labs(color = "Iter") +
  ylab(ts.obj.names) +
  theme_light()  +
  theme(
    strip.text.x = element_text(
      size = 12, color = "black"),
    strip.text.y = element_text(
      size = 12, color = "black"),
    strip.background = element_rect(
      color="grey75", fill="#FFFFFF", linewidth=1, linetype="solid")
  )
ggsave(obj.plot, filename = file.path(simdata_dir, paste0("om_",  prefix.plot, "_", ts.obj.names, ".png")),
       height=pheight, width=pwidth)



} # end function plot_ts_obj  ====



############# get_recr fn ====

get_recr <- function(y, slot.number)  {   
  tmp <- as.data.frame(y[slot.number])
  tmp.r <- tmp[,1]
  return(tmp.r)
}  #end function get_recr   ====



############# get_total fn ====

get_total <- function(y, slot.number)  {   
  tmp <- as.data.frame(y[slot.number])
  tmp.r <- apply(tmp, 1, sum)
  return(tmp.r)
}  #end function get_recr   ====



############# plot_SR fn ====

plot_SR <- function(sim_data=NULL, simdata_dir=NULL,
                    pheight=6, pwidth=6) {   

  if (!dir.exists(simdata_dir))  dir.create(simdata_dir)
  
  n_om <- length(sim_data)
  n_sim <- length(sim_data[[1]])
  om.nums <-  as.integer(sapply( names(sim_data), function(x)  substr(x, 4, nchar(x)) ) )
  sim.nums <-  as.integer(sapply( names(sim_data[[1]]), function(x)  substr(x, 4, nchar(x)) ) )
  prefix.plot <- paste0(om.nums, sep=".", collapse="")
  

  SSB.slot <- which((names(sim_data[[1]][[1]])== "SSB" ))
  NAA.slot <- which((names(sim_data[[1]][[1]])== "NAA" ))

  SSB.obj <-  lapply(sim_data, function(x) lapply(x, function(y)  y[SSB.slot] )   ) 
  SSB.tib <-  as_tibble(lapply(SSB.obj, function(x) as_tibble(unlist(x) )) )
  
  Rec.obj <-  lapply(sim_data, function(x) lapply(x, function(y)  get_recr(y, NAA.slot) )   ) 
  Rec.tib <-  as_tibble(lapply(Rec.obj, function(x) as_tibble(unlist(x) )) )
  
  
  len.SSB <- dim(as.data.frame(SSB.obj[[1]][1]) ) [1]
  years <- seq(1, len.SSB)
  len.Rec <- dim(as.data.frame(Rec.obj[[1]][1]) ) [1]
  years.rec <- years-1  #   #build in 1 year lag
  
  SSB.tib1 <- SSB.tib %>%
    mutate(sim = rep(rep(sim.nums, each=len.SSB)))  %>%
    mutate( Year = rep(years, n_sim)) %>%
    pivot_longer(cols=starts_with("om"), names_to="OM", values_to="Value")  %>%
    mutate( OM=as.factor(OM))  %>%
    mutate(SSB = unlist(Value))  %>%
    select(-Value)
  
  
  Rec.tib1 <- Rec.tib %>%
    mutate(sim = rep(rep(sim.nums, each=len.Rec)))  %>%
    mutate( Year = rep(years.rec, n_sim)) %>%
    pivot_longer(cols=starts_with("om"), names_to="OM", values_to="Value")  %>%
    mutate( OM=as.factor(OM))  %>%
    mutate(Recruits = unlist(Value)) %>%
    filter(Year>0)   %>%
    select(-Value)
  
  Rec.tib.all.yrs <- Rec.tib %>%
    mutate(sim = rep(rep(sim.nums, each=len.Rec)))  %>%
    mutate( Year = rep(years, n_sim)) %>%
    pivot_longer(cols=starts_with("om"), names_to="OM", values_to="Value")  %>%
    mutate( OM=as.factor(OM))  %>%
    mutate(Recruits = unlist(Value)) %>%
    select(-Value)
  
   SR.tib <- SSB.tib1 %>%
     left_join(Rec.tib1)
   
   SR.line.plot <- ggplot(SR.tib, aes(x=SSB, y=Recruits, color=as.factor(sim)))  +
     facet_wrap(~OM)  +
     geom_line(alpha=0.45) +
     theme_light()  +
     theme(axis.text.x = element_text(angle = 90))  +
     labs(color = "Iter") +
     theme(
       strip.text.x = element_text(
         size = 12, color = "black"),
       strip.text.y = element_text(
         size = 12, color = "black"),
       strip.background = element_rect(
         color="grey75", fill="#FFFFFF", linewidth=1, linetype="solid")
     )
   ggsave(SR.line.plot, filename = file.path(simdata_dir, paste0("om_", prefix.plot, "_SR.line.plot",  ".png")),
          height=pheight, width=pwidth)
   
   
   SR.point.plot <- ggplot(SR.tib, aes(x=SSB, y=Recruits, color=as.factor(sim)))  +
     facet_wrap(~OM)  +
     geom_point(alpha=0.45) +
     theme_light()  +
     theme(axis.text.x = element_text(angle = 90))  +
     labs(fill = "Iter") +
     theme(
       strip.text.x = element_text(
         size = 12, color = "black"),
       strip.text.y = element_text(
         size = 12, color = "black"),
       strip.background = element_rect(
         color="grey75", fill="#FFFFFF", linewidth=1, linetype="solid")
     )
   ggsave(SR.point.plot, filename = file.path(simdata_dir, paste0("om_", prefix.plot, "_SR.point.plot",  ".png")),
          height=pheight, width=pwidth)
   
   Recr.plot <- ggplot(Rec.tib.all.yrs, aes(x=Year, y=Recruits, color=as.factor(sim)))  +
     facet_wrap(~OM)  +
     geom_line(alpha=0.45) +
     theme_light()  +
     labs(color = "Iter") +
     theme(
       strip.text.x = element_text(
         size = 12, color = "black"),
       strip.text.y = element_text(
         size = 12, color = "black"),
       strip.background = element_rect(
         color="grey75", fill="#FFFFFF", linewidth=1, linetype="solid")
     )
   ggsave(Recr.plot, filename = file.path(simdata_dir, paste0("om_", prefix.plot, "_Recr.plot",  ".png")),
          height=pheight, width=pwidth)
   
  
}  # end function plot_SR ====

