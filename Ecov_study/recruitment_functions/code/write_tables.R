library(here)
library(xtable)

df.oms          <- readRDS(file.path(here(),"Ecov_study","recruitment_functions", "inputs", "df.oms.RDS"))
df.oms$Ecov_how <- as.integer(df.oms$Ecov_how)

df.ems          <- readRDS(file.path(here(),"Ecov_study","recruitment_functions", "inputs", "df.ems.RDS"))
df.ems$ecov_how <- as.integer(df.ems$ecov_how)

df.oms.tex <- xtable(df.oms)
print(df.oms.tex, file=file.path(here(), 'Ecov_study', 'recruitment_functions', 'inputs', 'df.oms.tex'))

df.ems.tex <- xtable(df.ems)
print(df.ems.tex, file=file.path(here(), 'Ecov_study', 'recruitment_functions', 'inputs', 'df.ems.tex'))

