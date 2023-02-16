# Mortality

In all there are 288 Operating models which were used to simulate 100 time series and observations. Any process errors configured in the operating model were simulated and, conditional on the simulated process errors, all aggregate, age composition and environmental covariates were simulated.

### Population

The population structure and demographic parameters are the same as those in the Project 0 simulation study.

There are 10 age classes: ages 1 to 10+.

Spawning was assumed to occur 1/4 of the way through the year.

Natural mortality rate was assumed 0.2 when it was constant and the mean of the time series process for operating models with M random effects.
maturity a50 = 2.89, slope = 0.88
  
Weight at age was generated with a LVB growth

$$L_a = L_{\infty}\left(1 - e^{-k(a - t_0)}\right) $$

where $t_0 = 0$, $L_\infty = 85$, and $k = 0.3$, and a L-W relationship such that 

$$W_a = \theta_1 L_a^{\theta_2}$$

where $\theta_1 = e^{-12.1}$ and $\theta_2 = 3.2$.

We assumed a Beverton-Holt stock recruit function with constant pre-recruit mortality parameters for all operating models. All post-recruit productivity components other than natural mortality are constant in the operating models. For operating models with no M process errors and no environmental covariate effects M would also be constant. We specified stock-recruit parameters by equating Fmsy and F40 = 0.348 under the model where M is constant and unfished recruitment was assumed $R_0 = e^{10}$. This equates to a steepness of 0.69 and $\alpha=0.60$ and $\beta = 2.4 \times 10^{-5}$ for the Beverton-Holt parameterization
$$N_{1,y} = \frac{\alpha \text{SSB}_{y-1}}{1 + \beta \text{SSB}_{y-1}}. $$


The magnitude of the overfishing assumptions is based on average estimates of overfishing for NE groundfish stocks from Wiedenmann et al. (20XX). Legault et al. (2023) also used similar approaches to defining fishing mortality histories for operating models.

### Fleets

We assumed a single fleet operating year round for catch observations with logistic selectivity for the fleet with $a_{50} = 5$ and slope = 1. This selectivity is was used to define $F_{\text{MSY}}$ for the Beverton-Holt stock recruitment parameters above. We assumed a logistic-normal distribution for the age-composition observations for the fleet.


Initial population was configured at the equilibrium distribution fishing at either $F = 2.5\times F_{\text{MSY}}$ or $F = F_{\text{MSY}}$ for the two alternative fishing histories. That is for a deterministic model, the age composition would not change over time when the fishing mortality was constant at the respective level. 

### Indices
Two time series of surveys are assumed and observed in numbers rather than biomass for the entire 40 year period with one occurring in the spring (0.25 way through the year) and one in the fall (0.75 way through the year). Catchability of both surveys are assumed to be 0.1.  We assumed logistic selectivity for both indices with $a_{50} = 5$ and slope = 1. We assumed a logistic-normal distribution for the age-composition observations.

## Factors defining operating models

### Process errors in numbers at age and natural mortality

We assumed one of three different configurations of process errors in numbers at age and/or natural mortality. In the first ("rec"), we assume iid deviation of log-recruitment from the stock-recruit function predicted log-recruitment with standard errors of 0.5. In the second ("rec+1") we assumed iid process errors in both recruitment (sd = 0.5) and transitions in log-abundance between years 
$$ log N_{y,a} \sim \text{N} \left(\log {\hat N} , sd = 0.3\right). $$
In the third configuration ("rec+M") we assumed iid process errors for recruitment (sd = 0.5) and iid annual process errors in log-natural mortality (sd = 0.3) constant at all ages.

### Environmental covariates

We assumed state-space models for covariates. The latent "true" covariate was assumed to arise as a first order autoregressive process with mean 0 and standard deviation of either 0.1 or 0.5. The correlation of the process was either 0 (iid) or 0.5. The error in the observations of this process have standard deviations of 0.1 or 0.5. Finally the linear effect on log-natural mortality was 0 (no effect), 0.25, or 0.5. All 24 combinations of these assumptions were used to defined operating models.

### Fishing histories
All operating models assumed one of two different fishing histories. 
One : Fishing mortality is equal to $F_{\text{MSY}}$ (0.348) for the whole 40 year period.
Two : Fishing mortality is 2.5 times $F_{\text{MSY}}$ for the first 20 years then changes to $F_{\text{MSY}}$ for the last 20 years.

### Observation Uncertainty 

Standard deviation for log-aggregate catch was 0.1. There were two levels of observation error variance for indices and age composition for both indices and fleet catch. A low uncertainty specification assumed standard deviation of both series of log-aggregate index observations was 0.1 and the standard deviation of the logistic-normal for age composition observations was 0.3 In the high uncertainty specification the standard deviation for log-aggregate indices was 0.4 and that for the age composition observations was 1.5. For all estimating models, standard deviation for log-aggregate observations was assumed known whereas that for the logitic-normal age composition observations was estimated.
