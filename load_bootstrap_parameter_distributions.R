# Simulation parameters
nruns=100000

# What are the varying parameters
# 1) Sigma
mean_sigma=SummaryTable_Efficacy_NeutRatio_SD_SEM$PooledSD[1]
sd_sigma=SummaryTable_Efficacy_NeutRatio_SD_SEM$SE_PooledSD[1]

# 2) Hill Coeff
# 3) IC50
cov_hill_IC50<-solve(FittedLogistic_RawEfficacy_MeanRept_SDPool$hessian)[9:10,9:10]
mean_hill_IC50=c(tail(FittedLogistic_RawEfficacy_MeanRept_SDPool$estimate,2))
sd_hill_IC50=sqrt(diag(cov_hill_IC50))

# 4) Efficacy -> Starting Neuts assume normal dist of starting log10 neut ratios
# Use pfizer as a baseline
#study = 'Pfizer'
#eff_init= SummaryTable_Efficacy_NeutRatio_SD_SEM$Efficacy[SummaryTable_Efficacy_NeutRatio_SD_SEM$Study==study]
#eff_lower95 = SummaryTable_Efficacy_NeutRatio_SD_SEM$Lower[SummaryTable_Efficacy_NeutRatio_SD_SEM$Study==study]/100
#eff_upper95 = SummaryTable_Efficacy_NeutRatio_SD_SEM$Upper[SummaryTable_Efficacy_NeutRatio_SD_SEM$Study==study]/100

#4) Starting neut values (by vaccine)
vaccines_to_use = c('Pfizer','Moderna','Astra','mRNA','Boosted')
vaccines = names(log10NeutRatios)[names(log10NeutRatios) %in% vaccines_to_use]
mean_log10neutR = log10NeutRatios[names(log10NeutRatios) %in% vaccines]
# USE SEM from SummaryTable_Efficacy_NeutRatio_SD_SEM as the SD from this distribution

nserums = length(mean_log10neutR)

#sd_log10neutR = estimate_sd_from_CIs(mean_log10neutR,log10neutR_low,log10neutR_high)

# 5) Drops for variants
# meanDelta = log10(1/3.9)
# sdDelta = .1
# meanOmicron = log10(1/9.7)
# sdOmicron = .1
omicron_m_l_u = log10(1/c(22,30,16))

#omicron_m_l_u = log10(1/c(9.7,17.1,5.5))
delta_m_l_u = log10(1/c(3.9,5.5,3.5))
mean_variants=c(delta_m_l_u[1],omicron_m_l_u[1])
low_variants =c(delta_m_l_u[2],omicron_m_l_u[2])
high_variants = c(delta_m_l_u[3],omicron_m_l_u[3])
sd_variants = estimate_sd_from_CIs(mean_variants,low_variants,high_variants)

# 6) Shift to severe
#mean_shift = 20/3
#sd_shift = .05*mean_shift
#Note that  think the sd ofo this shift is too small becuase the off diagonal elements for the hill coefficient are not considered. 
#Really we should look at a 3x3 matrix that includes the last 3 elements and accounts for hill, IC50 and shift
indicies<-c(-1,0)+length(FittedLogistic_RawEfficacy_MeanRept_SDPool_Diff_EC50_Only$estimate)
CovS<-solve(FittedLogistic_RawEfficacy_MeanRept_SDPool_Diff_EC50_Only$hessian)[indicies,indicies]
mean_shift=sum(FittedLogistic_RawEfficacy_MeanRept_SDPool_Diff_EC50_Only$estimate[indicies]*c(-1,1))
sd_shift=sqrt(sum(c(diag(CovS)^2,2*CovS[1,2]*CovS[2,1])))

# 7) hill, IC50 for severe
cov_hill_IC50_severe = solve(FittedLogistic_RawEfficacy_MeanRept_SDPool_Diff_EC50_Only$hessian)[15:16,15:16]
mean_hill_IC50_severe=c(tail(FittedLogistic_RawEfficacy_MeanRept_SDPool_Diff_EC50_Only$estimate,3))[1:2]
sd_hill_IC50_severe=sqrt(diag(cov_hill_IC50_severe))

# 8) Ab decay rate 
mean_decay = log(2)/108
sd_decay = estimate_sd_from_CIs(mean_decay,log(2)/159,log(2)/82)

#9) Neut ratios
mean_log10neutR = log10NeutRatios[names(log10NeutRatios) %in% vaccines]
sd_log10neutR = log10NeutRatios_SEM[names(log10NeutRatios_SEM) %in% vaccines]

select_parameters_from_distributions = function(nruns, paramsIn=NULL){
  
  paramnames = c('sigma','hill','IC50','decay',paste0('log10neutR_',vaccines,sep=""),'severe_shift','deltaDrop','omicronDrop','hill_severe','IC50_severe','variant','serum')
  means = c(mean_sigma, mean_hill_IC50, mean_decay,mean_log10neutR,mean_shift,mean_variants,mean_hill_IC50_severe,0,0)
  sds = c(sd_sigma, sd_hill_IC50, sd_decay, sd_log10neutR,sd_shift,sd_variants,sd_hill_IC50_severe,1,1)
  names(means)=paramnames
  names(sds)=paramnames
  paramseeds = randomLHS(nruns,length(means))
  params = data.frame(paramseeds)
  colnames(params)<-paramnames
  tstart=Sys.time()
  for (i in c(1:(ncol(paramseeds)-2))){
    params[,i] = qnorm(paramseeds[,i], mean=means[i],sd = sds[i])   
  }
  params[colnames(params) %in% c('hill','IC50')]= rmvnorm(nruns,mean=mean_hill_IC50, sigma = cov_hill_IC50)
  params[colnames(params) %in% c('hill_severe','IC50_severe')]= rmvnorm(nruns,mean=mean_hill_IC50_severe, sigma = cov_hill_IC50_severe)
  # Replace with input params if the same length
  if (!is.null(paramsIn)){
    if (nrow(paramsIn)==nruns){
      for (cname in colnames(params)){
        ind = which(colnames(paramsIn)==cname)
        if (length(ind)==1){
          params[,cname]=paramsIn[,ind]
        }
      }
    }
  }
  params  
}
    
select_base_parameters_from_distributions = function(withBoosting=F){
  paramnames = c('sigma','hill','IC50','decay',paste0('log10neutR_',vaccines,sep=""),'severe_shift','deltaDrop','omicronDrop','hill_severe','IC50_severe','variant','serum')
  means = c(mean_sigma, mean_hill_IC50, mean_decay,mean_log10neutR,mean_shift,mean_variants,mean_hill_IC50_severe,0,0)
  sds = c(sd_sigma, sd_hill_IC50, sd_decay, sd_log10neutR,sd_shift,sd_variants,sd_hill_IC50_severe,1,1)
  names(means)=paramnames
  names(sds)=paramnames
  paramseeds = rep(0.5,length(means))
  params = paramseeds
  names(params)<-paramnames
  for (i in c(1:(length(paramseeds)-2))){
    params[i] = qnorm(paramseeds[i], mean=means[i],sd = sds[i])   
  }
  
  if(withBoosting){
    nvax = length(vaccines)
  } else{
    nvax = sum(vaccines!='Boosted')
  }
  nvariants = length(variants_to_use)
  paramstruct = data.frame(matrix(rep(params,each=nvariants*nvax),nrow=nvariants*nvax))
  colnames(paramstruct) = paramnames
  paramstruct$variant = rep(c(0:(nvariants-1))/nvariants,each=nvax)
  paramstruct$serum = rep(c(0:(nvax-1))/length(vaccines),nvariants)
  paramstruct
}

