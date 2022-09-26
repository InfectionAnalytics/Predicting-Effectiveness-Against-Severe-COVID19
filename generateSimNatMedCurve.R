print('Inside generateSimNatMedCurve')
source('./load_bootstrap_parameter_distributions.R')
params = select_parameters_from_distributions(nruns)

generateSimNatMedCurve=T
if(generateSimNatMedCurve){
  # paramnames = c('sigma','hill','IC50','severe_shift','hill_severe','IC50_severe')
  # means = c(mean_sigma, mean_hill_IC50,mean_shift,mean_hill_IC50_severe)
  # sds = c(sd_sigma, sd_hill_IC50,sd_shift, sd_hill_IC50_severe)
  # names(means)=paramnames
  # names(sds)=paramnames
  # paramseeds = randomLHS(nruns,length(means))
  # params = data.frame(paramseeds)
  # colnames(params)<-paramnames
  # tstart=Sys.time()
  # for (i in c(1:ncol(paramseeds))){
  #   params[,i] = qnorm(paramseeds[,i], mean=means[i],sd = sds[i])   
  # }
  log10neuts = log10(2^seq(-7,4,.25))
  simulated_sympt_efficacy = matrix(ncol=length(log10neuts),nrow=nruns)
  simulated_severe_efficacy = matrix(ncol=length(log10neuts),nrow=nruns)
  #simulated_severe_efficacy_old_wrong = matrix(ncol=length(log10neuts),nrow=nruns)
  #for (p in c(1:nruns)){
  ind=1
  for (log10N in log10neuts){
    #paramvect = params[p,]
    #simulated_sympt_efficacy[p,] = 100*LogisticModel_PercentUninfected(log10neuts,paramvect$sigma,paramvect$hill,paramvect$IC50)
    simulated_sympt_efficacy[,ind] = 100*LogisticModel_PercentUninfected(log10N,params$sigma,params$hill,params$IC50)
    #simulated_severe_efficacy_old_wrong[p,] = 100*LogisticModel_PercentUninfected(log10neuts,paramvect$sigma,paramvect$hill,paramvect$IC50-paramvect$severe_shift)
    #simulated_severe_efficacy[p,] = 100*LogisticModel_PercentUninfected(log10neuts,paramvect$sigma,paramvect$hill_severe,paramvect$IC50_severe)
    simulated_severe_efficacy[,ind] = 100*LogisticModel_PercentUninfected(log10N,params$sigma,params$hill_severe,params$IC50_severe)
    #if ((p %% 10)==0){
      #print(paste0('Iter ',p))
      #save(simulated_sympt_efficacy,simulated_severe_efficacy,file=paste0('intermediate_simulated_efficacy_data',nruns,'.Rdata'))
    #}
    #if ((ind %% 5)==0){
      print(paste0('Iter ',ind,'/',length(log10neuts)))
      save(simulated_sympt_efficacy,simulated_severe_efficacy,file=paste0('intermediate_simulated_efficacy_data',nruns,'.Rdata'))
    #}
    ind=ind+1
  }
}
simulated_progression = 100*(1-(100-simulated_severe_efficacy)/(100-simulated_sympt_efficacy))

efficacy_summary_curve_sympt = data.frame(t(sapply(data.frame(simulated_sympt_efficacy),quantile,c(.5,.025,.975)))) %>% mutate(EndpointCat = 'Symptomatic')
colnames(efficacy_summary_curve_sympt)[1:3]=c('eff_median','eff_lower','eff_upper')
efficacy_summary_curve_sympt$eff_mean = 100*LogisticModel_PercentUninfected(log10neuts,mean_sigma,mean_hill_IC50[1],mean_hill_IC50[2])

efficacy_summary_curve_severe = data.frame(t(sapply(data.frame(simulated_severe_efficacy),quantile,c(.5,.025,.975)))) %>% mutate(EndpointCat = 'Severe')
colnames(efficacy_summary_curve_severe)[1:3]=c('eff_median','eff_lower','eff_upper')
efficacy_summary_curve_severe$eff_mean = 100*LogisticModel_PercentUninfected(log10neuts,mean_sigma,mean_hill_IC50_severe[1],mean_hill_IC50_severe[2])

#efficacy_summary_curve_severe_old_wrong = data.frame(t(sapply(data.frame(efficacy_summary_curve_severe_old_wrong),quantile,c(.5,.025,.975)))) %>% mutate(EndpointCat = 'SevereOldWrong')
#colnames(efficacy_summary_curve_severe_old_wrong)[1:3]=c('eff_median','eff_lower','eff_upper')

efficacy_summary_curve_prog = data.frame(t(sapply(data.frame(simulated_progression),quantile,c(.5,.025,.975)))) %>% mutate(EndpointCat = 'Progression')
colnames(efficacy_summary_curve_prog)[1:3]=c('prog_protect_median','prog_protect_low','prog_protect_high')
efficacy_summary_curve_prog$prog_protect_mean = 100*(1-(100-efficacy_summary_curve_severe$eff_mean)/(100-efficacy_summary_curve_sympt$eff_mean))
#symptsevereeff_summary = data.frame(SymptEff=c(simulated_sympt_efficacy),SevereEff=c(simulated_severe_efficacy)) %>%
#  mutate(SymptEffRounded = round(SymptEff,1)) %>%
#  group_by(SymptEffRounded) %>%
#  summarise(SevereEffLower = quantile(SevereEff,.025),SevereEffMid = quantile(SevereEff,.5),SevereEffUpper = quantile(SevereEff,.975))

efficacy_summary_curve = rbind(efficacy_summary_curve_sympt, efficacy_summary_curve_severe)
efficacy_summary_curve$log10neutR = c(log10neuts,log10neuts)
efficacy_summary_curve$neutR = 10^efficacy_summary_curve$log10neut
efficacy_summary_curve_prog$log10neutR = log10neuts
efficacy_summary_curve_prog$neutR = 10&log10neuts
efficacy_summary_curve$EndpointCat = factor(efficacy_summary_curve$EndpointCat, levels=c('Symptomatic','Severe'))
#median_eff = sapply(data.frame(simulated_efficacy),quantile,c(.5))
save(efficacy_summary_curve,efficacy_summary_curve_prog,file=paste0('simulated_natmed_efficacy_curve_n',nruns,'.RData'))
  