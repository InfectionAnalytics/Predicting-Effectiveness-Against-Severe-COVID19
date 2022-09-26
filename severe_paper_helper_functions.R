# make a rescaler
alpharescale = function(CIwidth, from=c(0,100)){
  # Make them between 0 and 100
  scaledvals = pmax(pmin(CIwidth,from[2]),from[1]) 
  #Do a basic rescale
  scaledvals = 1-scales::rescale(scaledvals, from=from)
  scaledvals = case_when(scaledvals<.75~scaledvals/3,
                         scaledvals<.9~scaledvals/2,
                         T~scaledvals)
  
}


# Make the param struct to use in simulations
make_param_struct_for_efficacy_simulation = function (params, withBoosting=F){
  ##### Now simulate the data based on parameter choices
  if (withBoosting){
    boosterDose = 'mRNA'
  } else {
    boosterDose = 'None'
  }
  param_struct <- params %>% mutate(
    Vaccine=case_when(serum<1/nserums~vaccines_to_use[1],
                      serum<2/nserums~vaccines_to_use[2],
                      serum<3/nserums~vaccines_to_use[3],
                      serum<4/nserums~vaccines_to_use[4],
                      serum<(nserums-3/4)/nserums~vaccines_to_use[1], # These are now the boosted ones
                      serum<(nserums-2/4)/nserums~vaccines_to_use[2],
                      serum<(nserums-1/4)/nserums~vaccines_to_use[3],
                      T~vaccines_to_use[4]),
    Booster=ifelse(serum<(nserums-1)/nserums,'None',boosterDose),
    boosted=ifelse(Booster=='None','Unboosted','Boosted'),
    VariantGroup = case_when(variant<1/3~'Delta',
                             variant<2/3~'Omicron',
                             T~'WT'),
    log10neutR=case_when(boosted=='Boosted' ~log10neutR_Boosted,#log10NeutRatios['Boosted']
                         Vaccine=='Pfizer' ~log10neutR_Pfizer, #log10NeutRatios['Pfizer']
                         Vaccine=='Moderna'~log10neutR_Moderna,#log10NeutRatios['Moderna']
                         Vaccine=='Astra' ~log10neutR_Astra,#log10NeutRatios['Astra']
                         Vaccine=='mRNA' ~log10neutR_mRNA)+case_when(VariantGroup=='Delta'~deltaDrop,
                                                                     VariantGroup=='Omicron'~omicronDrop,
                                                                     T~0),
    regime = paste0(Vaccine,Booster),
    renormFactor = exp(trialdays*decay/2),
    startingLog10neutR = log10neutR+log10(renormFactor)) 
  param_struct = param_struct %>% 
    mutate(
    sympt_efficacy = LogisticModel_PercentUninfected(startingLog10neutR,sigma,hill,IC50),
    #severe_efficacy = LogisticModel_PercentUninfected(log10neutR+severe_shift,sigma,hill,IC50),
    severe_efficacy = LogisticModel_PercentUninfected(startingLog10neutR,sigma,hill_severe,IC50_severe),
    progression_protection = 1-(1-severe_efficacy)/(1-sympt_efficacy),
    progression_ratio = 100*(1-progression_protection),
    boosted = factor(boosted,levels=c('Unboosted','Boosted')))

  param_struct$Vaccine = factor(param_struct$Vaccine,levels=vaccine_list[vaccine_list %in% unique(param_struct$Vaccine)])
  param_struct
}

#run the simulations
run_simulations_from_params = function (simulation_params,daygap=1,daystart=0,dayend=240){
   #Next simulate severe and symptomatic efficacy over time
   TimeCourse=seq(daystart, dayend,by=daygap)
   log10neutRs = matrix(ncol=length(TimeCourse),nrow=nrow(simulation_params))
   simulated_sympt_efficacy = log10neutRs
   simulated_severe_efficacy = log10neutRs
   
   for (dayind in c(1:length(TimeCourse))){
     day = TimeCourse[dayind]
     #paramvect = c(simulation_params[p,c('sigma','hill','IC50','severe_shift','log10neutR','decay', 'hill_severe','IC50_severe')])
     # Renormalise to get the efficacy over the whole trial
     # astra has median follow up of 2 months after second dose
     # pfizer has median followup of 68 days
     # moderna has median follow up of 63 days
     renormFactor = exp(trialdays*simulation_params$decay/2)
     
     startingLog10neutR = simulation_params$log10neutR+log10(renormFactor)
     log10neutRs[,dayind] = log10((10^startingLog10neutR)*exp(-day*abs(simulation_params$decay)))
     
     simulated_sympt_efficacy[,dayind] = LogisticModel_PercentUninfected(log10neutRs[,dayind],simulation_params$sigma,simulation_params$hill,simulation_params$IC50)  
     simulated_severe_efficacy[,dayind]= LogisticModel_PercentUninfected(log10neutRs[,dayind],simulation_params$sigma,simulation_params$hill_severe,simulation_params$IC50_severe)  
     
     if ((dayind %% 10)==0){
       print(paste0('Iter ',dayind,' Day: ',day))
       save(log10neutRs,simulated_severe_efficacy,simulated_sympt_efficacy,file=paste0('intermediate_timecourse_simulation_data_n',nruns,'.RData'))
     }
   }
   
   # Name the columns
   colnames(log10neutRs)<-paste0('log10neutR_',c(0:(ncol(log10neutRs)-1)))
   colnames(simulated_sympt_efficacy)<-paste0('sympt_eff_',c(0:(ncol(simulated_sympt_efficacy)-1)))
   colnames(simulated_severe_efficacy)<-paste0('severe_eff_',c(0:(ncol(simulated_severe_efficacy)-1)))
   
   # Calculate the progression data
   progression_protection_data = 1-(1-simulated_severe_efficacy)/(1-simulated_sympt_efficacy)
   colnames(progression_protection_data)=str_replace(colnames(progression_protection_data),'severe_eff','progression_protection')
   progression_ratio_data = 100*(1-progression_protection_data)
   colnames(progression_ratio_data)=str_replace(colnames(progression_ratio_data),'progression_protection','progression_ratio')
   
   # Make into a list
   simulation_results = list(param_struct=simulation_params,
                             log10neutR = log10neutRs,
                             sympt_efficacy = simulated_sympt_efficacy, 
                             severe_efficacy = simulated_severe_efficacy,
                             progression_protection = progression_protection_data, 
                             progression_ratio = progression_ratio_data,
                             time_course = TimeCourse)
}

simulation_list_into_melted_frame = function(simulation_list, summarise = T){
  i=1
  ncolsleft  = ncol(simulation_list$log10neutR)
  maxcols = ifelse(summarise,10,ncolsleft)
  startcol = 1
  daygap = diff(simulation_list$time_course)[1]
  while(ncolsleft > maxcols | i==1){
    thiscols = min(ncolsleft,maxcols)
    inds = c(startcol:(startcol+thiscols-1))
    this_simulated_data<-cbind(simulation_list$param_struct,
                               simulation_list$log10neutR[,inds],
                               simulation_list$sympt_efficacy[,inds],
                               simulation_list$severe_efficacy[,inds],
                               simulation_list$progression_protection[,inds],
                               simulation_list$progression_ratio[,inds])
    
    this_simulated_data_melted = melt(this_simulated_data, id.vars=colnames(simulation_list$param_struct))
    
    if(summarise){
      this_simulated_time_course_summary = distinct(this_simulated_data_melted  %>% group_by(Vaccine, Booster, boosted, VariantGroup, regime, variable) %>%  
                                                    summarise(lower = quantile(value, probs = 0.025, na.rm=T),
                                                              mid = quantile(value, probs = 0.5, na.rm=T),
                                                              upper = quantile(value, probs = 0.975, na.rm=T),
                                                              Day = daygap*as.numeric(str_sub(variable,stri_locate_last(variable,fixed='_')[,2]+1,100)),
                                                              variable = str_sub(variable,1,stri_locate_last(variable,fixed='_')[,2]-1)))
      if (i==1){
        simulated_time_course_summary=this_simulated_time_course_summary
      } else {
        simulated_time_course_summary=rbind(simulated_time_course_summary,this_simulated_time_course_summary)
      }
    } else {
      simulated_time_course_summary=distinct(this_simulated_data_melted  %>% group_by(Vaccine, Booster, boosted, VariantGroup, regime, variable) %>%  
                                               summarise(modelMean = value,
                                                         Day = daygap*as.numeric(str_sub(variable,stri_locate_last(variable,fixed='_')[,2]+1,100)),
                                                         variable = str_sub(variable,1,stri_locate_last(variable,fixed='_')[,2]-1)))
    }
      
    startcol=startcol+thiscols
    ncolsleft = ncolsleft-thiscols
    i=i+1
  }
  simulated_time_course_summary
}