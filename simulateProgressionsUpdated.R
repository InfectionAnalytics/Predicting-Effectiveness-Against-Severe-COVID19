# Set up the number of runs first
source('./load_bootstrap_parameter_distributions.R')


##### Combine all the changing parameters into a vector, and create a matrix of param values

# If you want to regenerate teh curves, set these to true, otherwise set to false.
# They have already been generated for nruns = 1e5
generateNatureMedCurve = F
generateParamStructs=F
generateTimeCourseData=F

if (generateNatureMedCurve){
  ##### Helper function to Estimate the standard deviation from the 95% CIs
  source('./generateSimNatMedCurve.R')
} else {
  load(paste0('./simulated_natmed_efficacy_curve_n',nruns,'.RData'))
}

if(generateParamStructs){
  params = select_parameters_from_distributions(nruns,params)
  base_params = select_base_parameters_from_distributions()
  param_struct_for_efficacy_simulation = make_param_struct_for_efficacy_simulation(params)
  base_param_struct_for_efficacies = make_param_struct_for_efficacy_simulation(base_params)
  save(params,base_params, param_struct_for_efficacy_simulation,base_param_struct_for_efficacies, nserums,file=paste0('non_timecourse_paramseter_structures_n',nruns,'.RData'))
} else {
  load(file=paste0('non_timecourse_paramseter_structures_n',nruns,'.RData'))
}

# Now run the time course simulations
if (generateTimeCourseData){
  simulated_data_list = run_simulations_from_params(param_struct_for_efficacy_simulation)    
  base_simulated_data_list = run_simulations_from_params(base_param_struct_for_efficacies)
  
  # Now make the summaries
  time_course_summary=simulation_list_into_melted_frame(simulated_data_list, summarise=T)
  base_time_course=simulation_list_into_melted_frame(base_simulated_data_list, summarise=F)
  
  # Add in the means
  time_course_summary = merge(time_course_summary,base_time_course, by=c("Vaccine","Booster","boosted","VariantGroup","regime", "variable","Day"))
  
  # Extract the efficacies
  time_course_efficacies_summary = time_course_summary %>% 
    filter(variable %in% c('sympt_eff','severe_eff')) %>%
    rename(Efficacy = modelMean, EfficacyMin=lower,EfficacyMax=upper,EfficacyMedian=mid) %>%
    mutate(plotjoin = regime, fitgroup=regime, Variant = VariantGroup, fillVar = Vaccine,
           EndpointCat = ifelse(variable=='sympt_eff','Symptomatic','Severe'),
           TimeMidpt = Day/30,
           Efficacy = 100*Efficacy,EfficacyMedian=100*EfficacyMedian, EfficacyMin = 100*EfficacyMin, EfficacyMax=100*EfficacyMax) %>% 
    ungroup %>% select(-c('regime','variable'))
  time_course_efficacies_summary$EndpointCat = factor(time_course_efficacies_summary$EndpointCat, levels = endpoint_categories)
  
  save(time_course_summary,time_course_efficacies_summary, base_time_course,file=paste0('time_course_summary_n',nruns,'.RData'))
  
} else {
  load(file=paste0('time_course_summary_n',nruns,'.RData'))
}

