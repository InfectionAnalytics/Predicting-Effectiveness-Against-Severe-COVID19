source('./plot_efficacies_revised.R')
source('./load_bootstrap_parameter_distributions.R')

# Set up the real world data for plotting and analysis
imported_efficacy_data=import_efficacy_data()
rw_unboosted_data = imported_efficacy_data %>% 
  filter(UseInAnalysis, AgeGp != 'Under 15', 
         !(Vaccine %in% c('Pfizer+Moderna','JnJ')), TimeMax>.5, 
        boosted=='Unboosted', VariantGroup%in% c('Omicron','Delta','Non-Delta')) %>%
  mutate(fitgroup = paste0(VariantGroup,Vaccine, EndpointCat),
    VariantGroup = factor(VariantGroup, levels = names(variant_colours)),
    TimeMidpt = as.numeric(TimeMidpt),
    studyname = paste0(FirstAuthor,' et. al.'),
    CI_width = pmin(EfficacyMax-EfficacyMin,100))

rw_unboosted_data$CI_width[is.na(rw_unboosted_data$CI_width)]=0
rw_unboosted_symptsev_data = filter(rw_unboosted_data, EndpointCat %in% c('Symptomatic','Severe'))


# Load the data we have simulated previously
loadsimulateddata=T
if (loadsimulateddata){
  # This is the time course summary
  load(file=paste0('time_course_summary_n',nruns,'.RData'))
  #time_course_summary,time_course_efficacies_summary, base_time_course
}
loadSimNatmedEffCurve=T
if (loadSimNatmedEffCurve){
  load(paste0('./simulated_natmed_efficacy_curve_n',nruns,'.RData'))
  efficacy_summary_curve$EndpointCat = factor(efficacy_summary_curve$EndpointCat, levels=endpoint_categories)
}
variants = c('Non-Delta','Delta','Omicron')

rw_and_sim_efficacy_data= bind_rows(rw_unboosted_symptsev_data,time_course_efficacies_summary)
rw_and_sim_efficacy_data$VariantGroup[rw_and_sim_efficacy_data$VariantGroup=='WT']='Non-Delta'

rw_and_sim_efficacy_data = rw_and_sim_efficacy_data %>% filter(VariantGroup %in% variants) %>%
  mutate(VariantGroup = factor(VariantGroup, levels=variants),
  CI_width=pmin(EfficacyMax-EfficacyMin,100))

#progression_base_data = filter(rw_unboosted_data, UseInAnalysis, Variant %in% c('Delta','Omicron'))

#progression_base_data$Efficacy = pmin(progression_base_data$Efficacy,100)
rw_unboosted_sev_data = filter(rw_unboosted_symptsev_data, EndpointCat=='Severe') %>% 
                  rename(SevereEff = Efficacy, SevereMin = EfficacyMin, SevereMax = EfficacyMax) %>% 
                  select(-c("Endpoint","EndpointCat","join","VariantVaccine","notboosted","fitgroup","plotjoin" ))
#symptomaticEfficacy = filter(progression_base_data, EndpointCat %in% c('Symptomatic','Infection')) %>% # This is only if we want to include infection
rw_unboosted_sympt_data = filter(rw_unboosted_symptsev_data, EndpointCat %in% c('Symptomatic')) %>% 
                  rename(SymptomaticEff = Efficacy,SymptomaticMin = EfficacyMin, SymptomaticMax = EfficacyMax) %>% 
                  select(-c("Endpoint","EndpointCat","join","VariantVaccine","notboosted","fitgroup","plotjoin" ))

rw_matched_symptsev_data = merge(rw_unboosted_sev_data,rw_unboosted_sympt_data, 
                                 by = c('Paper', 'FirstAuthor','Vaccine', 'Variant','VariantGroup', 'Age','AgeGp', 'Booster','boosted','TimeMidpt','TimeMin','TimeMax'))%>% 
  mutate(matcheffs = paste0(Paper, Vaccine, Variant, AgeGp, Booster,TimeMidpt,TimeMin,TimeMax), 
         progression_ratio =  100*(100-SevereEff)/(100-SymptomaticEff), 
         regime = paste0(Vaccine,Booster),
         studyname = paste0(FirstAuthor,' et. al.'))
rw_matched_symptsev_data$regime = factor(rw_matched_symptsev_data$regime, levels = c('JnJNone','AstraNone','PfizerNone','ModernaNone','AstraPfizer','AstraModerna','PfizerPfizer','PfizerModerna','ModernaPfizer','ModernaModerna'))
levels(rw_matched_symptsev_data$regime) =  c('AD26 only','ChAd-Ox only','BNT162b2 only','mRNA1273 only','ChAd-Ox + BNT162b2','ChAd-Ox + mRNA1272','BNT162b2 + BNT162b2','BNT162b2 + mRNA1273','mRNA1272 + BNT162b2','mRNA1272 + mRNA1272')
#rw_matched_symptsev_data = rw_matched_symptsev_data %>% mutate(overlap = (SevereMin>=SymptomaticMin & SevereMin<=SymptomaticMax)| SymptomaticMin>=SevereMin & SymptomaticMin<=SevereMax)
#rw_matched_symptsev_data = rw_matched_symptsev_data %>% mutate(smalloverlap = (SevereMin<SymptomaticEff |SymptomaticMax>SevereEff))


# Add in 95% CIs of progression
niter = 50000
rw_matched_symptsev_data$prog_protect_low=NA
rw_matched_symptsev_data$prog_protect_high=NA
for (p in c(1:nrow(rw_matched_symptsev_data))){
  if (sum(is.na(c(rw_matched_symptsev_data[p,'SymptomaticEff'],rw_matched_symptsev_data[p,'SymptomaticMin'],rw_matched_symptsev_data[p,'SymptomaticMax'])))==0){
    symptSD = estimate_sd_from_CIs(rw_matched_symptsev_data[p,'SymptomaticEff'],rw_matched_symptsev_data[p,'SymptomaticMin'],rw_matched_symptsev_data[p,'SymptomaticMax'])
  } else {symptSD=0}
  symptEffs = rnorm(niter,rw_matched_symptsev_data[p,'SymptomaticEff'],symptSD)
  if (sum(is.na(c(rw_matched_symptsev_data[p,'SevereEff'],rw_matched_symptsev_data[p,'SevereMin'],rw_matched_symptsev_data[p,'SevereMax'])))==0){
    severeSD = estimate_sd_from_CIs(rw_matched_symptsev_data[p,'SevereEff'],rw_matched_symptsev_data[p,'SevereMin'],rw_matched_symptsev_data[p,'SevereMax'])
  } else {severeSD=0}
  severeEffs = rnorm(niter,rw_matched_symptsev_data[p,'SevereEff'],severeSD)
  progression_protections = 100-100*(100-severeEffs)/(100-symptEffs)
  rw_matched_symptsev_data[p,c('prog_protect_low','prog_protect_high')]=quantile(progression_protections, probs = c(.025,.975),na.rm=T)
  #protect_CI_width
}
rw_matched_symptsev_data$prog_protect_CI_width=rw_matched_symptsev_data$prog_protect_high-rw_matched_symptsev_data$prog_protect_low
rw_matched_symptsev_data= rw_matched_symptsev_data %>% mutate(prog_protect_low=pmax(prog_protect_low,0),prog_protect_high=pmin(prog_protect_high,100),prog_protect_CI_width = pmin(prog_protect_CI_width,100))
rw_matched_symptsev_data=filter(rw_matched_symptsev_data, VariantGroup !='Mix', Vaccine!='JnJ')

