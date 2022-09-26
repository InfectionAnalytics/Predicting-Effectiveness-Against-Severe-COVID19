library(data.table)
# get into the right directory
thisdir = getwd()
currentdirectory<-dirname(rstudioapi::getActiveDocumentContext()$path)
if (currentdirectory!=""){setwd(currentdirectory)}

# Extract parameters that aren't going to change
load('./natmed_parameters.RData')
sig=SummaryTable_Efficacy_NeutRatio_SD_SEM$PooledSD[1]
logk = tail(FittedLogistic_RawEfficacy_MeanRept_SDPool$estimate,2)[1]
# Fitted IC50 values
C50=tail(FittedLogistic_RawEfficacy_MeanRept_SDPool$estimate,1)
C50_Severe=tail(FittedLogistic_RawEfficacy_MeanRept_SDPool_Severe$estimate,1)
# Guessed IC50 values
# C50_Death = -1.8
# C50_Infection = -.3
# C50_Transmission = 0

# Make a lookup table and set the key of the table to the log10 neut ratios
# get_efficacies needs to be set to T the first time in an R session or if efficacies are not already defined

make_efficacies = F
save_csv = F
read_efficacies = T

log10_neut_ratio_range =seq(log10(0.001),log10(20),by=0.00001)

# Calculate data for lookup matrices
if (make_efficacies){ 
  efficacies_transmission = LogisticModel_PercentUninfected(log10_neut_ratio_range,sig,logk,C50_Transmission)
  efficacies_infection = LogisticModel_PercentUninfected(log10_neut_ratio_range,sig,logk,C50_Infection)
  efficacies_symptoms = LogisticModel_PercentUninfected(log10_neut_ratio_range,sig,logk,C50)
  efficacies_severe = LogisticModel_PercentUninfected(log10_neut_ratio_range,sig,logk,C50_Severe)
  efficacies_death = LogisticModel_PercentUninfected(log10_neut_ratio_range,sig,logk,C50_Death)
  # Make lookup data table
  log10_neut_efficacy_lookup = data.table(log10_neutR = log10_neut_ratio_range, efficacy1 = efficacies_transmission, efficacy2=efficacies_infection, efficacy3 = efficacies_symptoms, efficacy4 = efficacies_severe, efficacy5=efficacies_death)
  saveRDS(log10_neut_efficacy_lookup,'natmed_efficacylookup.RDS') 
}
if (save_csv){ write.csv(log10_neut_efficacy_lookup, 'natmed_efficacylookup.csv')}
if (read_efficacies){ log10_neut_efficacy_lookup=readRDS('natmed_efficacylookup.RDS')}

# Main function
estimate_current_efficacy = function(VE0, decay_rate, time_since_vaccination, efficacy_type){
  init_neut = get_neut_from_efficacy(VE0, efficacy_type) 
  current_neut = get_neut_over_time(init_neut, decay_rate,time_since_vaccination) 
  currentVE = get_efficacy_from_neut(current_neut, efficacy_type)
  currentVE
}

# Helper functions
get_neut_from_efficacy<-function(VE, efficacy_type=3){
  setkeyv(log10_neut_efficacy_lookup,paste0('efficacy',efficacy_type))
  loweff=log10_neut_efficacy_lookup[J(VE), roll = 1]
  higheff=log10_neut_efficacy_lookup[J(VE), roll = -1]
  
  # Indices where we haven't been given an exact value from the lookup table
  if(loweff[[paste0('efficacy',efficacy_type)]]!=VE & !is.na(loweff[[paste0('efficacy',efficacy_type)]])){
    log10neut = approx(c(loweff[[paste0('efficacy',efficacy_type)]],higheff[[paste0('efficacy',efficacy_type)]]),c(loweff$log10_neutR,higheff$log10_neutR), xout=VE)
  } else if (is.na(loweff$log10_neutR)){
    log10neut = -1000 
  } else {
    log10neut = loweff$log10_neutR
  }
  neut = 10^log10neut
  neut
}

get_efficacy_from_neut<-function(neut, efficacy_type=3){
  setkey(log10_neut_efficacy_lookup,log10_neutR)
  loweff=log10_neut_efficacy_lookup[J(log10(neut)), roll = 1]
  higheff=log10_neut_efficacy_lookup[J(log10(neut)), roll = -1]
  
  # Indices where we haven't been given an exact value from the lookup table
  if(loweff$log10_neutR!=log10(neut) & !is.na(loweff$log10_neutR)){
    eff = approx(c(loweff$log10_neutR,higheff$log10_neutR),c(loweff[[paste0('efficacy',efficacy_type)]],higheff[[paste0('efficacy',efficacy_type)]]), xout=log10(neut))
  } else if (is.na(loweff[[paste0('efficacy',efficacy_type)]])){
    if (loweff$log10_neutR < 0) {eff = 0} 
    else {eff=1}
  } else {
    eff = loweff[[paste0('efficacy',efficacy_type)]]
  }
  eff
}

get_neut_over_time = function (starting_neut, decay_rate, time){
  if (decay_rate > 0) {decay_rate = -decay_rate}
  current_neut = starting_neut*exp(decay_rate*time)
}


## This is the code from Khoury et.al. Nature Med and Cromer et. al. Lancet Microbe
####Logistic Model
ProbRemainUninfected=function(logTitre,logk,C50){1/(1+exp(-exp(logk)*(logTitre-C50)))}

LogisticModel_PercentUninfected=function(mu_titre,sig_titre,logk,C50){
  NumInteration<-max(length(mu_titre),length(sig_titre),length(logk),length(C50))
  Output<-NULL
  
  if (length(C50)==1) {
    C50=rep(C50,NumInteration)
  }
  
  if (length(logk)==1) {
    logk=rep(logk,NumInteration)
  }
  
  if (length(sig_titre)==1) {
    sig_titre=rep(sig_titre,NumInteration)
  }
  
  if (length(mu_titre)==1) {
    mu_titre=rep(mu_titre,NumInteration)
  }
  
  for (i in 1:NumInteration) {
    Step=sig_titre[i]*0.001
    IntegralVector=seq(mu_titre[i]-5*sig_titre[i],mu_titre[i]+5*sig_titre[i],by=Step)
    Output[i]=sum(ProbRemainUninfected(IntegralVector,logk[i],C50[i])*dnorm(IntegralVector,mu_titre[i],sig_titre[i]))*Step
  }
  Output
}

setwd(thisdir)
