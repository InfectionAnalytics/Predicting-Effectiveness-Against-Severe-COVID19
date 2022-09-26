library(data.table)

#### Set directory to same directory as the r-script
currentdirectory<-dirname(rstudioapi::getActiveDocumentContext()$path)
if (!currentdirectory==""){ setwd(currentdirectory) }
source('./init-decay.R')

#source('ab_efficacy_correlations_DK.R')

# Extract parameters that aren't going to change
load(paste0(basedir,'natmed_parameters.RData'))
sig=SummaryTable_Efficacy_NeutRatio_SD_SEM$PooledSD[1]
logk = tail(FittedLogistic_RawEfficacy_MeanRept_SDPool$estimate,2)[1]
C50=tail(FittedLogistic_RawEfficacy_MeanRept_SDPool$estimate,1)

C50_Severe=tail(FittedLogistic_RawEfficacy_MeanRept_SDPool_Severe$estimate,1)

# Make a lookup table and set the key of the table to the log10 neut ratios
# get_efficacies needs to be set to T the first time in an R session or if efficacies are not already defined
get_efficacies = F
read_efficacies = T
log10_neut_ratio_range =seq(log10(0.001),log10(100),by=0.00001)
if (get_efficacies){ 
  efficacies = LogisticModel_PercentUninfected(log10_neut_ratio_range,sig,logk,C50)
  saveRDS(efficacies,'efficacies.RDS')
}
if(read_efficacies){
  efficacies=readRDS('efficacies.RDS')
}

LogisticModel_PercentUninfected_lookup = data.table(neut_key = log10_neut_ratio_range, log10_neutR = log10_neut_ratio_range, efficacy = efficacies)
setkey(LogisticModel_PercentUninfected_lookup, neut_key)

InitialTitre_func_lookup=data.table(eff_key = efficacies, log10_neutR = log10_neut_ratio_range, efficacy = efficacies)
setkey(InitialTitre_func_lookup, eff_key)

calc_eff <-function(log10_nr){
  loweff=LogisticModel_PercentUninfected_lookup[J(log10_nr), roll = 1]
  loweff$log10_neutR[is.na(loweff$efficacy)]=loweff$neut_key[is.na(loweff$efficacy)]
  loweff$efficacy[is.na(loweff$efficacy)]=0
  
  eff = loweff$efficacy
  
  higheff=LogisticModel_PercentUninfected_lookup[J(log10_nr), roll = -1]
  higheff$log10_neutR[is.na(higheff$efficacy)] = higheff$neut_key[is.na(higheff$efficacy)]
  higheff$efficacy[is.na(higheff$efficacy)]=1
  
  # Indices where we haven't been given an exant neut ration from the lookup table
  nonExactInds = loweff$log10_neutR!=log10_nr & !is.na(loweff$log10_neutR)
  if (sum(nonExactInds)>0){
    eff[nonExactInds] = approx(c(loweff$log10_neutR[nonExactInds],higheff$log10_neutR[nonExactInds]),c(loweff$efficacy[nonExactInds],higheff$efficacy[nonExactInds]), xout=log10_nr[nonExactInds])$y
  }
  eff[is.na(loweff$log10_neutR)] = 0
  eff
}
eff_from_neut_ratio <-function(nr){
  eff = nr
  if (!is.null(nrow(nr)) & !is.null(ncol(nr))){
    if(nrow(nr)>1 & ncol(nr)>1){
      option=1
    } else{
      option=2
    }
  }
  else{
    option=2
  }
  
  if (option==1){
    for (c in 1:ncol(nr)){
      #eff[,c]=LogisticModel_PercentUninfected(log10(nr[,c]),sig,logk,C50)
      # Do as a lookup instead (really should use apply)
      loweff=LogisticModel_PercentUninfected_lookup[J(log10(nr[,c])), roll = 1]
      # This is all just  to try to get to the bottom of errors which occured when decay rate is too fast and so nr drops too low.
      if (sum(is.na(loweff$log10_neutR))>0){
        #rint('Warning')
      }
      eff[,c] = tryCatch(calc_eff(log10(nr[,c])), 
                         error = function(cond){
                           print(paste0('Eff has ',nrow(eff),' rows and ',ncol(eff),' cols'))
                           print(paste0('nr has ',nrow(nr),' rows and ',ncol(nr),' cols'))
                           print(paste0('We are up to column ',c))
                           loweff=LogisticModel_PercentUninfected_lookup[J(log10(nr[,c])), roll = 1]
                           eff = loweff$efficacy
                           higheff=LogisticModel_PercentUninfected_lookup[J(log10(nr[,c])), roll = -1]
                           # Indices where we haven't been given an exant neut ration from the lookup table
                           nonExactInds = loweff$log10_neutR!=log10(nr[,c])
                           print('loweff$log10_neutR=')
                           print(loweff$log10_neutR)
                           #print(log10(nr[,c]))
                           #print(nonExactInds)
                           print('Efficacy=')
                           print(eff)
                           #print(paste0('There are ',sum(nonExactInds),' non exact inds'))
                           #a=calc_eff(log10(nr[,c]))
                           #print(paste0('Calculated eff has ',nrow(a),' rows and ',ncol(a),' cols'))        
                           message(cond)
                           
                         }
      )
    }
  } else {
    eff=calc_eff(log10(nr))
    if(!is.null(colnames(nr))){
      eff = as.matrix(eff)
      colnames(eff)=colnames(nr)
    }
  }
  eff
} 


eff_over_time_from_neutR<-function(start_neut_ratio, decay_rate_neut_ratio, tstart = 0, tend=180, as.column=T, tswitch=NULL, late_decay=NULL){
  start_eff = round(eff_from_neut_ratio(start_neut_ratio),2)*100
  total_length = length(start_neut_ratio)*length(decay_rate_neut_ratio)
  start_neut_ratio_new=rep(start_neut_ratio,length(decay_rate_neut_ratio))
  start_eff_new=rep(start_eff,length(decay_rate_neut_ratio))
  decay_rate_neut_ratio=rep(decay_rate_neut_ratio,each=length(start_neut_ratio))
  if(!is.null(late_decay)){
    late_decay=rep(late_decay,each=length(start_neut_ratio))
  }
  #decay_rate_neut_ratio[decay_rate_neut_ratio>0]=-decay_rate_neut_ratio[decay_rate_neut_ratio>0]
  start_neut_ratio = start_neut_ratio_new
  start_eff = start_eff_new
  rm(start_eff_new,start_neut_ratio_new)
  
  eff_over_time=NULL
  time = seq(tstart, tend, by =1)
  
  # requires matrix multiplication
  neut_ratios = matrix(rep(start_neut_ratio,each=length(time)),ncol=length(start_neut_ratio))*exp(tcrossprod(time,decay_rate_neut_ratio))
  
  # This is the code where a shift in decay rates is implemented.
  if(is.null(tswitch)){tswitch=83}
  if (is.null(late_decay)){late_decay=decay_rate_neut_ratio}
  
  neut_ratios1 = matrix(rep(start_neut_ratio,each=length(time[time<=tswitch])),ncol=length(start_neut_ratio))*exp(tcrossprod(time[time<=tswitch],decay_rate_neut_ratio))
  neut_ratios2 = matrix(rep(neut_ratios1[nrow(neut_ratios1),],each=length(time[time>tswitch])),ncol=length(start_neut_ratio))*exp(tcrossprod(time[time>tswitch]-tswitch,late_decay))
  neut_ratios = rbind(neut_ratios1,neut_ratios2)
  colnames(neut_ratios) = paste0('VE0_',start_eff,'_d_',decay_rate_neut_ratio)
  eff_over_time = data.frame(time=time)
  calc_eff = eff_from_neut_ratio(neut_ratios)
  eff_over_time = cbind(eff_over_time,calc_eff)
  
  #colnames(eff_over_time)[2:ncol(eff_over_time)] = paste('VE0',round(eff_over_time[1,2:7],2)*100,'d',-decay_rate_neut_ratio,sep='_')
  colnames(neut_ratios) = paste0('neutR_',colnames(neut_ratios))
  effneut_over_time = cbind(eff_over_time,neut_ratios)
  
  if (as.column){  
    #Melt the dataframe
    eff_over_time=melt(effneut_over_time, id.vars=c('time',colnames(effneut_over_time)[startsWith(colnames(effneut_over_time), 'neutR')]), value.name='efficacy')
    eff_over_time=rename(eff_over_time,c('eff0_decay'='variable'))
    eff_over_time$eff0 = substr(eff_over_time$eff0_decay,5,6)
    eff_over_time$decay = substr(eff_over_time$eff0_decay,10,20)
    neut_over_time=melt(effneut_over_time, id.vars=c('time',colnames(effneut_over_time)[startsWith(colnames(effneut_over_time), 'VE0')]), value.name='neutR')
    neut_over_time$eff0 = substr(neut_over_time$variable,11,12)
    neut_over_time$decay = substr(neut_over_time$variable,16,26)
    effneut_over_time = merge(select(eff_over_time,time,efficacy,eff0,decay),select(neut_over_time,time,neutR,eff0,decay), by=c('time','eff0','decay'))
  }
  effneut_over_time
}

neut_ratio_from_efficacy<-function(eff_for_lookup){
  neut_ratio = NULL
  if(eff_for_lookup[1]>1){eff_for_lookup=eff_for_lookup/100 }
  #Also do this as a lookup
  neut_ratio = 10^(InitialTitre_func_lookup[J(eff_for_lookup), roll = "nearest"]$log10_neutR)
  # for (i in 1:length(efficacy)) {
  #   log10_neutR=nlm(function(x){InitialTitre_func(x,efficacy[i])},0)
  #   neut_ratio[i] = 10^log10_neutR$estimate
  # }
  neut_ratio
}

eff_over_time_from_VE<-function(start_eff, decay_rate_neut_ratio, tstart = 0, tend=180,as.column=T, tswitch=NULL, late_decay=NULL){
  neutR<-neut_ratio_from_efficacy(start_eff)
  eff_over_time_from_neutR(neutR, decay_rate_neut_ratio,tstart, tend,as.column, tswitch, late_decay)
}

efficacy_from_fit_params_natmed = function(params ,tf=270, groupCol='fitgroup', variants=c('Non-Delta','Delta'), vaccines = c('Pfizer','Moderna','Astra')){
  renormalisingFactor = params[1]
  decay = params[2]
  variantShift = params[names(params)=='variantShift']
  severeShift = params[names(params)=='severeShift']
  infectionShift = params[names(params)=='infectionShift']
  tswitch = params[names(params)=='tswitch']
  late_decay=params[names(params)=='late_decay']
  
  if(length(variantShift)==0){variantShift=NA}
  if(length(severeShift)==0){severeShift=NA}
  if(length(infectionShift)==0){infectionShift=NA}
  if(length(tswitch)==0){tswitch=NULL}
  if(length(late_decay)==0){late_decay=NULL}
  
  log10NeutRatios = SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutRatio_cens
  names(log10NeutRatios) = SummaryTable_Efficacy_NeutRatio_SD_SEM$Study
  names(log10NeutRatios)[names(log10NeutRatios)=='AstraZeneca']='Astra'
  log10NeutRatios = log10NeutRatios[names(log10NeutRatios) %in% vaccines]+log10(renormalisingFactor)
  neutRatios= 10^log10NeutRatios
  eff0 = eff_from_neut_ratio(neutRatios)
  names(eff0) = names(neutRatios)
  
  # Add in the Severe / Infection Shifts
  shifts = c(infectionShift,severeShift)
  shift_names = c('Symptomatic','Infection','Severe')
  use_shifts = !is.na(shifts)
  shifts = shifts[use_shifts]
  shift_names = shift_names[c(T,use_shifts)]
  if (length(shifts)>0){
    eff0 = extractEff0UsingShift(eff0, shifts, shift_names, nameposition = 'end')
  }
  # not sure this is really necessesary
  eff0 = reorderEff0ByEndpoint(eff0, shift_names)
  # Decide if we need to add in the variant Shift
  if(!is.na(variantShift) & length(variants)>=2){
    usevariantShift = T
    effparams = c(eff0,variantShift,decay)
  } else {
    usevariantShift = F
    effparams=c(eff0,decay)
  }
  
  # Make sure none of the efficacies start above 1
  est_eff = estimate_efficacy(eff0,decay,c(eff0names,'decay'), tf, NULL, groupCol, variantShift, variants, tswitch=tswitch, late_decay=late_decay)
  
}

eff_error_natmed_neutratios = function(params ,tf=270, eff_data_to_fit=NULL, groupCol, variants, vaccines){
  renormalisingFactor = params[1]
  decay = params[2]
  variantShift = params[names(params)=='variantShift']
  severeShift = params[names(params)=='severeShift']
  infectionShift = params[names(params)=='infectionShift']
  tswitch = params[names(params)=='tswitch']
  late_decay=params[names(params)=='late_decay']
  
  if(length(variantShift)==0){variantShift=NA}
  if(length(severeShift)==0){severeShift=NA}
  if(length(infectionShift)==0){infectionShift=NA}
  if(length(tswitch)==0){tswitch=NULL}
  if(length(late_decay)==0){late_decay=NULL}
  
  log10NeutRatios = SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutRatio_cens
  names(log10NeutRatios) = SummaryTable_Efficacy_NeutRatio_SD_SEM$Study
  names(log10NeutRatios)[names(log10NeutRatios)=='AstraZeneca']='Astra'
  log10NeutRatios = log10NeutRatios[names(log10NeutRatios) %in% vaccines]+log10(renormalisingFactor)
  neutRatios= 10^log10NeutRatios
  eff0 = eff_from_neut_ratio(neutRatios)
  names(eff0) = names(neutRatios)
  
  # Add in the Severe / Infection Shifts
  shifts = c(infectionShift,severeShift)
  shift_names = c('Symptomatic','Infection','Severe')
  use_shifts = !is.na(shifts)
  shifts = shifts[use_shifts]
  shift_names = shift_names[c(T,use_shifts)]
  if (length(shifts)>0){
    eff0 = extractEff0UsingShift(eff0, shifts, shift_names, nameposition = 'end')
  }
  # not sure this is really necessesary
  eff0 = reorderEff0ByEndpoint(eff0, shift_names)
  # Decide if we need to add in the variant Shift
  if(!is.na(variantShift) & length(variants)>=2){
    usevariantShift = T
    effparams = c(eff0,variantShift,decay)
  } else {
    usevariantShift = F
    effparams=c(eff0,decay)
  }
    
  err = eff_error(effparams,weighted=T,eff0names=names(eff0), tf, eff_data_to_fit,groupCol, usevariantShift, variants, tswitch=tswitch,late_decay=late_decay )
  err
  
}

reorderEff0ByEndpoint = function(eff0, ordered_shift_names){
  eff0
}

eff_error<-function(eff0_decay,weighted=FALSE,eff0names=NULL, tf=270, eff_data_to_fit,groupCol = 'VariantVaccine', usevariantShift = F, variants='', tswitch=NULL,late_decay=NULL){
  eff0 = head(eff0_decay,-1)
  decay = tail(eff0_decay,1)
  if(usevariantShift){
    # At the moment this only works with one variant
    variantShift=tail(eff0,1)
    eff0 = head(eff0,-1)
    
  } else{
    variantShift=NA
  }
  
  # Make sure none of the efficacies start above 1
  above1 = eff0>=1
  eff0[above1]=.999
  est_eff = estimate_efficacy(eff0,decay,c(eff0names,'decay'), tf, eff_data_to_fit, groupCol, variantShift, variants,tswitch,late_decay)
  if (weighted){
    sigmas = (eff_data_to_fit$EfficacyMax-eff_data_to_fit$EfficacyMin)/100/3.92
    sigmas[is.na(sigmas)]=mean(sigmas[!is.na(sigmas)])
  } else {
    sigmas = rep(1,nrow(eff_data_to_fit)) 
  }
  penalty = (5/mean(sigmas))^2
  err_eff = sum(((eff_data_to_fit$Efficacy/100-est_eff)/sigmas)^2, na.rm=T)+penalty*sum(above1)
  
  #Print one value every 5 seconds
  t=round(as.numeric(Sys.time()))
  modVal = 20
  if (t %% modVal ==0){
    if (printed != t){
      #print(est_eff)
      print(paste0('Error = ',round(err_eff,2)))
      printed<<-t
    } 
  }
  #
  err_eff
}

extractEff0UsingShift = function(eff0, shifts, shift_names=c('Non-Delta','Delta'), nameposition = 'start'){
  # This means we are using a variant drop and we need to work on including it
  # Now accounts for multiple drops
  neutR0 = neut_ratio_from_efficacy(eff0)
  eff0_final = eff0
  neutR0_final = neutR0
  if (nameposition=='start'){ names(eff0_final) = paste0(shift_names[1],names(eff0))}
  else if (nameposition=='end') { names(eff0_final) = paste0(names(eff0),shift_names[1]) }
  for (i in c(1:length(shifts))){
    shift = shifts[i]
    log10shift = log10(shift)
    neutR0_new = 10^(log10(neutR0)+log10shift)
    eff0_new = eff_from_neut_ratio(neutR0_new)
    if (nameposition=='start'){ names(eff0_new) = paste0(shift_names[i+1],names(eff0))}
    else if (nameposition=='end') { names(eff0_new) = paste0(names(eff0),shift_names[i+1]) }
    eff0_final = c(eff0_final,eff0_new)
  }
  
  eff0=eff0_final
}
estimate_efficacy<-function(eff0,decay,paramnames, tf=270, eff_data_to_fit=NULL, groupCol = 'VariantVaccine', variantShift=NA,variants='', tswitch=NULL, late_decay=NULL){
  months_to_days_factor = 30
  est_eff = NULL
  if(!is.na(variantShift)){
    eff0 = extractEff0UsingShift(eff0, variantShift,variants)
  } else {
    names(eff0) = head(paramnames,-1)
  }
  eff=eff_over_time_from_VE(eff0, decay, tend=tf, as.column = F, tswitch=tswitch, late_decay=late_decay)
  calceff0=eff[1,startsWith(colnames(eff),'VE0')]
  #print(paste0('max(abs(calceff0-eff0))=',max(abs(calceff0-eff0))))
  if (is.na(max(abs(calceff0-eff0)))){
    print('max(abs(calceff0-eff0)) is NA')
  }
  if (max(abs(calceff0-eff0))>1e-4){
    nomatch=1
  }
  colnames(eff)[startsWith(colnames(eff),'VE0')]=names(eff0)
  if(!is.null(eff_data_to_fit)){
    for (i in 1:nrow(eff_data_to_fit)){
      effInds = which(eff$time>=eff_data_to_fit$TimeMin[i]*months_to_days_factor & eff$time<=eff_data_to_fit$TimeMax[i]*months_to_days_factor)
      colInd = which(colnames(eff)==eff_data_to_fit[i,groupCol])
      # if we can't find the column by variant / vaccine then find the one called efficacy
      if(length(colInd)==0){
        colInd = which(colnames(eff)==paramnames[1])
      }
      est_eff[i] = mean(eff[effInds,colInd])
    }
  } else {
    est_eff=eff
  }
    
  est_eff
}

addAIC<-function(fit, calcmin=NULL, nparams=NULL){
  fit$N = nrow(fit$data)
  if (is.null(calcmin)){minval = fit$minimum}
  if (is.null(nparams)){nparams = length(fit$estimate)}
  else {minval = calcmin}
  fit$AIC = fit$N*log(minval/fit$N)+2*nparams
  fit
}

perform_fit<-function(data, initParams, weighted=FALSE, endptcat='Infection',groupCol = 'VariantVaccine', usevariantShift=F, variants='',maxiter = 1000,nreps = 10){
  if (is.null(names(initParams))){
    if(length(initParams)>2){
      names(initParams) = c(paste0('Gp',c(1:length(initParams)-1)),'decay')
    } else {
      names(initParams) = c("efficacy","decay") 
    }
  }
  # need to start from other points
  fitfunction <-function(eff0_d){eff_error(eff0_d,weighted,names(initParams),tf=270, eff_data_to_fit=data, groupCol, usevariantShift, variants)}
  steptolval = 1e-8
  stepmaxval=1
  gradtolval =1e-9
  
  for(i in c(1:nreps)){
    #fit<-nlm(fitfunction,initParams, stepmax=stepmaxval, steptol=steptolval,gradtol = gradtolval, iterlim = maxiter)
    if(usevariantShift){
      low = c(rep(0,length(initParams)-2),-Inf,-1)
      high = c(rep(1,length(initParams)-2),Inf,0)
    } else {
      low = c(rep(0,length(initParams)-1),-1) 
      high = c(rep(1,length(initParams)-1),0)
    }
    fit<-optim(initParams, fitfunction,method="L-BFGS-B", lower = low, upper = high)
    saveRDS(fit,'./intermediateFit.RDS')
    if(is.null(fit$estimate)){fit$estimate = fit$par}
    if(is.null(fit$iterations)){fit$iterations = fit$counts[1]}
    if(is.null(fit$code)){fit$code = fit$convergence}
    if(is.null(fit$minimum)){fit$minimum = fit$value}
    initParams = replace(initParams, c(1:length(initParams)),fit$estimate)
    if (i==1){
      iter = fit$iterations
    } else {
      iter = iter + fit$iterations -1
    }
    print(paste0('Completed iter ',i,' of 1st fit. Used ',fit$iterations,' iterations. Code = ',fit$code,'. Error = ',round(fit$minimum,1)))
    saveRDS(fit,'./intermediateFit.RDS')
  }
  fit$iterations = iter
  
  # Now try better start params
  trybetter=F
  usefit = 'First'
  if (trybetter){
    if (length(initParams)-1 == nrow(distinct(select(data,groupCol))) ){
      starteffs = data %>% group_by(VariantVaccine) %>% summarise(starteff = max(Efficacy))
      initParams = replace(initParams,starteffs$VariantVaccine,starteffs$starteff)
    } else if (length(initParams)==2){
      initParams[1] = max(data$Efficacy)
    }
    
    print('Next try')
    for(i in c(1:nreps)){
      fit2<-nlm(fitfunction,initParams, stepmax=stepmaxval, steptol=steptolval,gradtol = gradtolval, iterlim = maxiter)
      initParams = replace(initParams, c(1:length(initParams)),fit2$estimate)
      print(paste0('Completed iter ',i,' of 2nd fit. Used ',fit2$iterations,' iterations. Code = ',fit2$code,'. Error = ',round(fit2$minimum,1)))
      if (i==1){
        iter = fit2$iterations
      } else {
        iter = iter + fit2$iterations -1
      }
    }
    fit2$iterations = iter
    
    if (fit2$minimum<fit$minimum){
      fit = fit2
      usefit = 'Second'
    } 
  }
  fit = addInfoToFit(fit,usevariantShift,data,variants, endptcat, groupCol)
  saveRDS(fit, './perform_fit_result.RDS')
  print(paste0('Using ',usefit,' fit. Took ',fit$iterations,' iterations, Error = ',round(fit$minimum,1),' AIC = ',fit$AIC))
  
  fit
}
addInfoToFit = function(fit,usevariantShift=NULL,data=NULL,variants=NULL,endptcat=NULL,initParams=NULL,vaccines=NULL,natmed_neutratios=F, includePlot=T, groupCol=NULL, nparams=NULL){
  if (!is.null(initParams)){names(fit$estimate) = names(initParams)}
  if(!is.null(fit$gradient)){names(fit$gradient) = names(estimate)}
  if(!is.null(data)){fit$data = data}
  if (!is.null(endptcat)){fit$endptcat = endptcat}
  if (!is.null(groupCol)){fit$groupCol = groupCol}
  if (is.null(nparams)){nparams = length(fit$estimate)}
  #fit$calcminval = eff_error_natmed_neutratios(fit$estimate,tf=270, fit$data, fit$groupCol, fit$variants, fit$vaccines)
  #fit$calcn = eff_error_natmed_neutratios(fit$estimate,tf=270, fit$data, fit$groupCol, fit$variants, fit$vaccines)
  #fit<-addAIC(fit, fit$calcminval, nparams)
  if(!is.null(usevariantShift)){
    fit$usevariantShift = usevariantShift
  } else{
    usevariantShift = fit$usevariantShift
  }
    
  fit$group_names = names(fit$estimate)
  if(!is.null(variants)){fit$variants = variants}
  if(!is.null(vaccines)){fit$vaccines = vaccines}
  
  if(usevariantShift){
    variantShiftdecay = tail(fit$estimate,2)
    fit$decay = tail(variantShiftdecay,1)
    fit$variantShift = head(variantShiftdecay,1)
    fit$eff0 = extractEff0UsingShift(head(fit$estimate,-2), fit$variantShift, fit$variants)
  } else {
    fit$variantShift=0
    fit$eff0 = head(fit$estimate,-1)
    fit$decay = tail(fit$estimate,1)
  }
  fit$neut0 = neut_ratio_from_efficacy(fit$eff0)
  names(fit$neut0)=names(fit$eff0)
  #fit$weighted = fitstruct$weighted
  #endpt = fit$endptcat
  #groupCol = fit$groupCol
  if (includePlot){
    fit = addplot(fit,natmed_neutratios)
  }
  fit
}

addplot = function(fit,natmed_neutratios=F){
  fit$plot = plot_eff_fit_from_struct(fit,natmed_neutratios)
  fit
}

perform_fit_natmed_neutratios<-function(data, initParams, endptcat='Infection',groupCol = 'VariantVaccine',variants=c('Non-Delta','Delta'), vaccines = c('Astra','Pfizer','Moderna'),maxiter = 1000,nreps = 10){
  #maxiter = 10
  nreps = 10
  nfits=20
  weighted=T
  useVariantShift = T
  # need to start from other points
  fitfunction <-function(params){eff_error_natmed_neutratios(params ,tf=270, eff_data_to_fit=data, groupCol, variants, vaccines)}
  steptolval = 1e-8
  stepmaxval=1
  gradtolval =1e-9
  
  low = c(0,-1,0,0,0)
  high = c(Inf,0,Inf,Inf,Inf)
  low = low[c(1:length(initParams))]
  high = high[c(1:length(initParams))]
  startParams = initParams
  usefit = NULL
  errs = NULL
  startParams = initParams
  fittedParams = NULL
  for (fitno in c(1:nfits)) { 
    if (fitno > 1) {
      initParams =replace(initParams, c(1:length(initParams)), c(rexp(1,.4),runif(1,-.1,0),rexp(3,.4)))
      startParams = rbind(allparams,initParams)
    }
    for(i in c(1:nreps)){
      
      fit<-optim(initParams, fitfunction,method="L-BFGS-B", lower = low, upper = high)
      saveRDS(fit,'./intermediateFit.RDS')
      if(is.null(fit$estimate)){fit$estimate = fit$par}
      if(is.null(fit$iterations)){fit$iterations = fit$counts[1]}
      if(is.null(fit$code)){fit$code = fit$convergence}
      if(is.null(fit$minimum)){fit$minimum = fit$value}
      initParams = replace(initParams, c(1:length(initParams)),fit$estimate)
      if (i==1){
        iter = fit$iterations
      } else {
        iter = iter + fit$iterations -1
      }
      print(paste0('Completed iter ',i,' of Fit ',fitno,'. Used ',fit$iterations,' iterations. Code = ',fit$code,'. Error = ',round(fit$minimum,1)))
      saveRDS(fit,'./intermediateFit.RDS')
      if (fitno == 1){
        fittedParams = fit$estimate
      } else {
        fittedParams = rbind(fittedParams,fit$estimate)
      }
      errs = c(errs,fit$minimum)
    }
    fit$iterations = iter
    if (fitno == 1){
      usefit = fit
    } else if (fit$minimum < usefit$minimum) {
      usefit = fit
    }
  }
  fit = usefit
  fit$fittedParams=fittedParams
  fit$errs=errs
  fit$startParams = startParams
  saveRDS(fit, './perform_fit_result.RDS')
  
  # Now try better start params
  # trybetter=5
  # usefit = 1
  # while (trybetter>0){
  #   fitno = fitno+1
  #   newInitParams = c(rexp(1,.4),runif(1,-.1,0),rexp(3,.4))
  #   newInitParams = newInitParams[c(1:length(initParams))]
  #   print(paste('Trying fit ',fitno,' with starting params of ',paste0(round(newInitParams,2),collapse=', ')))
  #   for(i in c(1:nreps)){
  #     fit2<-optim(newInitParams, fitfunction,method="L-BFGS-B", lower = low, upper = high)
  #     if(is.null(fit2$estimate)){fit2$estimate = fit2$par}
  #     if(is.null(fit2$iterations)){fit2$iterations = fit2$counts[1]}
  #     if(is.null(fit2$code)){fit2$code = fit2$convergence}
  #     if(is.null(fit2$minimum)){fit2$minimum = fit2$value}
  #     saveRDS(fit2,'./intermediateFit.RDS')
  #     newInitParams = replace(newInitParams, c(1:length(newInitParams)),fit2$estimate)
  #     print(paste0('Completed iter ',i,' of fit ',fitno,'. Used ',fit2$iterations,' iterations. Code = ',fit2$code,'. Error = ',round(fit2$minimum,1)))
  #     if (i==1){
  #       iter = fit2$iterations
  #     } else {
  #       iter = iter + fit2$iterations -1
  #     }
  #   }
  #   fit2$iterations = iter
  #   
  #   if (fit2$minimum<fit$minimum){
  #     fit = fit2
  #     usefit = fitno
  #   } 
  #   trybetter=trybetter-1
  # }
  #fit = readRDS('./intermediateFit.RDS')
  fit = addInfoToFit(fit,useVariantShift,data,variants, endptcat,vaccines =vaccines,natmed_neutratios=T)
  saveRDS(fit, './perform_fit_result.RDS')
  print(paste0('Using ',usefit,' fit. Took ',fit$iterations,' iterations, Error = ',round(fit$minimum,1),' AIC = ',fit$AIC))
  
  fit
}

# Estimate the standard deviation in the log10 neut ratio (normally distributed) from the 95% CIs in the efficacy data
estimate_sd_log10neutR<-function(eff_init, eff_lower95, eff_upper95,CIlim = .95){
  # Efficacies 
  effs = c(eff_lower95,eff_init,eff_upper95)
  # Estimate the mean and SD of the log10 neut ratios
  log10neutRs = log10(neut_ratio_from_efficacy(effs))
  mean_log10neutR = log10neutRs[2]
  sd_log10neutR = mean(diff(log10neutRs)/qnorm((1+CIlim)/2))
  
  # Now take this back and see if this is the best estimate
  sd_log10neutR_est = nlm(function(s,m){sum((calc_eff(m+s*c(qnorm((1-CIlim)/2),0,qnorm((1+CIlim)/2)))-effs)^2)},p=sd_log10neutR,m=mean_log10neutR)
}

eff_over_time_full <- function(params,ind=NA,tf=270){
  sigma = params[1]
  hill=params[2]
  IC50 = params[3]
  log10neutR = params[4]
  decay = params[5]
  eff_time = data.frame(time = c(0:tf))
  neutRs = (10^log10neutR)*exp(eff_time$time*decay)
  log10neutRs = log10(neutRs)
  eff_time$efficacy = LogisticModel_PercentUninfected(log10neutRs,sigma,hill,IC50)  
  eff_time$sigma = sigma
  eff_time$hill = hill
  eff_time$IC50 = IC50
  eff_time$log10neutR = log10neutR
  eff_time$decay = decay
  eff_time$ind=ind
  eff_time
}

  
  