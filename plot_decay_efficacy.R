effplot=(ggplot(efficacy_data[efficacy_data$HighFidelity.DC==1,], aes(x=TimeMidpt, y=Efficacy, colour=VariantVaccine, fill=VariantVaccine,shape = Paper, group=join))
         +geom_point(size=3)
         +geom_line(linetype='dashed')
         + geom_errorbar(aes(ymin=EfficacyMin, ymax=EfficacyMax), width=.03)
         + geom_errorbar(aes(xmin=TimeMin, xmax=TimeMax), width=.03)
         +scale_colour_manual(values = variantvaccine_colours)
         +scale_fill_manual(values = variantvaccine_colours)
         +scale_shape_manual(values = study_shapes_eff)
         +theme_bw()
         +ylim(0,100)
)
effplot_all = effplot
effplot_all$data = efficacy_data
  effplot_all$data[,c("Efficacy","EfficacyMin","EfficacyMax")][effplot_all$data[,c("Efficacy","EfficacyMin","EfficacyMax")]<0]=0

#neutplot=(ggplot(neut_data, aes(x=Day, y=Neut, colour=Vaccine, fill=Vaccine,shape = Paper, group=join))
#         +geom_point(size=3)
#         +geom_line(linetype='dashed')
#         + geom_errorbar(aes(ymin=NeutMin, ymax=NeutMax), width=.03)
#         +scale_colour_manual(values = vaccine_colours)
#         +scale_fill_manual(values = vaccine_colours)
#         +scale_shape_manual(values = study_shapes_neut)
#         +scale_y_log10()
#         +theme_bw()
#)

plot_eff_fit_from_struct = function(fitstruct,natmed_neutratios=T, makeTitle=T,withBoosting=F){
  effplot=plot_eff_fit(fitstruct, fitstruct$data, group_names=names(fitstruct$estimate),weighted=fitstruct$weighted,endpt = fitstruct$endptcat, groupCol = fitstruct$groupCol,usevariantShift=fitstruct$usevariantShift, variants = fitstruct$variants,vaccines = fitstruct$vaccines,natmed_neutratios=natmed_neutratios,makeTitle=makeTitle, withBoosting=withBoosting)
}


# Plot the estimate
plot_eff_fit<- function(eff_fit, data, group_names=NULL,origplot=NULL, weighted=FALSE,endpt = '', groupCol = 'VariantVaccine',usevariantShift=F, variants='', vaccines = '',natmed_neutratios=F, makeTitle = T, withBoosting=F) {
  if (natmed_neutratios){
    params=eff_fit$estimate
    #eff_estimate = efficacy_from_fit_params_natmed(eff_fit$estimate,groupCol='fitgroup', variants=eff_fit$variants, vaccines =eff_fit$vaccines)
    eff_estimate = efficacy_from_paramvect(eff_fit$estimate, groupCol='fitgroup', variants=eff_fit$variants, vaccines = eff_fit$vaccines)

    eff0 = eff_estimate[1,]
    eff0 = tail(eff0,-1)
    eff0 = eff0[!startsWith(names(eff0),'neutR')]
    renormalisingFactor = eff_fit$estimate[1]
    decay = eff_fit$estimate[2]
    variantShift = params[names(params)=='variantShift']
    severeShift = params[names(params)=='severeShift']
    infectionShift = params[names(params)=='infectionShift']
    
  } else { 
    renormalisingFactor=NA
    if(natmed_neutratios){
      renormalisingFactor = eff_fit$estimate[1]
      decay = eff_fit$estimate[2]
      if (length(eff_fit$estimate)>2){
        variantShift = eff_fit$estimate[3]
      } else {
        variantShift = 0
      }
      
      log10NeutRatios = SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutRatio_cens
      names(log10NeutRatios) = SummaryTable_Efficacy_NeutRatio_SD_SEM$Study
      names(log10NeutRatios)[names(log10NeutRatios)=='AstraZeneca']='Astra'
      log10NeutRatios = log10NeutRatios[names(log10NeutRatios) %in% vaccines]+log10(renormalisingFactor)
      neutRatios= 10^log10NeutRatios
      eff0 = eff_from_neut_ratio(neutRatios)
      names(eff0) = names(neutRatios)
      eff0=extractEff0UsingShift(eff0, variantShift, variants)
    } else if(usevariantShift){
      variantShiftdecay = tail(eff_fit$estimate,2)
      decay = tail(variantShiftdecay,1)
      variantShift = head(variantShiftdecay,1)
      eff0 = extractEff0UsingShift(head(eff_fit$estimate,-2), variantShift, variants)
    } else {
      variantShift=0
      eff0 = head(eff_fit$estimate,-1)
      decay = tail(eff_fit$estimate,1)
    }
    eff_estimate  = eff_over_time_from_VE(eff0,decay,tend = 270, as.column=F)
    colnames(eff_estimate)[startsWith(colnames(eff_estimate),'VE0')]=names(eff0)
  }
  
  # if (is.null(group_names)) {
  #   if(is.null(names(eff_fit$estimate))){
  #     group_names = c("efficacy", "neutR")
  #     names(group_names) = group_names
  #   } else {
  #     group_names = names(eff_fit$estimate)
  #     names(group_names) = group_names
  #   }
  # }
  
  #colnames(eff_estimate)[startsWith(colnames(eff_estimate),'VE0')]=head(group_names,-1)
  # The line below actually isnt needed any more
  eff_estimate = subset(eff_estimate, select=which(!startsWith(colnames(eff_estimate), 'neutR')&!startsWith(colnames(eff_estimate), 'variant')&!startsWith(colnames(eff_estimate), 'decay')))
  
  eff_estimate_melted = reshape2::melt(eff_estimate, id.var = 'time', value.name='efficacy')
  eff_estimate_melted$variable = as.character(eff_estimate_melted$variable)
  eff_estimate_melted$join=eff_estimate_melted$variable
  eff_estimate_melted$EndpointCat = case_when(endsWith(eff_estimate_melted$variable,'Symptomatic')~'Symptomatic',
                                           endsWith(eff_estimate_melted$variable,'Severe')~'Severe',
                                           endsWith(eff_estimate_melted$variable,'Infecton')~'Infection',
                                           T~'')
  eff_estimate_melted$EndpointCat = factor(eff_estimate_melted$EndpointCat, levels = c('Infection','Symptomatic','Severe',''))
  eff_estimate_melted$EndpointCat[eff_estimate_melted$EndpointCat=='']=unique(origplot$data$EndpointCat)[1]
  eff_estimate_melted$variable = str_remove(eff_estimate_melted$variable,as.character(eff_estimate_melted$Endpoint))
  myfun_vaccine = function(x){
    vaxStart = NA
    for (var in variant_list){
      varStart = str_locate(x,var)[2]
      if (!is.na(varStart)){
        break
      }
    }
    if(is.na(varStart)){varStart=0}
    substr(x,varStart+1,100)
  }
  myfun_variant = function(x){
    varStart = NA
    for (var in variant_list){
      varStart = str_locate(x,var)[2]
      if (!is.na(varStart)){
        break
      }
    }
    if (is.na(varStart)){variant='None'}
    else {variant=substr(x,1,varStart)}
    variant
  }
  eff_estimate_melted$Vaccine = sapply(eff_estimate_melted$variable, myfun_vaccine)
  eff_estimate_melted$Vaccine = factor(eff_estimate_melted$Vaccine, levels = vaccine_list)
  eff_estimate_melted$Variant = sapply(eff_estimate_melted$variable, myfun_variant)
  eff_estimate_melted$Variant = factor(eff_estimate_melted$Variant, levels = variant_list)
  eff_estimate_melted_starts = eff_estimate_melted%>% group_by(Vaccine, Variant, EndpointCat, join) %>% summarise(VE0 = max(efficacy))
  eff_estimate_melted$VariantGroup = eff_estimate_melted$Variant
  endptstr = ifelse(endpt=='Infection',endpt,paste(endpt,'Disease'))
  addFoldStr = function(str,fold,IC50){
    if (length(fold)==0){
      str = ''
    } else {
      str1 = paste0(ifelse(fold>1,paste0(round(fold,2),'-fold'),paste0('-',round(1/fold,1),'-fold')),' (IC50 = ',paste0(round(IC50*100,0),'%',collapse='/'),' of Conv against WT)')
      str = paste0(str,str1)
    }
    str
  }
  IC50=neut_ratio_from_efficacy(.5)
  renormalisingIC50 = IC50/renormalisingFactor
  variantIC50 = c(renormalisingIC50,renormalisingIC50/variantShift)
  severeIC50 = variantIC50/severeShift
  infectionIC50 = variantIC50/infectionShift
  
  variantShiftStr=addFoldStr(paste0('\n Drop from ',variants[1],' to ',variants[2],' = '),variantShift, variantIC50)
  renomalisingStr = addFoldStr('\n Normalisation Factor = ',renormalisingFactor,renormalisingIC50)
  severeStr = addFoldStr('\n Severe Shift Factor = ',severeShift,severeIC50)
  infectionStr = addFoldStr('\n Infection Shift Factor = ',infectionShift, infectionIC50)
  eff_fit = addAIC_fromfitanddata(eff_fit)
  if(natmed_neutratios){isweighted=T}
  if(makeTitle) {
    titleStr = paste0(ifelse(isweighted,'Weighted','Unweighted'),' fit to ',endptstr,'\n AIC=',round(eff_fit$AIC,2),
                    '\n d = ', round(decay,4),' (half-life =' ,round(-log(2)/decay,1),'days)',
                    renomalisingStr,
                    variantShiftStr,
                    severeStr,
                    infectionStr)
  }else {
    titleStr = ''
  }
    #Make the base plot if not supplied
  variantsInPlot = unique(eff_estimate_melted$Variant)
  variantsInData = unique(data$Variant)
  variantsNotInData = variantsInPlot[!(variantsInPlot %in% variantsInData)]
  # add dummy variants into data if they aren't there
  if (length(variantsNotInData)>0) { if (!is.na(variantsNotInData)){
    for (i in c(1:length(variantsNotInData))){
      data = rbind(data[1,],data)
      data[1,c("Efficacy","EfficacyMin","EfficacyMax")]=0
      oldVariant = data[1,c("Variant")]
      data[1,sapply(data[1,],is.character)]=sapply(data[1,sapply(data[1,],is.character)], str_replace_all, as.character(oldVariant), as.character(variantsNotInData[i]))
      data[1,'Variant'] = variantsNotInData[i]
    }
  }}
  if (is.null(origplot)){
    origplot = make_eff_plot_new(data, withBoosting)
  }
  effplot_fit = (origplot
                 + labs(x='Time since vaccination (months)', y = 'Vaccine Efficacy', title = titleStr)
                 +facet_wrap(~EndpointCat)
                 +scale_linetype_manual(values=c('solid','longdash'), guide='none')
                 +theme(plot.title = element_text(hjust = 0.5))
                 #+ggtitle(TeX(titleStr)) 
                 #+geom_text(data = eff_estimate_melted_starts, mapping = aes(x=-5,y=VE0, label = round(100*VE0,1), colour=Vaccine, alpha = Variant))
  )
  
  if(withBoosting){
    eff_estimate_melted = filter(eff_estimate_melted, Vaccine=='Pfizer')
    effplot_fit=effplot_fit+geom_line(data=eff_estimate_melted,mapping = aes(x=time/30, y=efficacy*100, linetype = Variant, alpha=Variant, group=join), size=1.5, inherit.aes = F, color='purple')
  }
  else {
    effplot_fit=effplot_fit+geom_line(data=eff_estimate_melted,mapping = aes(x=time/30, y=efficacy*100, colour=Vaccine, linetype = Variant, alpha=Variant, group=join), size=1.5, inherit.aes = F)
  }  
}

make_eff_plot<-function(data){
  these_colours=variantvaccine_colours[names(variantvaccine_colours) %in% distinct(select(data,c('VariantVaccine')))[,1]]
  these_shapes = study_shapes_eff[names(study_shapes_eff) %in% distinct(select(data,c('Paper')))[,1]]
  this_eff_plot=(ggplot(data, aes(x=TimeMidpt, y=Efficacy, colour=VariantVaccine, fill=VariantVaccine,shape = Paper, group=join))
           +geom_point(size=3)
           +geom_line(linetype='dashed')
           + geom_errorbar(aes(ymin=EfficacyMin, ymax=EfficacyMax), width=.03)
           + geom_errorbar(aes(xmin=TimeMin, xmax=TimeMax), width=.03)
           +scale_colour_manual(values = these_colours)
           +scale_fill_manual(values = these_colours)
           +scale_shape_manual(values = study_shapes_eff)
           +theme_bw()
           +ylim(0,100)
  )

}

make_eff_plot_new<-function(data, withBoosting=F){
  ymin=0
  shapes = c(21:25,11,13,15,17)
  #shapes = c(21:25,21:25)
  FirstAuthors = unique(data$FirstAuthor)
  names(shapes) = paste0(FirstAuthors,' et. al.')
  shapes = shapes[c(1:length(FirstAuthors))]
  xmax = 7
  #shapes_infection
  if(withBoosting){
    data$fillVar = data$Booster
    fillName = 'Booster'
  } else {
    data$fillVar = data$Vaccine
    fillName = 'Vaccine'
  }
  efficacy_plot = (ggplot(data, aes(x=TimeMidpt, y=Efficacy, colour=Vaccine,fill=fillVar,shape = paste0(FirstAuthor,' et. al.'),alpha=Variant, group=join))
                   +geom_point(size=3)
                   # This additional point is a fudge so that we can get the legend to be present
                   +geom_point(data = filter(data,FirstAuthor ==FirstAuthors[1]),shape=21, size=3)
                   +geom_line(linetype='dashed')
                   + geom_errorbar(aes(ymin=EfficacyMin, ymax=EfficacyMax), width=.03)
                   +scale_colour_manual(values = vaccine_colours[names(vaccine_colours) %in% unique(data$Vaccine)])
                   +scale_fill_manual(values = vaccine_colours[names(vaccine_colours) %in% unique(data$fillVar)], name=fillName)
                   +scale_shape_manual(values = shapes, name = 'Study')
                   +scale_alpha_discrete(range=c(1,.5))
                   +theme_bw()
                   +theme(text=element_text(size=14))
                   #+ facet_grid(AgeGp~Vaccine, switch='y')
                   #+ylim(ymin,100)
                   +xlim(0,xmax)
                   + scale_y_continuous(breaks=seq(ymin, 100, 20), limits=c(ymin, 100))
                   +labs(x='Time since vaccination (In Months)')
  )
  
}

# pfizer_moderna_neutplot=(ggplot(filter(neut_data, Paper %in% c('DECAY07', 'DECAY12'), Variant=='WT'), aes(x=Day, y=Neut, colour=Vaccine, fill=Vaccine,shape = Cohort, group=join))
#           +geom_point(size=3)
#           +geom_line(linetype='dashed')
#           + geom_errorbar(aes(ymin=NeutMin, ymax=NeutMax), width=.03)
#           +scale_colour_manual(values = c('Pfizer'='blue','Moderna'='red'), labels=c('Pfizer (FDA report)','Moderna (Doria-Rose NEJM)'))
#           +scale_fill_manual(values = c('Pfizer'='blue','Moderna'='red'), guide = 'none')
#           +scale_shape_manual(values = c('18-55'=19, '56-70'=2, 'over71'=6,'65-85'=11), labels = c('18-55', '56-70', 'over 71','65-85'))
#           +scale_y_log10()
#           +theme_bw()
# )
# pfizer_moderna_neutplot
