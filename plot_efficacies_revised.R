#### Set directory to same directory as the r-script
currentdirectory<-dirname(rstudioapi::getActiveDocumentContext()$path)
if (!currentdirectory==""){ 
  setwd(currentdirectory) 
  currentdirectory='.'
}#source('./init-decay.R')
source('./ab_efficacy_estimates.R')

########
# Plot the estimate. Takes the following as input
# fitstruct - a fitting structure that contains
#   estimate - a named parameter vector that goes into the efficacy function to esitmate the efficacy
#       renormalisingFactor = determines the renormalisation from nat med params (multiplicative)
#       variantShift = The shift (up or down) for a variant. If not included, or 1, no shift is made (multiplicative)
#       severeShift = The shift (up or down) for severe infection (vs symptomatic). If not included, or 1, no severe estimate is included (multiplicative)
#       infectionShift = The shift (up or down) for infection only (vs symptomatic). If not included, or 1, no infection only estimate is included (multiplicative)
#       decay = The (early) decay rate of neutralising antibodies
#       late_decay = The late decay rate of neutralising antibodies
#       tswitch = The time that the switch is made between early and late decay.
#   data - the data frame against wihchi the fit was made
#   variants - These are the variants for which the estimate was (and will be) calculated, They should match the variants in the associated data frame
#   vaccines - These are the vaccines for which the estimate was (and will be) calculated, They should match the vaccines in the associated data frame
#   data - the data frame that includes the data that is to be plotted onto the fit
#   weighted - boolean stating whether a weighted fit was performed or not.
#   groupCol - Which variable should we group the data by when plotting
#
# The plotting function also takes the following input:
# makeTitle - boolean determining whether or not to include a (detailed) title for the plot
# plotBoosterDataInfo - boolean determining whether we fill data points according to booster type?

plot_eff_fit_from_struct_revised = function(fitstruct,makeTitle=T, plotBoosterDataInfo=F){
  effplot=plot_eff_fit_revised(fitstruct$estimate, fitstruct$data, weighted=fitstruct$weighted, groupCol = fitstruct$groupCol, variants = fitstruct$variants,vaccines = fitstruct$vaccines,makeTitle=makeTitle, plotBoosterDataInfo=plotBoosterDataInfo, nparams=length(fitstruct$paramestimate))
}

plot_eff_fit_revised<- function(fit_params, data, weighted=FALSE, groupCol = 'fitgroup', variants = '',  vaccines = '', makeTitle = T , plotBoosterDataInfo=F, basePlot=NULL,nparams=NULL) {
  #####
  # Check if we are meant to be plotting variants or variant groups. Usually it will be variant group
  # if all the variants input are actually from the Variant variable in the data, then we are ok, 
  # but if they aren't we might need to use VariantGroup 
  if((sum(variants %in% data$Variant)!=length(variants))| (length(variants)!=length(unique(data$Variants)))){
    if(sum(variants %in% data$VariantGroup)==length(variants)){
      data$Variant = data$VariantGroup
    }
  }
  #####
  # Get the estimates to plot
  eff_estimate = efficacy_from_paramvect(fit_params, groupCol=groupCol, variants=variants, vaccines = vaccines)
  eff_estimate_melted = reshape2::melt(eff_estimate, id.var = 'time', value.name='efficacy')
  eff_estimate_melted$variable = as.character(eff_estimate_melted$variable)
  eff_estimate_melted$join=eff_estimate_melted$variable
  eff_estimate_melted$EndpointCat = case_when(endsWith(eff_estimate_melted$variable,'Symptomatic')~'Symptomatic',
                                              endsWith(eff_estimate_melted$variable,'Severe')~'Severe',
                                              endsWith(eff_estimate_melted$variable,'Infection')~'Infection',
                                              T~'')
  eff_estimate_melted$EndpointCat = factor(eff_estimate_melted$EndpointCat, levels = c('Infection','Symptomatic','Severe',''))
  eff_estimate_melted$EndpointCat[eff_estimate_melted$EndpointCat=='']=unique(data$EndpointCat)[1]
  eff_estimate_melted$variable = str_remove(eff_estimate_melted$variable,as.character(eff_estimate_melted$Endpoint))
  
  #### check nparams entered
  if (is.null(nparams)){
    nparams = length(fit_params)
  }
  ### Internal functions
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
  
  #### 
  # Fix up vaccine and variants
  eff_estimate_melted$Vaccine = sapply(eff_estimate_melted$variable, myfun_vaccine)
  eff_estimate_melted$Vaccine = factor(eff_estimate_melted$Vaccine, levels = vaccine_list)
  eff_estimate_melted$Variant = sapply(eff_estimate_melted$variable, myfun_variant)
  eff_estimate_melted$Variant = factor(eff_estimate_melted$Variant, levels = variant_list)
  eff_estimate_melted_starts = eff_estimate_melted%>% group_by(Vaccine, Variant, EndpointCat, join) %>% summarise(VE0 = max(efficacy))
  eff_estimate_melted$boosted = 'Unboosted'
  eff_estimate_melted$boosted[eff_estimate_melted$Vaccine=='Boosted'] = 'Boosted'
  eff_estimate_melted$boosted = factor(eff_estimate_melted$boosted, levels = c('Unboosted','Boosted'))
  
  
  #### 
  # These are for the title
  aic = estimateAIC_fromparamsanddata(fit_params, data, groupCol, variants, vaccines)
  if (makeTitle){
    titleStr=make_plot_title_revised(fit_params, variants,aic)
  } else {
    titleStr = ''
  }
  #Make the base plot if not supplied
  variantsInPlot = unique(eff_estimate_melted$Variant)
  variantsInData = unique(data$Variant)
  # If there was no variants in the plot, add them in
  if (is.na(variantsInPlot)){
    eff_estimate_melted$Variant = variantsInData[1]
  }
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
  if (is.null(basePlot)){
    basePlot = make_eff_plot_revised(data, plotBoosterDataInfo)
  }
  effplot_fit = (basePlot
                 + labs(x='Time since vaccination (months)', y = 'Vaccine Efficacy', title = titleStr)
                 +scale_linetype_manual(values=c('solid','longdash','dotdash'), guide='none')
                 +theme(plot.title = element_text(hjust = 0.5))
                 #+ggtitle(TeX(titleStr)) 
                 #+geom_text(data = eff_estimate_melted_starts, mapping = aes(x=-5,y=VE0, label = round(100*VE0,1), colour=Vaccine, alpha = Variant))
  )
  # if(length(unique(eff_estimate_melted$Variant))<=1){
  #   effplot_fit=effplot_fit+facet_wrap(~EndpointCat)
  # } else{
    effplot_fit=effplot_fit+facet_grid(Variant~boosted+EndpointCat, switch='y')
  #}
  # If we are using boosted data, then only plot one line - the actual shift should be the pfizer shift multiplied by the renorm factor
  #if(plotBoostedLine){
  #  eff_estimate_melted = filter(eff_estimate_melted, Vaccine=='Pfizer')
  #  effplot_fit=effplot_fit+geom_line(data=eff_estimate_melted,mapping = aes(x=time/30, y=efficacy*100, linetype = Variant, alpha=Variant, group=join), size=1.5, inherit.aes = F, color='purple')
  #}
  #else {
    effplot_fit=effplot_fit+geom_line(data=filter(eff_estimate_melted,!is.na(Vaccine)),mapping = aes(x=time/30, y=efficacy*100, colour=Vaccine, linetype = Variant, alpha=Variant, group=join), size=1.5, inherit.aes = F)
  #}
  
}

make_eff_plot_revised<-function(data, withBoosting=F){
  ymin=0
  #shapes = c(11,13,15,17,19,21:25,1,2,3,4,5)
  #shapes = c(21:25,21:25)
  #FirstAuthors = unique(data$FirstAuthor)
  #names(shapes) = paste0(FirstAuthors,' et. al.')
  #shapes = shapes[c(1:length(FirstAuthors))]
  xmax = 8
  data = data %>% mutate(colourVar = Vaccine, fillVar = Vaccine)
  colourName = 'Vaccine'
  fillName = 'Vaccine'
  #shapes_infection
  if(withBoosting){
    data = data %>% mutate(colourVar = boosted, fillVar = boosted)
    fillName = 'Booster'
    colourName = 'Boosting'
  } 
  #check if there is a specific join variable for plotting - if not, just make it the join one
  if(!('plotjoin' %in% names(data))){
    data$plotjoin = data$join
  }
  if (!'CI_width' %in% colnames(data)){
    data$CI_width=100
  }
  efficacy_plot = (ggplot(data, aes(x=TimeMidpt, y=Efficacy, colour=colourVar,fill=fillVar,shape = paste0(FirstAuthor,' et. al.'), group=plotjoin))+
                   geom_point(mapping=aes(alpha=CI_width),size=3)+
                   # This additional point is a fudge so that we can get the legend to be present
                   #geom_point(data = filter(data,FirstAuthor ==names(study_shapes_eff)[1]),shape=study_shapes_eff[1], size=3)
                   geom_line(linetype='dashed')+
                   geom_errorbar(aes(ymin=EfficacyMin, ymax=EfficacyMax,alpha=CI_width), width=.03)+
                     
                   scale_colour_manual(values = vaccine_colours[names(vaccine_colours) %in% unique(data$colourVar)], name = colourName)+
                   scale_fill_manual(values = vaccine_colours[names(vaccine_colours) %in% unique(data$fillVar)], name=fillName) +
                   scale_shape_manual(values = study_shapes_eff[names(study_shapes_eff) %in% unique(paste0(data$FirstAuthor,' et. al.'))], name = 'Study') +
                   scale_alpha_continuous(range = c(.05,1), rescaler=alpharescale)+
                   
                   theme_bw()+
                   theme(text=element_text(size=14))+
                   #+ facet_grid(AgeGp~Vaccine, switch='y')
                   #+ylim(ymin,100)
                   #xlim(0,xmax)+
                   scale_y_continuous(breaks=seq(ymin, 100, 20))+ #, limits=c(ymin, 100)
                   coord_cartesian(ylim=c(ymin,100),xlim=c(0,xmax))+
                   labs(x='Time since vaccination (In Months)')
  )
  
}

make_shaded_plot = function(plotdata, boosted=F){
  plotdata = plotdata %>%  mutate(colourVar = Vaccine, fillVar = Vaccine)
  if(boosted){
    plotdata = plotdata %>% mutate(colourVar = boosted, fillVar = boosted)
  } 
  plot = make_eff_plot_revised(filter(plotdata, !is.na(Endpoint)),boosted)
  plot = plot+
    geom_ribbon(data=filter(plotdata, is.na(Endpoint)),mapping=aes(x=TimeMidpt,ymin=EfficacyMin,ymax=EfficacyMax,fill=fillVar),alpha=0.15,col=NA)+
    geom_line(data=filter(plotdata, is.na(Endpoint)),mapping=aes(x=TimeMidpt,y=Efficacy),size=1)
  if(boosted){
    plot = plot+facet_grid(EndpointCat~VariantGroup, switch='y')
  }else{
    plot = plot+facet_grid(Vaccine~VariantGroup, switch='y')
  }
  plot
}


make_plot_title_revised<-function(fit_params, variants,aic){
  #    endptstr = ifelse(endpt=='Infection',endpt,paste(endpt,'Disease'))
  endptstr = ''
  addFoldStr = function(str,folds,IC50s){
    if (length(folds)==0){
      str = ''
    } else {
      str1=''
      for (fold in folds){
        str1 = paste0(str1,', ',ifelse(fold>1,round(fold,2),-1*round(1/fold,2)))
      }
      str1 = substr(str1,3,nchar(str1))
      str1 = paste0(str1,' fold (IC50 = ')
      str2=''
      for (IC50 in IC50s){
        str2 = paste0(str2,', ',round(IC50*100,0))
      }
      str2 = substr(str2,3,nchar(str2))
      str2 =paste0(str2,'%)')
      str = paste0(str,str1, str2)
    }
    str
  }
  IC50=neut_ratio_from_efficacy(.5)
  decay = fit_params[names(fit_params)=='decay']
  
  renormalisingFactor = fit_params[names(fit_params) %in% c('renormalisingFactor','renorm_factor')]
  #if (plotBoostedLine){
  #  renormalisingFactor = renormalisingFactor&10^(log10NeutRatios[names(log10NeutRatios)=='Pfizer'])
  #}
  renormalisingIC50 = ifelse(renormalisingFactor<0,-IC50/renormalisingFactor,IC50*renormalisingFactor)
  renormalisingIC50s = c(IC50,renormalisingIC50)
  useBaseIC50 = renormalisingIC50
  useBaseIC50 = IC50
  variantShift = fit_params[names(fit_params)=='variantShift']
  variantIC50 = ifelse(variantShift<0,-useBaseIC50/variantShift,useBaseIC50*variantShift)
  variantShiftB = fit_params[names(fit_params)=='variantShiftB']
  variantIC50B = ifelse(variantShiftB<0,-useBaseIC50/variantShiftB,useBaseIC50*variantShiftB)
  variantShiftC = fit_params[names(fit_params)=='variantShiftC']
  variantIC50C = ifelse(variantShiftC<0,-useBaseIC50/variantShiftC,useBaseIC50*variantShiftC)
  variantShiftD = fit_params[names(fit_params)=='variantShiftD']
  variantIC50D = ifelse(variantShiftD<0,-useBaseIC50/variantShiftD,useBaseIC50*variantShiftD)
  variantShifts = c(variantShift,variantShiftB,variantShiftC,variantShiftD)
  variantIC50s=c(variantIC50,variantIC50B,variantIC50C,variantIC50D)
  variantShifts = c(1,variantShifts[!is.na(variantShifts)])
  variantIC50s = c(useBaseIC50,variantIC50s[!is.na(variantIC50s)])
  names(variantShifts) = tail(variants,-1)
  names(variantIC50s) = tail(variants,-1)
  severeShift = fit_params[names(fit_params)=='severeShift']
  if(severeShift>0){severeIC50s=variantIC50s/severeShift} else {severeIC50s=-variantIC50s*severeShift}
  infectionShift = fit_params[names(fit_params)=='infectionShift']
  if(infectionShift>0){infectionIC50s = variantIC50s/infectionShift}else {infectionIC50s=-variantIC50s*infectionShift}
  
  varStartString = paste0('\n Drop from ',variants[1],' to ',paste0(tail(variants,-1),collapse=', '),' = ')
  variantShiftStr=addFoldStr(varStartString,variantShifts, variantIC50s)
  
  renomalisingStr = addFoldStr('\n Renormalisation = ',renormalisingFactor,renormalisingIC50s)
  severeStr = addFoldStr('\n Severe Shift Factor = ',severeShift,severeIC50s)
  infectionStr = addFoldStr('\n Infection Shift Factor = ',infectionShift, infectionIC50s)
  
  titleStr = paste0('AIC=',round(aic,2),
                    '\n d = ', round(decay,4),' (half-life =' ,round(-log(2)/decay,1),'days)',
                    renomalisingStr,
                    variantShiftStr,
                    severeStr,
                    infectionStr)
} 

save_fit_plot<-function(name='fit_plot',w1=5,w2=2.5,h1=3,h2=3, savename=NULL){
  # get the structure
  if(class(name)=='character'){
    eval(parse(text=paste0('fitstruct=',name)))
    if (is.null(savename)){savename=name}
  } else{
    fitstruct=name
    if (is.null(savename)){savename=paste0('fit_',as.character(timeDate::Sys.timeDate()))}
  }
  fitstruct$plot = plot_eff_fit_from_struct_revised(fitstruct)
  ngpings = length(fitstruct$plot$facet$params$rows)
  nrows = 1
  if (ngpings>0){
    for (gp_no in c(1:ngpings)){
      gp_name = quo_name(fitstruct$plot$facet$params$rows[[gp_no]])
      nrows = nrows*length(unique(fitstruct$plot$data[,gp_name]))
    }
  }
  
  ngpings = length(fitstruct$plot$facet$params$cols)
  ncols = 1
  if (ngpings>0){
    for (gp_no in c(1:ngpings)){
      gp_name = quo_name(fitstruct$plot$facet$params$cols[[gp_no]])
      ncols = ncols*length(unique(fitstruct$plot$data[,gp_name]))
    }
  }
  
  ggsave(paste0('./Figures/',savename,'.pdf'),fitstruct$plot, width=(w1+ncols*w2), height=(h1+nrows*h2))
         
}

make_shaded_plot_formatted = function(plotdata, boosted=F){
  if (boosted){
    letterLabels = fct_cross(plotdata$EndpointCat,plotdata$VariantGroup)
    plotdata$letterLabel = paste(plotdata$EndpointCat,plotdata$VariantGroup,sep=':') 
  } else {
    letterLabels = fct_cross(plotdata$Vaccine,plotdata$EndpointCat,plotdata$VariantGroup)
    plotdata$letterLabel = paste(plotdata$Vaccine,plotdata$EndpointCat,plotdata$VariantGroup,sep=':') 
  }
  plotdata$letterLabel = factor(plotdata$letterLabel,levels=levels(letterLabels),labels=strsplit(intToUtf8(64+c(1:length(levels(letterLabels)))),"")[[1]])
  #levels(plotdata$letterLabel)=strsplit(intToUtf8(64+c(1:length(levels(letterLabels)))),"")[1]]
  plot = make_shaded_plot(plotdata, boosted)+decaypaper_theme
  # if(boosted){
  #   plot = plot+facet_wrap(~VariantGroup+EndpointCat, nrow=3)
  # }else{
  #   plot = plot+facet_wrap(~EndpointCat+VariantGroup+Vaccine, nrow=3)
  # }
  plot = plot+facet_wrap(~letterLabel, nrow=3)
}

decaypaper_theme = theme_classic()+
                    theme(legend.spacing = unit(.1,'pt'), 
                         legend.text = element_text(size=8), 
                         legend.key.height = unit(5,'mm'), 
                         legend.title = element_text(size=8),
                         strip.placement = 'outside',
                         strip.text = element_text(size=12, hjust=0, face ='bold', vjust = 0),
                         strip.background = element_rect(size=0, colour='white'),
                         text=element_text(size=12))
