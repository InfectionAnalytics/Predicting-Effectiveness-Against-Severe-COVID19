#source('./plot_efficacies_revised.R')
#source('./load_bootstrap_parameter_distributions.R')
vaccine_comparisons <- list( c("Moderna", "Pfizer"), c("Pfizer", "Astra"), c("Moderna", "Astra") )
variant_comparisons <- list( c("Non-Delta", "Delta"), c("Delta", "Omicron"), c("Non-Delta", "Omicron") )

spacesep="            "
variantAlpha = .6
#testmethod='wilcox.test'
#scattercomparisons_plotdata=filter(rw_unboosted_data,Vaccine !='mRNA') 

progressionColour='darkorange'
progressionShape=19

######## DATA ########
makeData=T 
if (makeData){
# Set up the data we will need for plotting
# This is to plot predicted vs observed
source('./load_bootstrap_parameter_distributions.R')
deltaDrop = mean_variants[1]
omicronDrop = mean_variants[2]
#mean_decay = log(2)/108
renormFactor = exp(-trialdays*mean_decay/2)
# Matching real world data points with predicted efficacies
sympt_summary_curve=filter(efficacy_summary_curve, EndpointCat=='Symptomatic')
severe_summary_curve=filter(efficacy_summary_curve, EndpointCat=='Severe')

rw_unboosted_symptsev_data_with_predict = rw_unboosted_symptsev_data %>% mutate(log10neutR=case_when(
  Vaccine=='Pfizer' ~log10NeutRatios['Pfizer'],
  Vaccine=='Moderna'~log10NeutRatios['Moderna'],
  Vaccine=='Astra' ~log10NeutRatios['Astra'],
  Vaccine=='mRNA' ~log10NeutRatios['mRNA'])+
    case_when(VariantGroup=='Delta'~deltaDrop,
              VariantGroup=='Omicron'~omicronDrop,
              T~0)-log10(renormFactor)-
    abs(mean_decay)*TimeMidpt*30*log10(exp(1)),
  predictedEfficacy = case_when(
    EndpointCat=='Symptomatic' ~ spline(sympt_summary_curve$log10neutR,sympt_summary_curve$eff_mean,xout=log10neutR)$y,
    EndpointCat=='Severe' ~ spline(severe_summary_curve$log10neutR,severe_summary_curve$eff_mean,xout=log10neutR)$y,
    T~0),
  predictedEfficacyMin = case_when(
      EndpointCat=='Symptomatic' ~ spline(sympt_summary_curve$log10neutR,sympt_summary_curve$eff_lower,xout=log10neutR)$y,
      EndpointCat=='Severe' ~ spline(severe_summary_curve$log10neutR,severe_summary_curve$eff_lower,xout=log10neutR)$y,
      T~0),
  predictedEfficacyMax= case_when(
    EndpointCat=='Symptomatic' ~ spline(sympt_summary_curve$log10neutR,sympt_summary_curve$eff_upper,xout=log10neutR)$y,
    EndpointCat=='Severe' ~ spline(severe_summary_curve$log10neutR,severe_summary_curve$eff_upper,xout=log10neutR)$y,
    T~0),
  log10neutR_var = case_when(Vaccine=='Pfizer' ~sd_log10neutR['Pfizer']^2,
                             Vaccine=='Moderna'~sd_log10neutR['Moderna']^2,
                             Vaccine=='Astra' ~sd_log10neutR['Astra']^2,
                             Vaccine=='mRNA' ~sd_log10neutR['mRNA']^2)+
                    case_when(VariantGroup=='Delta'~ sd_variants[1]^2,
                             VariantGroup=='Omicron'~sd_variants[2]^2,
                             T~0)+(trialdays/2 + TimeMidpt*30*log10(exp(1)))*sd_decay^2, # add in the variance of the neuts
  log10neutR_sd = sqrt(log10neutR_var),
  Protection = Efficacy, lower=EfficacyMin, upper = EfficacyMax # this is for later plots
)

#matchedRWandPredictDataSymptomatic = rw_unboosted_symptsev_data_with_predict %>% 
#  filter(EndpointCat=='Symptomatic')
  
#matchedRWandPredictDataSevere = rw_unboosted_symptsev_data_with_predict %>% 
#  filter(EndpointCat=='Severe') #%>%
  #mutate(predictEffMedian = spline(severe_summary_curve$log10neutR,severe_summary_curve$eff_median,xout=log10neutR)$y,
  #       predictEffMin=spline(severe_summary_curve$log10neutR,severe_summary_curve$eff_lower,xout=log10neutR)$y,
   #      predictEffMax=spline(severe_summary_curve$log10neutR,severe_summary_curve$eff_upper,xout=log10neutR)$y  )

matchedRWandPredictData = rw_unboosted_symptsev_data_with_predict %>%
  mutate(dataMin = EfficacyMin, dataMax = EfficacyMax, dataEff = Efficacy, 
         modelMin = predictedEfficacyMin, modelMax = predictedEfficacyMax, modelEff = predictedEfficacy,
         CIsoverlap = (dataMin<modelMin&dataMax>modelMin)|(dataMax>modelMax&dataMin<modelMax)|(dataMin>modelMin&dataMax<modelMax)|(dataMax>modelMax&dataMin<modelMin), 
         dataCIsOverCurve = (dataMin<modelEff & dataMax>modelEff),
         modelCIsOverData = (modelMin<dataEff & modelMax>dataEff), 
         datamodeldiff=dataEff-modelEff, 
         CIwidth=dataMax-dataMin,
         dataOutsideCI = pmin(dataEff-modelMin, modelMax-dataEff)
  )
matchedRWandPredictDataSymptomatic = matchedRWandPredictData %>% 
  filter(EndpointCat=='Symptomatic')

matchedRWandPredictDataSevere = matchedRWandPredictData %>% 
  filter(EndpointCat=='Severe')


efficacy_summary_curve_for_logistic_plots = efficacy_summary_curve %>% 
  mutate(Protection = eff_median, lower = eff_lower, upper = eff_upper) %>%
  select(c(log10neutR, neutR, Protection, lower, upper, EndpointCat)) 

letterLabelValues = unique(as.character(rw_unboosted_symptsev_data_with_predict$EndpointCat))
letterLabelNames=letterLabelValues
letterLabelNames[1:2]=paste(letterLabelValues[1:2],'Efficacy')
#letterLabelNames[3]=paste(letterLabelValues[3],'to Severe')
spacesep=" "
rw_unboosted_symptsev_data_with_predict$letterLabel = factor(rw_unboosted_symptsev_data_with_predict$EndpointCat,levels=(letterLabelValues),paste(strsplit(intToUtf8(64+c(1:length(letterLabelValues))),"")[[1]],letterLabelNames,sep=spacesep))
efficacy_summary_curve_for_logistic_plots$letterLabel = factor(efficacy_summary_curve_for_logistic_plots$EndpointCat,levels=(letterLabelValues),paste(strsplit(intToUtf8(64+c(1:length(letterLabelValues))),"")[[1]],letterLabelNames,sep=spacesep))

sympt_summary_curve=filter(efficacy_summary_curve, EndpointCat=='Symptomatic')
severe_summary_curve=filter(efficacy_summary_curve, EndpointCat=='Severe')
effList = c(1:100)
symptsevere_summary_curve = data.frame(symptomaticEfficacy=effList) %>%
  mutate(log10neutR = spline(sympt_summary_curve$eff_median,sympt_summary_curve$log10neutR,xout=symptomaticEfficacy)$y, 
         severeEfficacy = spline(severe_summary_curve$log10neutR,severe_summary_curve$eff_median,xout=log10neutR)$y,
         severeEfficacyMin=spline(severe_summary_curve$log10neutR,severe_summary_curve$eff_lower,xout=log10neutR)$y,
         severeEfficacyMax=spline(severe_summary_curve$log10neutR,severe_summary_curve$eff_upper,xout=log10neutR)$y  )



}
######## LEGENDS ########
makeLegends=T
if (makeLegends){
# Legend creation
variantlegend = get_legend( ggplot(filter(rw_unboosted_data, EndpointCat%in% c('Symptomatic','Severe')), aes(x=VariantGroup, y=Efficacy))+
   geom_point(mapping = aes(fill=VariantGroup),shape=22,size=3.5,colour='white',alpha=.5)+
   geom_point(mapping = aes(shape = studyname, colour=Vaccine))+
   
   scale_colour_manual(values = vaccine_colours[names(vaccine_colours) %in% unique(rw_unboosted_data$Vaccine)],labels =  vaccine_labels[names(vaccine_labels) %in% unique(rw_unboosted_data$Vaccine)], name = 'Vaccine')+
   scale_fill_manual(values = variant_fills, labels = str_replace_all(names(variant_fills),'Non','pre'), name='Variant')+
   scale_shape_manual(values = study_shapes_eff[names(study_shapes_eff) %in% unique(rw_unboosted_data$studyname)], name = 'Study')+
     
   decaypaper_theme)
nonvariantlegend = get_legend( ggplot(filter(rw_unboosted_data, EndpointCat%in% c('Symptomatic','Severe')), aes(x=VariantGroup, y=Efficacy, alpha=CI_width))+
  geom_point(mapping = aes(fill=VariantGroup),shape=22,size=3.5,colour='white',alpha=.5)+
  geom_point(mapping = aes(shape = studyname, colour=Vaccine))+
  
  scale_colour_manual(values = vaccine_colours[names(vaccine_colours) %in% unique(rw_unboosted_data$Vaccine)],labels =  vaccine_labels[names(vaccine_labels) %in% unique(rw_unboosted_data$Vaccine)], name = 'Vaccine')+
  scale_fill_manual(values = variant_fills, labels = str_replace_all(names(variant_fills),'Non','pre'), name='Variant', guide='none')+
  scale_shape_manual(values = study_shapes_eff[names(study_shapes_eff) %in% unique(rw_unboosted_data$studyname)], name = 'Study')+
  scale_alpha_continuous(range = c(.05,1), rescaler=alpharescale, name = 'Width of confidence interval')+
  
  decaypaper_theme)
studylegend = get_legend( ggplot(filter(rw_unboosted_data, EndpointCat%in% c('Symptomatic','Severe')), aes(x=VariantGroup, y=Efficacy, alpha=CI_width))+
                                 geom_point(mapping = aes(fill=VariantGroup),shape=22,size=3.5,colour='white',alpha=.5)+
                                 geom_point(mapping = aes(shape = studyname, colour=Vaccine))+
                                 
                                 scale_colour_manual(values = vaccine_colours[names(vaccine_colours) %in% unique(rw_unboosted_data$Vaccine)],labels =  vaccine_labels[names(vaccine_labels) %in% unique(rw_unboosted_data$Vaccine)], name = 'Vaccine', guide='none')+
                                 scale_fill_manual(values = variant_fills, labels = str_replace_all(names(variant_fills),'Non','pre'), name='Variant', guide='none')+
                                 scale_shape_manual(values = study_shapes_eff[names(study_shapes_eff) %in% unique(rw_unboosted_data$studyname)], name = 'Study')+
                                 scale_alpha_continuous(range = c(.05,1), rescaler=alpharescale, name = 'Width of confidence interval', guide='none')+
                                 
                                 decaypaper_theme)

nonstudylegend = get_legend( ggplot(filter(rw_unboosted_data, EndpointCat%in% c('Symptomatic','Severe')), aes(x=VariantGroup, y=Efficacy, alpha=CI_width))+
                                 geom_point(mapping = aes(fill=VariantGroup),shape=22,size=3.5,colour='white',alpha=.5)+
                                 geom_point(mapping = aes( colour=Vaccine))+
                                 
                                 scale_colour_manual(values = vaccine_colours[names(vaccine_colours) %in% unique(rw_unboosted_data$Vaccine)],labels =  vaccine_labels[names(vaccine_labels) %in% unique(rw_unboosted_data$Vaccine)], name = 'Vaccine')+
                                 scale_fill_manual(values = variant_fills, labels = str_replace_all(names(variant_fills),'Non','pre'), name='Variant')+
                                 scale_alpha_continuous(range = c(.05,1), rescaler=alpharescale, name = 'Width of confidence interval')+
                                 
                                 decaypaper_theme)
}

######## REGRESSION #######
doRegression=T
if (doRegression){
  regressionData = filter(rw_unboosted_symptsev_data_with_predict, EndpointCat=='Severe') %>%
  select(c(Paper,Vaccine,VariantGroup,TimeMidpt, Efficacy, EfficacyMin, EfficacyMax))
  vaccine_variant_time_regression_severe=lm(Efficacy~Vaccine+VariantGroup+TimeMidpt+VariantGroup:TimeMidpt, weights = sqrt(EfficacyMax-EfficacyMin)/1.96, data=regressionData)
  severe_predicts = predict(vaccine_variant_time_regression_severe, interval = "confidence")
  #severe_predicts[severe_predicts>100]=100
  #severe_predicts[severe_predicts<0]=0

}
######## PLOTS OVER TIME (BLOCKS) ########
plotsOverTime=T
if (plotsOverTime){
make_eff_over_time_plot = function (plotdata){
  vaccines = unique(plotdata$Vaccine)
  variants = unique(plotdata$VariantGroup)
  dummytable = data.table(Vaccine = rep(vaccines,times=length(variants)),VariantGroup=rep(variants, each = length(vaccines)), EndpointCat = unique(plotdata$EndpointCat)[1], Endpoint = unique(plotdata$EndpointCat)[1])
  dummytable$Vaccine = factor(dummytable$Vaccine, levels = levels(plotdata$Vaccine))
  plotdata = bind_rows(plotdata,dummytable)
  letterLabels = fct_cross(plotdata$Vaccine,plotdata$EndpointCat,plotdata$VariantGroup)
  plotdata$letterLabel = paste(plotdata$Vaccine,plotdata$EndpointCat,plotdata$VariantGroup,sep=':') 
  plotdata$letterLabel = factor(plotdata$letterLabel,levels=levels(letterLabels),labels=strsplit(intToUtf8(64+c(1:length(levels(letterLabels)))),"")[[1]])

  plotdata = plotdata %>%  mutate(colourVar = Vaccine, fillVar = Vaccine)
  plot = make_eff_plot_revised(filter(plotdata, !is.na(Endpoint)),F)
  #plot = plot+facet_grid(Vaccine~VariantGroup, switch='y')
  #  geom_ribbon(data=filter(plotdata, is.na(Endpoint)),mapping=aes(x=TimeMidpt,ymin=EfficacyMin,ymax=EfficacyMax,fill=fillVar),alpha=0.15,col=NA)+
  #  geom_line(data=filter(plotdata, is.na(Endpoint)),mapping=aes(x=TimeMidpt,y=Efficacy),size=1)
  plot = plot+decaypaper_theme+facet_wrap(~letterLabel, nrow=3)
  plot
}

severe_data = filter(rw_unboosted_symptsev_data_with_predict, EndpointCat=='Severe') %>%
  select(c(Paper,FirstAuthor,Vaccine,VariantGroup,TimeMidpt, Endpoint,EndpointCat,Efficacy, EfficacyMin, EfficacyMax, plotjoin,join)) %>%
  cbind(severe_predicts)
severe_unboosted_plot = make_eff_over_time_plot(severe_data)+
  facet_wrap(~letterLabel)+
  geom_ribbon(mapping=aes(ymin = lwr, ymax = upr, group=paste0(VariantGroup,Vaccine)),col=NA, alpha=.15)+
  geom_line(mapping=aes(y=fit, group=paste0(VariantGroup,Vaccine)), size=1)
severe_unboosted_plot_formatted = plot_grid(severe_unboosted_plot+ theme(legend.position="none"),
                                  nonvariantlegend,
                                  rel_widths=c(4,1.2))
ggsave(paste0(figuresdir,'severe_unboosted_plot_formatted_line.pdf'),severe_unboosted_plot_formatted,width=8, height=6)
#ggsave(paste0(manuscriptfiguresdir, 'Fig2new_R.pdf'),severe_unboosted_plot_formatted,width=8.35, height=5.85)
ggsave(paste0(manuscriptfiguresdir, 'Fig2_R.pdf'),severe_unboosted_plot_formatted,width=8, height=6)


severe_unboosted_plot2 = make_shaded_plot_formatted(filter(rw_and_sim_efficacy_data, EndpointCat=='Severe'))#,boosted=='Unboosted'))
severe_unboosted_plot_formatted2 = plot_grid(severe_unboosted_plot2+ theme(legend.position="none"),
                                           nonvariantlegend,
                                           rel_widths=c(4,1.2))
ggsave(paste0(figuresdir,'severe_unboosted_plot_formatted.pdf'),severe_unboosted_plot_formatted2,width=8, height=6)
ggsave(paste0(manuscriptfiguresdir,'FigS4_R.pdf'),severe_unboosted_plot_formatted2,width=8, height=6)


sympt_unboosted_plot = make_shaded_plot_formatted(filter(rw_and_sim_efficacy_data, EndpointCat=='Symptomatic'))#,boosted=='Unboosted'))
sympt_unboosted_plot_formatted = plot_grid(sympt_unboosted_plot+ theme(legend.position="none"),
                                            nonvariantlegend,
                                            rel_widths=c(4,1.2))
ggsave(paste0(figuresdir,'sympt_unboosted_plot_formatted.pdf'),sympt_unboosted_plot_formatted,width=8, height=6)
ggsave(paste0(manuscriptfiguresdir,'FigS3_R.pdf'),sympt_unboosted_plot_formatted,width=8, height=6)

}
######## COMPARISON SIGMOID CURVES ########
comparisonSigmoid=T
if (comparisonSigmoid){
# Curves - all vaccines / variants on one plot
# Now do the all vaccine plot
rw_vs_model_predict_sympt_sev_plot = ggplot(filter(rw_unboosted_symptsev_data_with_predict, EndpointCat %in%c('Symptomatic','Severe')), aes(x=log10neutR, y=Protection))+
  #geom_line(data=efficacy_summary_curve,mapping=aes(x=log10neutR,y=eff_median), colour='darkred', size=1,alpha=1)+
  #geom_ribbon(data=efficacy_summary_curve,mapping=aes(x=log10neutR,ymin=eff_lower, ymax=eff_upper,y=eff_median), fill='darkred',alpha=.15,colour=NA)+
  geom_line(data=filter(efficacy_summary_curve_for_logistic_plots, EndpointCat %in%c('Symptomatic','Severe')), colour='darkred', size=1,alpha=1)+
  geom_ribbon(data=filter(efficacy_summary_curve_for_logistic_plots, EndpointCat %in%c('Symptomatic','Severe')),mapping = aes(ymin = lower, ymax = upper), fill='darkred',alpha=.15,colour=NA)+
  
  geom_point(mapping = aes(fill=VariantGroup,alpha=CI_width),shape=22,size=3.5, colour='white')+
  scale_alpha_continuous(range = c(.05,variantAlpha), rescaler=alpharescale, guide=F)+
  new_scale("alpha") +
  
  geom_point(mapping = aes(shape = studyname, colour=Vaccine,alpha=(CI_width)))+
  #geom_errorbar(aes(ymin=EfficacyMin, ymax=EfficacyMax,colour=Vaccine,alpha=(CI_width)))+
  #geom_errorbar(mapping=aes(xmin=log10neutR-log10neutR_sd*1.96, xmax=log10neutR+log10neutR_sd*1.96,colour=Vaccine,alpha=(CI_width)))+
  geom_errorbar(mapping=aes(alpha = CI_width, colour=Vaccine,ymin = lower, ymax = upper))+
  geom_errorbar(mapping=aes(xmin=log10neutR-log10neutR_sd*1.96, xmax=log10neutR+log10neutR_sd*1.96,colour=Vaccine,alpha=50))+
  
  scale_alpha_continuous(range = c(.05,1), rescaler=alpharescale, name='Width of CI')+
  
  scale_colour_manual(values = vaccine_colours[names(vaccine_colours) %in% unique(rw_unboosted_data$Vaccine)], name = 'Vaccine')+
  scale_fill_manual(values = variant_fills, name='Variant')+
  scale_shape_manual(values = study_shapes_eff[names(study_shapes_eff) %in% unique(rw_unboosted_data$studyname)], name = 'Study', guide='none')+
  scale_x_continuous(breaks=log10(4^seq(-3,3,1)), labels = MASS::fractions(4^seq(-3,3,1)))+
  scale_y_continuous(breaks=seq(0,100,25), labels = paste0(seq(0,100,25),'%'))+ #limits = c(0,100), 
  coord_cartesian(ylim=c(0,100))+
  
  stat_cor(method='spearman',cor.coef.name='rho',label.sep="\n", label.x=-1.2,label.y=5)+
  
  decaypaper_theme+
  facet_wrap(~letterLabel)+
  labs(x='Neutralising Antibody Levels (fold of Convalescent)', y='Vaccine Effectiveness')

ggsave(paste0(figuresdir,'rw_vs_model_predict_sympt_sev_plot.pdf'),rw_vs_model_predict_sympt_sev_plot,width=8, height=4)

multiplier = 1
#symptsevere_usedata = simulated_full_data_melted %>% 
#  select(c('sympt_efficacy', 'severe_efficacy'))%>%
#  mutate(sympt_efficacy = round(multiplier*sympt_efficacy,1)/multiplier) %>% 
#  group_by(sympt_efficacy) %>%
#  summarise(lower = quantile(severe_efficacy,.025), mid = quantile(severe_efficacy,.5),upper = quantile(severe_efficacy,.975))

#symptsevere_matchedplotdata = data.frame(t(sapply(symptsevere_usedata,quantile,c(.5,.025,.975)))) %>% mutate(EndpointCat = 'Symptomatic')


symptvssevereplot = ggplot(rw_matched_symptsev_data, aes(y=SevereEff, x=SymptomaticEff))+
  geom_point(mapping = aes(fill=VariantGroup),shape=22,size=3.5,colour='white',alpha=.5)+
  geom_point(mapping = aes(shape = studyname, colour=Vaccine))+
  geom_errorbar(mapping = aes(ymin=SevereMin, ymax=SevereMax, colour=Vaccine, alpha=(SevereMax-SevereMin)), width=0)+
  geom_errorbar(mapping = aes(xmin=SymptomaticMin, xmax=SymptomaticMax, colour=Vaccine, alpha = (SymptomaticMax-SymptomaticMin)), width=0)+
  #annotate(geom='line',x=seq(0,100,1),y=seq(0,100,1))+
  #stat_poly_eq(formula = y ~ x,aes(label = ..rr.label..), parse = TRUE)+
  
  #geom_line(data = symptsevereeff_summary, aes(x=SymptEffRounded, y=SevereEffMid), colour='darkred')+
  #geom_ribbon(data = symptsevereeff_summary, aes(x=SymptEffRounded, y=SevereEffMid,ymin=SevereEffLower, ymax=SevereEffUpper), fill='darkred', alpha=.15)+
  geom_line(data = symptsevere_summary_curve, aes(x=symptomaticEfficacy, y=severeEfficacy), colour='darkred')+
  geom_ribbon(data = symptsevere_summary_curve, aes(x=symptomaticEfficacy, y=severeEfficacy,ymin=severeEfficacyMin, ymax=severeEfficacyMax), fill='darkred', alpha=.15)+
  
  scale_alpha_continuous(range = c(.05,1), rescaler=alpharescale, name='Width of CI', guide='none')+
  scale_colour_manual(values = vaccine_colours[names(vaccine_colours) %in% unique(rw_unboosted_data$Vaccine)], name = 'Vaccine')+#, guide='none')+
  scale_fill_manual(values = variant_fills, name='Variant', guide='none')+
  scale_shape_manual(values = study_shapes_eff[names(study_shapes_eff) %in% unique(rw_unboosted_data$studyname)], name = 'Study')+#, guide='none')+
  
  stat_cor(method='spearman',cor.coef.name='rho',label.sep="\n", label.x=25,label.y=5)+
  
  decaypaper_theme+
  theme(plot.title = element_text(face = 'bold'))+
  labs(x='Symptomatic Protection', y='Severe Protection', title='C Severe vs Symptomatic Efficacy')


ggsave(paste0(figuresdir,'symptvssevereplot.pdf'),symptvssevereplot,width=6, height=4)


rw_vs_model_predict_ABC_plot_formatted = plot_grid(rw_vs_model_predict_sympt_sev_plot+ theme(legend.position="none"),
                                               symptvssevereplot+ theme(legend.position="none"),
                                               nonstudylegend,
                                               studylegend,
                                               #get_legend(predictedVsObservedPlot),
                                               rel_widths=c(4,2,1.2,1.2),
                                               ncol = 4)
ggsave(paste0(figuresdir,'rw_vs_model_predict_ABC_plot_formatted.pdf'),rw_vs_model_predict_ABC_plot_formatted,width=14, height=4)
ggsave(paste0(manuscriptfiguresdir,'Fig3_R.pdf'),rw_vs_model_predict_ABC_plot_formatted,width=14, height=4)

}
######## 1:1 plots ########
oneToOne = T
if (oneToOne){
# This is the 1:1 x-y plot
useFormula <-y ~ x
predictedVsObservedLine = ggplot(matchedRWandPredictData, aes(x=modelEff, y=dataEff))+
  geom_point(mapping = aes(fill=VariantGroup),shape=22,size=3.5,colour='white',alpha=.5)+
  geom_point(mapping = aes(shape = studyname, colour=Vaccine))+
  geom_ribbon(data=efficacy_summary_curve,aes(x=eff_mean,ymin=eff_lower, ymax=eff_upper, y=eff_mean),alpha=.15,fill='darkred', colour=NA, inherit.aes=F)+
  geom_errorbar(mapping=aes(xmin = modelMin, xmax = modelMax, alpha = (modelMax-modelMin), colour=Vaccine),width=0)+
  geom_errorbar(mapping=aes(ymin = dataMin, ymax = dataMax, alpha = (dataMax-dataMin), colour=Vaccine),width=0)+
  geom_line(data=efficacy_summary_curve,mapping=aes(x=eff_mean,y=eff_mean),colour='darkred')+
  stat_cor(method='pearson',cor.coef.name='R',label.sep="\n")+
  
  scale_colour_manual(values = vaccine_colours[names(vaccine_colours) %in% unique(rw_unboosted_data$Vaccine)], name = 'Vaccine')+
  scale_fill_manual(values = variant_fills, name='Variant')+
  scale_shape_manual(values = study_shapes_eff[names(study_shapes_eff) %in% unique(rw_unboosted_data$studyname)], name = 'Study')+
  scale_alpha_continuous(range = c(.05,1), rescaler=alpharescale, name='Width of CI')+
    
  decaypaper_theme+
  facet_wrap(~EndpointCat)

ggsave(paste0(figuresdir,'predictedVsObservedLine.pdf'),predictedVsObservedLine,width=8, height=4)
predictedVsObservedLine_formatted = plot_grid(predictedVsObservedLine+ theme(legend.position="none"),
                                           nonstudylegend,
                                           studylegend,
                                           rel_widths=c(5,1,1.2),
                                           ncol=3)
ggsave(paste0(manuscriptfiguresdir,'FigS2_R.pdf'),predictedVsObservedLine_formatted,width=10, height=4)
}
######## Text stuff ######
a=filter(time_course_summary, Vaccine=='Pfizer',VariantGroup=='Delta',variable == 'log10neutR',Day %in% c(45,75,105,135))
b=cbind(a$Day,a[,c('lower','mid','upper','modelMean')])
round(10^b[c(3,4,1,2),c(1,2,5,4,3)],2)

######## Supplementary Sigmoid curves #######
suppSigmoidFigs=T
if(suppSigmoidFigs){
suppVariantColours = c('hotpink','deepskyblue2','seagreen3')
suppStudyTypeColours = c('orange','royalblue2','gray57')

xlabpos = -.5
ylabpos=50
ylabpos2=30
rw_unboosted_symptsev_data_with_predict$n = 0
rw_vs_model_predict_sympt_sev_plot_variant_all = ggplot(filter(rw_unboosted_symptsev_data_with_predict, EndpointCat %in%c('Symptomatic','Severe')), aes(x=log10neutR, y=Protection))+
  #geom_line(data=efficacy_summary_curve,mapping=aes(x=log10neutR,y=eff_median), colour='darkred', size=1,alpha=1)+
  #geom_ribbon(data=efficacy_summary_curve,mapping=aes(x=log10neutR,ymin=eff_lower, ymax=eff_upper,y=eff_median), fill='darkred',alpha=.15,colour=NA)+
  geom_line(data=filter(efficacy_summary_curve_for_logistic_plots, EndpointCat %in%c('Symptomatic','Severe')), colour='darkred', size=1,alpha=1)+
  geom_ribbon(data=filter(efficacy_summary_curve_for_logistic_plots, EndpointCat %in%c('Symptomatic','Severe')),mapping = aes(ymin = lower, ymax = upper), fill='darkred',alpha=.15,colour=NA)+
  
  geom_point(mapping = aes(shape = studyname, colour=VariantGroup,alpha=(CI_width)))+
  geom_errorbar(mapping=aes(alpha = CI_width, colour=VariantGroup,ymin = lower, ymax = upper))+
  geom_errorbar(mapping=aes(xmin=log10neutR-log10neutR_sd*1.96, xmax=log10neutR+log10neutR_sd*1.96,colour=VariantGroup,alpha=50))+
  
  scale_alpha_continuous(range = c(.05,1), rescaler=alpharescale, name='Width of CI')+
  scale_colour_manual(values = suppVariantColours, name = 'Variant')+
  scale_shape_manual(values = study_shapes_eff[names(study_shapes_eff) %in% unique(rw_unboosted_data$studyname)], name = 'Study', guide='none')+
  scale_x_continuous(breaks=log10(4^seq(-3,3,1)), labels = MASS::fractions(4^seq(-3,3,1)))+
  scale_y_continuous(breaks=seq(0,100,25), limits = c(0,100), labels = paste0(seq(0,100,25),'%'))+
  
  stat_cor(method='spearman',cor.coef.name='rho',label.sep="\n", label.x=xlabpos,label.y=ylabpos)+
  geom_text(aes(label=paste0('n = ',..count..),x=n),y=ylabpos2,stat="count", inherit.aes=F)+
  decaypaper_theme+
  facet_grid(~letterLabel)+
  labs(x='Neutralising Antibody Levels (fold of Convalescent)', y='Vaccine Effectiveness')

ggsave(paste0(figuresdir,'rw_vs_model_predict_sympt_sev_plot_variant_all.pdf'),rw_vs_model_predict_sympt_sev_plot_variant_all,width=8, height=4)

rw_vs_model_predict_sympt_sev_plot_variant_split = rw_vs_model_predict_sympt_sev_plot_variant_all+
  facet_wrap(VariantGroup~letterLabel, ncol=2)
ggsave(paste0(figuresdir,'rw_vs_model_predict_sympt_sev_plot_variant_split.pdf'),rw_vs_model_predict_sympt_sev_plot_variant_split,width=8, height=12)

rw_vs_model_predict_sympt_sev_plot_studytype_all = ggplot(filter(rw_unboosted_symptsev_data_with_predict, EndpointCat %in%c('Symptomatic','Severe')), aes(x=log10neutR, y=Protection))+
  #geom_line(data=efficacy_summary_curve,mapping=aes(x=log10neutR,y=eff_median), colour='darkred', size=1,alpha=1)+
  #geom_ribbon(data=efficacy_summary_curve,mapping=aes(x=log10neutR,ymin=eff_lower, ymax=eff_upper,y=eff_median), fill='darkred',alpha=.15,colour=NA)+
  geom_line(data=filter(efficacy_summary_curve_for_logistic_plots, EndpointCat %in%c('Symptomatic','Severe')), colour='darkred', size=1,alpha=1)+
  geom_ribbon(data=filter(efficacy_summary_curve_for_logistic_plots, EndpointCat %in%c('Symptomatic','Severe')),mapping = aes(ymin = lower, ymax = upper), fill='darkred',alpha=.15,colour=NA)+
  
  geom_point(mapping = aes(shape = studyname, colour=Studytype,alpha=(CI_width)))+
  geom_errorbar(mapping=aes(alpha = CI_width, colour=Studytype,ymin = lower, ymax = upper))+
  geom_errorbar(mapping=aes(xmin=log10neutR-log10neutR_sd*1.96, xmax=log10neutR+log10neutR_sd*1.96,colour=Studytype,alpha=50))+
  
  scale_alpha_continuous(range = c(.05,1), rescaler=alpharescale, name='Width of CI')+
  scale_colour_manual(values = suppStudyTypeColours, name = 'StudyType')+
  scale_shape_manual(values = study_shapes_eff[names(study_shapes_eff) %in% unique(rw_unboosted_data$studyname)], name = 'Study', guide='none')+
  scale_x_continuous(breaks=log10(4^seq(-3,3,1)), labels = MASS::fractions(4^seq(-3,3,1)))+
  scale_y_continuous(breaks=seq(0,100,25), limits = c(0,100), labels = paste0(seq(0,100,25),'%'))+
  
  stat_cor(method='spearman',cor.coef.name='rho',label.sep="\n", label.x=xlabpos,label.y=ylabpos)+
  geom_text(aes(label=paste0('n = ',..count..),x=n),y=ylabpos2,stat="count", inherit.aes=F)+
  decaypaper_theme+
  facet_grid(~letterLabel)+
  labs(x='Neutralising Antibody Levels (fold of Convalescent)', y='Vaccine Effectiveness')

ggsave(paste0(figuresdir,'rw_vs_model_predict_sympt_sev_plot_studytype_all.pdf'),rw_vs_model_predict_sympt_sev_plot_studytype_all,width=8, height=4)

rw_vs_model_predict_sympt_sev_plot_studytype_split = rw_vs_model_predict_sympt_sev_plot_studytype_all+
  facet_wrap(Studytype~letterLabel, ncol=2)
  #stat_cor(method='spearman',cor.coef.name='rho',label.sep="\n")
ggsave(paste0(figuresdir,'rw_vs_model_predict_sympt_sev_plot_studytype_split.pdf'),rw_vs_model_predict_sympt_sev_plot_studytype_split,width=8, height=12)

theme_blank_axes = theme(legend.position="none")+theme(axis.title = element_blank(), 
                         axis.line.x = element_blank(),
                         axis.ticks.x = element_blank(),
                         axis.text.x = element_blank())
#Combine them all together
supp_sigmoids_combined =  plot_grid(
                          plot_grid(rw_vs_model_predict_sympt_sev_plot_variant_all+theme_blank_axes,
                                   rw_vs_model_predict_sympt_sev_plot_studytype_all+ theme(legend.position="none")+theme_blank_axes, 
                                   ncol=2),
                          get_legend(rw_vs_model_predict_sympt_sev_plot_variant_all),
                          plot_grid(rw_vs_model_predict_sympt_sev_plot_variant_split+ theme(legend.position="none"),
                                   rw_vs_model_predict_sympt_sev_plot_studytype_split+ theme(legend.position="none"),
                                   ncol=2),
                          get_legend(rw_vs_model_predict_sympt_sev_plot_studytype_all),
                          nrow=2,ncol=2,rel_heights=c(1,3), rel_widths=c(5,1))
ggsave(paste0(figuresdir,'supp_sigmoids_combined.pdf'),supp_sigmoids_combined, height=12, width=12)
ggsave(paste0(manuscriptfiguresdir,'FigS1_R.pdf'),supp_sigmoids_combined,width=12, height=12)

# nonstudylegend,
#                                                    studylegend,
#                                                    #get_legend(predictedVsObservedPlot),
#                                                    rel_widths=c(4,2,1.2,1.2),
#                                                    ncol = 4)
# ggsave(paste0(figuresdir,'rw_vs_model_predict_ABC_plot_formatted.pdf'),rw_vs_model_predict_ABC_plot_formatted,width=14, height=4)
# ggsave(paste0(manuscriptfiguresdir,'Fig3_R.pdf'),rw_vs_model_predict_ABC_plot_formatted,width=14, height=4)

}

########
# DONT NEED PAST HERE
makeExtraPlots = F
if(makeExtraPlots){

  # Progression over time plots
  progression_over_time_data_and_theory_unboosted = ggplot(simulated_progression_plot_data_use,aes(x=Day/30,y=100*mid,color=Vaccine))+ 
    geom_ribbon(aes(ymin=100*lower,ymax=100*upper, fill=Vaccine),alpha=0.15,col=NA)+
    geom_line()+
    facet_wrap(~letterLabel, nrow=3)+
    
    geom_point(data=unboosted_progression_data,aes(x=TimeMidpt, y = 100-progression_ratio, colour = Vaccine,fill=boosted, alpha=(prog_protect_CI_width), shape=studyname),size=3)+
    geom_errorbar(data=unboosted_progression_data,aes(x=TimeMidpt, ymin = prog_protect_low,ymax=prog_protect_high, alpha=(prog_protect_CI_width), colour = Vaccine), inherit.aes=F)+
    
    scale_alpha_continuous(range = c(.05,1), rescaler=alpharescale, name='Width of CI')+
    scale_colour_manual(values = vaccine_colours[names(vaccine_colours) %in% unique(simulated_progression_plot_data_use$Vaccine)],labels =  vaccine_labels[names(vaccine_labels) %in% unique(simulated_progression_plot_data_use$Vaccine)], name = 'Vaccine')+
    scale_fill_manual(values = vaccine_colours[names(vaccine_colours) %in% c(as.character(unique(simulated_progression_plot_data_use$Vaccine)),'Boosted')], name='Vaccine', guide='none')+
    scale_shape_manual(values = study_shapes_eff[names(study_shapes_eff) %in% unique(unboosted_progression_data$studyname)], name = 'Study')+
    
    decaypaper_theme+
    
    xlim(0,8)+
    scale_y_continuous(breaks=seq(0, 100, 20), limits=c(0, 100))+
    labs(x='Time since vaccination (In Months)', y='Conditional Severe Protection')
  
  progression_over_time_data_and_theory_unboosted_formatted = plot_grid(progression_over_time_data_and_theory_unboosted+ theme(legend.position="none"),
                                                                        nonvariantlegend,
                                                                        rel_widths=c(3,1.1))
  ggsave(paste0(figuresdir,'progression_over_time_data_and_theory_unboosted_shaped.pdf'),progression_over_time_data_and_theory_unboosted_formatted,width=5.5,height=4)
  
  # This is the 1:1 x-y plot
  predictedVsObserved_progressionplot_Nyberg = ggplot(new_progression_plot_data_matched, aes(y=100*mid, x=Efficacy))+
    geom_point(mapping = aes(fill=VariantGroup),shape=22,size=4,colour='white',alpha=.5)+
    geom_point(mapping = aes(alpha= CImin),shape = 21, size=2,fill=progressionColour, colour=progressionColour)+
    #geom_point(shape=21)+
    geom_errorbar(mapping = aes(ymin=100*lower, ymax=100*upper, alpha=CIsim), width=0, colour=progressionColour)+
    geom_errorbar(mapping = aes(xmin=EfficacyMin, xmax=EfficacyMax, alpha = CIdata), width=0, colour=progressionColour)+
    annotate(geom='line',x=seq(0,100,1),y=seq(0,100,1))+
    #stat_poly_eq(formula = y ~ x,aes(label = ..rr.label..), parse = TRUE)+
    
    scale_alpha_continuous(range = c(.05,1), rescaler=alpharescale, name='Width of CI')+#, guide='none')+
    #scale_colour_manual(values = vaccine_colours[names(vaccine_colours) %in% unique(new_progression_plot_data_matched$Vaccine)], name = 'Vaccine')+#, guide='none')+
    scale_fill_manual(values = variant_fills[names(variant_fills)%in% unique(new_progression_plot_data_matched$VariantGroup)], name='Variant')+#, guide='none')+
    #scale_shape_manual(values = study_shapes_eff[names(study_shapes_eff) %in% unique(new_progression_plot_data_matched$studyname)], name = 'Study')+#, guide='none')+
    
    decaypaper_theme+
    
    labs(x='Observed Conditional Protection', y='Predicted Conditional Protection', title='')
  
  ggsave(paste0(figuresdir,'predictedVsObserved_progressionplot_Nyberg.pdf'),predictedVsObserved_progressionplot_Nyberg,width=5, height=4)
  
  #progression_over_time_data_and_new = plot_grid(plot_grid(progression_over_time_data_and_theory_unboosted+theme(legend.position="none"),
  #                                                         predictedVsObservedPlot3,
  #                                                         rel_heights =c(1.6,1), ncol=1, rel_widths = c(1,1)),
  #                                               variantlegend,
  #                                               rel_widths=c(3,1.7),ncol = 2)
  #ggsave(paste0(figuresdir,'progression_over_time_data_and_new.pdf'),progression_over_time_data_and_new,width=5, height=7)
  
  
pfizer_progressions_2to4mo_simulated = simulated_full_data_melted %>% ungroup %>%
  mutate(group='simulated', studyname='simulation', TimeMidpt = as.numeric(str_sub(variable,stri_locate_last(variable,fixed='_')[,2]+1,100)), prog_protect=100*progression_protection) %>% 
  filter(Vaccine=='Pfizer', TimeMidpt>60, TimeMidpt<=120, VariantGroup %in% c('Delta','Omicron')) %>%
  select(c(studyname, group, Vaccine,Booster, boosted, VariantGroup, regime, TimeMidpt,prog_protect))

pfizer_progressions_2to4mo_simulated_ranges = pfizer_progressions_2to4mo_simulated %>% group_by(Vaccine,VariantGroup) %>% 
  summarise(low=quantile(prog_protect, .025),mid=quantile(prog_protect, .5),high=quantile(prog_protect, .975))

pfizer_progressions_2to4mo_data_narrow = unboosted_progression_data %>% filter(Vaccine=='Pfizer', TimeMin>=2&TimeMax<=4) %>% #TimeMidpt>=2, TimeMidpt<=4
  mutate(group='data_narrow', prog_protect=100-progression_ratio, xval=runif(n())) %>%
  select(c(studyname, group, Vaccine,Booster, boosted, VariantGroup, regime, TimeMidpt,prog_protect, prog_protect_low, prog_protect_high,xval, prog_protect_CI_width)) %>% 
  merge(pfizer_progressions_2to4mo_simulated_ranges, by=c('Vaccine','VariantGroup'))

pfizer_progressions_2to4mo_simulated_ranges = rbind(pfizer_progressions_2to4mo_simulated_ranges,pfizer_progressions_2to4mo_simulated_ranges)
pfizer_progressions_2to4mo_simulated_ranges$xval=c(0,0,1,1)

pfizer_progression_plot=ggplot(pfizer_progressions_2to4mo_data_narrow, aes(x=xval,y=prog_protect,ymin=prog_protect_low, ymax=prog_protect_high,shape=studyname,colour=Vaccine, alpha = prog_protect_CI_width))+
  #geom_ribbon(aes(ymin=low, ymax=high),alpha=.15,fill='blue', colour=NA)+
  #geom_line(aes(y=mid),colour='blue')+
  geom_ribbon(data=pfizer_progressions_2to4mo_simulated_ranges,aes(x=xval,ymin=low, ymax=high),alpha=.15,fill='blue', colour=NA, inherit.aes=F)+
  geom_line(data=pfizer_progressions_2to4mo_simulated_ranges,aes(x=xval,y=mid),colour='blue', inherit.aes=F,size=1)+
  
  geom_point()+
  geom_errorbar()+
  
  scale_colour_manual(values = vaccine_colours[names(vaccine_colours) %in% unique(simulated_progression_plot_data$Vaccine)], name = 'Vaccine')+
  scale_shape_manual(values = study_shapes_eff[names(study_shapes_eff) %in% unique(unboosted_progression_data$studyname)], name = 'Study')+
  scale_alpha_continuous(range = c(.05,1), rescaler=alpharescale, name='Width of CI')+
  
  scale_x_continuous(breaks=NULL)+
  labs(x='',y='Protection from Progression')+
  facet_wrap(~VariantGroup)+
  
  decaypaper_theme
ggsave(paste0(figuresdir,'pfizer_progression_plot.pdf'),pfizer_progression_plot,width=7,height=5)


modeldata = filter(all_plot_data, Variant =='WT')





# These are the scatter plot comparisons - not used anymore

scattercomparisons_plotdata$letterLabel = scattercomparisons_plotdata$EndpointCat
letterLabels = levels(scattercomparisons_plotdata$letterLabel)[levels(scattercomparisons_plotdata$letterLabel)%in% unique(scattercomparisons_plotdata$letterLabel)]
scattercomparisons_plotdata$letterLabel = factor(scattercomparisons_plotdata$letterLabel,levels=letterLabels,labels=paste(strsplit(intToUtf8(64+c(1:length(letterLabels))),"")[[1]],letterLabels,sep=spacesep))

# Plot with Sympt and Severe and Letters by Vaccine
symptseverevaccineplot = ggplot(filter(scattercomparisons_plotdata, EndpointCat %in% c('Symptomatic','Severe')), aes(x=Vaccine, y=Efficacy))+
  geom_boxplot(colour='grey')+
  geom_point(mapping = aes(fill=VariantGroup),shape=22,size=3.5,colour='white',alpha=.5)+
  geom_point(mapping = aes(shape = studyname, colour=Vaccine))+
  
  stat_compare_means(comparisons=vaccine_comparisons,size=3, method = testmethod)+
  scale_colour_manual(values = vaccine_colours[names(vaccine_colours) %in% unique(scattercomparisons_plotdata$Vaccine)], name = 'Vaccine')+
  scale_fill_manual(values = variant_fills, name='Variant')+
  scale_shape_manual(values = study_shapes_eff[names(study_shapes_eff) %in% unique(scattercomparisons_plotdata$studyname)], name = 'Study')+
  
  decaypaper_theme+
  facet_wrap(~letterLabel, nrow=2)

ggsave(paste0(figuresdir,'symptseverevaccineplot.pdf'),symptseverevaccineplot,w=6,h=6)

severevaccineplot = symptseverevaccineplot
severevaccineplot$data = filter(scattercomparisons_plotdata, EndpointCat=='Severe')
ggsave(paste0(figuresdir,'symptvaccineplot.pdf'),symptvaccineplot,w=6,h=5)

symptvaccineplot = symptseverevaccineplot
symptvaccineplot$data = filter(scattercomparisons_plotdata, EndpointCat=='Symptomatic')
ggsave(paste0(figuresdir,'symptvaccineplot.pdf'),symptvaccineplot,w=6,h=5)


# Plot with Sympt and Severe and Letters by Variant
scattercomparisons_plotdata$letterLabel = scattercomparisons_plotdata$EndpointCat
letterLabels = levels(scattercomparisons_plotdata$letterLabel)[levels(scattercomparisons_plotdata$letterLabel)%in% unique(scattercomparisons_plotdata$letterLabel)]
scattercomparisons_plotdata$letterLabel = factor(scattercomparisons_plotdata$letterLabel,levels=letterLabels,labels=paste(strsplit(intToUtf8(66+c(1:length(letterLabels))),"")[[1]],letterLabels,sep=spacesep))

symptseverevariantplot = ggplot(filter(scattercomparisons_plotdata, EndpointCat%in% c('Symptomatic','Severe')), aes(x=VariantGroup, y=Efficacy))+
  geom_boxplot(colour='grey')+
  geom_point(mapping = aes(fill=VariantGroup),shape=22,size=3.5,colour='white',alpha=.5)+
  geom_point(mapping = aes(shape = studyname, colour=Vaccine))+
  
  stat_compare_means(comparisons=variant_comparisons,size=3, method=testmethod)+
  scale_colour_manual(values = vaccine_colours[names(vaccine_colours) %in% unique(scattercomparisons_plotdata$Vaccine)], name = 'Vaccine')+
  scale_fill_manual(values = variant_fills, name='Variant')+
  scale_shape_manual(values = study_shapes_eff[names(study_shapes_eff) %in% unique(scattercomparisons_plotdata$studyname)], name = 'Study')+
  
  decaypaper_theme+
  facet_wrap(~letterLabel, nrow=2)

ggsave(paste0(figuresdir,'symptseverevariantplot.pdf'),symptseverevariantplot,w=6,h=6)

severevariantplot = symptseverevariantplot
severevariantplot$data = filter(scattercomparisons_plotdata, EndpointCat=='Severe')
ggsave(paste0(figuresdir,'severevariantplot.pdf'),symptvaccineplot,w=6,h=5)

symptvariantplot = symptseverevariantplot
symptvariantplot$data = filter(scattercomparisons_plotdata, EndpointCat=='Symptomatic')
ggsave(paste0(figuresdir,'symptvariantplot.pdf'),symptvaccineplot,w=6,h=5)

symptseverescatterplots = plot_grid(symptseverevaccineplot + theme(legend.position="none"),
                                    symptseverevariantplot + theme(legend.position="none"),
                                    variantlegend,
                                    ncol=3,
                                    rel_widths = c(2,2, 1.2))
ggsave(paste0(figuresdir,'symptseverescatterplots.pdf'),symptseverescatterplots,w=8,h=6)

# Stephen's suggestion - not used
predictedThreeCurvesNoData = ggplot(predictedvsobserved_data, aes(x=log10neutR, y=Protection, colour=letterLabel, fill=letterLabel))+
  
  geom_line(data=efficacy_summary_curve_for_logistic_plots, size=1,alpha=1)+
  geom_ribbon(data=efficacy_summary_curve_for_logistic_plots,mapping = aes(ymin = lower, ymax = upper),alpha=.15,colour=NA)+
  
  scale_colour_manual(values = c('blue','red','orange'), labels = c('Symptomatic','Severe (Acquision)','Severe (Progression)'),name = 'Type of Protection')+
  scale_fill_manual(values = c('blue','red','orange'), labels = c('Symptomatic','Severe (Acquision)','Severe (Progression)'), name = 'Type of Protection')+
  scale_x_continuous(breaks=log10(4^seq(-3,3,1)), labels = MASS::fractions(4^seq(-3,3,1)))+
  
  decaypaper_theme+
  
  labs(x='Neutralising Antibody Levels (fold of Convalescent)', y='Protection')

ggsave(paste0(figuresdir,'predictedThreeCurvesNoData.pdf'),predictedThreeCurvesNoData,width=6, height=4)

# Nyberg progression over time
progression_over_time_newdata = ggplot(simulated_progression_plot_data,aes(x=Day/30,y=100*mid,color=Vaccine))+ 
  geom_ribbon(aes(ymin=100*lower,ymax=100*upper, fill=Vaccine),alpha=0.15,col=NA)+
  geom_line()+
  facet_wrap(~letterLabel, nrow=2)+
  
  geom_point(data=new_progression_data,aes(x=TimeMidpt, y = Efficacy, colour = Vaccine,fill=boosted, alpha=(EfficacyMax-EfficacyMin), shape=Paper),size=3)+
  geom_errorbar(data=new_progression_data,aes(x=TimeMidpt, ymin = EfficacyMin,ymax=EfficacyMax, alpha=(EfficacyMax-EfficacyMin), colour = Vaccine), inherit.aes=F)+
  
  scale_alpha_continuous(range = c(.05,1), rescaler=alpharescale, name='Width of CI')+
  scale_colour_manual(values = vaccine_colours[names(vaccine_colours) %in% unique(simulated_progression_plot_data$Vaccine)], name = 'Vaccine')+
  scale_fill_manual(values = vaccine_colours[names(vaccine_colours) %in% c(as.character(unique(simulated_progression_plot_data$Vaccine)),'Boosted')], name='Vaccine', guide='none')+
  #scale_shape_manual(values = study_shapes_eff[names(study_shapes_eff) %in% unique(unboosted_progression_data$studyname)], name = 'Study')+
  
  decaypaper_theme+
  
  xlim(0,8)+
  scale_y_continuous(breaks=seq(0, 100, 20), limits=c(0, 100))+
  labs(x='Time since vaccination (In Months)', y='Conditional Protection against Severe Disease', title = 'Conditional Protection')

progression_over_time_newdata_formatted = plot_grid(progression_over_time_newdata+ theme(legend.position="none"),
                                                    nonvariantlegend,
                                                    rel_widths=c(4,1.2))
ggsave(paste0(figuresdir,'progression_over_time_newdata_formatted.pdf'),progression_over_time_newdata_formatted,width=8,height=5)

# Predicted vs Obseved conditional protection - not very good plot
predictedVsObserved_progressionplot = ggplot(progression_plot_data_matched, aes(y=sim_mid, x=data_mid))+
  geom_point(mapping = aes(fill=VariantGroup),shape=22,size=3.5,colour='white',alpha=.5)+
  geom_point(mapping = aes(shape = studyname, colour=Vaccine))+
  geom_errorbar(mapping = aes(ymin=sim_lower, ymax=sim_upper, colour=Vaccine, alpha=(sim_upper-sim_lower)), width=0)+
  geom_errorbar(mapping = aes(xmin=data_lower, xmax=data_upper, colour=Vaccine, alpha = (data_upper-data_lower)), width=0)+
  annotate(geom='line',x=seq(0,100,1),y=seq(0,100,1))+
  #stat_poly_eq(formula = y ~ x,aes(label = ..rr.label..), parse = TRUE)+
  
  scale_alpha_continuous(range = c(.05,1), rescaler=alpharescale, name='Width of CI', guide='none')+
  scale_colour_manual(values = vaccine_colours[names(vaccine_colours) %in% unique(progression_plot_data_matched$Vaccine)], name = 'Vaccine')+#, guide='none')+
  scale_fill_manual(values = variant_fills, name='Variant', guide='none')+
  scale_shape_manual(values = study_shapes_eff[names(study_shapes_eff) %in% unique(progression_plot_data_matched$studyname)], name = 'Study')+#, guide='none')+
  
  decaypaper_theme+
  
  labs(x='Observed Conditional Protection', y='Predicted Conditional Protection', title='G')
ggsave(paste0(figuresdir,'predictedVsObserved_progressionplot.pdf'),predictedVsObserved_progressionplot,width=6, height=4)


# single plot progression
predictedVsObservedProgression = ggplot(unboosted_progression_data, aes(x=log10neutR, y=prog_protect, ymin = prog_protect_low, ymax = prog_protect_high))+
  #geom_line(aes(x=log10neutR,y=predictedEfficacy), colour='darkred', size=2)+
  #annotate(geom='line',x=efficacy_summary_curve$log10neut,y=efficacy_summary_curve$eff_median, colour='darkred', size=1)+
  #annotate(geom='ribbon',x=efficacy_summary_curve$log10neut,ymin=efficacy_summary_curve$eff_lower,ymax=efficacy_summary_curve$eff_upper,alpha=0.15,fill='darkred', colour=NA)+
  #annotate(geom='line',x=efficacy_summary_curve$log10neut,y=efficacy_summary_curve$eff_severe_median, colour='darkred', size=1)+
  #annotate(geom='ribbon',x=efficacy_summary_curve$log10neut,ymin=efficacy_summary_curve$eff_severe_lower,ymax=efficacy_summary_curve$eff_severe_upper,alpha=0.15,fill='darkred', colour=NA)+
  
  geom_line(data=efficacy_summary_curve_prog, colour='darkred', size=1,alpha=1)+
  geom_ribbon(data=efficacy_summary_curve_prog, fill='darkred',alpha=.15,colour=NA)+
  
  geom_point(mapping = aes(colour=VariantGroup, alpha = prog_protect_CI_width),shape=15,size=3.5)+
  scale_colour_manual(values = variant_fills[names(variant_fills) %in% unique(unboosted_progression_data$VariantGroup)], name='Variant')+
  scale_alpha_continuous(range = c(.05,.6), rescaler=alpharescale, guide='none')+
  new_scale("alpha") +
  new_scale("color") +
  
  geom_point(mapping=aes(shape=studyname, alpha = prog_protect_CI_width, colour=Vaccine))+
  geom_errorbar(mapping=aes(alpha = prog_protect_CI_width, colour=Vaccine))+
  
  scale_colour_manual(values = vaccine_colours[names(vaccine_colours) %in% unique(unboosted_progression_data$Vaccine)], name = 'Vaccine')+
  scale_shape_manual(values = study_shapes_eff[names(study_shapes_eff) %in% unique(unboosted_progression_data$studyname)], name = 'Study', guide='none')+
  scale_x_continuous(breaks=log10(4^seq(-3,1,1)), labels = MASS::fractions(4^seq(-3,1,1)), limits = c(min(efficacy_summary_curve_prog$log10neutR),0.1))+
  scale_y_continuous(breaks=seq(0,100,25), limits = c(0,100), labels = paste0(seq(0,100,25),'%'))+
  
  scale_alpha_continuous(range = c(.05,1), rescaler=alpharescale, name='Width of CI')+
  
  decaypaper_theme+
  #facet_wrap(~EndpointCat)+
  labs(x='Neutralising Antibody Levels (fold of Convalescent)', y='Protection against Progression')

ggsave(paste0(figuresdir,'predictedVsObservedProgression.pdf'),predictedVsObservedProgression,width=5, height=4)

# With the Nyberg data

#change the letter labels for this plot only
new_progression_data$letterLabel = unique(predictedvsobserved_data$letterLabel)[3]
predictedVsObservedTotal = ggplot(predictedvsobserved_data, aes(x=log10neutR, y=Protection))+
  
  geom_line(data=efficacy_summary_curve_for_logistic_plots, colour='darkred', size=1,alpha=1)+
  geom_ribbon(data=efficacy_summary_curve_for_logistic_plots,mapping = aes(ymin = lower, ymax = upper), fill='darkred',alpha=.15,colour=NA)+
  
  geom_point(mapping = aes(colour=VariantGroup, alpha = CI_width),shape=15,size=3.5)+
  scale_colour_manual(values = variant_fills[names(variant_fills) %in% unique(predictedvsobserved_data$VariantGroup)], name='Variant', guide='none')+
  scale_alpha_continuous(range = c(.05,.6), rescaler=alpharescale, guide='none')+
  new_scale("alpha") +
  new_scale("color") +
  
  geom_point(mapping=aes(shape=studyname, alpha = CI_width, colour=Vaccine))+
  geom_errorbar(mapping=aes(alpha = CI_width, colour=Vaccine,ymin = lower, ymax = upper))+
  geom_errorbar(mapping=aes(xmin=log10neutR-log10neutR_sd*1.96, xmax=log10neutR+log10neutR_sd*1.96,colour=Vaccine,alpha=CI_width))+
  
  geom_point(data=new_progression_data, mapping=aes(alpha = CI_width), colour=progressionColour, size=3,shape=progressionShape)+
  geom_errorbar(data=new_progression_data, mapping=aes(alpha = CI_width,ymin = lower, ymax = upper), colour=progressionColour,shape=progressionShape)+
  geom_errorbar(data=new_progression_data, mapping=aes(xmin=log10neutR-log10neutR_sd*1.96, xmax=log10neutR+log10neutR_sd*1.96,alpha=CI_width),colour=progressionColour, width = 0)+
  
  scale_colour_manual(values = vaccine_colours[names(vaccine_colours) %in% unique(predictedvsobserved_data$Vaccine)], name = 'Vaccine', guide='none')+
  scale_shape_manual(values = study_shapes_eff[names(study_shapes_eff) %in% unique(predictedvsobserved_data$studyname)], name = 'Study')+
  scale_x_continuous(breaks=log10(4^seq(-3,3,1)), labels = MASS::fractions(4^seq(-3,3,1)))+
  scale_y_continuous(breaks=seq(0,100,25), limits = c(0,100), labels = paste0(seq(0,100,25),'%'))+
  
  scale_alpha_continuous(range = c(.05,1), rescaler=alpharescale, name='Width of CI', guide='none')+
  
  facet_wrap(~letterLabel)+
  decaypaper_theme+
  
  labs(x='Neutralising Antibody Levels (fold of Convalescent)', y='Protection')

predictedVsObservedTotal_formatted = plot_grid(predictedVsObservedTotal+ theme(legend.position="none"),
                                               nonstudylegend,
                                               get_legend(predictedVsObservedTotal),
                                               rel_widths=c(6,1,1.2),
                                               ncol = 3)
ggsave(paste0(figuresdir,'predictedVsObservedTotal.pdf'),predictedVsObservedTotal_formatted,width=12, height=4)

} # end if(makeExtraPlots) {}