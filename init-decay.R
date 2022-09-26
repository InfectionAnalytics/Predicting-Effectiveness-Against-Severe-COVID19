library(readxl)
library(stringr)
library(stringi)
#library(MASS)
library(dplyr) 
library(ggplot2)
library(ggnewscale)
library(pdftools)
library(data.table)
library(reshape2)
library(lhs)
library(latex2exp)
library(ggh4x)
library(cowplot)
library(forcats)
library(ggpubr)
library(ggpmisc)
library(mvtnorm)

username = str_split(str_remove(getwd(),'/Users/'),'/')[[1]][1]
location=str_c('/Users/',username,'/OneDrive - UNSW/')
currentdirectory<-dirname(rstudioapi::getActiveDocumentContext()$path)
if (!currentdirectory==""){ 
  setwd(currentdirectory) 
  currentdirectory='.'
}
currentdirectory='./'
basedir = str_c(location,"Documents - Vaccine Abs and Efficacy/05 AB AND EFFICACY DECAY/DECAY-Code/")
basedir = './'
manuscriptdir = str_c(location,"Documents - Vaccine Abs and Efficacy/05 AB AND EFFICACY DECAY/Manuscript/")
manuscriptfiguresdir = paste0(manuscriptdir,'Figures/')
figuresdir = paste0(basedir,'Figures/')

# for rerun

figuresdir = paste0(basedir,'Figures_Latest/')
manuscriptfiguresdir=figuresdir

source('./WP4_modelImplementationFunctions.R')
source('./severe_paper_helper_functions.R')
basedir='./'
manuscriptdir = './'
figuresdir = './Figures/'
manuscriptfiguresdir = './Figures/'

#setwd(basedir)
plot_type='pdf'
effdatasheet= 'EfficacyDataCombined'
effdatarange= 'A1:X1000'

neutdatasheet= 'NeutDecayDataCombined'
neutdatarange= 'A1:Q1000'

version = 1
source(str_c(basedir,'import_decay_data.R'))
dir = list(base = './',
           #excel = './ExcelFiles/',
           #results = './Results/',
           data = './',
           plots = './Figures/')

#dir = list(base = './',
#           data = './DECAY-Papers/',
#           plots = './Plots/')

files = list(data = str_c(dir$data,'DECAY_Master_ForPaper_github.xlsx')
             )

# Colours for plotting
variant_colours = c('WT'='black', 'Alpha' = 'darkgreen', 'Beta' = 'darkorange', 
                    'Non-Delta'='grey','preDelta'='red', 'upToDelta'='grey','Delta' = 'blue',
                    'Mix'='grey','OmicronDelta'='grey','Omicron' = 'purple')
variant_list = names(variant_colours)

variant_fills=c('Non-Delta'='pink','Delta'='lightblue','Omicron'='lightgreen')
variants_to_use = names(variant_fills)

variantvaccine_colours = c('AlphaPfizer' = 'deepskyblue', 
                           'BetaPfizer' = 'lavender', 
                           'preDeltaPfizer'='dodgerblue',
                           'DeltaPfizer' = 'blue3',
                           'Non-DeltaPfizer'='slategray2',
                           'DeltaPfizer+Moderna'='turquoise4',
                           'preDeltaModerna'='seagreen2',
                           'DeltaModerna'='seagreen4',
                           'preDeltaJnJ'='orchid1',
                           'DeltaJnJ'='darkmagenta',
                           'AlphaAZ'='orange', 
                           'DeltaAZ'='darkorange3',
                           'WTConvalescent'='lightgoldenrod1',
                           'AlphaConvalescent'= 'darkgoldenrod',
                           'BetaConvalescent'='wheat1')
vaccine_colours = c('Moderna' = 'red',
                    'Pfizer' = 'blue3',
                    'mRNA'='purple',
                    'AstraZeneca' = 'orchid1',
                    'Astra' = 'orchid1', 
                    'AZ' = 'orchid1',
                    'JnJ'='seagreen', 
                    'Convalescent'='orange', 
                    'PfizerModerna'='purple', 
                    'None'='grey', 
                    'Boosted'='seagreen'
                    )
vaccine_labels = c('Moderna' = 'mRNA1273',
                   'Pfizer' = 'BNT162b2',
                   'mRNA'='mRNA',
                   'AstraZeneca' = 'ChAdOx-nCoV-1',
                   'Astra' = 'ChAdOx-nCoV-1', 
                   'AZ' = 'ChAdOx-nCoV-1',
                   'JnJ'='JnJ', 
                   'Convalescent'='Conv', 
                   'PfizerModerna'='mRNA', 
                   'None'='NA', 
                   'Boosted'='Dose 3'
)
vaccine_list = names(vaccine_colours)

endpoint_categories = c('Infection', 'Symptomatic','Severe','All','Progression')
# Shapes for plotting
study_shapes_eff = c('DECAY01'=19, 
                     'DECAY03'=15,
                     'DECAY04'=4, 
                     'DECAY05'=13, 
                     'DECAY08'=18, 
                     'DECAY11'=25 , 
                     'DECAY13'=11, 
                     'DECAY14'=10, 
                     'DECAY16'=17, 
                     'DECAY23'=12)

studytypes=c('Retrospective','Retrospective','TNCC','Retrospective','TNCC','RCT','RCT','Retrospective','TNCC','TNCC','TNCC','Retrospective','TNCC','TNCC','TNCC','TNCC')
names(studytypes)=c('Tartof','Goldberg','Chemaitelly','Keehner','Andrews','Sahly','Thomas','Rosenberg','Andrews','Ferdinands','Bruxvoort','Poukka','Tseng','Skowronski _BC','Skowronski _Quebec','Thompson')
study_shapes_neut = c('DECAY06'=19, 'DECAY07'=15,'DECAY09'=16, 'DECAY10'=17, 'DECAY12'=18, 'DECAY15'=25 ,'DECAY17'=17, 'DECAY24'=11, 'DECAY26'=12)

efficacy_data = import_efficacy_data()
#neut_data = import_neut_data()

study_shapes_eff = c(11,12,13,14,15,16,17,18,0,1,2,5,6,7,8,9,10,35,37,94,64,36,38)
#shapes = c(21:25,21:25)
FirstAuthors = unique(efficacy_data$FirstAuthor)
names(study_shapes_eff) = paste0(FirstAuthors,' et. al.')
study_shapes_eff = study_shapes_eff[c(1:length(FirstAuthors))]

source(str_c(basedir,'plot_decay_efficacy.R'))

source(str_c(basedir,'natmed_modelImplementationFunctions.R'))
load(str_c(basedir,'natmed_parameters.RData'))
load(str_c(basedir,'natmed_parameters_sympt_and_severe.RData'))

sig=SummaryTable_Efficacy_NeutRatio_SD_SEM$PooledSD[1]
logk = tail(FittedLogistic_RawEfficacy_MeanRept_SDPool$estimate,2)[1]
# Fitted IC50 values
C50=tail(FittedLogistic_RawEfficacy_MeanRept_SDPool$estimate,1)
C50_Severe=tail(FittedLogistic_RawEfficacy_MeanRept_SDPool_Severe$estimate,1)

log10NeutRatios = SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutRatio_Reported
log10NeutRatios_SEM = SummaryTable_Efficacy_NeutRatio_SD_SEM$SEM

# Rename ratios 
names(log10NeutRatios) = SummaryTable_Efficacy_NeutRatio_SD_SEM$Study
names(log10NeutRatios)[names(log10NeutRatios)=='AstraZeneca']='Astra'
names(log10NeutRatios_SEM) = SummaryTable_Efficacy_NeutRatio_SD_SEM$Study
names(log10NeutRatios_SEM)[names(log10NeutRatios_SEM)=='AstraZeneca']='Astra'

# Boost astra zeneca to account for longer dosing schedule
astra_dose_shift = 1.59
log10NeutRatios[names(log10NeutRatios)=='Astra']=log10NeutRatios[names(log10NeutRatios)=='Astra']+log10(astra_dose_shift)


estimate_sd_from_CIs<-function(m, low95, high95,CIlim = .95){
  sdest=NULL
  for (i in c(1:length(m))){
    # Efficacies
    vals = c(low95[i],m[i],high95[i])
    # Estimate the mean and SD of the log10 neut ratios
    sd_init = mean((diff(vals))/qnorm((1+CIlim)/2))

    # Now take this back and see if this is the best estimate
    sdest[i] = nlm(function(s,m){sum((m+s*c(qnorm((1-CIlim)/2),0,qnorm((1+CIlim)/2))-vals)^2)},p=sd_init,m=m)$estimate
  }
  sdest
}


# Read in boosted information
BoostingTable_Convalescent<-filter(read.csv("./BoostingTableSummary_20210726.csv",fileEncoding="UTF-8-BOM"),PreviousImmune=="Convalescent")
BoostingTable_Convalescent$ReferenceNormalisation<-SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutRatio_Reported[match(BoostingTable_Convalescent$Reference,SummaryTable_Efficacy_NeutRatio_SD_SEM$Study)]
# There are some boosts that are Pfizer/Moderna, so need to work out the average of Pfizer and Moderna to get the reference level
ModernaLevelPreBoost=SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutRatio_Reported[SummaryTable_Efficacy_NeutRatio_SD_SEM$Study=="Moderna"]
ModernaLevelPreBoost_Pf=SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutRatio_Reported[SummaryTable_Efficacy_NeutRatio_SD_SEM$Study=="Pfizer"]
PfizerModernalAverageReferenceNormalisation=0.5*(ModernaLevelPreBoost_Pf+ModernaLevelPreBoost)
BoostingTable_Convalescent$ReferenceNormalisation[BoostingTable_Convalescent$Reference=="Pfizer/Moderna"]=PfizerModernalAverageReferenceNormalisation

# Now calculated boosted increase
BoostingTable_Convalescent$EstimatedLevelAfter<-log10((10^BoostingTable_Convalescent$ReferenceNormalisation)*BoostingTable_Convalescent$FoldChange)
Boosting_Average=mean(BoostingTable_Convalescent$EstimatedLevelAfter)
Boosting_SE=sd(BoostingTable_Convalescent$EstimatedLevelAfter)/sqrt(length(BoostingTable_Convalescent$EstimatedLevelAfter))

# Add in boosted and mRNA
# boosted_mean=12
# nboosted=5
# boosted_low = 6.1
# boosted_high=28.7
#boosted_SEM = estimate_sd_from_CIs(log10())
# boosted_mean=11.55
# boosted_SEM=1.27
#boosted_SEM = estimate_sd_from_CIs(log10(boosted_mean),log10(boosted_mean-1.96*boosted_SEM),log10(boosted_mean+1.96*boosted_SEM))
log10NeutRatios=c(log10NeutRatios,'Boosted'=Boosting_Average,'mRNA'=mean(log10NeutRatios[1:2]))
log10NeutRatios_SEM = c(log10NeutRatios_SEM,'Boosted'=Boosting_SE,'mRNA'=mean(log10NeutRatios_SEM[1:2]))

trialdays = 60




