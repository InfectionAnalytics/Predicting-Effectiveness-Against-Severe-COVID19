######## Printing #########
# Over time curves
realworlddata = rw_and_sim_efficacy_data %>% filter(!is.na(Paper))
realworlddata$modelEff=NA
realworlddata$modelMin=NA
realworlddata$modelMax=NA
# This estimates only from the time course curves adn I think is not as accurate as the other stuff
for (r in c(1:nrow(realworlddata))){
  predictData = rw_and_sim_efficacy_data %>% filter(is.na(Paper),
                                      Vaccine == realworlddata[r,'Vaccine'], 
                                      VariantGroup==realworlddata[r,'VariantGroup'],
                                      EndpointCat==realworlddata[r,'EndpointCat'])
  realworlddata$modelEff[r] = spline(predictData$TimeMidpt,predictData$Efficacy,xout=realworlddata$TimeMidpt[r])$y
  realworlddata$modelMin[r] = spline(predictData$TimeMidpt,predictData$EfficacyMin,xout=realworlddata$TimeMidpt[r])$y
  realworlddata$modelMax[r] = spline(predictData$TimeMidpt,predictData$EfficacyMax,xout=realworlddata$TimeMidpt[r])$y
}
realworlddata = realworlddata %>%
  mutate(dataEff = Efficacy,dataMin = EfficacyMin, dataMax = EfficacyMax,
         CIsoverlap = (dataMin<modelMin&dataMax>modelMin)|(dataMax>modelMax&dataMin<modelMax)|(dataMin>modelMin&dataMax<modelMax)|(dataMax>modelMax&dataMin<modelMin), 
         dataCIsOverCurve = (dataMin<modelEff & dataMax>modelEff),
         modelCIsOverData = (modelMin<dataEff & modelMax>dataEff), 
         datamodeldiff=dataEff-modelEff, 
         CIwidth=dataMax-dataMin,
         dataOutsideCI = pmin(dataEff-modelMin, modelMax-dataEff))

matchedRWandPredictDataSevere =filter(matchedRWandPredictData,!is.na(modelEff),EndpointCat=='Severe',!is.na(matchedRWandPredictData$dataEff))
realworlddataSevere =filter(realworlddata,!is.na(modelEff),EndpointCat=='Severe',!is.na(realworlddata$dataEff))
nSevereTimeCourse = nrow(realworlddataSevere)
nSevereTimeCourse_new = nrow(matchedRWandPredictDataSevere)
nSevereTimeCourseMatched = sum(realworlddataSevere$modelCIsOverData[!is.na(realworlddataSevere$modelMin)&!is.na(realworlddataSevere$modelMax)])
nSevereTimeCourseMatched_new = sum(matchedRWandPredictDataSevere$modelCIsOverData[!is.na(matchedRWandPredictDataSevere$modelMin)&!is.na(matchedRWandPredictDataSevere$modelMax)])
nSevereTimeCourse10pc =sum(realworlddataSevere$datamodeldiff<10)
nSevereTimeCourse10pc_new =sum(matchedRWandPredictDataSevere$datamodeldiff<10)
nSevereTimeCourse5pc =sum(realworlddataSevere$datamodeldiff<5)
nSevereTimeCourse5pc_new =sum(matchedRWandPredictDataSevere$datamodeldiff<5)

matchedRWandPredictDataSympt =filter(matchedRWandPredictData,!is.na(modelEff),EndpointCat=='Symptomatic',!is.na(matchedRWandPredictData$dataEff))
realworlddataSympt =filter(realworlddata,!is.na(modelEff),EndpointCat=='Symptomatic',!is.na(realworlddata$dataEff))
nSymptTimeCourse = nrow(realworlddataSympt)
nSymptTimeCourseMatched = sum(realworlddataSympt$modelCIsOverData[!is.na(realworlddataSympt$modelMin)&!is.na(realworlddataSympt$modelMax)])
nSymptTimeCourse10pc =sum(realworlddataSympt$datamodeldiff<10)
nSymptTimeCourse5pc =sum(realworlddataSympt$datamodeldiff<5)

print(paste0('For the time course data ',nSevereTimeCourseMatched,' / ',nSevereTimeCourse,' (',round(100*nSevereTimeCourseMatched/nSevereTimeCourse,1),
             '%) of the time course estimates for severe COVID-19 lie within the confidence intervals of the model.'))
print(paste0(nSevereTimeCourse10pc,' / ',nSevereTimeCourse,' (',round(100*nSevereTimeCourse10pc/nSevereTimeCourse,1),'%) lie within 10% of the predicted value.'))
print(paste0(nSevereTimeCourse5pc,' / ',nSevereTimeCourse,' (',round(100*nSevereTimeCourse5pc/nSevereTimeCourse,1),'%) lie within 5% of the predicted value.'))

print(paste0('For Symptomatic it is ',nSymptTimeCourseMatched,' / ',nSymptTimeCourse,' (',round(100*nSymptTimeCourseMatched/nSymptTimeCourse,1),
             '%) of the time course estimates for symptomatic COVID-19 lie within the confidence intervals of the model.'))
print(paste0(nSymptTimeCourse10pc,' / ',nSymptTimeCourse,' (',round(100*nSymptTimeCourse10pc/nSymptTimeCourse,1),'%) lie within 10% of the predicted value.'))
print(paste0(nSymptTimeCourse5pc,' / ',nSymptTimeCourse,' (',round(100*nSymptTimeCourse5pc/nSymptTimeCourse,1),'%) lie within 5% of the predicted value.'))


#Sigmoid curves
nSevereMatched =sum(matchedRWandPredictData$modelCIsOverData[!is.na(matchedRWandPredictData$modelMin)&matchedRWandPredictData$EndpointCat=='Severe'&!is.na(matchedRWandPredictData$dataEff)])
nSevere10pc =sum(matchedRWandPredictData$datamodeldiff[!is.na(matchedRWandPredictData$modelMin)&matchedRWandPredictData$EndpointCat=='Severe'&!is.na(matchedRWandPredictData$dataEff)]<10)
nSevere5pc =sum(matchedRWandPredictData$datamodeldiff[!is.na(matchedRWandPredictData$modelMin)&matchedRWandPredictData$EndpointCat=='Severe'&!is.na(matchedRWandPredictData$dataEff)]<5)
nSevere=length(matchedRWandPredictData$modelCIsOverData[!is.na(matchedRWandPredictData$modelMin)&matchedRWandPredictData$EndpointCat=='Severe'&!is.na(matchedRWandPredictData$dataEff)])

nSymptMatched =sum(matchedRWandPredictData$modelCIsOverData[!is.na(matchedRWandPredictData$modelMin)&matchedRWandPredictData$EndpointCat=='Symptomatic'&!is.na(matchedRWandPredictData$dataEff)])
nSympt10pc =sum(matchedRWandPredictData$datamodeldiff[!is.na(matchedRWandPredictData$modelMin)&matchedRWandPredictData$EndpointCat=='Symptomatic'&!is.na(matchedRWandPredictData$dataEff)]<10)
nSympt5pc =sum(matchedRWandPredictData$datamodeldiff[!is.na(matchedRWandPredictData$modelMin)&matchedRWandPredictData$EndpointCat=='Symptomatic'&!is.na(matchedRWandPredictData$dataEff)]<5)
nSympt=length(matchedRWandPredictData$modelCIsOverData[!is.na(matchedRWandPredictData$modelMin)&matchedRWandPredictData$EndpointCat=='Symptomatic'&!is.na(matchedRWandPredictData$dataEff)])

print("")
print(paste0('In fact ',nSevereMatched,' / ',nSevere,' (',round(100*nSevereMatched/nSevere,1),'%) real- world estimates for vaccine efficacy against severe COVID-19 lie within the confidence intervals of the model.'))
print(paste0(nSevere10pc,' / ',nSevere,' (',round(100*nSevere10pc/nSevere,1),'%) lie within 10% of the predicted value.'))
print(paste0(nSevere5pc,' / ',nSevere,' (',round(100*nSevere5pc/nSevere,1),'%) lie within 5% of the predicted value.'))

print(paste0('For Symptomatic it is ',nSymptMatched,' / ',nSympt,' (',round(100*nSymptMatched/nSympt,1),'%) real- world estimates for vaccine efficacy against symptomatic COVID-19 lie within the confidence intervals of the model.'))
print(paste0(nSympt10pc,' / ',nSympt,' (',round(100*nSympt10pc/nSympt,1),'%) lie within 10% of the predicted value.'))
print(paste0(nSympt5pc,' / ',nSympt,' (',round(100*nSympt5pc/nSympt,1),'%) lie within 5% of the predicted value.'))

matchCIs_symptseverematch = rw_matched_symptsev_data %>%
  #select(c(Paper, FirstAuthor,Vaccine, Variant,VariantGroup, Age,AgeGp,TimeMidpt,TimeMin,TimeMax)) %>% #,SymptomaticMin=EfficacyMin.x, SymptomaticEff=Efficacy.x, SymptomaticMax=EfficacyMax.x,SevereMin=EfficacyMin.y, SevereEff=Efficacy.y, SevereMax=EfficacyMax.y)) %>%
  mutate(predictSevereEff = spline(symptsevere_summary_curve$symptomaticEfficacy,symptsevere_summary_curve$severeEfficacy,xout=SymptomaticEff)$y,
         predictSevereEffMin = spline(symptsevere_summary_curve$symptomaticEfficacy,symptsevere_summary_curve$severeEfficacyMin,xout=SymptomaticEff)$y,
         predictSevereEffMax = spline(symptsevere_summary_curve$symptomaticEfficacy,symptsevere_summary_curve$severeEfficacyMax,xout=SymptomaticEff)$y,
         dataMin = SevereMin, dataMax = SevereMax, dataEff = SevereEff, 
         modelMin = predictSevereEffMin, modelMax = predictSevereEffMax, modelEff = predictSevereEff,
         CIsoverlap = (dataMin<modelMin&dataMax>modelMin)|(dataMax>modelMax&dataMin<modelMax)|(dataMin>modelMin&dataMax<modelMax)|(dataMax>modelMax&dataMin<modelMin), 
         dataCIsOverCurve = (dataMin<modelEff & dataMax>modelEff),
         modelCIsOverData = (modelMin<dataEff & modelMax>dataEff), 
         datamodeldiff=dataEff-modelEff, CIwidth=dataMax-dataMin,
        dataOutsideCI = pmin(dataEff-modelMin, modelMax-dataEff))

matchingdata_symptseverematch = matchCIs_symptseverematch %>% 
  filter(!is.na(modelMin)&!is.na(dataMin)) %>%
  summarise(datamodeloverlap = sum(CIsoverlap)/n(), modelCIsOverData = sum(modelCIsOverData)/n(), within10pc=sum(abs(datamodeldiff)<10)/n() )

print('Sympt Severe Match')
print(paste0('Data and Model CI overlap: ',round(matchingdata_symptseverematch$datamodeloverlap*100,1),'%'))
print(paste0('Model CIs overlap Data: ',round(matchingdata_symptseverematch$modelCIsOverData*100,1),'%'))
print(paste0('Data within 10% of Model: ',round(matchingdata_symptseverematch$within10pc*100,1),'%'))

nSymptSevere = nrow(matchCIs_symptseverematch)
nSymptSevereMatched = sum(matchCIs_symptseverematch$modelCIsOverData)
nSymptSevere10pc =sum(matchCIs_symptseverematch$datamodeldiff<10)
nSymptSevere5pc =sum(matchCIs_symptseverematch$datamodeldiff<5)

print("")
print(paste0('In addition, the relationship between symptomatic and severe protection in the real-world data set is consistent with the non-linear relationship predicted by the correlates model (',
             nSymptSevereMatched,'/',nSymptSevere,' (',round(100*nSymptSevereMatched/nSymptSevere,0),'%) datapoints lie within the confidence intervals of the model). '))

#sympt_summary_curve=filter(efficacy_summary_curve, EndpointCat=='Symptomatic')
#severe_summary_curve=filter(efficacy_summary_curve, EndpointCat=='Severe')
# This uses widerCIs - it is the correct option, but progressiondata only includes delta and omicron need pre delta too
#matchCIs_symptseverematch2 = progression_data %>% 
#  filter(Booster=='None') %>% 
#  mutate(log10neutR = spline(sympt_summary_curve$eff_median,sympt_summary_curve$log10neutR,xout=SymptomaticEff)$y, modelEff = spline(severe_summary_curve$log10neutR,severe_summary_curve$eff_median,xout=log10neutR)$y,modelMin=spline(severe_summary_curve$log10neutR,severe_summary_curve$eff_lower,xout=log10neutR)$y,modelMax=spline(severe_summary_curve$log10neutR,severe_summary_curve$eff_upper,xout=log10neutR)$y  ) %>%
#  select(c(Paper, FirstAuthor, Vaccine, VariantGroup,SymptomaticEff,dataMin=SevereMin,log10neutR,dataEff=SevereEff,dataMax=SevereMax,modelEff, modelMin, modelMax))%>%
#  mutate(CIsoverlap = (dataMin<modelMin&dataMax>modelMin)|(dataMax>modelMax&dataMin<modelMax)|(dataMin>modelMin&dataMax<modelMax)|(dataMax>modelMax&dataMin<modelMin), dataCIsOverCurve = (dataMin<modelEff & dataMax>modelEff),modelCIsOverData = (modelMin<dataEff & modelMax>dataEff), datamodeldiff=dataEff-modelEff, CIwidth=dataMax-dataMin)


print('')
print(paste0('SIGMOID CURVES'))
print(paste0('Symptomatic ',nSymptMatched,' / ',nSympt,' (',round(100*nSymptMatched/nSympt,1),'%) data points lie wihting the CIs of the model, and ',
             nSympt10pc,' / ',nSympt,' (',round(100*nSympt10pc/nSympt,1),'%) are within 10%. ',
             '\rho=',round((cor(matchedRWandPredictDataSympt$dataEff, matchedRWandPredictDataSympt$modelEff, method='spearman')),2),
             ' p= ',cor.test(matchedRWandPredictDataSympt$dataEff, matchedRWandPredictDataSympt$modelEff, method='spearman')$p.value,' (Spearman)'))

print(paste0('Severe ',nSevereMatched,' / ',nSevere,' (',round(100*nSevereMatched/nSevere,1),'%) data points lie within the CIs of the model, and ',
             nSevere10pc,' / ',nSevere,' (',round(100*nSevere10pc/nSevere,1),'%) are within 10%. ',
             '\rho=',round((cor(matchedRWandPredictDataSevere$dataEff, matchedRWandPredictDataSevere$modelEff, method='spearman')),2),
             ' p= ',cor.test(matchedRWandPredictDataSevere$dataEff, matchedRWandPredictDataSevere$modelEff, method='spearman')$p.value,' (Spearman)'))

print(paste0('Sympt vs Severe ',nSymptSevereMatched,' / ',nSymptSevere,' (',round(100*nSymptSevereMatched/nSymptSevere,1),'%) data points lie within the CIs of the model, and ',
             nSympt10pc,' / ',nSymptSevere,' (',round(100*nSymptSevere10pc/nSymptSevere,1),'%) are within 10%. ',
             '\rho=',round((cor(matchCIs_symptseverematch$SymptomaticEff, matchCIs_symptseverematch$SevereEff, method='spearman')),2),
             ' p= ',cor.test(matchCIs_symptseverematch$SymptomaticEff, matchCIs_symptseverematch$SevereEff, method='spearman')$p.value,' (Spearman)'))


