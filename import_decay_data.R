import_efficacy_data<-function(){
  efficacy_data =data.frame(read_xlsx(path = str_c(files$data),sheet = effdatasheet, range = effdatarange))
  efficacy_data = efficacy_data[!is.na(efficacy_data[,1]),]
  efficacy_data$TimeMax[efficacy_data$TimeMax==100]=8
  efficacy_data = efficacy_data %>% 
    mutate(join=paste0(Paper,Age,Variant, Vaccine, Booster), VariantVaccine = paste0(Variant,Vaccine),notboosted = is.na(Booster),plotjoin=paste0(Paper,Age,Variant, Vaccine, Booster, Endpoint))
  # Add in the newer efficacies where they have been published
  efficacy_data$Efficacy[!is.na(efficacy_data$Efficacy_New)]=efficacy_data$Efficacy_New[!is.na(efficacy_data$Efficacy_New)]
  efficacy_data$EfficacyMin[!is.na(efficacy_data$Efficacy_Min_New)]=efficacy_data$Efficacy_Min_New[!is.na(efficacy_data$Efficacy_Min_New)]
  efficacy_data$EfficacyMax[!is.na(efficacy_data$Efficacy_Max_New)]=efficacy_data$Efficacy_Max_New[!is.na(efficacy_data$Efficacy_Max_New)]
  
  # Change teh types of the variables
  efficacy_data$notboosted[efficacy_data$Booster == 'None'] = TRUE
  efficacy_data$Booster[is.na(efficacy_data$Booster)]='None'
  efficacy_data$boosted=factor(efficacy_data$notboosted, labels=c('Boosted','Unboosted'))
  efficacy_data$boosted=factor(efficacy_data$boosted, levels=c('Unboosted','Boosted'))
  efficacy_data$UseInAnalysis = as.logical(efficacy_data$UseInAnalysis)
  
  # Set up the variant and variant groups as factors
  efficacy_data$Variant = factor(efficacy_data$Variant,levels=variant_list[variant_list %in% unique(efficacy_data$Variant)])
  efficacy_data$VariantGroup = factor(efficacy_data$VariantGroup,levels=variant_list[variant_list %in% unique(efficacy_data$VariantGroup)])
  
  
  # Set up the Endpoint Cats as factor
  efficacy_data$EndpointCat = factor(efficacy_data$EndpointCat,levels=endpoint_categories[endpoint_categories %in% unique(efficacy_data$EndpointCat)])
  
  ymin=0
  efficacy_data$EfficacyMin[efficacy_data$EfficacyMin < ymin] = ymin
  efficacy_data$AgeGp = factor(efficacy_data$AgeGp, levels = unique(efficacy_data$AgeGp))
  
  # Set up the vaccine as a factor
  efficacy_data$Vaccine[efficacy_data$Vaccine=='AZ']='AstraZeneca'
  efficacy_data$Vaccine = factor(efficacy_data$Vaccine, levels =vaccine_list[vaccine_list %in% unique(efficacy_data$Vaccine)])
  efficacy_data$studyname = paste0(efficacy_data$FirstAuthor,' et. al.')
  
  #efficacy_data$Variant1=efficacy_data$Variant
  #efficacy_data$Variant = efficacy_data$VariantGroup
  #efficacy_data$VariantVaccine1 = efficacy_data$VariantVaccine
  efficacy_data$fitgroup=paste0(efficacy_data$VariantGroup,efficacy_data$Vaccine, efficacy_data$EndpointCat)
  
  efficacy_data
}

import_neut_data<-function(){
  neut_data =data.frame(read_xlsx(path = str_c(files$data),sheet = neutdatasheet, range = neutdatarange))
  neut_data = neut_data[!is.na(neut_data[,1]),]
  neut_data[neut_data=='LOD'] = NA
  
  numericalInds = which(startsWith(colnames(neut_data),'Neut'))
  neut_data[numericalInds]  <- lapply(neut_data[numericalInds], as.numeric)
  neut_data = neut_data %>% mutate(join=paste0(FirstAuthor,AgeGp,CohortDescription,PriorStatus,Variant, Vaccine,Group), VariantVaccine = paste0(Variant,Vaccine))
  neut_data
}

