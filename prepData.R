
#input: ct: ct data table
        #labnames: labs table, consists of labid and labname columns
        # gene: 'N' or 'E'
        #variant: 'Delta' or 'Omicron'
prepData <- function(ct,gene_,variant_,fig) {

  myformat="%d/%m/%Y"
  
  #set age category
  ct$age_category=''
  ct$age_category[which(ct$age>=12 & ct$age<=15)]='12-15'
  ct$age_category[which(ct$age>=16 & ct$age<=39)]='16-39'
  ct$age_category[which(ct$age>=40 & ct$age<=59)]='40-59'
  ct$age_category[which(ct$age>=60)]='60+'
  ct=ct[-which(ct$age<12),]
  if(length(which(is.na(ct$age)==TRUE))>0){ct=ct[-which(is.na(ct$age)==TRUE),]}
  ct$age_category=factor(ct$age_category,levels=c('16-39','12-15','40-59','60+'))
  
  if(fig==2){ct=ct[which(ct$age>=60),]}
  
  #set variant and gene
  if(variant_=='Delta'){ctt=ct[which(ct$new_gencodename==gene_ & ct$take4ratiodt>='2021-06-15' & ct$take4ratiodt<'2021-12-01' ),] }
  if(variant_=='Omicron'){ctt=ct[which(ct$new_gencodename==gene_ & ct$take4ratiodt>='2021-12-28'  & ct$take4ratiodt<='2022-01-29'),]}
  
  ctt=ctt[order(ctt$take4ratiodt),]
  
  
  ####################set cohorts 
  if(variant_=='Delta'){
    cohortname1='2-dose: 10-39 days'
    cohortname2='2-dose: 40-69 days'
    cohortname3='2-dose: 70+ days'
  }
  if(variant_=='Omicron'){
    cohortname1='2-dose: 10+ days'
    cohortname2='2-dose: 10+ days'
    cohortname3='2-dose: 10+ days'
  }
  
  ctt$cohort=NA
  ctt$cohort[which(ctt$recovered==0 & is.na(ctt$dff_vac1)==FALSE   & ctt$dff_vac1>=0 )]='1-dose'
  ctt$cohort[which(ctt$recovered==0 & is.na(ctt$dff_vac1)==FALSE   & ctt$dff_vac1<0)]='Unvaccinated'
  ctt$cohort[which(ctt$recovered==0 & is.na(ctt$dff_vac1)==TRUE )   ]='Unvaccinated'
  ctt$cohort[which(ctt$recovered==0 & is.na(ctt$dff_vac2)==FALSE  &  ctt$dff_vac2>=10  &  ctt$dff_vac2<40)]=cohortname1  
  ctt$cohort[which(ctt$recovered==0 & is.na(ctt$dff_vac2 )==FALSE  & ctt$dff_vac2>=40  &  ctt$dff_vac2<70)]=cohortname2
  ctt$cohort[which(ctt$recovered==0 & is.na(ctt$dff_vac2)==FALSE   & ctt$dff_vac2>=70    )]=cohortname3

  
  #remove if second dose was not given or was overdue 
  tempind=which( is.na(ctt$vacdt1)==FALSE &  is.na(ctt$vacdt2)==TRUE  & ctt$dff_vac1>27  )
  if(length(tempind)>0){ctt$cohort[tempind]=NA}
  tempind=which(is.na(ctt$vacdt1)==FALSE &   is.na(ctt$vacdt2)==FALSE &  ctt$vacdt2-ctt$vacdt1>27   & ctt$dff_vac1>27 )
  if(length(tempind)>0){ctt$cohort[tempind]=NA}
  
  
  #3-dose
  
  cohortname1='3-dose: 10-39 days'
  cohortname2='3-dose: 40-69 days'
  cohortname3='3-dose: 70+ days'
  if(fig==2){
    cohortname1='3-dose: 10+ days'
    cohortname2='3-dose: 10+ days'
    cohortname3='3-dose: 10+ days'
    
  }
  
  ctt$diff3=as.numeric(ctt$take4ratiodt-ctt$vacdt3)
  tempind=which(is.na(ctt$diff3)==FALSE & ctt$diff3>=0 &  ctt$diff3<10)
  if(length(tempind)>0){ctt$cohort[tempind]=NA}
  ctt$cohort[which(is.na(ctt$diff3)==FALSE & ctt$diff3>=10 & ctt$diff3<40)]=cohortname1
  ctt$cohort[which(is.na(ctt$diff3)==FALSE & ctt$diff3>=40 & ctt$diff3<70)]=cohortname2
  ctt$cohort[which(is.na(ctt$diff3)==FALSE & ctt$diff3>=70 )]=cohortname3
  
  #4-dose: 10+ days
  tempind=which(ctt$take4ratiodt>=ctt$vacdt4 & ctt$take4ratiodt<ctt$vacdt4+10 & is.na(ctt$vacdt4)==FALSE)
  if(length(tempind)>0){ctt=ctt[-tempind,]}
  tempind=which(ctt$take4ratiodt>=ctt$vacdt4+10 & is.na(ctt$vacdt4)==FALSE)
  if(length(tempind)>0){ctt$cohort[tempind]='4-dose: 10+ days'}
  
  #print(paste('tempind',length(tempind)))
  

  #recovered
  ctt$cohort[which(ctt$recovered==1)]='Recovered'
  
  
  #remove all recovered who received a dose before their 1st infection, or received their 
  #1st dose from 10 days before 2nd infection up to 2nd infection
  tempind=which(ctt$recovered==1 & is.na(ctt$vacdt1)==FALSE & ctt$vacdt1<=ctt$lastPCRfirstSequence)
  if(length(tempind)>0){ctt$cohort[tempind]=NA} 
  tempind=which(ctt$cohort=='Recovered' & is.na(ctt$vacdt1)==FALSE & ctt$vacdt1>ctt$take4ratiodt-10 & ctt$vacdt1<=ctt$take4ratiodt)
  if(length(tempind)>0){ctt$cohort[tempind]=NA}
  
  
  #Recovered+vaccine: recovered who received a single dose between two infection events (up to 10 days before 2nd infection)
  tempind=which(ctt$cohort=='Recovered' & is.na(ctt$vacdt1)==FALSE & ctt$vacdt1>=ctt$takedateFirst & ctt$vacdt1<=ctt$take4ratiodt-10)
  if(length(tempind)>0){ctt$cohort[tempind]='Recovered+vaccine'}
  

  #remove from Recovered+vaccine all those who received 2 or more doses between infection events
  tempind=which(ctt$cohort=='Recovered+vaccine' & is.na(ctt$vacdt2)==FALSE & ctt$vacdt2>=ctt$takedateFirst & ctt$vacdt2<=ctt$take4ratiodt)
  if(length(tempind)>0){ctt$cohort[tempind]=NA}
  tempind=which(ctt$cohort=='Recovered+vaccine' & is.na(ctt$vacdt3)==FALSE & ctt$vacdt3>=ctt$takedateFirst & ctt$vacdt3<=ctt$take4ratiodt)
  if(length(tempind)>0){ctt$cohort[tempind]=NA}
  tempind=which(ctt$cohort=='Recovered+vaccine' & is.na(ctt$vacdt4)==FALSE & ctt$vacdt4>=ctt$takedateFirst & ctt$vacdt4<=ctt$take4ratiodt)
  if(length(tempind)>0){ctt$cohort[tempind]=NA}
  
  
  if(fig==1){
  ctt=ctt[which(ctt$cohort=='Unvaccinated' | ctt$cohort=='2-dose: 10-39 days' |
                ctt$cohort=='2-dose: 40-69 days'  | ctt$cohort=='2-dose: 70+ days' | 
                ctt$cohort=='2-dose: 10+ days'    | ctt$cohort=='3-dose: 10-39 days' |
                ctt$cohort=='3-dose: 40-69 days'  | ctt$cohort=='3-dose: 70+ days'  | 
                ctt$cohort=='4-dose: 10+ days' |  ctt$cohort=='Recovered'  | ctt$cohort=='Recovered+vaccine'),]}
  
  if(fig==2){
    ctt=ctt[which(ctt$cohort=='Unvaccinated' |  ctt$cohort=='2-dose: 10+ days'    | 
                    ctt$cohort=='3-dose: 10+ days' |
                    ctt$cohort=='4-dose: 10+ days' |  ctt$cohort=='Recovered'  | ctt$cohort=='Recovered+vaccine'),]
  }
  
  if(fig==3){
    ctt=ctt[which(ctt$cohort=='Unvaccinated' |  ctt$cohort=='Recovered'),]
  }
  
  #partition calendar time to 1-week bins 
  dff_variant=7
  ctt$diffvariant=as.numeric(ctt$take4ratiodt-min(ctt$take4ratiodt)) 
  ctt$diffvariantbin=ctt$diffvariant%/%dff_variant
  
  return(ctt)
  
}
