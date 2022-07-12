
#Fig 3

#input: ctt: ct data table
#variant_: 'delta' or 'omicron'
prepFigTableRecovered <- function(ctt,variant_) {
  
  
  ctt$diffrecovered=0
  tempind=which(ctt$cohort=='Recovered')
  ctt$diffrecovered[tempind]=as.numeric(ctt$take4ratiodt[tempind] -ctt$takedateFirst[tempind]-90)%/%60+1 
  
  if(variant_=='Delta'){
    tempind=which( ctt$diffrecovered==1 | ctt$diffrecovered==2 )
    ctt$diffrecovered[tempind]='1-2'
    tempind=which( ctt$diffrecovered==7 | ctt$diffrecovered==8 |   ctt$diffrecovered==9 )
    ctt$diffrecovered[tempind]='7-9'
  }
  
  if(variant_=='Omicron'){
    tempind=which(ctt$diffrecovered==4 | ctt$diffrecovered==3 )
    ctt$diffrecovered[tempind]='3-4'
    tempind=which(ctt$diffrecovered==10 | ctt$diffrecovered==9)
    ctt$diffrecovered[tempind]='9-10'
  }
  
  
  
  ctreg=lm(new_genvalue~as.factor(diffvariantbin)+as.factor(diffrecovered)+labname+gender+age_category ,data=ctt)
  
  
  #intercepts
  
  vec1=as.vector(table(ctt$diffvariantbin)/sum(table(ctt$diffvariantbin)))
  vec1=vec1[2:length(vec1)]
  
  facs_=c()
  for(j_ in 1:length(vec1)){facs_=c(facs_,paste0('as.factor(diffvariantbin)',j_))}
  uu=match(facs_,rownames(as.data.frame(ctreg$coefficients)))
  vec2=as.vector(ctreg$coefficients[uu])
  intercept_variant_bin=sum(vec1*vec2)
  
  
  
  facs_=c('M' )
  tab1=as.data.frame(table(ctt$gender)/sum(table(ctt$gender)))
  uu=match(facs_,tab1$Var1)
  vec1=as.vector(tab1$Freq[uu])
  
  facs_=c('genderM')
  uu=match(facs_,rownames(as.data.frame(ctreg$coefficients)))
  vec2=as.vector(ctreg$coefficients[uu])
  intercept_variant_gender=sum(vec1*vec2)
  
  
  vec1=as.vector(table(ctt$labname)/sum(table(ctt$labname)))[2:4]
  facs_=c('labnameLab2','labnameLab3','labnameLab4')
  uu=match(facs_,rownames(as.data.frame(ctreg$coefficients)))
  vec2=as.vector(ctreg$coefficients[uu])
  intercept_variant_lab=sum(vec1*vec2)
  
  
  
  
  
  facs_=c('12-15','40-59','60+')
  tab1=as.data.frame(table(ctt$age_category)/sum(table(ctt$age_category)))
  uu=match(facs_,tab1$Var1)
  vec1=as.vector(tab1$Freq[uu])
  
  facs_=c('age_category12-15','age_category40-59','age_category60+')
  uu=match(facs_,rownames(as.data.frame(ctreg$coefficients)))
  vec2=as.vector(ctreg$coefficients[uu])
  intercept_variant_age=sum(vec1*vec2)
  
  
  labs_=unique(ctt$labname)
  labs_=c(labs_,'Combined')
  intercept=as.numeric(ctreg$coefficients[1])+intercept_variant_bin+intercept_variant_gender+intercept_variant_age
  vc=vcov(ctreg)
  
  
  
  
  
  dff3=as.character(unique(ctt$diffrecovered))
  
  
  regTable=c()
  for(i in 1:length(labs_)){
    print(i)
    
    for(i1 in 1:length(dff3)){
      
      if(labs_[i]=='Lab1' & dff3[i1]=='0'){
        regTable=rbind(regTable,c(labs_[i],dff3[i1],intercept,intercept,intercept  ))
        next
      }
      
      
      facs1=c()
      if(labs_[i]!='Lab1' & labs_[i]!='Combined'){facs1=c(facs1,paste0('labname',labs_[i]))}
      if(labs_[i]=='Combined'){facs1=c(facs1,'labnameLab2','labnameLab3','labnameLab4')}
      if(dff3[i1]!='0'){facs1=c(facs1,paste0('as.factor(diffrecovered)',dff3[i1]))}
      
      
      
      uu=match(colnames(vc),facs1)
      indices=which(is.na(uu)==FALSE)
      H=vc[indices,indices]
      fc=as.vector(ctreg$coefficients[indices])
      
      vector=as.vector(table(ctt$labname)/sum(table(ctt$labname)))
      vector=vector[2:4]
      
      
      if(labs_[i]=='Combined' & dff3[i1]!='0'){
        vector=c(1,vector)
        variance_1=vector%*%H%*%vector
        mean_1=intercept+sum(fc*vector)
      }
      
      if(labs_[i]=='Combined' & dff3[i1]=='0'){
        variance_1=vector%*%H%*%vector
        mean_1=intercept+sum(fc*vector)
          }
      
      if(labs_[i]!='Combined'  ){
        variance_1=sum(H)
        mean_1=intercept+sum(fc)
      }
      
      
      
      dff_=dff3[i1]
      vectorc=c()
      
      
      
      temp1=mean_1-1.96*sqrt(variance_1)
      temp2=mean_1+1.96*sqrt(variance_1)
      
      regTable=rbind(regTable,c(labs_[i],dff_,temp1,mean_1,temp2  ))
      
    }
    
  }
  
  
  regTable=as.data.frame(regTable)
  colnames(regTable)<-c('labname','cohort','qlow','mean','qhigh')
  
  
  regTable$mean=as.numeric(regTable$mean)
  regTable$qlow=as.numeric(regTable$qlow)
  regTable$qhigh=as.numeric(regTable$qhigh)
  regTable$variant=variant_
  
  
  if(variant_=='Delta'){
    regTable$cohort[which(regTable$cohort=='0')]= 'Unvaccinated'
    regTable$cohort[which(regTable$cohort=='1-2')]= '4-7'
    regTable$cohort[which(regTable$cohort=='3')]= '8-9'
    regTable$cohort[which(regTable$cohort=='4')]= '10-11'
    regTable$cohort[which(regTable$cohort=='5')]= '12-13'
    regTable$cohort[which(regTable$cohort=='6')]= '14-15'
    regTable$cohort[which(regTable$cohort=='7-9')]= '16-20'
  }
  
  
  if(variant_=='Omicron'){
    regTable$cohort[which(regTable$cohort=='0')]= 'Unvaccinated'
    regTable$cohort[which(regTable$cohort=='1')]= '4-5'
    regTable$cohort[which(regTable$cohort=='2')]= '6-7'
    regTable$cohort[which(regTable$cohort=='3-4')]= '8-11'
    regTable$cohort[which(regTable$cohort=='5')]= '12-13'
    regTable$cohort[which(regTable$cohort=='6')]= '14-15'
    regTable$cohort[which(regTable$cohort=='7')]= '16-17'
    regTable$cohort[which(regTable$cohort=='8')]= '18-19'
    regTable$cohort[which(regTable$cohort=='9-10')]= '20-23'
  }
  
  
  return(regTable)
}





