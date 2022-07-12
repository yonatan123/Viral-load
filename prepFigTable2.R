
#Fig. 2

#input: ct: ct data table
prepFigTable60 <- function(ctt) {
  
  ctt$cohort=relevel(as.factor(ctt$cohort),ref = "Unvaccinated")
  labs_=unique(ctt$labname)
  cohorts=unique(ctt$cohort)
  vbins=length(unique(ctt$diffvariantbin))-1
  
  print(paste('cohorts',cohorts))
  
  
  vec1=as.vector(table(ctt$diffvariantbin)/sum(table(ctt$diffvariantbin)))
  vec1=vec1[2:length(vec1)]
  
  
  ctreg=lm(new_genvalue~cohort+as.factor(diffvariantbin)+labname+gender ,data=ctt)
  
  
  #intercepts
  facs_=c()
  for(j_ in 1:length(vec1)){facs_=c(facs_,paste0('as.factor(diffvariantbin)',j_))}
  uu=match(facs_,rownames(as.data.frame(ctreg$coefficients)))
  vec2=as.vector(ctreg$coefficients[uu])
  intercept_variant_bin=sum(vec1*vec2)
  
  
  
  
  vec1=as.vector(table(ctt$labname)/sum(table(ctt$labname)))[2:4]
  facs_=c('labnameLab2','labnameLab3','labnameLab4')
  uu=match(facs_,rownames(as.data.frame(ctreg$coefficients)))
  vec2=as.vector(ctreg$coefficients[uu])
  intercept_variant_lab=sum(vec1*vec2)
  
  
  
  facs_=c('M' )
  tab1=as.data.frame(table(ctt$gender)/sum(table(ctt$gender)))
  uu=match(facs_,tab1$Var1)
  vec1=as.vector(tab1$Freq[uu])
  
  facs_=c('genderM')
  uu=match(facs_,rownames(as.data.frame(ctreg$coefficients)))
  vec2=as.vector(ctreg$coefficients[uu])
  intercept_variant_gender=sum(vec1*vec2)
  
  
  
  labs_=unique(ctt$labname)
  labs_=c(labs_,'Combined')
  intercept=as.numeric(ctreg$coefficients[1])+intercept_variant_bin+intercept_variant_gender
  vc=vcov(ctreg)
  
  
  
  regTable=c()

  
  for(i in 1:length(labs_)){
    print(i)
    
    for(i1 in 1:length(cohorts)){
      
      if(labs_[i]=='Lab1' & cohorts[i1]=='Unvaccinated'){
        regTable=rbind(regTable,c(labs_[i],'Unvaccinated',intercept,intercept,intercept  ))
        next
      }
      
      
      facs1=c()
      if(labs_[i]=='Combined'){facs1=c(facs1,'labnameLab2','labnameLab3','labnameLab4')}
      
      if(labs_[i]!='Lab1' & labs_[i]!='Combined' ){facs1=c(facs1,paste0('labname',labs_[i]))}
      if(cohorts[i1]!='Unvaccinated'){facs1=c(facs1,paste0('cohort',cohorts[i1]))}
      
      
      uu=match(colnames(vc),facs1)
      indices=which(is.na(uu)==FALSE)
      H=vc[indices,indices]
      fc=as.vector(ctreg$coefficients[indices])
      
      vector=as.vector(table(ctt$labname)/sum(table(ctt$labname)))
      vector=vector[2:4]
      
      if(labs_[i]=='Combined' & cohorts[i1]!='Unvaccinated'){
        vector=c(1,vector)
        variance_1=vector%*%H%*%vector
        mean_1=intercept+sum(fc*vector)
      }
      
      if(labs_[i]=='Combined' & cohorts[i1]=='Unvaccinated'){
        variance_1=vector%*%H%*%vector
        mean_1=intercept+sum(fc*vector)
      }
      
      if(labs_[i]!='Combined'  ){
        variance_1=sum(H)
        mean_1=intercept+sum(fc)
      }
      
      
      
      
      vectorc=c()
      
      
      temp1=mean_1-1.96*sqrt(variance_1)
      temp2=mean_1+1.96*sqrt(variance_1)
      
      regTable=rbind(regTable,c(labs_[i],as.character(cohorts[i1]),temp1,mean_1,temp2  ))
      
      
      
    }
    
  }
  
  
  
  regTable=as.data.frame(regTable)
  colnames(regTable)<-c('labname','cohort','qlow','mean','qhigh')
  
  
  regTable$mean=as.numeric(regTable$mean)
  regTable$qlow=as.numeric(regTable$qlow)
  regTable$qhigh=as.numeric(regTable$qhigh)

  return(regTable)
}









