
#This code accompanies the paper: Viral load dynamics of SARS-CoV-2 Delta and Omicron variants
#following multiple vaccine doses and previous infection
#This files contains code that produces the main figures and regression tables

#
# ctTable.RData: Ct dataset to be loaded, obtained from preprocessing. Not available in our folder.
#Contains 20 columns:
# "age"
#"new_gencodename": gene name
#"new_genvalue": Ct value
#"take4ratiodt": PCR sampling date
# "vacdt1": 1st dose date
# "vacdt2": 2nd dose date
# "vacdt3": 3rd dose date
# "vacdt4": 4th dose date
# "dff_vac1": take4ratiodt-vacdt1 (difference between)
# "dff_vac2": take4ratiodt-vacdt2 (difference between dates)
# "recovered": binary, recovered or not
# "labname"
# "gender"
# "takedateFirst": patient's first positive PCR 
# "lastPCRfirstSequence": latest PCR of patient's first infection event
# "patient_ID": patient ID
#"new_labresultid": Test ID
#"ct_date": Ct creation date
# "earlyPCR": patient's earliest positive PCR date 
# "dffEarlyPCR": take4ratiodt-earlyPCR (difference between dates)



#load package
#source('libraries.R')

#Fig 1

rm(list = ls())
load('ctTable.RData')


source('prepData.R') #load functions
source('prepFigTable1.R')

ctt=prepData(ct,'N','Delta',1)
regtable1=prepFigTable(ctt,'Delta')

ctt=prepData(ct,'N','Omicron',1)
regtable2=prepFigTable(ctt,'Omicron')

U=rbind(regtable1,regtable2)
U=U[-which(U$cohort=='Recovered+vaccine' | U$cohort=='4-dose: 10+ days'),]
U=U[which(U$labname=='Combined'),]

level_order2=c('Unvaccinated','2-dose: 10-39 days','2-dose: 40-69 days','2-dose: 70+ days','2-dose: 10+ days',
               '3-dose: 10-39 days','3-dose: 40-69 days','3-dose: 70+ days','Recovered')
cbPalette=c('gray0','gray90','gray90','gray90','gray90','gray60','gray60','gray60','gray30')
minc_=24
maxc_=28
dodge <- position_dodge(width=0.9)


p <- ggplot(U, aes(  y=as.numeric(mean),x = factor(cohort, level=level_order2),fill = factor(cohort, level=level_order2) ))+ 
  theme(text = element_text(size=12))  +labs(title=paste(''),x='',y ='Ct-value')+
  geom_col(position = dodge,col='white') +geom_errorbar(aes(ymin=qlow,ymax=qhigh), position = dodge, width = 0.25)+ 
  scale_y_continuous(limits=c(minc_,maxc_),oob = rescale_none)+ 
  theme(legend.title = element_blank())+scale_fill_manual(values=cbPalette) +theme(plot.title = element_text(hjust = 0.5))+ facet_wrap( ~variant,scales = "free_x") + 
  theme(legend.position="none")+ theme( strip.text = element_text(size=14,face = "bold"))
p<-p+scale_x_discrete(guide = guide_axis(angle = 45))+ theme(axis.text.x = element_text(face="bold"))
p<-p+ theme(panel.spacing = unit(4, "lines"))
p+theme(aspect.ratio = .5)







#Fig 2

rm(list = ls())
load('ctTable.RData')
source('prepData.R')
source('prepFigTable2.R')


ctt=prepData(ct,'N','Omicron',2)
U=prepFigTable60(ctt)
U=U[-which(U$cohort=='Recovered+vaccine' ),]
U=U[which(U$labname=='Combined'),]
level_order2=c('Unvaccinated','2-dose: 10+ days',
               '3-dose: 10+ days','4-dose: 10+ days','Recovered')
cbPalette=c('gray0','gray90','gray60','gray50','gray30')
minc_=24
maxc_=28
dodge <- position_dodge(width=0.9)

p <- ggplot(U, aes(  y=as.numeric(mean),x = factor(cohort, level=level_order2),fill = factor(cohort, level=level_order2) ))+ 
  theme(text = element_text(size=12))  +labs(title=paste(''),x='',y ='Ct-value')+
  geom_col(position = dodge,col='white') +geom_errorbar(aes(ymin=qlow,ymax=qhigh), position = dodge, width = 0.25)+ 
  scale_y_continuous(limits=c(minc_,maxc_),oob = rescale_none)+ 
  theme(legend.title = element_blank())+scale_fill_manual(values=cbPalette) +theme(plot.title = element_text(hjust = 0.5)) + 
  theme(legend.position="none")+ theme( strip.text = element_text(size=14,face = "bold"))
p<-p+scale_x_discrete(guide = guide_axis(angle = 45))+ theme(axis.text.x = element_text(face="bold"))
p<-p+ theme(panel.spacing = unit(4, "lines"))
p+theme(aspect.ratio = .5)



#Fig 3
rm(list = ls())
load('ctTable.RData')
source('prepData.R')
source('prepFigTable3.R')


ctt=prepData(ct,'N','Delta',3)
U1=prepFigTableRecovered(ctt,'Delta')
ctt=prepData(ct,'N','Omicron',3)
U2=prepFigTableRecovered(ctt,'Omicron')
U=rbind(U1,U2)
U=U[which(U$labname=='Combined'),]

minc_=24.25
maxc_=28.75
dodge <- position_dodge(width=0.9)

cbPalette=c('gray0','gray30','gray30','gray30','gray30','gray30','gray30','gray30','gray30','gray30','gray30'
            ,'gray30','gray30','gray30','gray30')





U=rbind(U,c('',' ',0,0,0,'Delta'))
U=rbind(U,c('',' ',0,0,0,'Omicron'))
level_order2=c('Unvaccinated',' ','4-5','6-7','4-7','8-9',
               '10-11','8-11',
               '12-13','14-15','16-17',
               '18-19','20-23','16-20')


U$qlow=as.numeric(U$qlow)
U$qhigh=as.numeric(U$qhigh)
U$mean=as.numeric(U$mean)


minc_=24.25
maxc_=28.75
U$cohort=factor(U$cohort,level=level_order2)

p <- ggplot(U, aes(  y=as.numeric(mean),x = cohort,fill = cohort ))+ 
  theme(text = element_text(size=12))  +labs(title=paste(''),x='',y ='Ct value')+
  geom_col(position = dodge,col='white') +geom_errorbar(aes(ymin=qlow,ymax=qhigh), position = dodge, width = 0.25)+ 
  scale_y_continuous(limits=c(minc_,maxc_),oob = rescale_none)+ 
  theme(legend.title = element_blank())+scale_fill_manual(values=cbPalette) +theme(plot.title = element_text(hjust = 0.5))+ 
  theme(legend.position="none")+ theme( strip.text = element_text(size=14,face = "bold"))
p<-p+scale_x_discrete(guide = guide_axis(angle = 0))+ theme(axis.text.x = element_text(face="bold"))
p<-p+ theme(panel.spacing = unit(1, "lines"))
p<-p+ geom_rect(data = data.frame(variant='Delta'), aes(xmin = 2, xmax = 8.8, ymin = 0, ymax = Inf), alpha = 0.3, fill="grey", inherit.aes = FALSE)
p<-p+ geom_rect(data = data.frame(variant='Omicron'), aes(xmin = 2, xmax = 10.5, ymin = 0, ymax = Inf), alpha = 0.3, fill="grey", inherit.aes = FALSE)
p<-p+ facet_wrap( ~variant,scales = "free_x")


dat_text <- data.frame(
  label = c( "Recovered (months from first infection)","Recovered (months from first infection)"),
  variant=c("Delta","Omicron"),
  cohort     = c( "10-11","8-11"),
  mean     = c(28.5,28.5)
)

p<- p + geom_text(data    = dat_text,mapping= aes(x = cohort, y = mean, label = label),size=3.5, position = position_dodge(width = 1),
                  hjust = 0.35)


p+theme(aspect.ratio = .5)







#Table 2A, gene='N', variant='Delta'
rm(list = ls())
load('ctTable.RData')
source('prepData1.R')
ctt=prepData(ct,'N','Delta',1)



ctt$cohort=relevel(as.factor(ctt$cohort),ref = "Unvaccinated")
ctreg=lm(new_genvalue~cohort+gender+age_category+labname+factor(diffvariantbin) ,data=ctt)
print(summary(ctreg))
ctreg=rq(new_genvalue~cohort+gender+age_category+labname+factor(diffvariantbin) ,data=ctt,tau=.5)
ctreg=rq(new_genvalue~cohort+gender+age_category+labname+factor(diffvariantbin) ,data=ctt,tau=.2)




# Table 4A - Delta
rm(list = ls())
load('ctTable.RData')
source('prepData1.R')
ctt=prepData(ct,'N','Delta',3)


ctt$diffrecovered=0
tempind=which(ctt$cohort=='Recovered')
ctt$diffrecovered[tempind]=(as.numeric(ctt$take4ratiodt[tempind] -ctt$takedateFirst[tempind])-90)%/%60+1 

tempind=which( ctt$diffrecovered==1 | ctt$diffrecovered==2 )
ctt$diffrecovered[tempind]='1to2'
tempind=which( ctt$diffrecovered==7 | ctt$diffrecovered==8 |   ctt$diffrecovered==9 )
ctt$diffrecovered[tempind]='7to9'
ctt$diffrecovered[which(ctt$diffrecovered==0)]="Unvaccinated"




ctt$cohort=relevel(as.factor(ctt$cohort),ref = "Unvaccinated")
ctreg=lm(new_genvalue~cohort+gender+age_category+labname+factor(diffvariantbin) ,data=ctt)








# Table 4A - Omicron
rm(list = ls())
load('ctTable.RData')
source('prepData1.R')
ctt=prepData(ct,'N','Omicron',3)


ctt$diffrecovered=0
tempind=which(ctt$cohort=='Recovered')
ctt$diffrecovered[tempind]=(as.numeric(ctt$take4ratiodt[tempind] -ctt$takedateFirst[tempind])-90)%/%60+1 

tempind=which(ctt$diffrecovered==4 | ctt$diffrecovered==3 )
ctt$diffrecovered[tempind]='3to4'
tempind=which(ctt$diffrecovered==10 | ctt$diffrecovered==9)
ctt$diffrecovered[tempind]='9to10'
ctt$diffrecovered[which(ctt$diffrecovered==0)]="Unvaccinated"

ctt$cohort=relevel(as.factor(ctt$cohort),ref = "Unvaccinated")
ctreg=lm(new_genvalue~cohort+gender+age_category+labname+factor(diffvariantbin) ,data=ctt)



