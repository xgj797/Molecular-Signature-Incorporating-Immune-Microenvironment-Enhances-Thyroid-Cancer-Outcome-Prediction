
# Data cleaning function-----------------------------------------------------------
data.cleaning<-function(d0,clic.var,survival.var,predictor.var,poor.var){
  
  d0<-within(d0,{
    
    TERT.mutation<-factor(TERT.mutation,levels=c(0,1),labels=c('No','Yes'))
    BRAF.mutation<-factor(BRAF.mutation,levels=c(0,1),labels=c('No','Yes'))
    TP53.mutation<-factor(TP53.mutation,levels=c(0,1),labels=c('No','Yes'))
    PIK3CA.mutation<-factor(PIK3CA.mutation,levels=c(0,1),labels=c('No','Yes'))
    RAS.mutation<-factor(RAS.mutation,levels=c(0,1),labels=c('No','Yes'))
    TERT.TP53.PIK3CA.mut<-factor(TERT.TP53.PIK3CA.mut,levels=c(0,1),labels=c('No','Yes'))
    
    Cancer.associated.fibroblast_EPIC_high<-factor(Cancer.associated.fibroblast_EPIC_high,levels=c(0,1),labels=c('Low','High'))
    Cancer.associated.fibroblast_MCPCOUNTER_high<-factor(Cancer.associated.fibroblast_MCPCOUNTER_high,levels=c(0,1),labels=c('Low','High'))
    
    BRS.PN<-ifelse(is.na(BRS),NA,ifelse(BRS>=0,'RAS-Like','BRAF-Like'))
    BRS.PN<-factor(BRS.PN,levels=c('BRAF-Like','RAS-Like'))
    
    ATS.PN<-ifelse(is.na(ATS),NA,ifelse(ATS>=0,'Positive','Negative'))
    ATS.PN<-factor(ATS.PN,levels=c('Negative','Positive'))

    label(Diagnosis)<-'Thyroid lesion subtype'
    label(Age.resection)<-'Age of patient when resection made'
    label(Sex)<-'Gender'
    label(Prior.cancer)<-'Prior cancer'
    
    label(Radiation.history)<-'History of radiation prior to primary tumor diagnosis'
    label(Fam.hist.thyroid.cancer)<-'Family history of thyroid cancer'
    label(TERT.mutation)<-'Mutation in the TERT promoter'
    label(BRAF.mutation)<-'Mutation in the BRAF gene'
    label(TP53.mutation)<-'Mutation in the TP53 gene'
    label(PIK3CA.mutation)<-'Mutation in the PIK3CA gene'
 
    label(RAS.mutation)<-'RAS mutation'
    label(TERT.TP53.PIK3CA.mut)<-'Mutation in TERT, TP53, PIK3CA'
    label(BRS)<-'BRS'
   
    label(BRS.PN)<-'BRS'
    label(ATS)<-'ATS'
    label(ATS.PN)<-'ATS'
  
    label(BRAFPoorOutcome_526)<-'BRAFPoorOutcome 526'
     
    label(Procedure.type)<-'Type of procedue done'
    label(Vascular.invasion)<-'Vascular invasion'
    label(Capsular.invasion)<-'Capsular invasion'
    label(Extrathyroidal.extension )<-'Extrathyroidal extension '
    label(Submitted.nodes)<-'Number of nodes submitted'
    label(Positive.nodes )<-'Number of nodes positive for cancer'
    label(T.stage)<-'Tumor stage'
    label(n.stage)<-'Lymph node stage'
    label(m.stage.at.sample.date)<-'Metstatic stage at diagnosis'
    label(Primary.tumor.size.cm)<-'Primary tumor size (cm)'
    label(p.AJCC.8th.ed.Stage)<-'8th edition (new) stage assignment'
    label(Local.met.size.mm)<-'Size of metastasis (mm)'
    label(Local.met.location)<-'Up-to-date local met location'
    
    label(Dist.met.location)<-'Up-to-date distant met location'
    label(Distant.met.date)<-'Date of discovery of distant met' #
    label(Recurrence.after.remission)<-'Did thyroid cancer recur after remission?'
    label(Recurrence.after.remission)<-'Location where recurrence was found'
    label(Number.of.surgeries)<-'Total number of thyroid surgeries the patient has had'
    label(TG.presurgery)<-'Thyroglobulin before surgery of primary (within 1 year)'
    label(TSH.presurgery)<-'Thyroid stimulating hormone before surgery of primary (within 1 year)'
    label(TG.6.weeks.post.surg)<-'Thyroglobulin 6 weeks after surgery of primary'
    label(TSH.6.weeks.post.surg)<-'Thyroid stimulating hormone 6 weeks after surgery of primary'
    label(Iodine.avid)<-'Radioactive iodine avid residual disease detected'
    
    label(WBC.presurg.k)<-'White blood count before surgery'
    label(Weight.presurg.kg)<-'Patient weight before primary surgery (kg)'
    label(Height.cm)<-'Patient height (cm)'
    label(BMI)<-'Body mass index'
    label(Ethnicity)<-'Ethnicity'
    label(Date.last.followup)<-'Date patient had last followup'
    label(Date.last.checked)<-'Date last checked patient charts'
    label(Death.from.disease)<-'Died from thyroid cancer'
    label(Date.of.death)<-'Date that patient died from thyroid cancer'
    
    #No event, then censored at date.last.followup
    Initial.therapy.complete.date<-as.Date(Initial.therapy.complete.date,"%Y-%m-%d") #start date 1
    Best.response.date<-as.Date(Best.response.date,"%Y-%m-%d")  #start date 2
    #  #progression-free survival ends at Recur.or.progress.date
    Recur.or.progress.date<-as.Date(Recur.or.progress.date,"%Y-%m-%d")
    Date.last.followup<-as.Date(Date.last.followup,"%Y-%m-%d")
   
    Date.of.death<-as.Date(Date.of.death,"%Y-%m-%d")
 
    label(Cancer.associated.fibroblast_EPIC_high)<-'EPIC CAF'
    label(Cancer.associated.fibroblast_MCPCOUNTER_high)<-'MCPCOUNTER CAF'
    
    label(Cancer.associated.fibroblast_EPIC)<-'EPIC CAF'
    label(Cancer.associated.fibroblast_MCPCOUNTER)<-'MCPCOUNTER CAF'
    
    
    Number.of.surgeries1<-ifelse(Number.of.surgeries>=2,'more than one',Number.of.surgeries)
    #ethnicity (W = white, B = black, A = asian, HL = hispanic latino, WHL = white hispanic latino, AHL = asian hispanic latino, NR = not recorded, WA=?, NH=non-hispanic,
    
    stage<-ifelse(p.AJCC.8th.ed.Stage%in%c('I','II'),'I and II','III and IV')
    
  })
  
  
 
  d0<-within(d0,{
    #PFS
    event.PFS<-ifelse(is.na(as.character(Recur.or.progress.date)),0,1)
    temp.time<-ifelse(is.na(as.character(Recur.or.progress.date)),as.character(Date.last.followup),as.character(Recur.or.progress.date))
    temp.time<-as.Date(temp.time,"%Y-%m-%d")
    PFS.from.best.response<-as.numeric(temp.time - Best.response.date)/30.44
    PFS.from.ini.therapy.complete<-as.numeric(temp.time - Initial.therapy.complete.date)/30.44
    label(event.PFS)<-'Recurr or progress'
    label(PFS.from.best.response)<-'Recurr or progress or censored time from best response date'
    label(PFS.from.ini.therapy.complete)<-'Recurr or progress or censored time from initial therapy complete'
 #OS
    event.OS<-ifelse(is.na(as.character(Date.of.death)),0,1)
    temp.time<-ifelse(is.na(as.character(Date.of.death)),as.character(Date.last.followup),as.character(Date.of.death))
    temp.time<-as.Date(temp.time,"%Y-%m-%d")
    OS.from.ini.therapy.complete<-as.numeric(temp.time - Initial.therapy.complete.date)/30.44
    label(event.OS)<-'Death'
    label(OS.from.ini.therapy.complete)<-'Death or censored time from initial therapy complete'
    
    
  })
  
  
  d1<-subset(d0,select=c(clic.var,
                         survival.var,
                         predictor.var,  
                     #    'Diagnosis.with.iFVPTC.and.eFVPTC', #add on 0509 for BRAFPoorOutcome_526_PTCandIFVPTConly, we want to subset BRAFPoorOutcome_526_PTCandIFVPTConly %in% c('PTC','IFVPTC')
                         'Initial.therapy.complete.date','Recur.or.progress.date','Date.last.followup','Date.of.death',
                         poor.var
                         ))
  
  d1
}


# Data --------------------------------------------------------------------
pacman::p_load(xtable,knitr,summarytools, survival, survminer, formula.tools,
               car,plyr,reshape,ggplot2,Hmisc,openxlsx,rms)
#library(plyr) #ddply
library(dplyr) #top_n
PC=FALSE
if (PC==TRUE){
  path.in<-'\\\\cqsresearch\\Research drive\\scchen\\2022\\project\\WeissVivian\\XuGeorge\\GitHub\\data\\'
} else{
  path.in<-'/Users/sheauchiannchen/Documents/chen/Sheau-Chiann_VU/2022/project/WeissVivian/XuGeorge/GitHub/data/'
}



d.data.name<-data.frame(dataset=c('A','B','C','D','E'),
                        names=c('malignant','malignant well-diff',
                                'PTC/IFVPTC','malignant well-diff before aggression','PTC/IFVPTC before aggression'))

file.name<-'VUMC.cohort.GX_9-15-22_for-Fei.xlsx' #data set a, b, f, g, h
clic.var<-c("Arbitrary.patient.identifier","IP",
            "Age.resection","Sex","Prior.cancer","Radiation.history",
            "Fam.hist.thyroid.cancer",
            "Primary.tumor.size.cm","p.AJCC.8th.ed.Stage","stage",
            "Number.of.surgeries", "Number.of.surgeries1","Iodine.avid",
            "Weight.presurg.kg","Height.cm",
            "Ethnicity")
survival.var<-c('event.OS','event.PFS','PFS.from.ini.therapy.complete','OS.from.ini.therapy.complete','Recur.or.progress.date',
               'Date.last.followup','Date.of.death')


predictor.var<-c('TERT.mutation', 'TP53.mutation', 'PIK3CA.mutation',
                 'TERT.TP53.PIK3CA.mut', 'Cancer.associated.fibroblast_EPIC_high', 'Cancer.associated.fibroblast_MCPCOUNTER_high',
                 'BRS.PN','ATS.PN'  
                 )

#poor variable-->aggressive.disease
poor.var<-c('Aggressive.disease')

# All local malignant (excel sheet A) ------------------------------------
d0.a<-openxlsx::read.xlsx((paste(path.in,file.name,sep='')),detectDates=TRUE,
                        sheet='a - local malignant',na.strings = c("NA","NaN"," "))


d0<-d0.a
t.names<-names(d0)
duplicate.var<-names(which(table(t.names)>1))
which.no<-which(t.names==duplicate.var) #duplicate name column
d0<-d0[,-which.no[1]]

d0.a.clean<-data.cleaning(d0=d0,clic.var=clic.var,survival.var=survival.var,predictor.var=predictor.var,poor.var=poor.var)



# local malignant well-diff (excel sheet B) -----------
d0.b<-openxlsx::read.xlsx((paste(path.in,file.name,sep='')),detectDates=TRUE,
                          sheet='b - well diff',na.strings = c("NA","NaN"," "))

d0<-d0.b
t.names<-names(d0)
duplicate.var<-names(which(table(t.names)>1))
which.no<-which(t.names==duplicate.var) #duplicate name column
d0<-d0[,-which.no[1]]
d0.b.clean<-data.cleaning(d0=d0,clic.var=clic.var,survival.var=survival.var,predictor.var=predictor.var,poor.var=poor.var)



# PTC/IFVPTC (excel sheet C) ---------------------------------------------------
 d0.c<-openxlsx::read.xlsx((paste(path.in,file.name,sep='')),detectDates=TRUE,
                          sheet='c - IFVPTC-PTC',na.strings = c("NA","NaN"," "))
d0<-d0.c
t.names<-names(d0)
duplicate.var<-names(which(table(t.names)>1))
which.no<-which(t.names==duplicate.var) #duplicate name column
d0<-d0[,-which.no[1]]
d0.c.clean<-data.cleaning(d0=d0,clic.var=clic.var,survival.var=survival.var,predictor.var=predictor.var,poor.var=poor.var)


#before aggression (excel sheet D) ----------------------
d0.d<-openxlsx::read.xlsx((paste(path.in,file.name,sep='')),detectDates=TRUE,
                          sheet='d - future aggression',na.strings = c("NA","NaN"," "))
d0<-d0.d
t.names<-names(d0)
duplicate.var<-names(which(table(t.names)>1))
which.no<-which(t.names==duplicate.var) #duplicate name column
d0<-d0[,-which.no[1]]
d0.d.clean<-data.cleaning(d0=d0,clic.var=clic.var,survival.var=survival.var,predictor.var=predictor.var,poor.var=poor.var)


#PTC/IFVPTC before aggression (excel sheet E) --------------------------------------
d0.e<-openxlsx::read.xlsx((paste(path.in,file.name,sep='')),detectDates=TRUE,
                          sheet='e - futureagression IFVPTC-PTC',na.strings = c("NA","NaN"," "))
d0<-d0.e
t.names<-names(d0)
duplicate.var<-names(which(table(t.names)>1))
which.no<-which(t.names==duplicate.var) #duplicate name column
d0<-d0[,-which.no[1]]
d0.e.clean<-data.cleaning(d0=d0,clic.var=clic.var,survival.var=survival.var,predictor.var=predictor.var,poor.var=poor.var)

save(d0.a.clean,d0.b.clean,d0.c.clean,d0.d.clean,d0.e.clean,
    clic.var, survival.var,predictor.var,
    d.data.name,
     file=paste0(path.in,'Thyroid.data.2022.Rdata'))

