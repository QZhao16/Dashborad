setwd("C:/zhao088/Desk20210304/Manuscript-SIR model simu/R code/0.analysis/figure plotting/Responses time/function loading2")
# load function Part 1
source("12 common EWs-1.R")
source("13.vest_functions-1.R")
source("0.surrogates5.R")
source("estimateRe function.R")
source("function.EWS.computation2.R")
source("function.loading R package.R")
source("response.time-function.EWS.computation.R")
source("find consectiveNumber.R")
source("other.functions.R")
source("modifys.R")


leng_window=10
# Part 1: ADIO data
setwd("C:/zhao088/Desk20210304/Manuscript-SIR model simu/R code/AIDO data(download)")
#disease_data<-read.csv("AIDO.disease-Anthrax.all.csv")

#for (leng_window in 7:21) { 
  
finals<-NULL

all_diseases<-c("Anthrax","Brucellosis","Campylobacteriosis","Chikungunya",
                "Cholera","Dengue","Ebola",
                #"Foot And Mouth Disease", #NO
                "Gastroenteritis",
                #"Japanese Encephalitis",
                "Lassa Fever","Leptospirosis",
                "Malaria","Marburg","Measles",
                #"Meningococcal Disease", #NO
                "Middle East Respiratory Syndrome (MERS-CoV)","Monkeypox",
                #"Mumps",
                "Nipah","Norovirus","Novel Influenza A","Pertussis","Plague",
                "Polio",
                #"Porcine Epidemic Diarrhea Virus (PEDV)", #NO
                "Q Fever","Rift Valley Fever",
                "Rubella","Salmonellosis",
                #"Severe Acute Respiratory Syndrome (SARS)",
                "Shiga Toxin-Producing E. Coli (STEC)",
                "Shigellosis","Tularemia","West Nile Virus","Yellow Fever",
                "Zika") #change 1
length(all_diseases)
#
ID_disease=1
for (ID_disease in 1:length(all_diseases)) {
#within each disease
the_diseases=all_diseases[ID_disease]
disease_data<-read.csv(file=paste("AIDO.disease-",the_diseases,".all.csv", sep=""))

the_length<-length( unique(disease_data$full_name) )
unique(disease_data$full_name)
#
i_disease=1#7
for (i_disease in 1:the_length){

the_countrys<-unique(disease_data$full_name)[i_disease]
the_countrys
subs_data=disease_data[which(disease_data$full_name %in% unique(disease_data$full_name)[i_disease]),]
#plot(subs_data$value, type="l")
subs_data$Incidence<- subs_data$value #Incidence
subs_data$date<- as.Date(subs_data$start, "%Y-%m-%d") #transfer date
subs_data=subs_data[!is.na(subs_data$date), ] #remove NA date
library(tidyverse)
library(dplyr)
subs_data=subs_data[!duplicated(subs_data$date), ] # remove duplicate
subs_data=subs_data[order(subs_data$date),] #sort from old to new date

if (nrow(subs_data)>leng_window+3 ){
  The_frequece<-The_frequecyss(subs_data)
#if (The_frequece=="daily"){
 
  # compute the Re (effective reproduction numebr)
  subs_data<-Re(subs_data) 
  # select columns that we use
subs_data<-subs_data %>% dplyr::select(date, Rt, CI_Rt_uppper, CI_Rt_lower,
         url, name, full_name, admin_level_name,
         longitude, latitude, disease_name, Incidence)
colnames(subs_data)<-c(colnames(subs_data)[1], "Re", "Up.confidence.interval.Re", "Low.confidence.interval.Re",colnames(subs_data)[5:ncol(subs_data)])
# total casese for a time series
subs_data$total_cases<-sum(subs_data$Incidence,na.rm = TRUE)
# data resolution (daily, weekly, or monthly)
subs_data$Resolution<-The_frequece
# determine the outbreak point
date_position=subs_data$date[floor(subs_data$Low.confidence.interval.Re)>=1]
# outbreak start and ends
start_end<-start_end_outbreak(date_position, subs_data)
#Gaussian Smoothing
subs_data$Incidence<-Gaussian_Smoothing(subs_data$Incidence)
#
# tau time series for Best early warning signals (Index of dispersion)
taus<-tau_time_series(start_end, subs_data)
if (length(taus)>0) {
# merge into original raw time series dataset.
colnames(subs_data)[1]<-'date2'; subs_data$date2<-as.Date(subs_data$date2,"%Y-%m-%d"); taus$date2<-as.Date(taus$date2,"%Y-%m-%d"); 
#
subs_data_2<-as.data.frame(left_join(subs_data, taus, by = "date2"))
if ("Re" %in% colnames(subs_data_2)) { subs_data_2<-subs_data_2 %>% dplyr::select(date2,Re,Up.confidence.interval.Re,Low.confidence.interval.Re,full_name, admin_level_name, longitude, latitude, disease_name, Incidence, total_cases, Resolution,seq_length ,Ews_timeSeries, tau_value, p_value, outbreak_posit1, outbreak_posit2,date1,EWs,length_time_series,outbreaks, number, leng_window, Yes_No,  start_day,    end_day)} else {
                                       subs_data_2<-subs_data_2 %>%  dplyr::select(date2,Re.x,Up.confidence.interval.Re.x,Low.confidence.interval.Re.x,full_name.x, admin_level_name.x, longitude.x, latitude.x, disease_name.x, Incidence.x, total_cases.x, Resolution.x,seq_length ,Ews_timeSeries,tau_value, p_value, outbreak_posit1, outbreak_posit2,date1,EWs,length_time_series, outbreaks, number, leng_window, Yes_No,  start_day,    end_day)
                                       colnames(subs_data_2)<-c("date2", "Re", "Up.confidence.interval.Re", "Low.confidence.interval.Re", "full_name", "admin_level_name", "longitude", "latitude", "disease_name", "Incidence", "total_cases", "Resolution", "seq_length", "Ews_timeSeries", "tau_value", "p_value", "outbreak_posit1", "outbreak_posit2", "date1", "EWs", "length_time_series", "outbreaks", "number", "leng_window", "Yes_No", "start_day", "end_day")}
}

finals<-rbind(finals,subs_data_2)

}# if( nrow(subs_data)>=5 )
#}#if (The_frequece=="daily")
}#if (length(date_position)>1)
}# the  for (i_disease in 1:the_length){

#}

dim(finals) #210  16
setwd("C:/zhao088/Desk20210304/Manuscript-SIR model simu/0.only.Tau(AUC)/minimumCumul = 1, remove 1/To Sam(Dashboard)")


write.csv(finals,"dataset_22.csv")
#write.csv(sumary_finals,paste("sumary_finals.",leng_window,"4.csv", sep=""))
#setwd("C:/zhao088/Desk20210304/Manuscript-SIR model simu/R code/AIDO data(download)")

#} # for (leng_window in 5:30) { 


