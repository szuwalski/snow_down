#==loking at NBS by size
#++combine the NBS and EBS data from 2017-2019
#==have a multipanel plot with columns of different slices of the population
library(dplyr)
library(reshape2)
library(tidyr)
library(ggplot2)
library(autoimage)
library(akima)  #install.packages("akima")
library(viridis)
library(broom)
library(maps)
library("rnaturalearth")
library(interp)
library(RColorBrewer)
library(reshape2) # for melt
library(mgcv)  
library(PBSmapping)
library(mapdata)    #some additional hires data
library(maptools)   #useful tools such as reading shapefiles
library(mapproj)
library(metR)
library(ggplot2)
#install.packages("rnaturalearthdata")
library(reshape)
library(dplyr)
library(ggplot2)
library(ggridges)
library(gridExtra)
library(png)
library(grid)
library(PBSmodelling)
library(patchwork)
library(cowplot)
#=================================
# EBS data
#==============================
survDAT<-read.csv("data/EBSCrab_Haul/EBSCrab_Haul.csv",header=T,skip=5)

drvYear<-as.numeric(substr(survDAT$CRUISE,1,4))
SurvYR<-unique(drvYear)
AllStation<-unique(survDAT$GIS_STATION)

#==plot GIS stations
AllStnLoc<-matrix(ncol=2,nrow=length(AllStation))
for(w in 1:length(AllStation))
{
  temp<-survDAT[survDAT$GIS_STATION==AllStation[w],]
  AllStnLoc[w,1]<-temp$MID_LATITUDE[1]
  AllStnLoc[w,2]<-temp$MID_LONGITUDE[1]
}

#================================
#==find survey densities (total)
#================================
nmiSurv<-140350
DensityM_45_55<-matrix(nrow=length(SurvYR),ncol=length(AllStation))
DensityF_45_55<-matrix(nrow=length(SurvYR),ncol=length(AllStation))
DensityM_55_65<-matrix(nrow=length(SurvYR),ncol=length(AllStation))
DensityM_65_75<-matrix(nrow=length(SurvYR),ncol=length(AllStation))
DensityM_45_85<-matrix(nrow=length(SurvYR),ncol=length(AllStation))
DensityM_78_100<-matrix(nrow=length(SurvYR),ncol=length(AllStation))
DensityM101<-matrix(nrow=length(SurvYR),ncol=length(AllStation))
DensityMge45<-matrix(nrow=length(SurvYR),ncol=length(AllStation))

num_M_45_85<-matrix(nrow=length(SurvYR),ncol=length(AllStation))
num_M_101<-matrix(nrow=length(SurvYR),ncol=length(AllStation))

StationYr<-matrix(nrow=length(SurvYR),ncol=length(AllStation))
station_bot_temp<-matrix(nrow=length(SurvYR),ncol=length(AllStation))
tot_stations<-rep(0,length(SurvYR))

for(y in 1:length(SurvYR))
{
  yrDAT<-survDAT[drvYear==SurvYR[y],]
  fileyr<-SurvYR[y]
  stationsUNQ<-(unique(yrDAT$GIS_STATION))
  #==density at station
  for(j in 1:length(stationsUNQ))
  {
    stationALL<-yrDAT[yrDAT$GIS_STATION==stationsUNQ[j],]
    StationYr[y,j]<-as.character(stationsUNQ[j])
    Hauls<-(unique(stationALL$HAUL))
    female_45_55<-stationALL[stationALL$SEX==2 &  stationALL$WIDTH>44 & stationALL$WIDTH<56,]
    male_45_55<-stationALL[stationALL$SEX==1 &  stationALL$WIDTH>44 & stationALL$WIDTH<56,]
    male_55_65<-stationALL[stationALL$SEX==1 &  stationALL$WIDTH>55 & stationALL$WIDTH<66,]
    male_65_75<-stationALL[stationALL$SEX==1 &  stationALL$WIDTH>65 & stationALL$WIDTH<76,]
    male_45_85<-stationALL[stationALL$SEX==1 &  stationALL$WIDTH>44 & stationALL$WIDTH<86,]
    male_78_100<-stationALL[stationALL$SEX==1 &  stationALL$WIDTH>78 & stationALL$WIDTH<101,]
    Males101<-stationALL[stationALL$SEX==1& stationALL$WIDTH>101,]
    male_ge45<-stationALL[stationALL$SEX==1& stationALL$WIDTH>44,]
    
    if(nchar(stationsUNQ[j])==4 & sum(stationALL$SAMPLING_FACTOR,na.rm=T)>0)
      tot_stations[y]<-tot_stations[y]+1
    
    #==densities across hauls in crabs per km^2
    tempDensM45<-NULL
    tempDensF45<-NULL
    tempDensM55<-NULL
    tempDensM65<-NULL
    tempDensM85<-NULL
    tempDensM75101<-NULL
    tempDensM101<-NULL
    tempDensF<-NULL
    temp_num_45<-0
    temp_num_101<-0
    tempDensMge45<-NULL
    for(k in 1:length(Hauls))
    {
      SampFactM<-female_45_55$SAMPLING_FACTOR[which(female_45_55$HAUL==Hauls[k])[1]] 
      AreaSweptM<-female_45_55$AREA_SWEPT[which(female_45_55$HAUL==Hauls[k])[1]]
      tempDensF45<-length(female_45_55$HAUL==Hauls[k])*SampFactM/AreaSweptM
      
      SampFactM<-male_45_55$SAMPLING_FACTOR[which(male_45_55$HAUL==Hauls[k])[1]] 
      AreaSweptM<-male_45_55$AREA_SWEPT[which(male_45_55$HAUL==Hauls[k])[1]]
      tempDensM45<-length(male_45_55$HAUL==Hauls[k])*SampFactM/AreaSweptM
      
      SampFactM<-male_55_65$SAMPLING_FACTOR[which(male_55_65$HAUL==Hauls[k])[1]] 
      AreaSweptM<-male_55_65$AREA_SWEPT[which(male_55_65$HAUL==Hauls[k])[1]]
      tempDensM55<-length(male_55_65$HAUL==Hauls[k])*SampFactM/AreaSweptM
      
      SampFactM<-male_65_75$SAMPLING_FACTOR[which(male_65_75$HAUL==Hauls[k])[1]] 
      AreaSweptM<-male_65_75$AREA_SWEPT[which(male_65_75$HAUL==Hauls[k])[1]]
      tempDensM65<-length(male_65_75$HAUL==Hauls[k])*SampFactM/AreaSweptM
      
      SampFactM<-male_45_85$SAMPLING_FACTOR[which(male_45_85$HAUL==Hauls[k])[1]] 
      AreaSweptM<-male_45_85$AREA_SWEPT[which(male_45_85$HAUL==Hauls[k])[1]]
      tempDensM85<-length(male_45_85$HAUL==Hauls[k])*SampFactM/AreaSweptM
      temp_num_45<-temp_num_45+length(Males101$HAUL==Hauls[k])*SampFactM
      
      SampFactM<-male_78_100$SAMPLING_FACTOR[which(male_78_100$HAUL==Hauls[k])[1]] 
      AreaSweptM<-male_78_100$AREA_SWEPT[which(male_78_100$HAUL==Hauls[k])[1]]
      tempDensM78101<-length(male_78_100$HAUL==Hauls[k])*SampFactM/AreaSweptM
      
      SampFactM<-male_ge45$SAMPLING_FACTOR[which(male_ge45$HAUL==Hauls[k])[1]] 
      AreaSweptM<-male_ge45$AREA_SWEPT[which(male_ge45$HAUL==Hauls[k])[1]]
      tempDensMge45<-length(male_ge45$HAUL==Hauls[k])*SampFactM/AreaSweptM
      
      SampFactM<-Males101$SAMPLING_FACTOR[which(Males101$HAUL==Hauls[k])[1]] 
      AreaSweptM<-Males101$AREA_SWEPT[which(Males101$HAUL==Hauls[k])[1]]
      tempDensM101<-length(Males101$HAUL==Hauls[k])*SampFactM/AreaSweptM
      temp_num_101<-temp_num_101+length(Males101$HAUL==Hauls[k])*SampFactM
    }
    DensityF_45_55[y,j]<-mean(tempDensF45)
    DensityM_45_55[y,j]<-mean(tempDensM45)
    DensityM_55_65[y,j]<-mean(tempDensM55)
    DensityM_65_75[y,j]<-mean(tempDensM65)
    DensityM_45_85[y,j]<-mean(tempDensM85)
    DensityM_78_100[y,j]<-mean(tempDensM78101)
    DensityMge45[y,j]<-mean(tempDensMge45)
    DensityM101[y,j]<-mean(tempDensM101)
    num_M_45_85[y,j]<-temp_num_45    
    num_M_101[y,j]<-temp_num_101
    
    station_bot_temp[y,j]<-stationALL$GEAR_TEMPERATURE[1]
  }
}
plot(apply(num_M_101,1,sum,na.rm=T)[-1]~SurvYR[-1],type='l',ylab='total number of observed crab >101mm')

#================================
# Northern Bering See
#================================
nbs_survDAT<-read.csv("C:/data/NBS crab data/nbs_opilio_crabhaul.csv",header=T)
nbs_drvYear<-as.numeric(substr(nbs_survDAT$CRUISE,1,4))
nbs_survDAT$AKFIN_SURVEY_YEAR<-nbs_drvYear
nbs_SurvYR<-unique(nbs_drvYear)
nbs_AllStation<-unique(nbs_survDAT$GIS_STATION)

#==plot GIS stations
nbs_AllStnLoc<-matrix(ncol=2,nrow=length(nbs_AllStation))
for(w in 1:length(nbs_AllStation))
{
  temp<-nbs_survDAT[nbs_survDAT$GIS_STATION==nbs_AllStation[w],]
  nbs_AllStnLoc[w,1]<-temp$MID_LATITUDE[1]
  nbs_AllStnLoc[w,2]<-temp$MID_LONGITUDE[1]
}

nbs_num_M_45_85<-matrix(nrow=length(nbs_SurvYR),ncol=length(nbs_AllStation))
nbs_num_M_101<-matrix(nrow=length(nbs_SurvYR),ncol=length(nbs_AllStation))

nbs_DensityF_45_55<-matrix(nrow=length(nbs_SurvYR),ncol=length(nbs_AllStation))
nbs_DensityM_45_55<-matrix(nrow=length(nbs_SurvYR),ncol=length(nbs_AllStation))
nbs_DensityM_55_65<-matrix(nrow=length(nbs_SurvYR),ncol=length(nbs_AllStation))
nbs_DensityM_65_75<-matrix(nrow=length(nbs_SurvYR),ncol=length(nbs_AllStation))
nbs_DensityM_45_85<-matrix(nrow=length(nbs_SurvYR),ncol=length(nbs_AllStation))
nbs_DensityM_78_100<-matrix(nrow=length(nbs_SurvYR),ncol=length(nbs_AllStation))
nbs_DensityM101<-matrix(nrow=length(nbs_SurvYR),ncol=length(nbs_AllStation))
nbs_StationYr<-matrix(nrow=length(nbs_SurvYR),ncol=length(nbs_AllStation))
nbs_station_bot_temp<-matrix(nrow=length(nbs_SurvYR),ncol=length(nbs_AllStation))
nbs_DensityM_ge45<-matrix(nrow=length(nbs_SurvYR),ncol=length(nbs_AllStation))
for(y in 1:length(nbs_SurvYR))
{
  yrDAT<-nbs_survDAT[nbs_drvYear==nbs_SurvYR[y],]
  fileyr<-nbs_SurvYR[y]
  stationsUNQ<-(unique(yrDAT$GIS_STATION))
  #==density at station
  for(j in 1:length(stationsUNQ))
  {
    stationALL<-yrDAT[yrDAT$GIS_STATION==stationsUNQ[j],]
    nbs_StationYr[y,j]<-as.character(stationsUNQ[j])
    Hauls<-(unique(stationALL$HAUL))
    female_45_55<-stationALL[stationALL$SEX==2 &  stationALL$WIDTH>44 & stationALL$WIDTH<56,]
    male_45_55<-stationALL[stationALL$SEX==1 &  stationALL$WIDTH>44 & stationALL$WIDTH<56,]
    male_55_65<-stationALL[stationALL$SEX==1 &  stationALL$WIDTH>55 & stationALL$WIDTH<66,]
    male_65_75<-stationALL[stationALL$SEX==1 &  stationALL$WIDTH>65 & stationALL$WIDTH<76,]
    male_45_85<-stationALL[stationALL$SEX==1 &  stationALL$WIDTH>44 & stationALL$WIDTH<86,]
    male_78_100<-stationALL[stationALL$SEX==1 &  stationALL$WIDTH>78 & stationALL$WIDTH<101,]
    Males101<-stationALL[stationALL$SEX==1& stationALL$WIDTH>101,]
    male_ge45<-stationALL[stationALL$SEX==1& stationALL$WIDTH>45,]
    
    #==densities across hauls in crabs per km^2
    tempDensF45<-NULL
    tempDensM45<-NULL
    tempDensM55<-NULL
    tempDensM65<-NULL
    tempDensM85<-NULL
    tempDensM75101<-NULL
    tempDensM101<-NULL
    tempDensF<-NULL
    temp_num_45<-0
    temp_num_101<-0
    tempDensMge45<-NULL
    
    for(k in 1:length(Hauls))
    {
      SampFactM<-female_45_55$SAMPLING_FACTOR[which(female_45_55$HAUL==Hauls[k])[1]] 
      AreaSweptM<-female_45_55$AREA_SWEPT[which(female_45_55$HAUL==Hauls[k])[1]]
      tempDensF45<-length(female_45_55$HAUL==Hauls[k])*SampFactM/AreaSweptM
      
      SampFactM<-male_45_55$SAMPLING_FACTOR[which(male_45_55$HAUL==Hauls[k])[1]] 
      AreaSweptM<-male_45_55$AREA_SWEPT[which(male_45_55$HAUL==Hauls[k])[1]]
      tempDensM45<-length(male_45_55$HAUL==Hauls[k])*SampFactM/AreaSweptM
      
      SampFactM<-male_55_65$SAMPLING_FACTOR[which(male_55_65$HAUL==Hauls[k])[1]] 
      AreaSweptM<-male_55_65$AREA_SWEPT[which(male_55_65$HAUL==Hauls[k])[1]]
      tempDensM55<-length(male_55_65$HAUL==Hauls[k])*SampFactM/AreaSweptM
      
      SampFactM<-male_65_75$SAMPLING_FACTOR[which(male_65_75$HAUL==Hauls[k])[1]] 
      AreaSweptM<-male_65_75$AREA_SWEPT[which(male_65_75$HAUL==Hauls[k])[1]]
      tempDensM65<-length(male_65_75$HAUL==Hauls[k])*SampFactM/AreaSweptM
      
      SampFactM<-male_45_85$SAMPLING_FACTOR[which(male_45_85$HAUL==Hauls[k])[1]] 
      AreaSweptM<-male_45_85$AREA_SWEPT[which(male_45_85$HAUL==Hauls[k])[1]]
      tempDensM85<-length(male_45_85$HAUL==Hauls[k])*SampFactM/AreaSweptM
      temp_num_45<-temp_num_45+length(Males101$HAUL==Hauls[k])*SampFactM
      
      SampFactM<-male_78_100$SAMPLING_FACTOR[which(male_78_100$HAUL==Hauls[k])[1]] 
      AreaSweptM<-male_78_100$AREA_SWEPT[which(male_78_100$HAUL==Hauls[k])[1]]
      tempDensM78101<-length(male_78_100$HAUL==Hauls[k])*SampFactM/AreaSweptM
      
      SampFactM<-male_ge45$SAMPLING_FACTOR[which(male_ge45$HAUL==Hauls[k])[1]] 
      AreaSweptM<-male_ge45$AREA_SWEPT[which(male_ge45$HAUL==Hauls[k])[1]]
      tempDensMge45<-length(male_ge45$HAUL==Hauls[k])*SampFactM/AreaSweptM
      
      SampFactM<-Males101$SAMPLING_FACTOR[which(Males101$HAUL==Hauls[k])[1]] 
      AreaSweptM<-Males101$AREA_SWEPT[which(Males101$HAUL==Hauls[k])[1]]
      tempDensM101<-length(Males101$HAUL==Hauls[k])*SampFactM/AreaSweptM
      temp_num_101<-temp_num_101+length(Males101$HAUL==Hauls[k])*SampFactM
    }
    nbs_DensityF_45_55[y,j]<-mean(tempDensF45)
    nbs_DensityM_45_55[y,j]<-mean(tempDensM45)
    nbs_DensityM_55_65[y,j]<-mean(tempDensM55)
    nbs_DensityM_65_75[y,j]<-mean(tempDensM65)
    nbs_DensityM_45_85[y,j]<-mean(tempDensM85)
    nbs_DensityM_78_100[y,j]<-mean(tempDensM78101)
    nbs_DensityM101[y,j]<-mean(tempDensM101)
    nbs_station_bot_temp[y,j]<-stationALL$GEAR_TEMPERATURE[1]
    nbs_num_M_45_85[y,j]<-temp_num_45    
    nbs_num_M_101[y,j]<-temp_num_101
    nbs_DensityM_ge45[y,j]<-mean(tempDensMge45)
    
  }
}

#==================================
# map
#===============================
ebs_dat<-data.frame(log_abund_101=log(c(t(DensityM101[-1,]))),
                    recruits=c(t(DensityM_78_100[-nrow(DensityM_78_100),])),
                    bot_temp=c(t(station_bot_temp[-nrow(station_bot_temp),])),
                    station=c(t(StationYr[-1,])),
                    lat=rep(AllStnLoc[,1],(nrow(DensityM_78_100)-1)),
                    lon=rep(AllStnLoc[,2],(nrow(DensityM_78_100)-1)),
                    year=rep(SurvYR[-1],each=ncol(DensityM101)),
                    log_abund_101=log(c(t(DensityM101[-1,]))),
                    log_abund_45=log(c(t(DensityM_45_55[-1,]))),
                    log_abund_45_f=log(c(t(DensityF_45_55[-1,]))),
                    log_abund_55=log(c(t(DensityM_55_65[-1,]))),
                    log_abund_65=log(c(t(DensityM_65_75[-1,]))),
                    log_abund_85=log(c(t(DensityM_45_85[-1,]))),
                    log_abund_ge45=log(c(t(DensityMge45[-1,]))),
                    loc="EBS",
                    obs_num_101=c(t(num_M_101[-1,])),
                    obs_num_45_85=c(t(num_M_45_85[-1,])))

for(x in 1:length(AllStation))
{
  ebs_dat$lat[ebs_dat$station==AllStation[x]]<-AllStnLoc[x,1]
  ebs_dat$lon[ebs_dat$station==AllStation[x]]<-AllStnLoc[x,2]
}

nbs_dat<-data.frame(log_abund_101=log(c(t(nbs_DensityM101))),
                    recruits=c(t(nbs_DensityM_78_100)),
                    bot_temp=c(t(nbs_station_bot_temp)),
                    station=c(t(nbs_StationYr)),
                    lat=rep(nbs_AllStnLoc[,1],(nrow(nbs_DensityM_78_100))),
                    lon=rep(nbs_AllStnLoc[,2],(nrow(nbs_DensityM_78_100))),
                    year=rep(nbs_SurvYR,each=ncol(nbs_DensityM101)),
                    log_abund_101=log(c(t(nbs_DensityM101))),
                    log_abund_45=log(c(t(nbs_DensityM_45_55))),
                    log_abund_45_f=log(c(t(nbs_DensityF_45_55))),
                    log_abund_55=log(c(t(nbs_DensityM_55_65))),
                    log_abund_65=log(c(t(nbs_DensityM_65_75))),
                    log_abund_85=log(c(t(nbs_DensityM_45_85))),
                    log_abund_ge45=log(c(t(nbs_DensityM_ge45))),
                    loc='NBS',
                    obs_num_101=c(t(nbs_num_M_101)),
                    obs_num_45_85=c(t(nbs_num_M_45_85)))

for(x in 1:length(nbs_AllStation))
{
  nbs_dat$lat[nbs_dat$station==nbs_AllStation[x]]<-nbs_AllStnLoc[x,1]
  nbs_dat$lon[nbs_dat$station==nbs_AllStation[x]]<-nbs_AllStnLoc[x,2]
}

#=========================
# plot highest year vs. lowest year
#=========================

library(dplyr)
in_dat<-rbind(filter(ebs_dat,year==2018 | year==2021),filter(nbs_dat,year==2018 | year==2021))
world <- ne_countries(scale = "medium", returnclass = "sf")
lon_1<-min(in_dat$lon,na.rm=T)
lon_2<-max(in_dat$lon,na.rm=T)*.99
lat_1<-min(in_dat$lat,na.rm=T)
lat_2<-max(in_dat$lat,na.rm=T)
in_gam_dat<-in_dat
in_gam_dat_2<-in_gam_dat[,c(5,6,7,14)]
in_fin<-melt(in_gam_dat_2,id=c('lat','lon','year'))
in_fin<-in_fin[complete.cases(in_fin),]

#==plot maps smaller crab
nbs_bound<-data.frame(lat=c(62.1,62.1,61.1,61.1,60.5,60.5),
                      lon=c(-177,-172.5,-172.5,-171,-171,-165))


op_map<-ggplot() + 
  geom_tile(data=in_fin, aes(x = lon, y = lat, fill = value),width=.5,height=.25) +
  #scale_fill_distiller(palette="Spectral", na.value="grey") +
  scale_fill_distiller(palette="RdYlBu", na.value="grey") +
  facet_wrap(~year,ncol=1) +
  geom_sf(data=world) +
  coord_sf(xlim = c(lon_1,lon_2), ylim = c(lat_1,lat_2), expand = FALSE) +
  theme_bw()+
  theme(axis.text.x = element_text(size = 10))+
  geom_line(data=nbs_bound,aes(x=lon,y=lat),col='red')+
  theme(legend.position=c(.88,.82),
        legend.background = element_rect(fill='transparent',color=NA),
        legend.box.background = element_rect(fill='transparent',color=NA))+
  labs(fill="ln(crab/nm^2)")

#########################################
# map for figure 3; shows density in 2018 vs. average
######################################
plot_st<-data.frame(Year=SurvYR,Stations=tot_stations)
stations_plot<-ggplot(data=filter(plot_st,Year>1988),aes(x=Year,y=Stations))+
  geom_line(data=filter(plot_st,Year>1988),aes(x=Year,y=Stations),lty=2)+
  theme_bw()+expand_limits(y=0)+
  geom_point(data=filter(plot_st,Year>1988),aes(x=Year,y=Stations))+
  geom_line(data=filter(plot_st,Year>1988&Year<2020),aes(x=Year,y=Stations))

in_dat<-rbind(filter(ebs_dat,year==2018),filter(nbs_dat,year==2018 ))
in_dat2<-rbind(filter(ebs_dat,year==1995),filter(nbs_dat,year==1995 ))
in_dat<-rbind(filter(ebs_dat,year==2018&nchar(station)==4),filter(nbs_dat,year==2018 &nchar(station)==4))
in_dat2<-rbind(filter(ebs_dat,year==1995&nchar(station)==4),filter(nbs_dat,year==1995&nchar(station)==4))
in_dat3<-rbind(filter(ebs_dat,year==1999&nchar(station)==4),filter(nbs_dat,year==1999&nchar(station)==4))
in_gam_dat<-in_dat[,c(5,6,7,14)]
in_gam_dat_2<-in_dat2[,c(5,6,7,14)]
in_gam_dat_3<-in_dat3[,c(5,6,7,14)]
in_fin<-melt(in_gam_dat,id=c('lat','lon','year'))
in_fin<-in_fin[complete.cases(in_fin),]
in_fin2<-melt(in_gam_dat_2,id=c('lat','lon','year'))
in_fin2<-in_fin2[complete.cases(in_fin2),]
in_fin3<-melt(in_gam_dat_3,id=c('lat','lon','year'))
in_fin3<-in_fin3[complete.cases(in_fin3),]

unique(nchar(ebs_dat$station))

dens_map<-ggplot() + 
  geom_tile(data=in_fin3, aes(x = lon, y = lat),width=.5,height=.25,fill='light grey') +
  geom_tile(data=in_fin, aes(x = lon, y = lat, fill = value),width=.5,height=.25) +
  scale_fill_distiller(palette="RdYlBu", na.value="grey") +
  geom_sf(data=world) +
  coord_sf(xlim = c(lon_1,lon_2), ylim = c(lat_1,lat_2), expand = FALSE) +
  theme_bw()+
  theme(axis.text.x = element_text(size = 10))+
  geom_line(data=nbs_bound,aes(x=lon,y=lat),col='red')+
  theme(legend.position=c(.8,.725),
        legend.background = element_rect(fill='transparent',color=NA),
        legend.box.background = element_rect(fill='transparent',color=NA))+labs(fill="ln(crab/nm^2)")


tmp<-stations_plot/dens_map + plot_layout(heights=c(.3,1)) + plot_annotation(tag_levels='a')

png('plots/figure_1_addon.png',height=10,width=7,res=350,units='in')
print(tmp)
dev.off()

###########################
# just survey stations
#########################

in_dat<-data.frame(rbind(nbs_AllStnLoc,AllStnLoc))
colnames(in_dat)<-c("Latitude","Longitude")

world <- ne_countries(scale = "medium", returnclass = "sf")
lon_1<-min(in_dat$Longitude,na.rm=T)
lon_2<-max(in_dat$Longitude,na.rm=T)*.99
lat_1<-min(in_dat$Latitude,na.rm=T)
lat_2<-max(in_dat$Latitude,na.rm=T)
#==plot maps smaller crab
nbs_bound<-data.frame(lat=c(62.1,62.1,61.1,61.1,60.5,60.5),
                      lon=c(-177,-172.5,-172.5,-171,-171,-165))


####################################
# other figures
####################################
annotation_custom2 <-   function (grob, xmin = -Inf, xmax = Inf, ymin = -Inf,ymax = Inf, data){ layer(data = data, 
                                                                                                      stat = StatIdentity, 
                                                                                                      position = PositionIdentity,
                                                                                                      geom = ggplot2:::GeomCustomAnn,
                                                                                                      inherit.aes = TRUE, 
                                                                                                      params = list(grob = grob,xmin = xmin, xmax = xmax,
                                                                                                                    ymin = ymin, ymax = ymax))}
opie_ts<-read.csv("data/EBSCrab_AB_Sizegroup.csv")
in_png_opie<-readPNG('data/snowcrab.png')

in_snow_ts<-filter(opie_ts,SIZE_GROUP=="MALE_TOTAL")
in_snow_ts$up_ci<-in_snow_ts$ABUNDANCE * exp(1.96*sqrt(log(1+in_snow_ts$ABUNDANCE_CV)))
in_snow_ts$dn_ci<-in_snow_ts$ABUNDANCE / exp(1.96*sqrt(log(1+in_snow_ts$ABUNDANCE_CV)))
colnames(in_snow_ts)[1]<-"Year"
div_n<-1000000000
time_snow<-ggplot(in_snow_ts)+
  geom_ribbon(aes(x=Year,ymin=dn_ci/div_n,ymax=up_ci/div_n),fill='light grey')+
  geom_line(aes(x=Year,y=ABUNDANCE/div_n))+
  geom_point(aes(x=Year,y=ABUNDANCE/div_n))+
  theme_bw()+
  expand_limits(y=0)+
   annotation_custom2(rasterGrob(in_png_opie, interpolate=TRUE), 
                      xmin=1995, xmax=2016, ymin=5, ymax=13,data=in_snow_ts[1:1280,])+
  xlab("")+
  scale_y_continuous(name="Abundance 
(billions)",position='right')+
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.title=element_text(size=8))


ice_extent<-read.csv("data/IceExtent.csv",skip=4)
cold_pool<-read.csv("data/ColdPool.csv",skip=4)
ice_extent$Index[which(ice_extent$Index=="null")]<-NA
cold_pool$Value[which(cold_pool$Value=="null")]<-NA

ice_plot<-data.frame(Year=rep(seq(1982,2021),2),
                     value=as.numeric(c(cold_pool$Value[which(cold_pool$Year>1981)],ice_extent$Index[which(ice_extent$Year>1981)])),
                     variable=c(rep("Cold pool",length(seq(1982,2021))),rep("Ice extent",length(seq(1982,2021)))))

ice_bckgrnd<-ggplot(ice_plot)+
  geom_line(aes(x=Year,y=value/10000,col=variable,group=variable),lwd=1.2)+
  geom_point(aes(x=Year,y=value/10000,col=variable,group=variable))+
  theme_bw()+  scale_colour_manual(values=c("dodgerblue","royalblue3"))+
  expand_limits(y=0)+
  theme(legend.position=c(.41,1),
        legend.title = element_blank(),
        legend.background = element_rect(fill = "transparent"),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.title=element_text(size=8),
        legend.direction='horizontal')+
  scale_y_continuous(name='Area (10,000 km^2)',position='right')

#========================================
# numbers at size figure
#========================================
kod_dat<-read.csv("data/EBSCrab_Abundance_Biomass_1mm.csv",skip=7)
kod_dat_1<-filter(kod_dat,SEX=='MALE')
kod_dat_m<-kod_dat_1 %>%
  group_by(SURVEY_YEAR,SIZE_CLASS_MM) %>%
  summarize(abund=sum(ABUNDANCE))

#==to make size composition
in_break<-seq(0, max(kod_dat_m$SIZE_CLASS_MM,na.rm=T), 5)
mid_pts<-rep(0,length(in_break)-1)
for(x in 1:length(mid_pts))
  mid_pts[x] <- (in_break[x]+in_break[x+1])/2

temp_kod<-kod_dat_m %>%
  group_by(SURVEY_YEAR ,group = cut(SIZE_CLASS_MM , breaks =in_break)) %>%
  summarise(numbers = sum(abund,na.rm=T))
temp_kod$mid_pt<-mid_pts[temp_kod$group]

# all data
size_yr <- ggplot(dat=kod_dat_m) 
size_yr <- size_yr + geom_density_ridges(aes(x=SIZE_CLASS_MM, y=SURVEY_YEAR, height = abund,
                                 group = SURVEY_YEAR, 
                                 alpha=.9999),stat = "identity",scale=5,fill='royalblue3') +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(legend.position = "none") +
  labs(x="Carapace width (mm)") +
  xlim(25,125)+
  scale_y_continuous(name='Year',position='right')+
  geom_vline(aes(xintercept=101),lty=2)

#========================================
# Temperature occupied
#========================================
in_break<-seq(0, max(survDAT$WIDTH,na.rm=T), 5)
mid_pts<-rep(0,length(in_break)-1)
for(x in 1:length(mid_pts))
  mid_pts[x] <- (in_break[x]+in_break[x+1])/2

#==to make size composition
temp_tot<-survDAT %>%
  group_by(AKFIN_SURVEY_YEAR ,group = cut(WIDTH, breaks =in_break),SEX) %>%
  summarise(temp_occ = weighted.mean(x=GEAR_TEMPERATURE,w=SAMPLING_FACTOR,na.rm=T))
temp_tot$mid_pt<-mid_pts[temp_tot$group]

#==temperature occupied by size class
in_dat<-filter(temp_tot,mid_pt==42.5 | mid_pt==102.5,SEX<2)
in_dat$mid_pt[in_dat$mid_pt==42.5]<-"42.5 mm"
in_dat$mid_pt[in_dat$mid_pt==102.5]<-"102.5 mm"
temp_occ<-ggplot(in_dat)+
  geom_line(aes(x=AKFIN_SURVEY_YEAR,y=temp_occ,col=as.factor(mid_pt),group=mid_pt),
            lwd=1)+
  geom_point(aes(x=AKFIN_SURVEY_YEAR,y=temp_occ,col=as.factor(mid_pt),group=mid_pt),
            lwd=1)+
  theme_bw()+
  scale_y_continuous(name='Temperature 
  occupied (C)',position='right')+
  scale_colour_manual(values=c("honeydew3","honeydew4"))+
  theme(legend.position=c(.45,.9),
        legend.key.height=unit(.15,'cm'),
        legend.title=element_blank(),
        legend.background = element_rect(fill = "transparent"),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.title=element_text(size=8),legend.direction='horizontal')+
  geom_hline(yintercept=0,lty=2)+xlim(1982,2021)+
  xlab("")

  patchwork <- ice_bckgrnd / plot_spacer() /  time_snow 
  patchwork[[1]] = patchwork[[1]] + theme(axis.text.x = element_blank(),
                                          axis.ticks.x = element_blank(),
                                          axis.title.x = element_blank() )
  tmp2<-(patchwork/plot_spacer()/size_yr) + plot_layout(heights=c(.3,-.11,.3,-.06,1))
  tmp2<-tmp2+plot_annotation(tag_levels='A')
  
  tmp3<-op_map + tmp2+plot_layout(widths=c(1,.5))
  
  png('plots/figure_1.png',height=10,width=8,res=350,units='in')
  tmp3+plot_annotation(tag_levels='A')
  dev.off()
  
  jpeg('plots/figure_1.jpeg',height=10,width=8,res=350,units='in')
  tmp3+plot_annotation(tag_levels='A')
  dev.off() 
  
  pdf('plots/figure_1.pdf',height=10,width=8)
  tmp3+plot_annotation(tag_levels='A')
  dev.off()