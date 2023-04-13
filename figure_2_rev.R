#===================================================================
#==ANALYSIS OF THE COLLAPSE OF SNOW CRAB IN THE EASTERN BERING SEA==
#===================================================================
#==Gameplan===============================
#==READ IN POP DY MODEL OUTPUT
#==CHECK POP DY MODEL FITS
#==READ IN COVARIATES
#==MAKE COVARIATE DATA FRAME
#==BUILD GAM TO EXPLAIN ESTIMATED M AND Q
#==DO ROBUSTNESS CHECKS OF GAM
#=========================================
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
library(GGally)
library(itsadug)
library(gratia)
library(scam)

annotation_custom2 <-   function (grob, xmin = -Inf, xmax = Inf, ymin = -Inf,ymax = Inf, data){ layer(data = data, 
                                     stat = StatIdentity, 
                                     position = PositionIdentity,
                                     geom = ggplot2:::GeomCustomAnn,
                                     inherit.aes = TRUE, 
                                     params = list(grob = grob,xmin = xmin, xmax = xmax,
                                                   ymin = ymin, ymax = ymax))}


#==============================
# READ IN POP DY MODEL OUTPUT
#==============================
#==========
# Input
#===============
input<-"C:/Users/cody.szuwalski/Work/snow_down_pub/models/model_base/snow_down.DAT"
styr <- scan(input,skip=2,n=1,quiet=T)
endyr <- scan(input,skip=4,n=1,quiet=T)
year_n <- scan(input,skip=6,n=1,quiet=T)
years <- scan(input,skip=8,n=year_n,quiet=T)
size_n <- scan(input,skip=10,n=1,quiet=T)
sizes <- scan(input,skip=12,n=size_n,quiet=T)

imm_n_at_size_obs <- matrix(scan(input,skip=18,n=year_n*size_n,quiet=T),ncol=size_n,byrow=T)
mat_n_at_size_obs <- matrix(scan(input,skip=52,n=year_n*size_n,quiet=T),ncol=size_n,byrow=T)

prob_term_molt <- matrix(scan(input,skip=86,n=year_n*size_n,quiet=T),ncol=size_n,byrow=T)
size_trans <- (matrix(scan(input,skip=120,n=size_n*size_n,quiet=T),ncol=size_n,byrow=T))

survey_select<-scan(input,skip=142,n=size_n,quiet=T)

imm_cv <- scan(input,skip=134,n=year_n,quiet=T)
mat_cv <- scan(input,skip=136,n=year_n,quiet=T)

#===============
# Output
#===============
rep_files<-c("models/model_base/snow_down.rep",
             "models/model_vary_m/snow_down.rep",
             "models/model_vary_q/snow_down.rep",
             "models/model_vary_q_m/snow_down.rep")
take_ind<-2
outs<-list(list())

for(x in 1:length(rep_files))
  outs[[x]]<-readList(rep_files[x])

obs_imm_comp<-imm_n_at_size_obs
obs_mat_comp<-mat_n_at_size_obs
for(x in 1:nrow(imm_n_at_size_obs))
{
  obs_imm_comp[x,]<-imm_n_at_size_obs[x,]/sum(imm_n_at_size_obs[x,])
  obs_mat_comp[x,]<-mat_n_at_size_obs[x,]/sum(mat_n_at_size_obs[x,]) 
}

out_like<-NULL
for(x in 1:length(outs))
 out_like<-c(out_like,sum(outs[[x]]$likelihoods))

used_par<-c(94,160,162,230)

AIC<-2*used_par+2*out_like

#===plot outputs building complexity====
plot_build_up<-data.frame(values=c(outs[[1]]$`natural mortality`[,1],
                    outs[[1]]$`mature natural mortality`[,1],
                    outs[[1]]$`survey selectivity`[,1],
                    outs[[1]]$`mature survey selectivity`[,1]),
           process=c(rep("Immature mortality",length(outs[[1]]$`natural mortality`[,1])),
                     rep("Mature mortality",length(outs[[1]]$`natural mortality`[,1])),
                     rep("Immature catchability",length(outs[[1]]$`natural mortality`[,1])),
                     rep("Mature catchability",length(outs[[1]]$`natural mortality`[,1]))),
           model="Time-invariant",year=seq(1989,2021))
in_mod<-c("Invariant","Vary M","Vary q","Vary both")
for(x in 2:4)
{  temp<-data.frame(values=c(outs[[x]]$`natural mortality`[,1],
                                     outs[[x]]$`mature natural mortality`[,1],
                                     outs[[x]]$`survey selectivity`[,1],
                                     outs[[x]]$`mature survey selectivity`[,1]),
                            process=c(rep("Immature mortality",length(outs[[x]]$`natural mortality`[,1])),
                                      rep("Mature mortality",length(outs[[x]]$`natural mortality`[,1])),
                                      rep("Immature catchability",length(outs[[x]]$`natural mortality`[,1])),
                                      rep("Mature catchability",length(outs[[x]]$`natural mortality`[,1]))),
                            model=in_mod[x],year=seq(1989,2021))
plot_build_up<-rbind(plot_build_up,temp)
} 

png("plots/mod_comp.png",height=8,width=8,res=400,units='in')
ggplot(filter(plot_build_up,year<2020))+
  geom_line(aes(x=year,y=values,col=model),lwd=2)+
  facet_wrap(~process,scales='free_y')+
  theme_bw()
dev.off()

#==uncertainty in estimates of imm and mat N
ret <- list()
parfile <- as.numeric(scan("models/model_vary_m/snow_down.par", what = '', n=16, quiet=TRUE)[c(6,11,16)])
ret$nopar <- as.integer(parfile[1])
ret$nlogl <- parfile[2]
ret$maxgrad <- parfile[3]
file <- "models/model_vary_m/snow_down.cor"
if (file.exists(file))
{
  lin <- readLines(file)
  ret$npar <- length(lin)-2
  ret$logDetHess <- as.numeric(strsplit(lin[1], '=')[[1]][2])
  sublin <- lapply(strsplit(lin[1:ret$npar+2], ' '),function(x)x[x!=''])
  ret$names <- unlist(lapply(sublin,function(x)x[2]))
  ret$est <- as.numeric(unlist(lapply(sublin,function(x)x[3])))
  ret$std <- as.numeric(unlist(lapply(sublin,function(x)x[4])))
  ret$cor <- matrix(NA, ret$npar, ret$npar)
  corvec <- unlist(sapply(1:length(sublin), function(i)sublin[[i]][5:(4+i)]))
  ret$cor[upper.tri(ret$cor, diag=TRUE)] <- as.numeric(corvec)
  ret$cor[lower.tri(ret$cor)]  <-  t(ret$cor)[lower.tri(ret$cor)]
  ret$cov <- ret$cor*(ret$std%o%ret$std)
}

imm_sd<-ret$std[which(ret$names=="imm_n_sd")]
mat_sd<-ret$std[which(ret$names=="mat_n_sd")]


#==find r square
div_n<-1000000000
in_ind<-1
plot_mmb<-outs[[in_ind]]$imm_n_obs/div_n
plot_mmb[plot_mmb==0]<-NA

preds<-outs[[in_ind]]$imm_numbers_pred/div_n
mod<-gam(plot_mmb~(preds))
mod<-lm(plot_mmb~(preds))
summary(mod)

plot_mmb<-outs[[in_ind]]$mat_n_obs/div_n
plot_mmb[plot_mmb==0]<-NA

preds<-outs[[in_ind]]$mat_numbers_pred/div_n
mod<-gam(plot_mmb~(preds))
summary(mod)



#==immature size comps
rownames(obs_imm_comp)<-years
colnames(obs_imm_comp)<-sizes
df_1<-melt(obs_imm_comp)
colnames(df_1)<-c("Year","Size","Proportion")
df_1$quant<-'Observed'

df2<-data.frame(pred=apply(outs[[take_ind]]$'immature numbers at size',2,median),
                Size=(sizes))

allLevels <- levels(factor(c(df_1$Size,df2$Size)))
df_1$Size <- factor(df_1$Size,levels=(allLevels))
df2$Size <- factor(df2$Size,levels=(allLevels))

imm_size<-ggplot(data=df2,aes(x = factor(Size), y =pred)) + 
  geom_boxplot(data=df_1,aes(x = Size, y =Proportion ),fill='grey') +
  stat_summary(fun.y=mean, geom="line", aes(group=1),lwd=1.5,col='blue',alpha=.8)  + 
  ylab("Size composition")+theme(axis.title=element_text(size=11))+
  theme_void()

#==mature size comps
rownames(obs_mat_comp)<-years
colnames(obs_mat_comp)<-sizes
df_1<-melt(obs_mat_comp)
colnames(df_1)<-c("Year","Size","Proportion")
df_1$quant<-'Observed'
df2<-data.frame(pred=apply(outs[[take_ind]]$'mature numbers at size',2,median),
                Size=(sizes))

allLevels <- levels(factor(c(df_1$Size,df2$Size)))
df_1$Size <- factor(df_1$Size,levels=(allLevels))
df2$Size <- factor(df2$Size,levels=(allLevels))


mat_size<-ggplot(data=df2,aes(x = factor(Size), y =pred)) + 
  geom_boxplot(data=df_1,aes(x = Size, y =Proportion ),fill='grey') +
  stat_summary(fun=mean, geom="line", aes(group=1),lwd=1.5,col='blue',alpha=.8)  + 
  ylab("")+
  xlab("Carapace width (mm)")+
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.background=element_blank(),
        panel.grid=element_blank(),
        axis.line.x=element_line(colour="black"))

#==fits index
# imm
df_1<-data.frame(pred=outs[[take_ind]]$imm_numbers_pred/div_n,
                 obs=outs[[1]]$imm_n_obs/div_n,
                 year=years,
                 ci_dn=(outs[[1]]$imm_n_obs/div_n) /  exp(1.96*sqrt(log(1+imm_cv^2))),
                 ci_up=(outs[[1]]$imm_n_obs/div_n) *  exp(1.96*sqrt(log(1+imm_cv^2))),
                 est_ci_up=outs[[take_ind]]$imm_numbers_pred/div_n + 1.96*imm_sd/div_n,
                 est_ci_dn=outs[[take_ind]]$imm_numbers_pred/div_n - 1.96*imm_sd/div_n,
                 recruits=(outs[[take_ind]]$recruits)/div_n)

df_1<-df_1[-32,]

imm_abnd<-ggplot(data=df_1)+
  geom_segment(aes(x=year,xend=year,y=ci_dn,yend=ci_up))+
  geom_point(aes(x=year,y=obs))+
  geom_line(aes(x=year,y=pred),lwd=1,col='blue',alpha=.8)+
  geom_ribbon(aes(x=year,ymin=est_ci_dn,ymax=est_ci_up),fill='blue',alpha=0.2)+
  #geom_line(aes(x=year,y=recruits),col='purple',lwd=.7,lty=1,alpha=0.8)+
  scale_y_continuous(breaks=c(0,2.5,5.0),limits=c(-0.1,6))+
  ylab("IMMATURE
Abundance (billions)")+
  theme(axis.line.y = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),)

# mature
df_1<-data.frame(pred=outs[[take_ind]]$mat_numbers_pred/div_n,
                 obs=outs[[1]]$mat_n_obs/div_n,
                 year=years,
                 ci_dn=(outs[[1]]$mat_n_obs/div_n) /  exp(1.96*sqrt(log(1+mat_cv^2))),
                 ci_up=(outs[[1]]$mat_n_obs/div_n) *  exp(1.96*sqrt(log(1+mat_cv^2))),
                 est_ci_up=outs[[take_ind]]$mat_numbers_pred/div_n + 1.96*mat_sd/div_n,
                 est_ci_dn=outs[[take_ind]]$mat_numbers_pred/div_n - 1.96*mat_sd/div_n)
df_1<-df_1[-32,]

mat_abnd<-ggplot(data=df_1)+
  geom_segment(aes(x=year,xend=year,y=ci_dn,yend=ci_up))+
  geom_point(aes(x=year,y=obs))+
  geom_line(aes(x=year,y=pred),lwd=1,col='blue',alpha=.8)+
  geom_ribbon(aes(x=year,ymin=est_ci_dn,ymax=est_ci_up),fill='blue',alpha=0.2)+
  theme_bw()+
  scale_y_continuous(position = "left",limits=c(0,2.1))+
  ylab("MATURE
Abundance (billions)")+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 


#==============================
# READ IN COVARIATES
#==============================
survDAT<-read.csv("data/EBSCrab_Haul/EBSCrab_Haul.csv",header=T,skip=5)
drvYear<-as.numeric(substr(survDAT$CRUISE,1,4))
SurvYR<-unique(drvYear)

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

unq_ebs_station<-unique(AllStation)


#========================
#==set up for map making
#========================
nmiSurv<-140350
DensityM_45_55<-matrix(nrow=length(SurvYR),ncol=length(AllStation))
DensityF_45_55<-matrix(nrow=length(SurvYR),ncol=length(AllStation))
DensityM_55_65<-matrix(nrow=length(SurvYR),ncol=length(AllStation))
DensityM_65_75<-matrix(nrow=length(SurvYR),ncol=length(AllStation))
DensityM_45_85<-matrix(nrow=length(SurvYR),ncol=length(AllStation))
DensityM_78_100<-matrix(nrow=length(SurvYR),ncol=length(AllStation))
DensityM101<-matrix(nrow=length(SurvYR),ncol=length(AllStation))

bcd_density<-matrix(nrow=length(SurvYR),ncol=length(AllStation))

num_M_45_85<-matrix(nrow=length(SurvYR),ncol=length(AllStation))
num_M_101<-matrix(nrow=length(SurvYR),ncol=length(AllStation))

StationYr<-matrix(nrow=length(SurvYR),ncol=length(AllStation))
station_bot_temp<-matrix(nrow=length(SurvYR),ncol=length(AllStation))
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
    
    DensityM101[y,j]<-mean(tempDensM101)
    num_M_45_85[y,j]<-temp_num_45    
    num_M_101[y,j]<-temp_num_101
    
    station_bot_temp[y,j]<-stationALL$GEAR_TEMPERATURE[1]
    bcd_density[y,j]<-sum(stationALL$SAMPLING_FACTOR[which(stationALL$DISEASE_CODE==2)])/sum(stationALL$SAMPLING_FACTOR)
  }
}


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
                    loc="EBS",
                    obs_num_101=c(t(num_M_101[-1,])),
                    obs_num_45_85=c(t(num_M_45_85[-1,])),
                    bitter_crab=(c(t(bcd_density[-1,]))))

for(x in 1:length(AllStation))
{
  ebs_dat$lat[ebs_dat$station==AllStation[x]]<-AllStnLoc[x,1]
  ebs_dat$lon[ebs_dat$station==AllStation[x]]<-AllStnLoc[x,2]
}




size_covar_1<-30
size_covar_2<-95

#==make another survDAT that makes uses a cutoff for mature and immature
pred_n_size_imm<-sweep(outs[[take_ind]]$'immature numbers at size', 1, outs[[take_ind]]$imm_numbers_pred, FUN="*")
rownames(pred_n_size_imm)<-years
colnames(pred_n_size_imm)<-sizes
pred_n_size_mat<-sweep(outs[[take_ind]]$'mature numbers at size', 1, outs[[take_ind]]$mat_numbers_pred, FUN="*")
rownames(pred_n_size_mat)<-years
colnames(pred_n_size_mat)<-sizes


ugh<-pred_n_size_imm/(pred_n_size_imm+pred_n_size_mat)
mid_size<-rep(0,nrow(ugh))
for(x in 1:nrow(ugh))
  mid_size[x]<-as.numeric(names(which(abs(ugh[x,]-0.5)==min(abs(ugh[x,]-0.5)))))

larp<-data.frame(Size=mid_size[-length(mid_size)],Year=c(seq(1988,2019)))

png("plots/mat_midpt.png",height=5,width=8,res=400,units='in')
ggplot(data=larp,aes(x=Year,y=Size))+
  geom_line()+theme_bw()+geom_point()
dev.off()

survDAT_mat<-NULL
survDAT_imm<-NULL
for(x in 1:length(years))
{
  if(years[x]>1988)
  {
    temp_d<-survDAT[survDAT$AKFIN_SURVEY_YEAR==years[x],]
    survDAT_mat<-rbind(survDAT_mat,temp_d[nchar(temp_d$GIS_STATION)==4 & temp_d$WIDTH>mid_size[x],])
    survDAT_imm<-rbind(survDAT_imm,temp_d[nchar(temp_d$GIS_STATION)==4 & temp_d$WIDTH<=mid_size[x],])
  }
}


#========================================
# Disease
#========================================
#  median disease prevalence

#========================================
library(dplyr)
in_break<-seq(0, max(survDAT$WIDTH,na.rm=T), 5)
mid_pts<-rep(0,length(in_break)-1)
for(x in 1:length(mid_pts))
  mid_pts[x] <- (in_break[x]+in_break[x+1])/2

#==mature disease
tmp_tot_mat<-survDAT_mat %>%
  group_by(AKFIN_SURVEY_YEAR ,group = cut(WIDTH, breaks =in_break),SEX) %>%
  summarise(bin_n1 = sum(SAMPLING_FACTOR ))
tmp_tot_mat$mid_pt<-mid_pts[tmp_tot_mat$group]

tmp_dis_mat<-survDAT_mat %>%
  group_by(AKFIN_SURVEY_YEAR ,group = cut(WIDTH, breaks =in_break),SEX) %>%
  summarise(bin_n1 = sum(DISEASE_CODE==2 * SAMPLING_FACTOR,na.rm=T ))
tmp_dis_mat$mid_pt<-mid_pts[tmp_dis_mat$group]

tot_dis_mat<-tmp_tot_mat
tot_dis_mat$disease<-tmp_dis_mat$bin_n1
tot_dis_mat$prop<-tot_dis_mat$disease/tot_dis_mat$bin_n1
tot_dis_mat$mid_pt<-mid_pts[tot_dis_mat$group]

use_disease_mat<-filter(tot_dis_mat,SEX==1 & AKFIN_SURVEY_YEAR>1988 & mid_pt>size_covar_1 & mid_pt<size_covar_2 & bin_n1>1)%>%
  group_by(AKFIN_SURVEY_YEAR)%>%
  summarize(med_dis=median(prop))

#==immature disease
tmp_tot_imm<-survDAT_imm %>%
  group_by(AKFIN_SURVEY_YEAR ,group = cut(WIDTH, breaks =in_break),SEX) %>%
  summarise(bin_n1 = sum(SAMPLING_FACTOR ))
tmp_tot_imm$mid_pt<-mid_pts[tmp_tot_imm$group]

tmp_dis_imm<-survDAT_imm %>%
  group_by(AKFIN_SURVEY_YEAR ,group = cut(WIDTH, breaks =in_break),SEX) %>%
  summarise(bin_n1 = sum(DISEASE_CODE==2 * SAMPLING_FACTOR,na.rm=T ))
tmp_dis_imm$mid_pt<-mid_pts[tmp_dis_imm$group]

tot_dis_imm<-tmp_tot_imm
tot_dis_imm$disease<-tmp_dis_imm$bin_n1
tot_dis_imm$prop<-tot_dis_imm$disease/tot_dis_imm$bin_n1
tot_dis_imm$mid_pt<-mid_pts[tot_dis_imm$group]

use_disease_imm<-filter(tot_dis_imm,SEX==1 & AKFIN_SURVEY_YEAR>1988 & mid_pt>size_covar_1 & mid_pt<size_covar_2 & bin_n1>1)%>%
  group_by(AKFIN_SURVEY_YEAR)%>%
  summarize(med_dis=median(prop))

#========================================
# Temperature occupied
#========================================
in_break<-seq(0, max(survDAT$WIDTH,na.rm=T), 5)
mid_pts<-rep(0,length(in_break)-1)
for(x in 1:length(mid_pts))
  mid_pts[x] <- (in_break[x]+in_break[x+1])/2

#==temperature index by 'maturity state'
temp_tot_ind_imm<-filter(survDAT_imm,  SEX<2) %>%
  group_by(AKFIN_SURVEY_YEAR ,SEX) %>%
  summarise(temp_occ = weighted.mean(x=GEAR_TEMPERATURE,w=SAMPLING_FACTOR,na.rm=T))
temp_tot_ind_imm$Maturity<-"Immature"
temp_tot_ind_mat<-filter(survDAT_mat, SEX<2) %>%
  group_by(AKFIN_SURVEY_YEAR ,SEX) %>%
  summarise(temp_occ = weighted.mean(x=GEAR_TEMPERATURE,w=SAMPLING_FACTOR,na.rm=T))
temp_tot_ind_mat$Maturity<-"Mature"

#==to make size composition
temp_tot<-survDAT %>%
  group_by(AKFIN_SURVEY_YEAR ,group = cut(WIDTH, breaks =in_break),SEX) %>%
  summarise(temp_occ = weighted.mean(x=GEAR_TEMPERATURE,w=SAMPLING_FACTOR,na.rm=T))
temp_tot$mid_pt<-mid_pts[temp_tot$group]

wide_tmp_occ<-dcast(filter(temp_tot,SEX==1), AKFIN_SURVEY_YEAR~mid_pt,value.var='temp_occ')
rownames(wide_tmp_occ)<-wide_tmp_occ[,1]
wide_tmp_occ<-wide_tmp_occ[,which(as.numeric(colnames(wide_tmp_occ))>size_covar_1&as.numeric(colnames(wide_tmp_occ))<size_covar_2)]
wide_tmp_occ<-wide_tmp_occ[rownames(wide_tmp_occ)>1988,]

wtd_avg_temp_occ_mat<-rep(0,nrow(wide_tmp_occ))
wtd_avg_temp_occ_imm<-rep(0,nrow(wide_tmp_occ))
for(x in 1:(nrow(wide_tmp_occ)))
{
  wtd_avg_temp_occ_imm[x]<-weighted.mean(x=wide_tmp_occ[x,],w=pred_n_size_imm[x,]) 
  wtd_avg_temp_occ_mat[x]<-weighted.mean(x=wide_tmp_occ[x,],w=pred_n_size_mat[x,])   
}
names(wtd_avg_temp_occ_imm)<-rownames(wide_tmp_occ)
names(wtd_avg_temp_occ_mat)<-rownames(wide_tmp_occ)
#=============================
# cannibalism index
#=============================
cannibals<-read.csv('data/cannibal_index.csv',header=T)

cannibal_ind_imm<-filter(cannibals,X2<52)%>%
  group_by(X1)%>%
  summarize(med_cannib=mean(value))
colnames(cannibal_ind_imm)[1]<-"year"


#============================
# fishery data (numbers at size--used in previous versions)
#============================
fish_dat<-readList("data/fisheries_data.txt")
all_bycatch<-apply(fish_dat$trawl_num_size_obs,1,sum)

ret_catch<-fish_dat$retained_num_size_obs
colnames(ret_catch)<-seq(27.5,132.5,5)
rownames(ret_catch)<-seq(1982,2020)

disc_cat<-fish_dat$discard_num_size_obs
colnames(disc_cat)<-seq(27.5,132.5,5)
rownames(disc_cat)<-seq(1992,2020)

bycatch<-fish_dat$trawl_num_size_obs
colnames(bycatch)<-seq(27.5,132.5,5)
rownames(bycatch)<-seq(1991,2020)

#============================
# fishery data (aggregate--used in previous versions)
#============================
fish_dat2<-read.csv("data/fisheries_dat.csv")

#=============================
# make covariate dataframe
#=============================
covars_imm<-data.frame(year=seq(1989,2020),
                   temperature=NA,
                   disease=NA,
                   cannibalism=NA,
                   bycatch=NA,
                   discard=NA,
                   predation=NA, 
                   imm_pop = NA,
                   mat_pop = NA)

covars_mat<-data.frame(year=seq(1989,2020),
                   temperature=NA,
                   disease=NA,
                   cannibalism=NA,
                   bycatch=NA,
                   discard=NA,
                   predation=NA, 
                   imm_pop = NA,
                   mat_pop = NA)


#===========================================================================
# Total population numbers (instead of observed)
# important for calculating impact of given catch/predation on population
#=========================================================================
use_mat_sel<-outs[[take_ind]]$pred_mat_pop_num
use_imm_sel<-outs[[take_ind]]$pred_imm_pop_num
tot_imm_pop<-apply(use_imm_sel,1,sum)
tot_mat_pop<-apply(use_mat_sel,1,sum)

covars_imm$temperature[which(!is.na(match(covars_imm$year,names(wtd_avg_temp_occ_imm))))]<-scale(wtd_avg_temp_occ_imm[which(!is.na(match(names(wtd_avg_temp_occ_imm),covars_imm$year)))])
covars_imm$disease[which(!is.na(match(covars_imm$year,use_disease_imm$AKFIN_SURVEY_YEAR)))]<-scale(use_disease_imm$med_dis[which(!is.na(match(use_disease_imm$AKFIN_SURVEY_YEAR,covars_imm$year)))])
covars_imm$cannibalism[which(!is.na(match(covars_imm$year,cannibal_ind_imm$year)))]<-scale(cannibal_ind_imm$med_cannib[which(!is.na(match(cannibal_ind_imm$year,covars_imm$year)))])
covars_imm$imm_pop[which(!is.na(match(covars_imm$year,years)))]<-scale(tot_imm_pop[which(!is.na(match(years,covars_imm$year)))])
covars_imm$mat_pop[which(!is.na(match(covars_imm$year,years)))]<-scale(tot_mat_pop[which(!is.na(match(years,covars_imm$year)))])

covars_mat$temperature[which(!is.na(match(covars_mat$year,names(wtd_avg_temp_occ_mat))))]<-scale(wtd_avg_temp_occ_mat[which(!is.na(match(names(wtd_avg_temp_occ_mat),covars_mat$year)))])
covars_mat$disease[which(!is.na(match(covars_mat$year,use_disease_mat$AKFIN_SURVEY_YEAR)))]<-scale(use_disease_mat$med_dis[which(!is.na(match(use_disease_mat$AKFIN_SURVEY_YEAR,covars_mat$year)))])
covars_mat$imm_pop[which(!is.na(match(covars_mat$year,years)))]<-scale(tot_imm_pop[which(!is.na(match(years,covars_mat$year)))])
covars_mat$mat_pop[which(!is.na(match(covars_mat$year,years)))]<-scale(tot_mat_pop[which(!is.na(match(years,covars_mat$year)))])

#========================
#==predation, bycatch, and discards need to be divide observed covar by the amount of 
#==est crab in water to make it comparable to a mortality rate

#==bycatch
common_yrs<-intersect(fish_dat2$Year,rownames(pred_n_size_imm))
use_bycatch<-fish_dat2$bycatch[which(!is.na(match(fish_dat2$Year,common_yrs)))]
bycatch_mort<-scale(use_bycatch/(tot_mat_pop[-length(tot_mat_pop)]+tot_imm_pop[-length(tot_mat_pop)]))

covars_imm$bycatch<-bycatch_mort
covars_mat$bycatch<-bycatch_mort

#==discard
common_yrs<-intersect(fish_dat2$Year,rownames(pred_n_size_imm))
use_discard<-fish_dat2$discard[which(!is.na(match(fish_dat2$Year,common_yrs)))]
discard_mort<-scale(use_discard/(tot_mat_pop[-length(tot_mat_pop)]))

covars_imm$discard<-NA
covars_mat$discard<-discard_mort


#==predation
alt_pred<-read.csv("data/cod_consump.csv")
wt_size<-read.csv("data/wt_at_size.csv")
use_wts<-filter(wt_size,Size%in%sizes)

common_yrs<-intersect(alt_pred$Year,rownames(pred_n_size_imm))
last_bin<-ncol(use_imm_sel)

tot_pred<-alt_pred$use_this[match(common_yrs,alt_pred$Year)]
cod_yrs<-alt_pred$Year[match(common_yrs,alt_pred$Year)]

#==uncomment this line of code and rerun analysis to explore impacts of different indices of predation
# tot_pred<-alt_pred$steve_50cm[match(common_yrs,alt_pred$Year)]
# tot_pred<-alt_pred$kir_nbs_ebs[match(common_yrs,alt_pred$Year)]
tot_crab_imm<-apply(outs[[1]]$`pred_imm_pop_num`[,1:13],1,sum)[-c(32,33)]
tot_crab_mat<-apply(outs[[1]]$`pred_mat_pop_num`[,1:13],1,sum)[-c(32,33)]
tot_crab_all<-(tot_crab_imm+tot_crab_mat)

covars_imm$predation[which(!is.na(match(covars_imm$year,alt_pred$Year)))]<-scale(tot_pred/tot_crab_all)
covars_mat$predation[which(!is.na(match(covars_mat$year,alt_pred$Year)))]<-scale(tot_pred/tot_crab_all)

#==show the comparison of different potential predation indices
all_cod<-alt_pred
colnames(all_cod)[4:5]<-c("Stomach contents","Abundance >50 cm")
tmp_pred<-all_cod$"Stomach contents"[match(common_yrs,all_cod$Year)]
tmp1<-data.frame(value=scale(tmp_pred/tot_crab_all),
           Year=seq(1989,2019),
           variable="Consumption ratio")
tmp_pred<-all_cod$"Abundance >50 cm"[match(common_yrs,all_cod$Year)]
tmp2<-data.frame(value=scale(tmp_pred/tot_crab_all),
                 Year=seq(1989,2019),
                 variable="Abundance ratio")

cod_fig<-melt(all_cod,id.vars="Year")
cod_fig_dat<-rbind(cod_fig,tmp1,tmp2)

raw_cod<-ggplot(filter(cod_fig_dat,variable=="Stomach contents" | variable =="Abundance >50 cm"))+
  geom_line(aes(x=Year,y=value))+
  facet_wrap(~variable,scale='free_y',ncol=1)+
  theme_bw()+
  expand_limits(y=0)+
  ylab('')

ratio_cod<-ggplot(filter(cod_fig_dat,variable=="Consumption ratio" | variable =="Abundance ratio"))+
  geom_line(aes(x=Year,y=value))+
  facet_wrap(~variable,scale='free_y',ncol=1)+
  theme_bw()+
  expand_limits(y=0)+
  ylab('') +
  scale_y_continuous(position='right')+
  geom_hline(data=data.frame(yint=c(0,0)),aes(yintercept=yint),lty=2)
 
png(paste("plots/cod_sensitivity.png",sep=''),height=6,width=8,res=350,units='in') 
raw_cod | ratio_cod
dev.off()

#==================================
# plot all covars
#==================================

plot_imm<-melt(covars_imm,id=c("year"))
plot_mat<-melt(covars_mat,id=c("year"))

#====pairs plots
 png(paste("plots/pairs_plot_mat.png",sep=''),height=8,width=8,res=350,units='in')
 in_pairs<- covars_mat[which(complete.cases(covars_mat[,c(2,3,5,6,7,8,9)])),c(2,3,5,6,7,8,9)]
 in_pairs$bycatch<-as.numeric(in_pairs$bycatch)
 in_pairs$discard<-as.numeric(in_pairs$discard)
 ggpairs(data=(in_pairs))
 dev.off()
 
 png(paste("plots/pairs_plot_imm.png",sep=''),height=8,width=8,res=350,units='in') 
 bleh<- covars_imm[which(complete.cases(covars_imm[,c(2:5,7,8,9)])),c(2:5,7,8,9)]
 bleh$bycatch<-as.numeric(bleh$bycatch)
 ggpairs(bleh)
 dev.off()

 #==time series plots 
 plot_mat$maturity<-"Mature covariates"
 plot_imm$maturity<-"Immature covariates"
 all_covars<-rbind(plot_mat,plot_imm)
 var_on_var_all<-ggplot()+
   geom_line(data=all_covars,aes(x=year,y=value,color=variable),lwd=1.25)+
   #geom_line(data=plt_1,aes(x=Year,y=value),lwd=2)+
   theme_bw()+
   theme(panel.grid.minor = element_blank(),
         panel.grid.major = element_blank())+
   #,       legend.position=c(.1,.85))+
   geom_hline(yintercept=0,lty=2)+facet_wrap(~maturity)

 var_on_var_all2<-ggplot()+
   geom_line(data=all_covars,aes(x=year,y=value,color=variable),lwd=1.25)+
   #geom_line(data=plt_1,aes(x=Year,y=value),lwd=2)+
   theme_bw()+
   theme(panel.grid.minor = element_blank(),
         panel.grid.major = element_blank())+
   #,       legend.position=c(.1,.85))+
   geom_hline(yintercept=0,lty=2)+facet_grid(vars(variable),vars(maturity))+
   xlim(1992,2019)
 
 png("plots/grid_covariates.png",height=8,width=5,res=350,units='in') 
 print(var_on_var_all2)
 dev.off()
 
 #========================================================
 # simulate mortality trajectories using covariance matrix
 # to be used in 
 #========================================================
 
 getADMBHessian <- function(direc){
   ## Thank you, Cole Monnahan
   filename <- file(paste(direc,"admodel.hes",sep=""), "rb")
   on.exit(close(filename))
   num.pars <- readBin(filename, "integer", 1)
   hes.vec <- readBin(filename, "numeric", num.pars^2)
   hes <- matrix(hes.vec, ncol=num.pars, nrow=num.pars)
   hybrid_bounded_flag <- readBin(filename, "integer", 1)
   scale <- readBin(filename, "numeric", num.pars)
   result <- list(num.pars=num.pars, hes=hes,
                  hybrid_bounded_flag=hybrid_bounded_flag, scale=scale)
   return(result)
 }
 
 admb_outs<-getADMBHessian(direc="C:/Users/cody.szuwalski/Work/snow_down/models/model_vary_m/")
 cov <- solve(admb_outs$hes)
 
 #==pull the rows of interest for covariance matrix
 imm_m_ind<-seq(27,59)
 mat_m_ind<-seq(60,92)
 
 imm_m_cov<-cov[imm_m_ind,imm_m_ind]
 mat_m_cov<-cov[mat_m_ind,mat_m_ind] 
 
 #==pull the estimated parameters
 styr<-1989
 endyr<-2021
 years<-seq(styr,endyr)
 par_file<-"models/model_vary_m/snow_down.par"
 imm_devs <- scan(par_file,skip=6,n=length(years),quiet=T)
 mat_devs <- scan(par_file,skip=8,n=length(years),quiet=T)
 imm_m_mu <- scan(par_file,skip=18,n=2,quiet=T)[1]
 mat_m_mu <- scan(par_file,skip=18,n=2,quiet=T)[2]
 
 #==little matrix multiplication for simulated 
 library(MASS)
 sim_imm_m<-mvrnorm(n=100,mu=imm_devs,Sigma=imm_m_cov)
 sim_mat_m<-mvrnorm(n=100,mu=mat_devs,Sigma=mat_m_cov) 
 
 trans_sim_imm <- exp(imm_m_mu + sim_imm_m)
 trans_sim_mat <- exp(mat_m_mu + sim_mat_m)
 
 mle_imm<-exp(imm_m_mu + imm_devs)
 mle_mat<- exp(mat_m_mu + mat_devs)
 
 colnames(trans_sim_imm)<-years
 colnames(trans_sim_mat)<-years
 
 melted<-melt(t(trans_sim_imm))
 melted2<-melt(t(trans_sim_mat))
 melted$state<-"Immature"
 melted2$state<-"Mature"
 in_plot<-rbind(melted,melted2)
 
 mle_est<-data.frame(X1=seq(1989,2021),value=c(mle_imm,mle_mat),state=c(rep("Immature",33),rep("Mature",33)))
 
 png("plots/covariance_sims.png",height=8,width=5,res=400,units='in')
 ggplot(filter(in_plot,X1<2020))+
   geom_line(aes(x=X1,y=value,group=X2),alpha=0.1,color='grey')+
   geom_line(data=filter(mle_est,X1<2020),aes(x=X1,y=value))+
   theme_bw()+
   facet_wrap(~state,scales='free_y',ncol=1)+
   expand_limits(y=0)+
   ylab("Estimated mortality")+
   xlab("Year")
 dev.off()
 
#=============================================
# build models for mortality 
#=============================================
##############################################
est_vars<-data.frame(year=years,
                     est_mat_m=(outs[[take_ind]]$'mature natural mortality')[,1],
                     est_imm_m=(outs[[take_ind]]$'natural mortality')[,1])

#===================
#==mature mortality
#==================
#==GAM==============
 input_dat<-merge(est_vars[,c(1,2)],covars_mat)[,-c(5)]
 #write.csv(input_dat,"mat_mort.csv")
 input_dat<-input_dat[complete.cases(input_dat),]
 colnames(input_dat)[2]<-"dep_var"
 input_dat<-input_dat[-c(1),] #==chop first year of data because of initialization conditions
 #ggpairs(data=input_dat)
 
 mod_base<-gam(data=input_dat,(dep_var)~s(temperature,k=4)+s(disease,k=3)+s(discard,k=3)+s(bycatch,k=4)+s(mat_pop,k=3)+s(predation,k=3)) 
 mod_base1<-gam(data=input_dat,(dep_var)~s(temperature,k=4)+s(disease,k=3)+s(discard,k=3)+s(bycatch,k=4)+s(mat_pop,k=3)) 
 summary(mod_base)
 plot(mod_base,page=1)
 summary(mod_base1)
 plot(mod_base1,page=1)
 
 mod_base_temp<-gam(data=input_dat,(dep_var)~s(temperature,k=4)) 
 summary(mod_base_temp)
 plot(mod_base_temp,page=1)
 
 AIC(mod_base1)
 AIC(mod_base)

 mod_base_t<-gam(data=input_dat,(dep_var)~s(temperature,k=4)+s(mat_pop,k=3) )
 summary(mod_base_t)
 plot(mod_base_t,page=1)
 
#==try shape constrained additive models for mature population--monotonic increasing 
b <- scam(data=input_dat,dep_var~s(temperature,k=3)+s(mat_pop,bs='mpi')+s(discard,k=3)+s(predation,k=3)+s(bycatch,k=3) )
b <- scam(data=input_dat,dep_var~s(disease,k=3)+s(temperature,k=3)+s(mat_pop,bs='mpi')+s(discard,k=3)+s(predation,k=3)+s(bycatch,k=3) )
b <- scam(data=input_dat,dep_var~disease+s(temperature,k=3)+s(mat_pop,bs='mpi')+discard+predation+bycatch )
summary(b)
plot(b,pages=1)

#==try betar to see if assumptions about distribution impact the 
beta_dat<-input_dat
beta_dat$dep_var<-1-exp(-beta_dat$dep_var)
bm <- gam((dep_var)~s(temperature,k=4)+s(disease,k=3)+s(discard,k=3)+s(bycatch,k=4)+s(mat_pop,k=3)+s(predation,k=3),family=betar(link="logit"),data=beta_dat)
summary(bm)
plot(bm,pages=1)

#==prep GAM tables and diagnostics for supplementary materials
gamtabs(mod_base,caption=paste("GAM output for full model predicting mature mortality. Deviance explained = ",round(100*summary(mod_base)$dev,2),"%"))
gamtabs(mod_base_t,caption=paste("GAM output for trimmed model predicting mature mortality. Deviance explained = ",round(100*summary(mod_base_t)$dev,2),"%"))

png("plots/gam_smooths_mat.png",height=5,width=8,res=350,units='in') 
draw(mod_base,scales='fixed')&theme_bw()
dev.off()

png("plots/gam_check_full_mat.png",height=8,width=8,res=350,units='in') 
par(mfrow=c(2,2))
gam.check(mod_base)
dev.off()

png("plots/gam_check_trim_mat.png",height=8,width=8,res=350,units='in') 
par(mfrow=c(2,2))
gam.check(mod_base_t)
dev.off()


preds   <- predict(mod_base, se.fit = TRUE)
my_data <- data.frame(input_dat,
                            mu   = (preds$fit),
                            low  = (preds$fit - 1.96 * preds$se.fit),
                            high = (preds$fit + 1.96 * preds$se.fit))
dev_expl<-summary(mod_base)$dev

#==cross validation============
cv_plot_mat_m<-list()
use_dat<-input_dat
dev_expl_cv<-rep(0,nrow(use_dat))
predicted_cv<-matrix(nrow=nrow(use_dat),ncol=nrow(use_dat)-1)
cv_yrs<-matrix(nrow=nrow(use_dat),ncol=nrow(use_dat)-1)
smooth_table<-list()
p_table<-list()
for(x in 1:nrow(use_dat))
{
  cv_dat<-use_dat[-c(x),]
  #mod_cv<-gam(data=cv_dat,dep_var~s(temperature,k=4)+s(disease,k=3)+s(discard,k=3)+s(bycatch,k=4)+s(mat_pop,k=3) )
  mod_cv<-gam(data=cv_dat,(dep_var)~s(temperature,k=4)+s(disease,k=3)+s(discard,k=3)+s(bycatch,k=4)+s(mat_pop,k=3)+s(predation,k=3)) 
  dev_expl_cv[x]<-summary(mod_cv)$dev.expl
  cv_plot_mat_m[[x]]<-plot(mod_cv,page=1)
  predicted_cv[x,]<-predict(mod_cv)
  cv_yrs[x,]<-cv_dat$year
  smooth_table[[x]]<-summary(mod_cv)$s.table
  p_table[[x]]<-summary(mod_cv)$p.table
  print(x)
}

#==massage data for ggplots
mat_cv_M_outs<-data.frame(pred=c(predicted_cv),
                          year=c(cv_yrs),
                          replic=as.factor(rep(1:nrow(predicted_cv),times=nrow(predicted_cv)-1)),
                          fits="Mature mortality")

mat_cv_M_dev<-data.frame(deviance_explained=dev_expl_cv,
                         fits="Mature mortality")

in_names<-c("Temperature","Disease","Discard","Bycatch","Mature population","Predation")
smooth_cv_outs<-matrix(ncol=length(smooth_table),nrow=nrow(smooth_table[[1]]))
for(x in 1:length(smooth_table))
  try(smooth_cv_outs[,x]<-smooth_table[[x]][,4],silent=TRUE)
rownames(smooth_cv_outs)<-in_names
colnames(smooth_cv_outs)<-seq(1,ncol(smooth_cv_outs))
mat_cv_M_pval<-melt(smooth_cv_outs)
colnames(mat_cv_M_pval)<-c("Process","Rep","Value")
mat_cv_M_pval$fits<-"Mature mortality"

plot_pval_mat<-ggplot(data=mat_cv_M_pval)+
  geom_histogram(aes(x=Value,group=Process,fill=Process),bins=50,alpha=.8)+
  theme_bw()

#==test impact of uncertainty in estimates============
unc_plot_mat_m<-list()
use_dat<-input_dat
dev_expl_unc<-rep(0,nrow(trans_sim_mat))
predicted_unc<-matrix(nrow=nrow(trans_sim_mat),ncol=nrow(input_dat))
unc_yrs<-matrix(nrow=nrow(trans_sim_mat),ncol=nrow(input_dat))
smooth_table_unc<-list()

for(x in 1:nrow(trans_sim_mat))
{
  unc_dat<-use_dat
  use_dat$dep_var<-trans_sim_mat[x,  which(!is.na(match(colnames(trans_sim_mat),use_dat$year)))]
  mod_unc<-gam(data=use_dat,dep_var~s(temperature,k=4)+s(disease,k=3)+s(discard,k=3)+s(bycatch,k=4)+s(mat_pop,k=3)+s(predation,k=3) )
  
  dev_expl_unc[x]<-summary(mod_unc)$dev.expl
  predicted_unc[x,]<-predict(mod_unc)
  unc_yrs[x,]<-use_dat$year
  smooth_table_unc[[x]]<-summary(mod_unc)$s.table
}

#==massage data for ggplots
colnames(predicted_unc)<-unc_yrs[1,]
rownames(predicted_unc)<-seq(1,nrow(predicted_unc))
mat_unc<-melt(predicted_unc)
mat_unc$fits<-"Mature mortality"

mat_unc_M_dev<-data.frame(deviance_explained=dev_expl_unc,
                          fits="Mature mortality")

in_names<-c("Temperature","Disease","Discard","Bycatch","Mature population","Predation")
smooth_unc_outs<-matrix(ncol=length(smooth_table_unc),nrow=nrow(smooth_table_unc[[1]]))

for(x in 1:length(smooth_table_unc))
  smooth_unc_outs[,x]<-smooth_table_unc[[x]][,4]
rownames(smooth_unc_outs)<-in_names
colnames(smooth_unc_outs)<-seq(1,ncol(smooth_unc_outs))

mat_unc_M_pval<-melt(smooth_unc_outs)
colnames(mat_unc_M_pval)<-c("Process","Rep","Value")
mat_unc_M_pval$fits<-"Mature mortality"

plot_pval_mat<-ggplot(data=mat_unc_M_pval)+
  geom_histogram(aes(x=Value,group=Process,fill=Process),bins=50,alpha=.8)+
  theme_bw()


#==prediction==================

my_data_mat_m_pred<-NULL
for(y in 0:2)
{
train_dat<-input_dat[-seq(nrow(input_dat)-y,nrow(input_dat)),]
mod_pred<-gam(data=train_dat,dep_var~s(temperature,k=4)+s(mat_pop,k=3))
preds<-predict(mod_pred, se.fit = TRUE,newdata=input_dat)
my_data_mat_m_pred <- rbind(my_data_mat_m_pred,data.frame(input_dat,
                            mu   = (preds$fit),
                            low  = (preds$fit - 1.96 * preds$se.fit),
                            high = (preds$fit + 1.96 * preds$se.fit),
                            peel= rep(y+1,length(preds$fit))))

}

#==plot full model, but with 
my_data$peel<-0
cc <- scales::seq_gradient_pal('blue','light blue',"Lab")(seq(0,1,.33))
big_dat<-rbind(my_data,my_data_mat_m_pred)
mat_fits_plot<-ggplot( )+
  geom_line(data=mat_cv_M_outs,aes(x=year,y=pred,group=replic),col='grey')+
  geom_line(data=big_dat,aes(x=year, y = mu,col=as.factor(peel)),lwd=1.2,alpha=.5)+
  scale_color_manual(values=cc,name="Years projected")+
  geom_line(data=my_data, aes(x = year, y = mu),lwd=2,col='blue') +
  geom_line(data=my_data, aes(x = year, y = dep_var),lwd=1.2,col='red') + 
  geom_point(data=my_data, aes(x = year, y = dep_var),col='red') + 

  theme_bw()+
  theme(panel.border=element_blank(),
        panel.background=element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),legend.position='none')+
  scale_x_continuous(breaks=c(1995,2005,2015),name='Year')+
  ylab("Mortality")


mat_fits_plot_alt<-ggplot( )+
   geom_line(data=big_dat,aes(x=year, y = mu,col=as.factor(peel)),lwd=1.2,alpha=.5)+
  scale_color_manual(values=cc,name="Years projected")+
  geom_line(data=my_data, aes(x = year, y = mu),lwd=2,col='blue') +
  geom_line(data=my_data, aes(x = year, y = dep_var),lwd=1.2,col='red') + 
  geom_point(data=my_data, aes(x = year, y = dep_var),col='red') + 
  
  theme_bw()+
  theme(panel.border=element_blank(),
        panel.background=element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),legend.position='none')+
  scale_x_continuous(breaks=c(1995,2005,2015),name='Year')+
  ylab("Mortality")

#==randomization===============

n_sims<-1000

in_dat_sim<-input_dat
pred_out<-matrix(ncol=n_sims,nrow=nrow(in_dat_sim))
s_table<-list(list())
dev_expl_sim<-rep(0,n_sims)
for(y in 1:n_sims)
{
  for(z in c(3,4,5,6,7))
    in_dat_sim[,z]<-sample(input_dat[,z],length(input_dat[,z]),replace=TRUE)

  mod_sim<-gam(data=in_dat_sim,dep_var~s(temperature,k=3)+s(mat_pop,k=3))
  s_table[[y]]<-summary(mod_sim)$s.table
  dev_expl_sim[y]<-summary(mod_sim)$dev.expl
}

png("plots/rando_mat.png",height=5,width=8,res=350,units='in') 
hist(dev_expl_sim,breaks=seq(0,1,0.01),las=1,xlab="Deviance explained",main="")
abline(v=sort(dev_expl_sim)[.95*n_sims],lty=2,col=2)
abline(v=dev_expl,lwd=2,col=4)
legend('topleft',bty='n',col=c(2,4),lty=c(2,1),legend=c("Randomized 95th percentile","Observed data"))
dev.off()

#=======================================
#==immature mortality
#=======================================
#===GAM=====
input_dat<-merge(est_vars[,c(1,3)],covars_imm)[,-7]
#write.csv(input_dat,"imm_mort.csv")
input_dat<-input_dat[complete.cases(input_dat),]
colnames(input_dat)[2]<-"dep_var"
input_dat<-input_dat[-c(1),]
#ggpairs(data=input_dat)

mod_base_t<-gam(data=input_dat,dep_var~s(temperature,k=3)+s(mat_pop,k=3))
mod_base<-gam(data=input_dat,dep_var~s(disease,k=3)+s(temperature,k=3)+s(mat_pop,k=3)+s(imm_pop,k=3)+s(predation,k=3)+s(bycatch,k=3) +s(cannibalism,k=3))

summary(mod_base)
plot(mod_base,page=1)

summary(mod_base_t)
plot(mod_base_t,page=1)

mod_base_temp<-gam(data=input_dat,(dep_var)~s(temperature,k=4)) 
summary(mod_base_temp)
plot(mod_base_temp,page=1)


#==check a two dimensional smooth
mod_base_int<-gam(data=input_dat,dep_var~s(mat_pop,temperature,k=8))
summary(mod_base_int)
plot(mod_base_int,page=1,too.far=0,cex=3)

#==try shape constrained additive models
b <- scam(data=input_dat,dep_var~disease+s(temperature,k=3)+s(mat_pop,bs='mpi')+imm_pop+predation+bycatch +cannibalism)
summary(b)
png("plots/scam_output.png",height=5,width=8,res=350,units='in') 
plot(b,pages=1)
dev.off()

b_t <- scam(data=input_dat,dep_var~s(temperature,k=3)+s(mat_pop,bs='mpi'))
summary(b_t)
plot(b_t,pages=1)

dev_expl<-summary(mod_base)$dev


#==check a beta distribution
beta_dat<-input_dat
beta_dat$dep_var<-1-exp(-beta_dat$dep_var)
bm <- gam((dep_var)~s(disease,k=3)+s(temperature,k=3)+s(mat_pop,k=3)+s(predation,k=3)+s(bycatch,k=3) +s(cannibalism,k=3)+s(imm_pop,k=3),family=betar(link="logit"),data=beta_dat)
bm_trim <- gam((dep_var)~s(temperature,k=4)+s(mat_pop,k=3),family=betar(link="logit"),data=beta_dat)
summary(bm)
plot(bm,pages=1)

gamtabs(mod_base,caption=paste("GAM output for full model predicting immature mortality. Deviance explained = ",round(100*summary(mod_base)$dev,2),"%"))
gamtabs(mod_base_t,caption=paste("GAM output for trimmed model predicting immature mortality. Deviance explained = ",round(100*summary(mod_base_t)$dev,2),"%"))
gamtabs(b,caption=paste("GAM output for trimmed model predicting immature mortality. Deviance explained = ",round(100*summary(mod_base_t)$dev,2),"%"))

png("plots/gam_smooths_imm.png",height=5,width=8,res=350,units='in') 
draw(mod_base,scales='fixed')&theme_bw()
dev.off()

png("plots/gam_check_full_imm.png",height=8,width=8,res=350,units='in') 
par(mfrow=c(2,2))
gam.check(mod_base)
dev.off()

png("plots/gam_check_trim_imm.png",height=8,width=8,res=350,units='in') 
par(mfrow=c(2,2))
gam.check(mod_base_t)
dev.off()

preds   <- predict(mod_base, se.fit = TRUE)
my_data <- data.frame(input_dat,
                            mu   = (preds$fit),
                            low  = (preds$fit - 1.96 * preds$se.fit),
                            high = (preds$fit + 1.96 * preds$se.fit))

#==cross validation============
cv_plot_imm_m<-list()
use_dat<-input_dat
dev_expl_cv<-rep(0,nrow(use_dat))
predicted_cv<-matrix(nrow=nrow(use_dat),ncol=nrow(use_dat)-1)
cv_yrs<-matrix(nrow=nrow(use_dat),ncol=nrow(use_dat)-1)
smooth_table<-list()
for(x in 1:nrow(use_dat))
{
  cv_dat<-use_dat[-c(x),]
  #mod_cv<-gam(data=cv_dat,dep_var~s(disease,k=3)+s(temperature,k=3)+s(mat_pop,k=3)+s(predation,k=3)+s(bycatch,k=3) +s(cannibalism,k=3))
  mod_cv<-gam(data=cv_dat,dep_var~s(disease,k=3)+s(temperature,k=3)+s(mat_pop,k=3)+s(imm_pop,k=3)+s(predation,k=3)+s(bycatch,k=3) +s(cannibalism,k=3))
 # mod_cv <- scam(data=cv_dat,dep_var~disease+s(temperature,k=3)+s(mat_pop,bs='mpi')+s(imm_pop,k=3)+s(predation,k=3)+s(bycatch,k=3) +s(cannibalism,k=3))
  
  dev_expl_cv[x]<-summary(mod_cv)$dev.expl
  cv_plot_imm_m[[x]]<-plot(mod_cv,page=1)
  predicted_cv[x,]<-predict(mod_cv)
  cv_yrs[x,]<-cv_dat$year
  smooth_table[[x]]<-summary(mod_cv)$s.table
}

#==massage data for ggplots
imm_cv_M_outs<-data.frame(pred=c(predicted_cv),
                          year=c(cv_yrs),
                          replic=as.factor(rep(1:nrow(predicted_cv),times=nrow(predicted_cv)-1)),
                          fits="Immature mortality")

imm_cv_M_dev<-data.frame(deviance_explained=dev_expl_cv,
                         fits="Immature mortality")

in_names<-c("Disease","Temperature","Mature population","Immature population","Predation","Bycatch","Cannibalism")

smooth_cv_outs<-matrix(ncol=length(smooth_table),nrow=nrow(smooth_table[[1]]))
for(x in 1:length(smooth_table))
  smooth_cv_outs[,x]<-smooth_table[[x]][,4]
rownames(smooth_cv_outs)<-in_names
colnames(smooth_cv_outs)<-seq(1,ncol(smooth_cv_outs))
imm_cv_M_pval<-melt(smooth_cv_outs)
colnames(imm_cv_M_pval)<-c("Process","Rep","Value")
imm_cv_M_pval$fits<-"Immature mortality"

ggplot(data=imm_cv_M_pval)+
  geom_histogram(aes(x=Value,group=Process,fill=Process),bins=50,alpha=.8)+
  theme_bw()


#==test impact of uncertainty in mortality estimates============
unc_plot_imm_m<-list()
use_dat<-input_dat
dev_expl_unc<-rep(0,nrow(trans_sim_imm))
predicted_unc<-matrix(nrow=nrow(trans_sim_imm),ncol=nrow(input_dat))
unc_yrs<-matrix(nrow=nrow(trans_sim_imm),ncol=nrow(input_dat))
smooth_table_unc<-list()

for(x in 1:nrow(trans_sim_imm))
{
  unc_dat<-use_dat
  use_dat$dep_var<-trans_sim_imm[x,  which(!is.na(match(colnames(trans_sim_imm),use_dat$year)))]
  mod_unc<-gam(data=use_dat,dep_var~s(disease,k=3)+s(temperature,k=3)+s(mat_pop,k=3)+s(imm_pop,k=3)+s(predation,k=3)+s(bycatch,k=3) +s(cannibalism,k=3))
  
  dev_expl_unc[x]<-summary(mod_unc)$dev.expl
  predicted_unc[x,]<-predict(mod_unc)
  unc_yrs[x,]<-use_dat$year
  smooth_table_unc[[x]]<-summary(mod_unc)$s.table
}

#==massage data for ggplots
colnames(predicted_unc)<-unc_yrs[1,]
rownames(predicted_unc)<-seq(1,nrow(predicted_unc))
imm_unc<-melt(predicted_unc)
imm_unc$fits<-"Immature mortality"

imm_unc_M_dev<-data.frame(deviance_explained=dev_expl_unc,
                          fits="Immature mortality")

in_names<-c("Disease","Temperature","Mature population","Immature population","Predation","Bycatch","Cannibalism")
smooth_unc_outs<-matrix(ncol=length(smooth_table_unc),nrow=nrow(smooth_table_unc[[1]]))

for(x in 1:length(smooth_table_unc))
  smooth_unc_outs[,x]<-smooth_table_unc[[x]][,4]
rownames(smooth_unc_outs)<-in_names
colnames(smooth_unc_outs)<-seq(1,ncol(smooth_unc_outs))

imm_unc_M_pval<-melt(smooth_unc_outs)
colnames(imm_unc_M_pval)<-c("Process","Rep","Value")
imm_unc_M_pval$fits<-"Immature mortality"

ggplot(data=imm_unc_M_pval)+
  geom_histogram(aes(x=Value,group=Process,fill=Process),bins=50,alpha=.8)+
  theme_bw()


unc_plot<-rbind(imm_unc_M_pval,mat_unc_M_pval)
unc_devs<-rbind(imm_unc_M_dev,mat_unc_M_dev)

unc_p_pval<-ggplot(data=unc_plot)+
  geom_histogram(aes(x=Value,group=Process,fill=Process),bins=50,alpha=.8)+
  theme_bw()+
  facet_wrap(~fits)+
  theme(legend.position=c(0.8,0.7))+
  xlab("P-value")

unc_p_dev<-ggplot(data=unc_devs)+
  geom_histogram(aes(x=deviance_explained),bins=50,alpha=.8)+
  theme_bw()+
  facet_wrap(~fits,ncol=1)+
  xlim(0,1)+xlab("Deviance explained")

png("plots/unc_pvals.png",height=6,width=8,res=350,units='in') 
(unc_p_pval |  unc_p_dev) + plot_layout(widths=c(3,1))
dev.off()

#==prediction==================

my_data_imm_m_pred<-NULL
for(y in 0:2)
{
  train_dat<-input_dat[-seq(nrow(input_dat)-y,nrow(input_dat)),]
  mod_pred<-gam(data=train_dat,dep_var~s(temperature,k=3)+s(mat_pop,k=3))
  #mod_pred <- scam(data=train_dat,dep_var~s(temperature,k=3)+s(mat_pop,bs='mpi')+s(imm_pop,k=3))
  
  preds<-predict(mod_pred, se.fit = TRUE,newdata=input_dat)
  my_data_imm_m_pred <- rbind(my_data_imm_m_pred,data.frame(input_dat,
                                                            mu   = (preds$fit),
                                                            low  = (preds$fit - 1.96 * preds$se.fit),
                                                            high = (preds$fit + 1.96 * preds$se.fit),
                                                            peel= rep(y+1,length(preds$fit))))
  
}

#==plot full model, but with 
my_data$peel<-0
cc <- scales::seq_gradient_pal('blue','light blue',"Lab")(seq(0,1,.33))
big_dat<-rbind(my_data,my_data_imm_m_pred)
imm_fits_plot<-ggplot( )+
  geom_line(data=imm_cv_M_outs,aes(x=year,y=pred,group=replic),col='grey')+
  geom_line(data=big_dat,aes(x=year, y = mu,col=as.factor(peel)),lwd=2)+
  scale_color_manual(values=cc,name="Years projected")+
  geom_line(data=my_data, aes(x = year, y = mu),lwd=2,col='blue') +
  geom_line(data=my_data, aes(x = year, y = dep_var),lwd=1.2,col='red') + 
  geom_point(data=my_data, aes(x = year, y = dep_var),col='red') +   
  theme_bw()+
  theme(panel.border=element_blank(),
        panel.background=element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position=c(.5,.85),
        axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.title=element_text(size=8), 
        legend.text=element_text(size=7),
        legend.key.size = unit(1,"line"))+
  scale_y_continuous(position='left',name="Mortality")

imm_fits_plot_alt<-ggplot( )+
  geom_line(data=big_dat,aes(x=year, y = mu,col=as.factor(peel)),lwd=2)+
  scale_color_manual(values=cc,name="Years projected")+
  geom_line(data=my_data, aes(x = year, y = mu),lwd=2,col='blue') +
  geom_line(data=my_data, aes(x = year, y = dep_var),lwd=1.2,col='red') + 
  geom_point(data=my_data, aes(x = year, y = dep_var),col='red') +   
  theme_bw()+
  theme(panel.border=element_blank(),
        panel.background=element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position=c(.5,.85),
        axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.title=element_text(size=8), 
        legend.text=element_text(size=7),
        legend.key.size = unit(1,"line"))+
  scale_y_continuous(position='left',name="Mortality")


#==randomization===============

n_sims<-1000

in_dat_sim<-input_dat
pred_out<-matrix(ncol=n_sims,nrow=nrow(in_dat_sim))
s_table<-list(list())
dev_expl_sim<-rep(0,n_sims)
for(y in 1:n_sims)
{
  for(z in c(3,4,5,6,7))
    in_dat_sim[,z]<-sample(input_dat[,z],length(input_dat[,z]),replace=TRUE)
  
  mod_sim<-gam(data=in_dat_sim,dep_var~s(temperature,k=3)+s(mat_pop,k=3))

  s_table[[y]]<-summary(mod_sim)$s.table
  dev_expl_sim[y]<-summary(mod_sim)$dev.expl
}

png("plots/rando_imm.png",height=5,width=8,res=350,units='in') 
hist(dev_expl_sim,breaks=seq(0,1,0.01),las=1,xlab="Deviance explained",main="")
abline(v=sort(dev_expl_sim)[.95*n_sims],lty=2,col=2)
abline(v=dev_expl,lwd=2,col=4)
legend('topleft',bty='n',col=c(2,4),lty=c(2,1),legend=c("Randomized 95th percentile","Observed data"))
dev.off()


######################################
# make big crossvalidation figure
#####################################
big_cv_outs<-rbind(mat_cv_M_outs,imm_cv_M_outs)
big_cv_dev<-rbind(mat_cv_M_dev,imm_cv_M_dev)
big_cv_pval<-rbind(mat_cv_M_pval,imm_cv_M_pval)
sm_cv_pval<-rbind(mat_cv_M_pval,imm_cv_M_pval)

big_cv_dev_in<-big_cv_dev
big_cv_dev_in$deviance_explained<-100*big_cv_dev_in$deviance_explained
plot_dev<-ggplot(data=big_cv_dev_in)+
  geom_histogram(aes(x=deviance_explained),bins=50,fill='blue',alpha=0.5)+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  facet_wrap(~fits,ncol=1,scales='free_y')+
  theme( strip.background = element_blank(),
         strip.text.x = element_blank(),
         panel.border = element_blank(),
         panel.background = element_blank(),
         axis.title.y=element_blank(),
         axis.text.y=element_blank(),
         axis.ticks.y=element_blank()
  )+
  scale_x_continuous(breaks=c(0,50,100),limits=c(0,100))+
  xlab("Deviance explained")

plot_dev_mat<-ggplot(data=filter(big_cv_dev_in,fits=="Mature mortality"))+
  geom_histogram(aes(x=deviance_explained),bins=20,fill='blue',alpha=0.5)+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  theme( strip.background = element_blank(),
         strip.text.x = element_blank(),
         panel.border = element_blank(),
         panel.background = element_blank(),
         axis.title.y=element_blank(),
         axis.text.y=element_blank(),
         axis.ticks.y=element_blank()
  )+
  scale_x_continuous(breaks=c(0,50,100),limits=c(0,100))+
  xlab("Deviance 
explained")

plot_dev_imm<-ggplot(data=filter(big_cv_dev_in,fits=="Immature mortality"))+
  geom_histogram(aes(x=deviance_explained),bins=20,fill='blue',alpha=0.5)+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  theme( strip.background = element_blank(),
         strip.text.x = element_blank(),
         panel.border = element_blank(),
         panel.background = element_blank(),
         axis.title.y=element_blank(),
         axis.text.y=element_blank(),
         axis.ticks.y=element_blank(),
         axis.title.x=element_blank(),
         axis.text.x=element_blank(),
         axis.ticks.x=element_blank()
  )+
  scale_x_continuous(breaks=c(0,50,100),limits=c(0,100))

temp_cv<-big_cv_pval
temp_cv$Value[temp_cv$fits=="Mature mortality"]<-NA
plot_pval<-ggplot(data=big_cv_pval)+
  geom_histogram(aes(x=Value,group=Process,fill=Process),bins=50,alpha=.8)+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  xlim(-0.1,1)+
  facet_wrap(~fits,ncol=1,scales='free_y')+
  theme( strip.background = element_blank(),
         strip.text.x = element_blank(),
         panel.border = element_blank(),
         panel.background = element_blank(),
         axis.title.y=element_blank(),
         axis.text.y=element_blank(),
         axis.ticks.y=element_blank(),
         legend.position=c(.45,.85),
         legend.title = element_blank(),
         legend.key=element_blank())+
  geom_vline(xintercept=0.06,lty=2)+
  xlab("p-value")

temp_cv<-big_cv_pval
temp_cv$Value[temp_cv$fits=="Mature mortality"]<-1
temp_cv$fits[temp_cv$fits=="Mature mortality"]<-"Immature mortality"
plot_pval_imm<-ggplot(data=temp_cv)+
  geom_histogram(aes(x=Value,group=Process,fill=Process),bins=50,alpha=.8)+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  xlim(-0.1,1)+
  theme( strip.background = element_blank(),
         strip.text.x = element_blank(),
         panel.border = element_blank(),
         panel.background = element_blank(),
         # axis.title.y=element_blank(),
         # axis.text.y=element_blank(),
         # axis.ticks.y=element_blank(),
         axis.title.x=element_blank(),
         axis.text.x=element_blank(),
         axis.ticks.x=element_blank(),
         legend.position=c(.68,.6),
         legend.title = element_blank(),
         legend.key=element_blank(),
         legend.background = element_rect(fill = "transparent"),
         legend.text=element_text(size=6.8),
         legend.key.size = unit(.2, 'cm'))+
  geom_vline(xintercept=0.06,lty=2)+
  xlab("p-value")+
  scale_y_continuous(position="right",limits=c(0,60))

temp_cv<-big_cv_pval
temp_cv$Value[temp_cv$fits=="Immature mortality"]<-1
temp_cv$fits[temp_cv$fits=="Immature mortality"]<-"Mature mortality"
plot_pval_mat<-ggplot(data=temp_cv)+
  geom_histogram(aes(x=Value,group=Process,fill=Process),bins=50,alpha=.8)+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  xlim(-0.1,1)+
  theme( strip.background = element_blank(),
         strip.text.x = element_blank(),
         panel.border = element_blank(),
         panel.background = element_blank(),
         legend.position='none',
         legend.title = element_blank(),
         legend.key=element_blank())+
  geom_vline(xintercept=0.06,lty=2)+
  xlab("p-value")+
  scale_y_continuous(position="right",limits=c(0,50))

alf<-imm_abnd+mat_abnd+imm_size+mat_size+
  imm_fits_plot+mat_fits_plot+plot_dev_imm+plot_dev_mat +
  plot_pval_imm+plot_pval_mat+plot_layout(byrow=FALSE,ncol=5,nrow=2,widths=c(1,1,1,.3,.7))+
  plot_layout(tag_level = 'new') +
  plot_annotation(tag_levels = list(c('a','','b','','c','','d','','e','')))+
  plot_annotation(title="Population Dynamics                                                                      Generalized Additive Models")


png('plots/figure_2.png',height=5,width=10,res=350,units='in')
print(alf)
dev.off()

png('plots/predict_big.png',height=8,width=8,res=350,units='in')
imm_fits_plot_alt / mat_fits_plot_alt
dev.off()