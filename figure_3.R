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
library(gratia)
library(itsadug)

survDAT<-read.csv("data/EBSCrab_Haul/EBSCrab_Haul.csv",header=T,skip=5)
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
#==the trends are a lot different by size...
in_dat<-filter(temp_tot,SEX==1,mid_pt>27.5&mid_pt<135)
temp_occ<-ggplot(in_dat)+
  geom_line(aes(x=AKFIN_SURVEY_YEAR,y=temp_occ,col=as.factor(mid_pt),group=mid_pt),
            lwd=1,alpha=.7)+
  geom_point(aes(x=AKFIN_SURVEY_YEAR,y=temp_occ,col=as.factor(mid_pt),group=mid_pt),
             lwd=1)+
  theme_bw()+
  scale_y_continuous(name='Temperature 
  occupied (C)',position='right')+
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

############################################
# temperature mediated weight at size
##########################################
in_wt_m<-filter(survDAT,!is.na(WEIGHT)&SEX==1,GEAR_TEMPERATURE<5,AKFIN_SURVEY_YEAR>2001,AKFIN_SURVEY_YEAR!=2016&SHELL_CONDITION<=2,AKFIN_SURVEY_YEAR!=2021)
mod_m<-gam(data=in_wt_m,WEIGHT~s(WIDTH)+s(GEAR_TEMPERATURE,k=4)+as.factor(AKFIN_SURVEY_YEAR))
gamtabs(mod_m,caption='GAM output for full model predicting immature mortality. Deviance explained = 72.2%')

png('plots/wt_size_smooth.png',height=6,width=9,res=350,units='in')
draw(mod_m)&theme_bw()
dev.off()

sizes<-seq(35,135,1)
pred_17<-predict(mod_m,newdata=data.frame(WIDTH=sizes,GEAR_TEMPERATURE=2,AKFIN_SURVEY_YEAR=2017))
pred_18<-predict(mod_m,newdata=data.frame(WIDTH=sizes,GEAR_TEMPERATURE=2,AKFIN_SURVEY_YEAR=2018))

preds<-data.frame(pred=c(pred_17,pred_18),
                  width=rep(sizes,2),
                  temp=c(rep(2,length(pred_17)),rep(2,length(pred_17))),
                  year=c(rep(2017,length(pred_17)),rep(2018,length(pred_17))))

wt_wid<-ggplot()+
  geom_point(data=in_wt_m,
             aes(x=WIDTH,y=(WEIGHT),col=(GEAR_TEMPERATURE)),alpha=.5)+
  scale_colour_distiller(palette="Spectral")+
  geom_line(data=preds,aes(x=width,y=pred,group=year))+
  theme_bw()+
  theme(legend.position=c(.3,.8))+
  labs(colour="Bottom temperature")+
  xlab("Carapace width (mm)")+
  ylab("Weight (g)")+
  annotate(geom='text',x=sizes[length(sizes)]+10,y=pred_17[length(pred_17)]+10,
           label="2017")+
  annotate(geom='text',x=sizes[length(sizes)]+10,y=pred_18[length(pred_18)]-5,
           label="2018")+
  xlim(20,150)+ plot_annotation(tag_levels='a')

wt_at_s<-read.csv("data/wt_at_size.csv")  

png('plots/wt_size_year_temp.png',height=6,width=9,res=350,units='in')
ggplot()+
  geom_point(data=in_wt_m,
             aes(x=WIDTH,y=(WEIGHT),col=(GEAR_TEMPERATURE)),alpha=.5)+
  facet_wrap(~AKFIN_SURVEY_YEAR)+
  scale_colour_distiller(palette="Spectral")+
  labs(colour="Temp")+
  geom_line(data=wt_at_s,aes(x=Size,y=Weight*1000))+
  theme_bw()+  theme(legend.position=c(.83,.2))+
  xlab("Carapace width (mm)")+
  ylab("Weight (g)")
dev.off()



all_wt_m<-filter(survDAT,!is.na(WEIGHT)&SEX==1,GEAR_TEMPERATURE<5,AKFIN_SURVEY_YEAR>1982,SHELL_CONDITION<=2)
ggplot()+
  geom_point(data=all_wt_m,
             aes(x=WIDTH,y=(WEIGHT),col=(GEAR_TEMPERATURE)),alpha=.5)+
  facet_wrap(~AKFIN_SURVEY_YEAR)+
  scale_colour_distiller(palette="Spectral")+
  labs(colour="Temp")+
  geom_line(data=wt_at_s,aes(x=Size,y=Weight*1000))+
  theme_bw()+  theme(legend.position=c(.83,.2))+
  xlab("Carapace width (mm)")+
  ylab("Weight (g)")


#=================================
# calculate calories required
#=================================
# Kleiber's law: metabolic rate = mass^0.75
# crab from foyle were 85-95 mm (250-370g)
## 'cal_intake' comes from Foyle et al. figure 7

keep_cals2<-matrix(nrow=length(mid_pts),ncol=41)
use_tk_size<-seq(27.5,132.5,5)
for(x in 1:length(mid_pts))
{
  take_size<-mid_pts[x]
  in_temp<-filter(temp_tot,SEX==1,mid_pt==take_size,AKFIN_SURVEY_YEAR>1979)
  in_temp[(in_temp=="NaN")]<-NA
  cal_intake<-2.2*exp((-(in_temp$temp_occ-5.2)^2)/30.7)
  const<-cal_intake/(300^0.75)
  adj_cal_intake<-const*wt_at_s[match(take_size,wt_at_s[,2]),2]^0.75
  in_abund<-filter(temp_kod,mid_pt==take_size)
  keep_cals2[x,]<-in_abund$numbers*adj_cal_intake
}

in_temps<-seq(0,18,3)
foyle_dat<-2.2*exp((-(in_temps-5.2)^2)/30.7)

bar_dat<-data.frame(temperature=in_temps,
           kCal = foyle_dat)

kcal_dat<-ggplot(bar_dat,aes(x=as.factor(temperature),y=kCal))+
  geom_bar(stat='identity')+
  theme_bw()+
  labs(y=expression ("kCal kg "~crab^-1 ~day^-1))+
  xlab("Temperature (degrees C)")

rownames(keep_cals2)<-mid_pts
colnames(keep_cals2)<-seq(1980,2021)[-41]

gg_dat<-melt(keep_cals2)
colnames(gg_dat)<-c("Size","Year","kCal")
in_dat<-filter(gg_dat,Size>32&Size<95)

#==need colorramppalette here

cols <- colorRampPalette(brewer.pal(12, "RdYlBu"))
my_pal <- rev(cols(length(unique(in_dat$Size))))

kcal_plot<-ggplot(in_dat,aes(x=Year,y=kCal,fill=as.factor(Size)))+
  geom_area()+
  theme_bw()+
  theme(legend.position=c(.3,.75),
        legend.direction="vertical",legend.key.width = unit(.5, 'cm'),)+
  scale_fill_manual(values=my_pal,name = "Carapace width (mm)",)+
  ylab(expression ("kCal" ~day^-1))+
  guides(fill=guide_legend(ncol=3))
  
plot_1<-kcal_dat / kcal_plot+ plot_layout(heights=c(.3,1)) + plot_annotation(tag_levels='a')

#==have to run figure 1 first to get 'tmp'

png('plots/figure_3.png',height=6,width=13.5,res=350,units='in')
(plot_1 | wt_wid | tmp) + plot_annotation(tag_levels='a')
dev.off()

##+=======================
## plot the different kCal curves by size
#============
keep_cals3<-matrix(nrow=length(use_tk_size),ncol=6)
in_temps<-c(0,3,6,9,12,15)
for(x in 1:length(in_temps))
{
  cal_intake<-2.2*exp((-(in_temps[x]-5.2)^2)/30.7)
  const<-cal_intake/(300^0.75)
  keep_cals3[,x]<-const*wt_at_s[,2]^0.75
}
colnames(keep_cals3)<-in_temps
rownames(keep_cals3)<-wt_at_s[,2]

ugh<-melt(keep_cals3)
colnames(ugh)<-c("Size","Temperature","kCal")

png('plots/kcal_reqs_temp.png',height=6,width=6,res=350,units='in')

ggplot(data=ugh,aes(x=Size,y=kCal,group=Temperature,col=as.factor(Temperature)))+
  geom_line(lwd=2)+
  theme_bw()+
  theme(legend.position=c(.2,.8))+labs(col="Temperature")+xlab("Size (mm carapace width)")
dev.off()