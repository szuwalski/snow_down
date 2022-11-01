#==========
# Input
#===============
input<-"C:/snow_down/models/simulation_testing/estimation_model/3_m_vary_big/snow_down.DAT"
styr <- scan(input,skip=2,n=1,quiet=T)
endyr <- scan(input,skip=4,n=1,quiet=T)
year_n <- scan(input,skip=6,n=1,quiet=T)
years <- scan(input,skip=8,n=year_n,quiet=T)
size_n <- scan(input,skip=10,n=1,quiet=T)
sizes <- scan(input,skip=12,n=size_n,quiet=T)

imm_n_at_size_obs <- matrix(scan(input,skip=14,n=year_n*size_n,quiet=T),ncol=size_n,byrow=T)
mat_n_at_size_obs <- matrix(scan(input,skip=48,n=year_n*size_n,quiet=T),ncol=size_n,byrow=T)

prob_term_molt <- matrix(scan(input,skip=82,n=year_n*size_n,quiet=T),ncol=size_n,byrow=T)
size_trans <- (matrix(scan(input,skip=116,n=size_n*size_n,quiet=T),ncol=size_n,byrow=T))

imm_cv <- scan(input,skip=132,n=year_n,quiet=T)
mat_cv <- scan(input,skip=134,n=year_n,quiet=T)

#==========
# Output
#===============
rep_file<-"C:/snow_down/models/simulation_testing/estimation_model/3_m_vary_big/snow_down.REP"

log_recruits<-scan(rep_file,skip=2,n=year_n,quiet=T)
nat_m_est <- matrix(scan(rep_file,skip=4,n=year_n*size_n,quiet=T),ncol=size_n,byrow=T)
nat_m_mat_est <- matrix(scan(rep_file,skip=116,n=year_n*size_n,quiet=T),ncol=size_n,byrow=T)
imm_n_at_size_pred <- matrix(scan(rep_file,skip=38,n=year_n*size_n,quiet=T),ncol=size_n,byrow=T)
mat_n_at_size_pred <- matrix(scan(rep_file,skip=72,n=year_n*size_n,quiet=T),ncol=size_n,byrow=T)

tot_imm_pred<-scan(rep_file,skip=106,n=year_n,quiet=T)
tot_mat_pred<-scan(rep_file,skip=108,n=year_n,quiet=T)
tot_imm_obs<-scan(rep_file,skip=110,n=year_n,quiet=T)
tot_mat_obs<-scan(rep_file,skip=112,n=year_n,quiet=T)

all_imm<-apply(imm_n_at_size_obs,1,sum)
all_mat<-apply(mat_n_at_size_obs,1,sum)

survey_sel<-matrix(scan(rep_file,skip=114,n=size_n,quiet=T),ncol=size_n,byrow=T)

par_file<-"C:/snow_down/models/simulation_testing/estimation_model/3_m_vary_big/snow_down.PAR"
est_n_imm<-scan(par_file,skip=2,n=size_n,quiet=T)
est_n_mat<-scan(par_file,skip=4,n=size_n,quiet=T)
est_prop_rec<-scan(par_file,skip=80,n=year_n,quiet=T)
est_m_dev <- matrix(scan(par_file,skip=6,n=year_n*size_n,quiet=T),ncol=size_n,byrow=T)
est_m_mat_dev <- matrix(scan(par_file,skip=40,n=year_n*size_n,quiet=T),ncol=size_n,byrow=T)
hist(est_m_mat_dev)
#plot(all_imm+all_mat)
obs_imm_comp<-imm_n_at_size_obs
obs_mat_comp<-mat_n_at_size_obs
for(x in 1:nrow(imm_n_at_size_obs))
{
  obs_imm_comp[x,]<-imm_n_at_size_obs[x,]/tot_imm_obs[x]
  obs_mat_comp[x,]<-mat_n_at_size_obs[x,]/tot_mat_obs[x] 
}

#=========================================================
# model fits and estimated parameters
#+========================================================
plot(exp(params$log_recruits),type='b')
lines(exp(log_recruits),lty=2)

plot(params$log_n_imm,type='l')
lines(est_n_imm)

plot(params$log_n_mat,type='l')
lines(est_n_mat)

plot(params$prop_rec,type='l')
lines(est_prop_rec)

plot(survey_sel[1,],type='l',ylim=c(0,1))
for(x in 2:nrow(survey_sel))
  lines(survey_sel[x,])

#=============================================
# plot size compositions
#=============================================
library(reshape2)
library(ggplot2)
library(ggridges)
library(GGally)
par(mfrow=c(6,6),mar=c(.1,.1,.1,.1),oma=c(4,4,1,1))
for(x in 1:nrow(imm_n_at_size_obs))
{
  if(tot_imm_obs[x]>0)
  {
  plot(imm_n_at_size_obs[x,]/tot_imm_obs[x],type='l',ylim=c(0,max(imm_n_at_size_obs[x,]/tot_imm_obs[x],survey_sel*imm_n_at_size_pred[x,])),
       yaxt='n',xaxt='n',xlim=c(1,ncol(imm_n_at_size_obs)))
  lines(imm_n_at_size_pred[x,],lty=2,col=2)
  legend('topright',bty='n',legend=years[x])
  }
}

par(mfrow=c(6,6),mar=c(.1,.1,.1,.1),oma=c(4,4,1,1))
for(x in 1:nrow(imm_n_at_size_obs))
{ 
  if(tot_mat_obs[x]>0)
  {
  plot(mat_n_at_size_obs[x,]/tot_mat_obs[x],type='l',ylim=c(0,max(mat_n_at_size_obs[x,]/tot_mat_obs[x],survey_sel*mat_n_at_size_pred[x,])),
       yaxt='n',xaxt='n',xlim=c(3,ncol(imm_n_at_size_obs)))
  lines(mat_n_at_size_pred[x,],lty=2,col=2)
  legend('topright',bty='n',legend=years[x])
  }
}


#=============================================
# plot immature and mature numbers by year
#=============================================

#png(paste("plots/figure_2_rough.png",sep=''),height=5,width=8,res=350,units='in') 
div_n<-1000000000
use_col<-"#3366FF"
par(mfrow=c(2,2),mar=c(.1,.1,.3,.1),oma=c(4,4,1,4))
plot_mmb<-tot_mat_obs/div_n
plot_mmb[plot_mmb==0]<-NA

inylim<-c(0,max(plot_mmb*  exp(1.96*sqrt(log(1+mat_cv^2)))))
plot(plot_mmb~years,
     pch=20,ylab="Mature male abundance",
     xlab="Year",las=1,ylim=inylim,xaxt='n')

for(j in 1:length(years))
{
  segments(x0=years[j],x1=years[j],
           y0=plot_mmb[j] /  exp(1.96*sqrt(log(1+mat_cv[j]^2))),
           y1=plot_mmb[j] *  exp(1.96*sqrt(log(1+mat_cv[j]^2))))
}

lines(tot_mat_pred/div_n ~ years,lwd=2,col=use_col)
legend('topleft',bty='n',"Mature males")

boxplot(obs_mat_comp,yaxt='n',xaxt='n')
axis(side=4,las=1)
lines(apply(mat_n_at_size_pred,2,median),lwd=2,col=use_col)

plot_mmb<-tot_imm_obs/div_n
plot_mmb[plot_mmb==0]<-NA
inylim<-c(0,max(plot_mmb*  exp(1.96*sqrt(log(1+mat_cv^2)))))
plot(plot_mmb~years,
     pch=20,ylab="Mature male abundance",
     xlab="Year",las=1,ylim=inylim)

for(j in 1:length(years))
{
  segments(x0=years[j],x1=years[j],
           y0=plot_mmb[j] /  exp(1.96*sqrt(log(1+imm_cv[j]^2))),
           y1=plot_mmb[j] *  exp(1.96*sqrt(log(1+imm_cv[j]^2))))
}
lines(tot_imm_pred/div_n ~ years,lwd=2,col=use_col)
legend('topleft',bty='n',"Immature males")
mtext(side=2,line=2.5,"Abundance (billions)",outer=T)
mtext(side=1,line=2.5,"Year")
boxplot(obs_imm_comp,yaxt='n',names=sizes)
axis(side=4,las=1)
lines(apply(imm_n_at_size_pred, 2,median),lwd=2,col=use_col)
mtext(side=1,line=2.5,"Carapace width (mm)")
mtext(side=4,line=2.5,outer=T,"Proportion at size")
#dev.off()


#################################
# compare esetimated natural mortality to input#
################################################
nat_m_in<-exp(params$nat_m)
nat_m_mat_in<-exp(params$nat_m_mat)

diff_m<-(nat_m_est-nat_m_in)/nat_m_in
diff_m_mat<-(nat_m_mat_est-nat_m_mat_in)/nat_m_mat_in

par(mfrow=c(2,2))
boxplot(diff_m,ylim=c(-2,2),xaxt='n')
abline(h=0,lty=2)
boxplot(diff_m_mat,ylim=c(-2,2),yaxt='n',xaxt='n')
abline(h=0,lty=2)
boxplot(t(diff_m),ylim=c(-2,2))
abline(h=0,lty=2)
boxplot(t(diff_m_mat),ylim=c(-2,2),yaxt='n')
abline(h=0,lty=2)

plot(log(nat_m_est)~log(nat_m_in))
cor.test(log(nat_m_est),log(nat_m_in))


par(mfrow=c(6,6),mar=c(.1,.1,.1,.1),oma=c(4,4,1,1))
for(x in 1:nrow(nat_m_in))
{ 
  plot(nat_m_in[x,],xaxt='n',yaxt='n',type='b',ylim=c(0,1))
  lines(nat_m_est[x,],lty=2)
  legend('top',bty='n',legend=x)
}

par(mfrow=c(6,6),mar=c(.1,.1,.1,.1),oma=c(4,4,1,1))
for(x in 1:nrow(nat_m_in))
{ 
  plot(nat_m_mat_in[x,],xaxt='n',yaxt='n',type='b',ylim=c(0,1))
  lines(nat_m_mat_est[x,],lty=2)
}

#==plot size comps
# imm
rownames(obs_imm_comp)<-years
colnames(obs_imm_comp)<-sizes
df_1<-melt(obs_imm_comp)
colnames(df_1)<-c("Year","Size","Proportion")

df2<-data.frame(pred=apply(imm_n_at_size_pred,2,median),
                    Size=(sizes))

allLevels <- levels(factor(c(df_1$Size,df2$Size)))
df_1$Size <- factor(df_1$Size,levels=(allLevels))
df2$Size <- factor(df2$Size,levels=(allLevels))


imm_size<-ggplot(data=df2,aes(x = factor(Size), y =pred)) + 
  geom_boxplot(data=df_1,aes(x = Size, y =Proportion ),fill='grey') +
  stat_summary(fun.y=mean, geom="line", aes(group=1),lwd=1.5,col='blue',alpha=.8)  + 
  theme_bw()+
  ylab("Proportion")+theme(axis.title=element_text(size=11))+
  xlab("Carapace width (mm)")+
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 

#==mat
rownames(obs_mat_comp)<-years
colnames(obs_mat_comp)<-sizes
df_1<-melt(obs_mat_comp)
colnames(df_1)<-c("Year","Size","Proportion")

df2<-data.frame(pred=apply(mat_n_at_size_pred,2,median),
                Size=(sizes))

allLevels <- levels(factor(c(df_1$Size,df2$Size)))
df_1$Size <- factor(df_1$Size,levels=(allLevels))
df2$Size <- factor(df2$Size,levels=(allLevels))


mat_size<-ggplot(data=df2,aes(x = factor(Size), y =pred)) + 
  geom_boxplot(data=df_1,aes(x = Size, y =Proportion ),fill='grey') +
  stat_summary(fun.y=mean, geom="line", aes(group=1),lwd=1.5,col='blue',alpha=.8)  + 
  theme_bw()+
  ylab("")+
  scale_y_continuous(position = "right",limits=c(0,0.3))+
  xlab("Carapace width (mm)")+
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 

#==fits index
# imm
df_1<-data.frame(pred=tot_imm_pred/div_n,
                 obs=tot_imm_obs/div_n,
                 year=years,
                 ci_dn=(tot_imm_obs/div_n) /  exp(1.96*sqrt(log(1+imm_cv^2))),
                 ci_up=(tot_imm_obs/div_n) *  exp(1.96*sqrt(log(1+imm_cv^2))))

df_1<-df_1[-32,]

imm_abnd<-ggplot(data=df_1)+
  geom_segment(aes(x=year,xend=year,y=ci_dn,yend=ci_up))+
  geom_point(aes(x=year,y=obs))+
  geom_line(aes(x=year,y=pred),lwd=1.5,col='blue',alpha=.8)+
  theme_bw()+
  scale_x_continuous(position = "top",name='IMMATURE')+
  ylab("Abundance (billions)")+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
#+  annotate('text',label="Immature", x=2005,y=5.5,size=7)

# mature
df_1<-data.frame(pred=tot_mat_pred/div_n,
                 obs=tot_mat_obs/div_n,
                 year=years,
                 ci_dn=(tot_mat_obs/div_n) /  exp(1.96*sqrt(log(1+mat_cv^2))),
                 ci_up=(tot_mat_obs/div_n) *  exp(1.96*sqrt(log(1+mat_cv^2))))
df_1<-df_1[-32,]

mat_abnd<-ggplot(data=df_1)+
  geom_segment(aes(x=year,xend=year,y=ci_dn,yend=ci_up))+
  geom_point(aes(x=year,y=obs))+
  geom_line(aes(x=year,y=pred),lwd=1.5,col='blue',alpha=.8)+
  theme_bw()+
  scale_x_continuous(position = "top",name='MATURE')+
  scale_y_continuous(position = "right",limits=c(0,0.5))+
  ylab("")+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 
  

#=======================================
# plot natural mortality
#=======================================
plot_m<-(nat_m_est)
plot_m<-exp(params$nat_m)
colnames(plot_m)<-sizes
rownames(plot_m)<-years
in_dat<-melt(plot_m)
colnames(in_dat)<-c("Year","Size","nat_m")

nat_m_plot<- ggplot(dat=in_dat) 
nat_m_plot <- nat_m_plot + geom_density_ridges(aes(x=Size, y=Year, height = nat_m,
                                               group = Year, 
                                               fill=stat(y),alpha=.9999),stat = "identity",scale=2) +
  # scale_fill_viridis_c()+
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(legend.position = "none",
        axis.title.y=element_blank()) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

#png(paste("plots/figure_3_rough_cv.png",sep=''),height=5,width=3,res=350,units='in') 
#print(nat_m_plot)
#dev.off()


plot_m_mat<-(nat_m_mat_est)
#plot_m_mat<-exp(params$nat_m_mat)
colnames(plot_m_mat)<-sizes
rownames(plot_m_mat)<-years
in_dat_mat<-melt(plot_m_mat)
colnames(in_dat_mat)<-c("Year","Size","nat_m")

nat_m_mat_plot<- ggplot(dat=in_dat_mat) 
nat_m_mat_plot <- nat_m_mat_plot + geom_density_ridges(aes(x=Size, y=Year, height = nat_m,
                                                   group = Year, 
                                                   fill=stat(y),alpha=.9999),stat = "identity",scale=1.5) +
  # scale_fill_viridis_c()+
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(legend.position = "none",
        axis.title.y=element_blank()) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  scale_y_continuous(position = "right")

# why is there not a lot of correlation between mature and immature mortality?
#plot(nat_m_mat_est~nat_m_est,ylim=c(0,2.5),xlim=c(0,2.5))

#png(paste("plots/figure_3_rough_cv.png",sep=''),height=5,width=3,res=350,units='in') 
#print(nat_m_mat_plot)
#dev.off()

layout_mat<-matrix(c(1,2,
              3,4,
              3,4,
              5,6),ncol=2,byrow=T)

#png(paste("plots/figure_2_rough.png",sep=''),height=8,width=5,res=350,units='in') 
grid.arrange(imm_abnd,mat_abnd,
             nat_m_plot,nat_m_mat_plot,
             imm_size,mat_size,layout_matrix=layout_mat)
#dev.off()




png(paste("plots/figure_4_rough.png",sep=''),height=4,width=7,res=350,units='in') 
par(mfrow=c(2,2),mar=c(.1,.1,.1,.1),oma=c(4,4,1,1))
inylim<-c(0,max(plot_m,plot_m_mat))
inylim<-c(0,2.5)
boxplot(plot_m_mat[,1:ncol(plot_m)],las=1,ylim=inylim,xaxt='n')
boxplot(t(plot_m_mat[,1:ncol(plot_m)]),yaxt='n',ylim=inylim,xaxt='n')
boxplot(plot_m[,1:ncol(plot_m)],las=1,ylim=inylim)
mtext(side=2,line=2.5,"Natural mortality")
mtext(side=1,line=2.5,"Carapace width (mm)")
boxplot(t(plot_m[,1:ncol(plot_m)]),yaxt='n',ylim=inylim)
mtext(side=1,line=2.5,"Year")
dev.off()
#lines(apply(t(plot_m),2,median),lty=2,col=2,lwd=3)
#abline(h=0.31,lty=2)
#abline(h=median(plot_m),lty=2)
#abline(h=median(plot_m_mat),lty=2)

############################
#===========================
# build a GAM for immature
#==========================
############################
library(mgcv)
#=======================================
# plot immature natural mortality against covars
#=======================================
which_mat<-'imm'
melt_m<-melt(log(plot_m))
if(which_mat=='mat')
 melt_m<-melt(plot_m_mat)

colnames(melt_m)<-c("Year","Size","Variable")
melt_m$Process<-"Mortality"
big_input<-rbind(filter(gam_covars,Size<95),melt_m)
ugh<-dcast(big_input, Year+Size~Process,value.var='Variable')
ugh$Predation[is.na(ugh$Predation)]<-0
ugh<-ugh[complete.cases(ugh),]

#=========================================================
# make predation and fishery effects relative to abundance
# need to consider the fact that fishery mortality are both mature and immature
#=========================================================
pred_n_size_imm<-sweep(imm_n_at_size_pred, 1, tot_imm_pred, FUN="*")[4:31,2:13]
pred_n_size_mat<-sweep(mat_n_at_size_pred, 1, tot_mat_pred, FUN="*")[4:31,2:13]

ugh_yrs<-unique(ugh$Year)
for(x in 1:length(ugh_yrs))
{
  if(which_mat=='imm')
  {
  mult_ind<-which(ugh$Year==ugh_yrs[x]) 
  ugh$Bycatch[mult_ind]<-scale(ugh$Bycatch[mult_ind]/pred_n_size_imm[x,])
  ugh$Discards[mult_ind]<-scale(ugh$Discards[mult_ind]/pred_n_size_imm[x,])
  ugh$Predation[mult_ind]<-scale(ugh$Predation[mult_ind]/pred_n_size_imm[x,])
  }
  
  if(which_mat=='mat')
  {
  mult_ind<-which(ugh$Year==ugh_yrs[x]) 
  ugh$Bycatch[mult_ind]<-scale(ugh$Bycatch[mult_ind]/pred_n_size_mat[x,])
  ugh$Discards[mult_ind]<-scale(ugh$Discards[mult_ind]/pred_n_size_mat[x,])
  ugh$Predation[mult_ind]<-scale(ugh$Predation[mult_ind]/pred_n_size_mat[x,] )
  }
}
png(paste("plots/pairs_plot.png",sep=''),height=8,width=8,res=350,units='in') 
ggpairs(ugh)+theme_bw()
dev.off()
png(paste("plots/pairs_plot_size.png",sep=''),height=8,width=8,res=350,units='in') 
ggpairs(ugh,aes(colour=as.factor(Size),alpha=.4),
        upper = list(continuous = wrap("cor", size = 2.5)))+theme_bw()
dev.off()


#================================
# predict just to 2018 because 2019 cannot be 'trusted' because of no data in 2020
# also do this so we could see if we could have predicted the decline in 2019?...no lags, so not sure
# how that would work...well, the mortality that is estimated is the mortality after the survey...so it could work
# then try to predict mortality in 2019 and 2020?
#================================
input_dat<-ugh
mod_lm<-lm(data=input_dat,Mortality~Disease+Temperature+Predation+Discards+Bycatch+Cannibalism)
summary(mod_lm)
mod_lm<-lm(data=input_dat,Mortality~(Disease+Temperature+Predation+Discards+Bycatch+Cannibalism)^2)
summary(mod_lm)
mod_base<-gam(data=input_dat,Mortality~s(Disease)+s(Temperature)+s(Predation)+s(Discards)+s(Bycatch) +s(Cannibalism))
mod_base2<-gam(data=input_dat,Mortality~Disease+Temperature+Predation+Discards+Bycatch+Cannibalism)
mod<-gam(data=input_dat,Mortality~Disease+Temperature+Predation+Discards+Bycatch+Cannibalism+
            s(Size,Disease) +s(Size,Temperature) +
           s(Size,Bycatch)+  s(Size,Cannibalism)+s(Temperature,Disease)  +
           s(Temperature,Bycatch)+s(Temperature,Cannibalism)+s(Temperature,Discards)+
           s(Temperature, Predation)+s(Disease,Predation),method='REML' )
mod2<-gam(data=input_dat,Mortality~s(Disease,k=3)+s(Temperature,k=3)+s(Predation,k=3)+s(Discards,k=3)+s(Bycatch,k=3) +s(Cannibalism,k=3)+
           s(Size,Disease) +s(Size,Temperature) +
           s(Size,Bycatch)+  s(Size,Cannibalism)+s(Temperature,Disease)  +
           s(Temperature,Bycatch)+s(Temperature,Cannibalism)+s(Temperature,Discards)+
           s(Temperature, Predation),method='REML' )
par(mfrow=c(2,2))
gam.check(mod)
summary(mod_base)
summary(mod_base2)
summary(mod)
summary(mod2)
plot(mod_base,pages=1)
AIC(mod)
plot(mod,pages=1,seWithMean=TRUE,scheme=2,too.far=0,shade=TRUE,shade.col='red',
     all.terms=TRUE,cex=3)
plot(mod2,pages=1,seWithMean=TRUE,scheme=2,too.far=0,shade=TRUE,shade.col='red',
     all.terms=TRUE,cex=3)
#==organize plots
library(RColorBrewer)
in_cont_col<-rev(colorRampPalette(brewer.pal(11, 'RdBu'))(100))
#in_cont_col<-terrain.colors(50)
inscheme<-2
par(mfrow=c(3,2),mar=c(.3,.1,.1,.1),oma=c(4,4,1,4))
plot(mod,select=1,seWithMean=TRUE,scheme=inscheme,too.far=0,shade=TRUE,
     all.terms=TRUE,cex=3,xaxt='n',las=1,hcolors=in_cont_col)
mtext(side=2,line=2.5,"Disease")
plot(mod,select=2,seWithMean=TRUE,scheme=inscheme,too.far=0,shade=TRUE,
     all.terms=TRUE,cex=3,xaxt='n',yaxt='n,las=1',hcolors=in_cont_col)
axis(side=4,las=1)
mtext(side=4,line=2.5,"Disease")
plot(mod,select=5,seWithMean=TRUE,scheme=inscheme,too.far=0,shade=TRUE,
     all.terms=TRUE,cex=3,xaxt='n',las=1,hcolors=in_cont_col)
mtext(side=2,line=2.5,"Discards")
plot(mod,select=3,seWithMean=TRUE,scheme=inscheme,too.far=0,shade=TRUE,
     all.terms=TRUE,cex=3,xaxt='n',yaxt='n',las=1,hcolors=in_cont_col)
mtext(side=4,line=2.5,"Temperature")
axis(side=4,las=1)
plot(mod,select=6,seWithMean=TRUE,scheme=inscheme,too.far=0,shade=TRUE,
     all.terms=TRUE,cex=3,las=1,hcolors=in_cont_col)
mtext(side=2,line=2.5,"Bycatch")
mtext(side=1,line=2.5,"Temperature")
plot(mod,select=4,seWithMean=TRUE,scheme=inscheme,too.far=0,shade=TRUE,
     all.terms=TRUE,cex=3,yaxt='n',las=1,hcolors=in_cont_col)
mtext(side=4,line=2.5,"Bycatch")
axis(side=4,las=1)
mtext(side=1,line=2.5,"Carapace width (mm)")

preds<-matrix(predict(mod),ncol=length(unique(input_dat$Size)),byrow=T)

# plot predictions from the GAM against observed
png(paste("plots/figure_6_rough.png",sep=''),height=8,width=8,res=350,units='in') 
use_obs<-log(plot_m)
par(mfrow=c(2,2),mar=c(.1,.1,.1,.1),oma=c(4,4,1,1))
boxplot(use_obs,las=1,ylim=c(min(use_obs,preds),max(use_obs,preds)),xaxt='n')
legend("topleft",bty='n',"Observed")
boxplot(t(use_obs)[,4:31],yaxt='n',ylim=c(min(use_obs,preds),max(use_obs,preds)),xaxt='n')
boxplot(preds,las=1,ylim=c(min(use_obs,preds),max(use_obs,preds)),names=unique(input_dat$Size))
legend("topleft",bty='n',"Predicted")
mtext(side=1,line=2.5,"Carapace width (mm)")
boxplot(t(preds),yaxt='n',ylim=c(min(use_obs,preds),max(use_obs,preds)),names=seq(1992,2019),las=2)
mtext(side=1,line=2.75,"Year")
mtext(side=2,line=2.5,"Natural mortality",outer=T)
dev.off()


nerp<-data.frame(size=rep(unique(input_dat$Size),2),Mortality=c(apply(use_obs[4:31,],2,median),apply(preds,2,median)),
           ci_90=c(apply(use_obs[4:31,],2,quantile,0.9),apply(preds,2,quantile,0.9)),
           ci_10=c(apply(use_obs[4:31,],2,quantile,0.1),apply(preds,2,quantile,0.1)),
           data=c(rep('Observed',length(unique(input_dat$Size))),rep("Predicted",length(unique(input_dat$Size)))))
                   

size_m<-ggplot()+
  geom_ribbon(data=nerp,aes(x=size,ymin=ci_10,ymax=ci_90,fill=data),alpha=.3)+
  geom_line(data=nerp,aes(x=size,y=Mortality,col=data),lwd=2)+
  theme_bw()+ theme(legend.position=c(.7,.75),legend.title=element_blank())+
  ylim(0,2.5)

derp<-data.frame(year=rep(seq(1992,2019),2),medians=c(apply(t(use_obs[4:31,]),2,median),apply(t(preds),2,median)),
                 ci_90=c(apply(t(use_obs[4:31,]),2,quantile,0.9),apply(t(preds),2,quantile,0.9)),
                 ci_10=c(apply(t(use_obs[4:31,]),2,quantile,0.1),apply(t(preds),2,quantile,0.1)),
                 data=c(rep('Observed',length(seq(1992,2019))),rep("Predicted",length(seq(1992,2019)))) )


year_m<-ggplot()+
  geom_ribbon(data=derp,aes(x=year,ymin=ci_10,ymax=ci_90,fill=data),alpha=.3)+
  geom_line(data=derp,aes(x=year,y=medians,col=data),lwd=2)+
  theme_bw()+
  theme(legend.position='none',
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+
  ylim(0,2.52)

png(paste("plots/figure_7_rough.png",sep=''),height=4,width=8,res=350,units='in') 
grid.arrange(size_m,year_m,ncol=2)
dev.off()
cor.test(plot_m[4:31,],preds)
      plot(plot_m[4:31,],preds)   
      abline(a=0,b=1)
#===========================
# plot predicted vs. obs
#===========================
in_plot_m<-plot_m[-c(1,2,3,nrow(plot_m)),]
in_yr<-seq(1992,2019)
plot(-10,xlim=c(1,10),ylim=c(1992,2020),las=1,ylab='')
counter<-1992
for(x in 1:length(in_yr))
{
  input<-in_plot_m[x,]/2+counter
  polygon(x=c(seq(1,10),seq(10,1)),y=c(rep(counter,length(input)),input),col='#aadd2222')
  input2<-preds[x,]/2+counter 
  polygon(x=c(seq(1,10),seq(10,1)),y=c(rep(counter,length(input2)),input2),col='#3366FF88',border=NA)
  counter<-counter+1
  
}

#========================================================
# build up model to show how much each variable explains
#========================================================

mod_1<-gam(data=ugh,Mortality~ s(Size))
summary(mod_1)

mod_2<-gam(data=ugh,Mortality~ s(Size) + s(Size,Disease))
summary(mod_2)

mod_3<-gam(data=ugh,Mortality~ s(Size) + s(Size,Temperature))
summary(mod_3)

mod_4<-gam(data=ugh,Mortality~ s(Size) + s(Size,Predation))
summary(mod_4)

mod_5<-gam(data=ugh,Mortality~ s(Size) + s(Size,Bycatch))
summary(mod_5)

mod_6<-gam(data=ugh,Mortality~ s(Size) + s(Size,Discards))
summary(mod_6)


#===========================
# build a randomforest
#==========================
library(randomForest)
rf_input<-ugh[complete.cases(ugh),-1]
model_1 = randomForest(Mortality~., data = rf_input, importance = TRUE,mtry=4)
varImpPlot(model_1)
