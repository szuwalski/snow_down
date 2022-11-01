#==sensitivity tests for the priors and weightings
#==axes for sensitivity include: size comp weight, prior on M, sigmas on M, smoothness penalties, estimate mean m vs. prior
library(PBSmodelling)
#==define axes
comp_wt<-c(25,50,100)
log_mu_m<-c(-1.6,-1.2,-.8)
sigma_m<-c(0.01,0.1,0.2)
smooth_pen<-c(0.001,0.01,0.1,0.5,1)
combos<-expand.grid(comp_wt,log_mu_m,sigma_m,smooth_pen)

#==read in DAT file
orig_wd<-getwd()
dat_file<-readLines("models/model_vary_m/snow_down.dat")
for(y in 1:nrow(combos))
{
  setwd(orig_wd)
  use_dat<-dat_file
  wrk_drv<-paste("sensitivity/sensitivities/",paste(c(combos[y,]),collapse="_"),sep="") 
  dir.create(wrk_drv)
  
  #==modify DAT file
  use_dat[grep("mat_eff_samp",use_dat)+1]<-combos[y,1]
  use_dat[grep("imm_eff_samp",use_dat)+1]<-combos[y,1]
  use_dat[grep("log_m_mu_prior",use_dat)+1]<-paste(rep(combos[y,2],2),collapse=" ")
  use_dat[grep("sigma_m_mu",use_dat)+1]<-paste(rep(combos[y,3],2),collapse=" ")  
  use_dat[grep("smooth_q_weight",use_dat)+1]<-combos[y,4]  
  use_dat[grep("smooth_m_weight",use_dat)+1]<-combos[y,4]  
  
  write(use_dat,file=paste(wrk_drv,"/snow_down.dat",sep=''))
  file.copy(from="models/model_vary_m/snow_down.exe",to=wrk_drv)
  file.copy(from="models/model_vary_m/snow_down.pin",to=wrk_drv)
  setwd(wrk_drv)
  shell("snow_down.exe")
}
#==this is annoying..
setwd(orig_wd)

outs<-list(list())
#==read output, check for .COR file, check gradient, etc.
for(y in 1:nrow(combos))
{
  use_dat<-dat_file
  wrk_drv<-paste("sensitivity/sensitivities/",paste(c(combos[y,]),collapse="_"),sep="") 
  rep_file_name<-paste(wrk_drv,"/snow_down.rep",sep='')
  outs[[y]]<-readList(rep_file_name)
  rep_file_name<-paste(wrk_drv,"/snow_down.rep",sep='')
  outs[[y]]$cor_file<-file.exists(paste(wrk_drv,"/snow_down.cor",sep=''))
}
#==read in likelihoods for data sources (that are comparable)
out_like<-matrix(ncol=length(outs),nrow=13)
for(x in 1:length(outs))
  out_like[,x]<-outs[[x]]$likelihoods

rownames(out_like)<-c("Total","Imm_num","Mat_num","Imm_size","Mat_size","Imm_mort","Mat_mort","Imm_mort_mu","Mat_mort_mu",
                      "smooth_q","smooth_m",'','')
hist(out_like[2,],breaks=seq(0,1000,1))
hist(out_like[3,],breaks=seq(0,100,.5))


#==identify replicates that did not converge or don't fit the data?
good_mat<-which(out_like[2,]<10 & out_like[3,]<2)
good_mat<-seq(1,nrow(combos))
in_col<-as.factor(combos[good_mat,4])
yr_in<-seq(1989,2021)
#==plot estimated M and q
png("plots/sens_m_q_2.png",height=8,width=8,res=350,units='in') 

par(mfrow=c(2,2),mar=c(.1,.1,.3,.1),oma=c(4,4,1,4))
plot(outs[[1]]$`natural mortality`[,1]~yr_in,type='l',col='grey',ylim=c(0,2),xaxt='n',las=1)
for(x in 1:length(outs))
{
  if(outs[[x]]$cor_file & !is.na(match(x,good_mat)))
    lines(outs[[x]]$`natural mortality`[,1]~yr_in,col=in_col[x],)
}
mtext(side=2,line=2.5,"Mortality")
mtext(side=3,line=0.05,"Immature")
legend('topleft',bty='n',col=seq(1,4),pch=15,legend=smooth_pen,title="Smoothness on M")

plot(outs[[1]]$`mature natural mortality`[,1]~yr_in,type='l',col='grey',ylim=c(0,2),yaxt='n',xaxt='n')
for(x in 1:length(outs))
{
  if(outs[[x]]$cor_file & !is.na(match(x,good_mat)))
  lines(outs[[x]]$`mature natural mortality`[,1]~yr_in,col=in_col[x])
}
mtext(side=3,line=0.05,"Mature")

plot(outs[[1]]$`survey selectivity`[,4]~yr_in,type='l',col='grey',ylim=c(0,1),las=1)
for(x in 1:length(outs))
{
  if(outs[[x]]$cor_file & !is.na(match(x,good_mat)))
    lines(outs[[x]]$`survey selectivity`[,4]~yr_in,col=in_col[x])
}
mtext(side=2,line=2.5,"Catchability")
plot(outs[[1]]$`mature survey selectivity`[,4]~yr_in,type='l',col='grey',ylim=c(0,1),yaxt='n')
for(x in 1:length(outs))
{
  if(outs[[x]]$cor_file & !is.na(match(x,good_mat)))
    lines(outs[[x]]$`mature survey selectivity`[,4]~yr_in,col=in_col[x])
}
mtext(outer=T,side=1,"Year",line=2.5)
dev.off()

#=============================
# plot fits to the data
#=============================
png("plots/sens_mod_fits.png",height=8,width=8,res=350,units='in') 

div_n<-1000000000
use_col<-"#3366FF"
par(mfrow=c(2,2),mar=c(.1,.1,.3,.1),oma=c(4,4,1,4))
plot_mmb<-outs[[1]]$mat_numbers_obs/div_n
plot_mmb<-outs[[1]]$mat_n_obs/div_n
plot_mmb[plot_mmb==0]<-NA

plot(plot_mmb~years,
     pch=20,ylab="Mature male abundance",
     xlab="Year",las=1,ylim=c(0,2),xaxt='n')

for(j in 1:length(years))
{
  segments(x0=years[j],x1=years[j],
           y0=plot_mmb[j] /  exp(1.96*sqrt(log(1+mat_cv[j]^2))),
           y1=plot_mmb[j] *  exp(1.96*sqrt(log(1+mat_cv[j]^2))))
}

for(x in 1:length(outs))
{
  if(outs[[x]]$cor_file & !is.na(match(x,good_mat)))
  lines(outs[[x]]$mat_numbers_pred/div_n ~ years,lwd=2,col=in_col[x])
}
legend('topleft',bty='n',"Mature males")

input_box_obs<-sweep(mat_n_at_size_obs, 1, outs[[1]]$mat_numbers_obs, FUN="/")[4:31,]
input_box_obs<-sweep(mat_n_at_size_obs, 1, outs[[1]]$mat_n_obs, FUN="/")[4:31,]

boxplot(input_box_obs,yaxt='n',xaxt='n')
axis(side=4,las=1)
for(x in 1:length(outs))
{
  if(outs[[x]]$cor_file & !is.na(match(x,good_mat)))
  lines(apply(outs[[x]]$'mature numbers at size',2,median),lwd=2,col=in_col[x])
}

plot_mmb<-outs[[1]]$imm_numbers_obs/div_n
plot_mmb<-outs[[1]]$imm_n_obs/div_n
plot_mmb[plot_mmb==0]<-NA

plot(plot_mmb~years,
     pch=20,ylab="Mature male abundance",
     xlab="Year",las=1,ylim=c(0,5))

for(j in 1:length(years))
{
  segments(x0=years[j],x1=years[j],
           y0=plot_mmb[j] /  exp(1.96*sqrt(log(1+imm_cv[j]^2))),
           y1=plot_mmb[j] *  exp(1.96*sqrt(log(1+imm_cv[j]^2))))
}
for(x in 1:length(outs))
{
  if(outs[[x]]$cor_file & !is.na(match(x,good_mat)))
  lines(outs[[x]]$imm_numbers_pred/div_n ~ years,lwd=2,col=in_col[x])
}
legend('topleft',bty='n',"Immature males")
mtext(side=2,line=2.5,"Abundance (billions)",outer=T)
mtext(side=1,line=2.5,"Year")

input_box_obs<-sweep(imm_n_at_size_obs, 1, outs[[1]]$imm_numbers_obs, FUN="/")[4:31,]
input_box_obs<-sweep(imm_n_at_size_obs, 1, outs[[1]]$imm_n_obs, FUN="/")[4:31,]
boxplot(input_box_obs,yaxt='n',names=sizes)
legend('topright',bty='n',col=seq(1,4),pch=15,legend=smooth_pen,title="Smoothness on M")

axis(side=4,las=1)
for(x in 1:length(outs))
{
  if(outs[[x]]$cor_file & !is.na(match(x,good_mat)))
  lines(apply(outs[[x]]$'immature numbers at size', 2,median),lwd=2,col=in_col[x])
}
mtext(side=1,line=2.5,"Carapace width (mm)")
mtext(side=4,line=2.5,outer=T,"Proportion at size")

dev.off()
#++++++++++++++++++++++++++++++++++++++
# do only over smoothness
#++++++++++++++++++++++++++++++++++++++
#==process: 
#==define axes
smooth_pen<-c(0.0001,0.001,0.01,0.1,0.25,0.5,1,5,10,1000)

#==read in DAT file
orig_wd<-getwd()
dat_file<-readLines("models/model_vary_m/snow_down.dat")
for(y in 1:length(smooth_pen))
{
  setwd(orig_wd)
  use_dat<-dat_file
  wrk_drv<-paste("sensitivity/sensitivities/smooth_pen_",smooth_pen[y],sep="") 
  dir.create(wrk_drv)
  
  #==modify DAT file
  use_dat[grep("smooth_q_weight",use_dat)+1]<-smooth_pen[y] 
  use_dat[grep("smooth_m_weight",use_dat)+1]<-smooth_pen[y]
  
  write(use_dat,file=paste(wrk_drv,"/snow_down.dat",sep=''))
  file.copy(from="models/model_vary_m/snow_down.exe",to=wrk_drv)
  file.copy(from="models/model_vary_m/snow_down.pin",to=wrk_drv)
  setwd(wrk_drv)
  shell("snow_down.exe")
}
#==this is annoying..
setwd(orig_wd)

outs2<-list(list())
#==read output, check for .COR file, check gradient, etc.
for(y in 1:length(smooth_pen))
{
  use_dat<-dat_file
  wrk_drv<-paste("sensitivity/sensitivities_2/smooth_pen_",smooth_pen[y],sep="")
  rep_file_name<-paste(wrk_drv,"/snow_down.rep",sep='')
  outs2[[y]]<-readList(rep_file_name)
  rep_file_name<-paste(wrk_drv,"/snow_down.rep",sep='')
  outs2[[y]]$cor_file<-file.exists(paste(wrk_drv,"/snow_down.cor",sep=''))
}
#==read in likelihoods for data sources (that are comparable)
out_like2<-matrix(ncol=length(outs2),nrow=13)
for(x in 1:length(outs2))
  out_like2[,x]<-outs2[[x]]$likelihoods

rownames(out_like2)<-c("Total","Imm_num","Mat_num","Imm_size","Mat_size","Imm_mort","Mat_mort","Imm_mort_mu","Mat_mort_mu",
                      "smooth_q","smooth_m","","")

tot_like<-apply(out_like2[2:8,],2,sum)
plot(tot_like~log(smooth_pen))


