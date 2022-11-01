# =operating model to test estimation ability of the model

library(reshape2)
library(ggplot2)
library(ggridges)
library(doParallel)
library(parallel)
library(foreach)


# Detect the number of available cores and create cluster
 cl <- parallel::makeCluster(detectCores())
 doParallel::registerDoParallel(cl)


pop_dy<-function(params,vary_m=FALSE)
{
 nat_m_in<- rep(params$log_m_mu[1],length(params$yrs))
 nat_m_mat_in<- rep(params$log_m_mu[2],length(params$yrs)) 
 if(vary_m==TRUE)
 {
   nat_m_in<- params$log_m_mu[1]+params$nat_m
   nat_m_mat_in<-params$log_m_mu[2]+params$nat_m_mat
 }
 
 imm_n_size<-matrix(ncol=length(params$sizes),nrow=length(params$yrs))  
 mat_n_size<-matrix(ncol=length(params$sizes),nrow=length(params$yrs))  
 
 imm_n_size[1,]<-exp(params$log_n_imm)
 mat_n_size[1,]<-exp(params$log_n_mat) 
 
 for(x in 1:(length(params$yrs)-1))
 {
  # mortality
  temp_imm<-imm_n_size[x,]*exp(-0.5*exp(nat_m_in[x]) )
  temp_mat<-mat_n_size[x,]*exp(-0.5*exp(nat_m_mat_in[x]))   
  
  # growth
  trans_imm <- params$size_trans%*% temp_imm
  trans_imm[1] <- trans_imm[1] + exp(params$log_recruits[x])*params$prop_rec[x]
  trans_imm[2] <- trans_imm[2] + exp(params$log_recruits[x])*(1-params$prop_rec[x]) 

  imm_n_size[x+1,] <- (trans_imm * (1-params$prop_term_molt[x,])) * exp(-0.5*exp(nat_m_in[x]) )
  mat_n_size[x+1,] <- (trans_imm * (params$prop_term_molt[x,]) + temp_mat) * exp(-0.5*exp(nat_m_mat_in[x]) )   
 }
 
  list(imm_n_size,mat_n_size)
}


#======================================================
# plug parameters from the fitted model in, run it back
#======================================================
styr<-1989
endyr<-2021
params<-list()
params$yrs<-seq(styr,endyr)
params$sizes<-seq(32.5,92.5,5)

#==test against the most complicated model for OM
par_file<-"models/model_vary_q_m/snow_down.par"
dat_file<-"models/model_vary_q_m/snow_down.dat"

params$log_n_imm <- scan(par_file,skip=2,n=length(params$sizes),quiet=T)
params$log_n_mat <- scan(par_file,skip=4,n=length(params$sizes),quiet=T)
params$nat_m <- scan(par_file,skip=6,n=length(params$yrs),quiet=T)
params$nat_m_mat <- scan(par_file,skip=8,n=length(params$yrs),quiet=T)
params$q <- scan(par_file,skip=10,n=length(params$yrs),quiet=T)
params$q_mat <- scan(par_file,skip=12,n=length(params$yrs),quiet=T)
params$log_recruits <- scan(par_file,skip=14,n=length(params$yrs),quiet=T)
params$log_m_mu <- scan(par_file,skip=18,n=2,quiet=T)
params$prop_rec <- scan(par_file,skip=20,n=length(params$yrs),quiet=T)

params$size_trans <- matrix(scan(dat_file,skip=120,n=length(params$sizes)*length(params$sizes),quiet=T),ncol=length(params$sizes),byrow=T)
params$prop_term_molt <- matrix(scan(dat_file,skip=86,n=length(params$sizes)*length(params$yrs),quiet=T),ncol=length(params$sizes),byrow=T)
params$selectivity <- scan(dat_file,skip=142,n=length(params$sizes),quiet=T)

params$imm_cv <- scan(dat_file,skip=134,n=length(params$yrs),quiet=T)
params$mat_cv <- scan(dat_file,skip=136,n=length(params$yrs),quiet=T)

params$mat_cv[params$mat_cv==0]<-mean(params$mat_cv)
params$imm_cv[params$imm_cv==0]<-mean(params$imm_cv)

in_vary_m<-TRUE
outs<-pop_dy(params,vary_m=in_vary_m)


#===========================
# plot the simulated numbers
#===========================
plot_stuff<-FALSE

if(plot_stuff)
{
plot_m<-(outs[[1]])
colnames(plot_m)<-params$sizes
rownames(plot_m)<-params$yrs
in_dat<-melt(plot_m)
colnames(in_dat)<-c("Year","Size","N")

nat_m_plot<- ggplot(dat=in_dat) 
nat_m_plot <- nat_m_plot + geom_density_ridges(aes(x=Size, y=Year, height = N,
                                                   group = Year, 
                                                   fill=stat(y),alpha=.9999),stat = "identity",scale=3) +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(legend.position = "none",
        axis.title.y=element_blank()) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
print(nat_m_plot)
}

###################################################
# simulate over different combinations of OM and EM
###################################################
OM<-c("m_vary_q_vary")
EM<-c("m_vary_q_vary","m_vary")
nsims<-10
in_sds<-c(0.01,0.1,0.3)
orig_drv<-getwd()
for(f in 1:length(in_sds))
for(z in 1:length(OM))
{
#======================================
# Simulate time- and size-varying catchability
#======================================
varying_q<-matrix(0,ncol=length(params$sizes),nrow=length(params$yrs))
varying_q_mat<-matrix(0,ncol=length(params$sizes),nrow=length(params$yrs))
if(OM[z] == "m_vary_q_vary" | OM[z] == "q_vary")
{
   for(x in 1:nrow(varying_q))
   {
     varying_q[x,]<-params$q[x]
     varying_q_mat[x,]<-params$q_mat[x]
   }
}

if(plot_stuff)
{
plot_m<-(varying_q)
colnames(plot_m)<-params$sizes
rownames(plot_m)<-params$yrs
in_dat<-melt(plot_m)
colnames(in_dat)<-c("Year","Size","N")

nat_m_plot<- ggplot(dat=in_dat) 
nat_m_plot <- nat_m_plot + geom_ridgeline(aes(x=Size, y=Year, height = N,
                                                   group = Year, 
                                                   fill=stat(y),alpha=.9999),stat = "identity",scale=3,min_height = -5) +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(legend.position = "none",
        axis.title.y=element_blank()) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
print(nat_m_plot)

}

#============================
# finalize selectivity matrix
#============================
in_selectivity<-sweep(varying_q,2,params$selectivity,FUN="+")
in_selectivity_mat<-sweep(varying_q_mat,2,params$selectivity,FUN="+")

if(plot_stuff)
{
plot(in_selectivity[1,],type='l',ylim=c(0,1))
for(x in 2:nrow(in_selectivity))
  lines(in_selectivity[x,])
}
#=========================================
# extract data & write dat file
#=========================================
for(w in 1:length(EM))
{
  setwd(orig_drv)
  in_sd<-in_sds[f]
  dir.create(paste("simulation/OM_",OM[z],"_EM_",EM[w],"_cv_",in_sd,sep="") )
 
  # set the mod_driv to pull from
  if(EM[w]=='no_vary')
  {
   mod_drv_exe<- paste(orig_drv,"/models/model_base/snow_down.exe",sep='')
   mod_drv_pin<- paste(orig_drv,"/models/model_base/snow_down.pin",sep='')
  }
  
  if(EM[w]=='m_vary')
  {
    mod_drv_exe<- paste(orig_drv,"/models/model_vary_m/snow_down.exe",sep='')
    mod_drv_pin<- paste(orig_drv,"/models/model_vary_m/snow_down.pin",sep='')
  }
  if(EM[w]=='q_vary')
  {
    mod_drv_exe<- paste(orig_drv,"/models/model_vary_q/snow_down.exe",sep='')
    mod_drv_pin<- paste(orig_drv,"/models/model_vary_q/snow_down.pin",sep='')
  }
  if(EM[w]=='m_vary_q_vary')
  {
    mod_drv_exe<- paste(orig_drv,"/models/model_vary_q_m/snow_down.exe",sep='')
    mod_drv_pin<- paste(orig_drv,"/models/model_vary_q_m/snow_down.pin",sep='')  
  }

foreach(y = 1:nsims)%dopar%
 {
  # for(y in 1:nsims)
  # {
  wrk_drv<-paste(orig_drv,"/simulation/OM_",OM[z],"_EM_",EM[w],"_cv_",in_sd,"/",y,sep="") 
  dir.create(wrk_drv)
  setwd(wrk_drv)
  dat_file<-"snow_down.DAT"
  file.create(dat_file)
 
cat("# snow crab decline analysis",file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)
cat("# start year",file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)
cat(min(params$yrs),file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)
cat("# endyr",file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)
cat(max(params$yrs),file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)

cat("# number of years",file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)
cat(length(params$yrs),file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)

cat("# years",file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)
cat(params$yrs,file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)

cat("# number of sizes",file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)
cat(length(params$sizes),file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)

cat("# sizes",file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)
cat(params$sizes,file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)

cat("# immature numbers",file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)
tmp3<-apply(outs[[1]],1,sum)*exp(rnorm(nrow(outs[[1]]),0,in_sd))
cat(tmp3,file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)

cat("# mature numbers",file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)
tmp3<-apply(outs[[2]],1,sum)*exp(rnorm(nrow(outs[[2]]),0,in_sd))
cat(tmp3,file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)

cat("# immature numbers at carapace width",file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)
tmp1<-as.matrix(outs[[1]])*in_selectivity
for(x in 1:nrow(tmp1))
{
  tot_num<-sum(tmp1[x,])
  size_comp<-hist(sample(params$sizes,size=10000,prob=tmp1[x,],replace=TRUE),plot=FALSE,breaks=seq(30,95,5))$density
  cat(tot_num*size_comp,file=dat_file,append=TRUE)
  cat("\n",file=dat_file,append=TRUE)
}

cat("# mature numbers at carapace width",file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)
tmp1<-as.matrix(outs[[2]])*in_selectivity_mat
for(x in 1:nrow(tmp1))
{
  tot_num<-sum(tmp1[x,])
  size_comp<-hist(sample(params$sizes,size=10000,prob=tmp1[x,],replace=TRUE),plot=FALSE,breaks=seq(30,95,5))$density
  cat(tot_num*size_comp,file=dat_file,append=TRUE)
  cat("\n",file=dat_file,append=TRUE)
}

cat("# proportion undergoing terminal molt at carapace width",file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)
for(x in 1:nrow(params$prop_term_molt))
{
  cat(params$prop_term_molt[x,],file=dat_file,append=TRUE)
  cat("\n",file=dat_file,append=TRUE)
}

cat("# size transition matrix",file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)
tmp<-(params$size_trans)
for(x in 1:nrow(tmp))
{
  cat(unlist(tmp[x,]),file=dat_file,append=TRUE)
  cat("\n",file=dat_file,append=TRUE)
}

 cat("# sigma_numbers_imm ",file=dat_file,append=TRUE)
 cat("\n",file=dat_file,append=TRUE)
 cat(rep(in_sd,nrow(outs[[1]])),file=dat_file,append=TRUE)
 cat("\n",file=dat_file,append=TRUE)

 cat("# sigma_numbers_mat ",file=dat_file,append=TRUE)
 cat("\n",file=dat_file,append=TRUE)
 cat(rep(in_sd,nrow(outs[[1]])),file=dat_file,append=TRUE)
 cat("\n",file=dat_file,append=TRUE)

cat("# mat_eff_samp ",file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)
cat(40,file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)

cat("# imm_eff_samp ",file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)
cat(40,file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)

cat("# survey_sel ",file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)
cat(params$selectivity ,file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)

cat("# log_mu_m",file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)
cat(-1.2,-1.2,file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)

cat("# est_m_devs ",file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)
if(EM[w] == 'm_vary' | EM[w]=='m_vary_q_vary')
 cat(1,file=dat_file,append=TRUE)
if(EM[w] == 'no_vary' | EM[w]=='q_vary')
cat(-1,file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)

cat("# est_q_devs ",file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)
if(EM[w] == 'q_vary' | EM[w]=='m_vary_q_vary')
  cat(1,file=dat_file,append=TRUE)
if(EM[w] == 'no_vary' | EM[w]=='m_vary')
  cat(-1,file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)

cat("# est_m_mat_devs ",file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)
if(EM[w] == 'm_vary' | EM[w]=='m_vary_q_vary')
  cat(1,file=dat_file,append=TRUE)
if(EM[w] == 'no_vary' | EM[w]=='q_vary')
  cat(-1,file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)

cat("# est_q_mat_devs ",file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)
if(EM[w] == 'q_vary' | EM[w]=='m_vary_q_vary')
  cat(1,file=dat_file,append=TRUE)
if(EM[w] == 'no_vary' | EM[w]=='m_vary')
  cat(-1,file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)

cat("# est_sigma_m ",file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)
cat(-1,file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)

cat("# sigma_m_mu",file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)
cat(c(0.01,0.01),file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)

cat("# smooth_q_weight ",file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)
cat(.001,file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)

cat("# smooth_m_weight ",file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)
cat(.001,file=dat_file,append=TRUE)
cat("\n",file=dat_file,append=TRUE)

  file.copy(from=mod_drv_exe,to=getwd())
  file.copy(from=mod_drv_pin,to=getwd())

shell("snow_down.exe")
}
}
}
setwd(orig_drv)


