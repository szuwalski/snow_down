#==pull in the simulated stuff, plot to see changes
#==run after running 'operating model.R'
library(PBSmodelling)
library(gridExtra)
library(ggplot2)
library(data.table)
OM<-c("OM_m_vary_q_vary")
EM<-c("EM_m_vary_q_vary","EM_m_vary")
cvs<-c(0.01,0.1,0.3)
model_drv<-expand.grid(OM,EM,cvs)

nEM<-c("M & q vary","M vary")
ncvs<-c('cv = 0.01','cv = 0.1','cv = 0.3')
model_names<-expand.grid(nEM,ncvs)

nsims<-10
all_output<-list(list())
relative_errors_m<-list(list())
relative_errors_m_mat<-list(list())
for(x in 1:nrow(model_drv))
{
  scenario_output<-list(list())
  for(y in 1:(nsims))
  {

input<-paste("simulation/",paste(model_drv[x,1],"_",model_drv[x,2],"_cv_",model_drv[x,3],sep=''),"/",y,"/snow_down.DAT",sep='')
temp<-list()
temp$styr <- scan(input,skip=2,n=1,quiet=T)
temp$endyr <- scan(input,skip=4,n=1,quiet=T)
temp$year_n <- scan(input,skip=6,n=1,quiet=T)
temp$years <- scan(input,skip=8,n=temp$year_n,quiet=T)
temp$size_n <- scan(input,skip=10,n=1,quiet=T)
temp$sizes <- scan(input,skip=12,n=temp$size_n,quiet=T)

temp$imm_n_at_size_obs <- matrix(scan(input,skip=18,n=temp$year_n*temp$size_n,quiet=T),ncol=temp$size_n,byrow=T)
temp$mat_n_at_size_obs <- matrix(scan(input,skip=52,n=temp$year_n*temp$size_n,quiet=T),ncol=temp$size_n,byrow=T)

temp$prob_term_molt <- matrix(scan(input,skip=86,n=temp$year_n*temp$size_n,quiet=T),ncol=temp$size_n,byrow=T)
temp$size_trans <- (matrix(scan(input,skip=120,n=temp$size_n*temp$size_n,quiet=T),ncol=temp$size_n,byrow=T))

temp$imm_cv <- scan(input,skip=134,n=temp$year_n,quiet=T)
temp$mat_cv <- scan(input,skip=136,n=temp$year_n,quiet=T)

#==========
# REP Output
#===============
rep_file<-paste("simulation/",paste(model_drv[x,1],"_",model_drv[x,2],"_cv_",model_drv[x,3],sep=''),"/",y,"/snow_down.REP",sep='')
temp2<-readList(rep_file)
temp<-c(temp,temp2)

#==================
# Relative errors in m an q
###################
nat_m_in<-exp(params$log_m_mu[1]+params$nat_m)
nat_m_mat_in<-exp(params$log_m_mu[2]+params$nat_m_mat)

#plot(temp2$'mature natural mortality'[,1],type='l',ylim=c(0,3))
#lines(nat_m_mat_in,lty=2)

temp$re_est_m<-(temp$'natural mortality'[,1]-nat_m_in)/nat_m_in
temp$re_est_m_mat<-(temp$'mature natural mortality'[,1]-nat_m_mat_in)/nat_m_mat_in

scenario_output[[y]]<-temp

  }
  all_output[[x]]<-scenario_output
}


sim_out_dat<-NULL
for(x in 1:length(all_output))
{
  temp_df<-NULL
  for(y in 1:length(all_output[[1]]))
  {
   temp_df1<-data.frame(year=all_output[[x]][[y]]$years,
                       imm_mort=all_output[[x]][[y]]$"natural mortality"[,1],
                       mat_mort=all_output[[x]][[y]]$"mature natural mortality"[,1],
                       imm_sel=all_output[[x]][[y]]$"survey selectivity"[,5],
                       mat_sel=all_output[[x]][[y]]$"mature survey selectivity"[,5],
                       recruits=exp(all_output[[x]][[y]]$"log recruits"),
                       sim_n=y,
                       sim_name=x,
                       true_imm_mort=exp(params$log_m_mu[1]+params$nat_m),
                       true_mat_mort=exp(params$log_m_mu[2]+params$nat_m_mat),
                       true_imm_sel=params$selectivity[5]+params$q,
                       true_mat_sel=params$selectivity[5]+params$q_mat,
                       true_recruit=exp(params$log_recruits),
                       obs_error=model_names[x,2],
                       est_model=model_names[x,1])
   temp_df<-rbind(temp_df,temp_df1)
   
  }
  sim_out_dat<-rbind(sim_out_dat,temp_df)
}


im<-ggplot(sim_out_dat)+
  geom_line(aes(x=year,y=imm_mort,group=sim_n),alpha=.2)+
  geom_line(aes(x=year,y=true_imm_mort),col=2,lwd=1.15,alpha=.8)+
  facet_grid(rows=vars(obs_error),cols=vars(est_model))+
  theme_bw()+ylab("Immature mortality")

mm<-ggplot(sim_out_dat)+
  geom_line(aes(x=year,y=mat_mort,group=sim_n),alpha=.2)+
  geom_line(aes(x=year,y=true_mat_mort),col=2,lwd=1.15,alpha=.8)+
  facet_grid(rows=vars(obs_error),cols=vars(est_model))+
  theme_bw()+ylab("Mature mortality")

png("plots/sim_mortality.png",height=6,width=8,res=350,units='in') 
grid.arrange(im,mm,ncol=2)
dev.off()

is<-ggplot(sim_out_dat)+
  geom_line(aes(x=year,y=imm_sel,group=sim_n),alpha=.2)+
  geom_line(aes(x=year,y=true_imm_sel),col=2,lwd=1.15,alpha=.8)+
  facet_grid(rows=vars(obs_error),cols=vars(est_model))+
  theme_bw()+ylab("Immature catchability")

ms<-ggplot(sim_out_dat)+
  geom_line(aes(x=year,y=mat_sel,group=sim_n),alpha=.2)+
  geom_line(aes(x=year,y=true_mat_sel),col=2,lwd=1.15,alpha=.8)+
  facet_grid(rows=vars(obs_error),cols=vars(est_model))+
  theme_bw()+ylab("Mature catchability")

png("plots/sim_catchability.png",height=6,width=8,res=350,units='in') 
grid.arrange(is,ms,ncol=2)
dev.off()

png("plots/sim_recruitment.png",height=6,width=8,res=350,units='in') 
ggplot(sim_out_dat)+
  geom_line(aes(x=year,y=recruits,group=sim_n),alpha=.2)+
  geom_line(aes(x=year,y=true_recruit),col=2,lwd=1.15,alpha=.8)+
  facet_grid(rows=vars(obs_error),cols=vars(est_model))+
  theme_bw()+ylab("Recruitment")
dev.off()


dt <- data.table(sim_out_dat)
dtCor_imm <- dt[, .(mCor = cor(imm_mort,true_imm_mort)), by=sim_name]
dtCor_mat <- dt[, .(mCor = cor(mat_mort,true_mat_mort)), by=sim_name]
dtCor_imm_s <- dt[, .(mCor = cor(imm_sel,true_imm_sel)), by=sim_name]
dtCor_mat_s <- dt[, .(mCor = cor(mat_sel,true_mat_sel)), by=sim_name]
dtCor_rec <- dt[, .(mCor = cor(recruits,true_recruit)), by=sim_name]

cor_out<-round(cbind(dtCor_imm[,2],dtCor_mat[,2],dtCor_imm_s[,2],dtCor_mat_s[,2],dtCor_rec[,2]),2)
colnames(cor_out)<-c("Immature mortality","Mature mortality","Immature catchability","Mature catchability","Recruitment")


