#ifdef DEBUG
  #ifndef __SUNPRO_C
    #include <cfenv>
    #include <cstdlib>
  #endif
#endif
#ifdef DEBUG
  #include <chrono>
#endif
#include <admodel.h>
#ifdef USE_ADMB_CONTRIBS
#include <contrib.h>

#endif
  extern "C"  {
    void ad_boundf(int i);
  }
#include <snow_down.htp>

model_data::model_data(int argc,char * argv[]) : ad_comm(argc,argv)
{
  adstring tmpstring;
  tmpstring=adprogram_name + adstring(".dat");
  if (argc > 1)
  {
    int on=0;
    if ( (on=option_match(argc,argv,"-ind"))>-1)
    {
      if (on>argc-2 || argv[on+1][0] == '-')
      {
        cerr << "Invalid input data command line option"
                " -- ignored" << endl;
      }
      else
      {
        tmpstring = adstring(argv[on+1]);
      }
    }
  }
  global_datafile = new cifstream(tmpstring);
  if (!global_datafile)
  {
    cerr << "Error: Unable to allocate global_datafile in model_data constructor.";
    ad_exit(1);
  }
  if (!(*global_datafile))
  {
    delete global_datafile;
    global_datafile=NULL;
  }
  styr.allocate("styr");
  endyr.allocate("endyr");
  dat_yr.allocate("dat_yr");
  years.allocate(1,dat_yr,"years");
  size_n.allocate("size_n");
  sizes.allocate(1,size_n,"sizes");
  imm_n_obs.allocate(styr,endyr,"imm_n_obs");
  mat_n_obs.allocate(styr,endyr,"mat_n_obs");
  imm_n_size_obs.allocate(styr,endyr,1,size_n,"imm_n_size_obs");
  mat_n_size_obs.allocate(styr,endyr,1,size_n,"mat_n_size_obs");
  prop_term_molt.allocate(styr,endyr,1,size_n,"prop_term_molt");
  size_trans.allocate(1,size_n,1,size_n,"size_trans");
  sigma_numbers_imm.allocate(styr,endyr,"sigma_numbers_imm");
  sigma_numbers_mat.allocate(styr,endyr,"sigma_numbers_mat");
  mat_eff_samp.allocate("mat_eff_samp");
  imm_eff_samp.allocate("imm_eff_samp");
  survey_sel.allocate(1,size_n,"survey_sel");
  log_mu_m_prior.allocate(1,2,"log_mu_m_prior");
  est_m_devs.allocate("est_m_devs");
  est_q_devs.allocate("est_q_devs");
  est_m_mat_devs.allocate("est_m_mat_devs");
  est_q_mat_devs.allocate("est_q_mat_devs");
  est_sigma_m.allocate("est_sigma_m");
  sigma_m_mu.allocate(1,2,"sigma_m_mu");
  smooth_q_weight.allocate("smooth_q_weight");
  smooth_m_weight.allocate("smooth_m_weight");
  est_log_m_mu.allocate("est_log_m_mu");
  est_sigma_q.allocate("est_sigma_q");
cout<<survey_sel<<endl;
cout<<size_trans<<endl;
cout<<sigma_m_mu<<endl;
  if (global_datafile)
  {
    delete global_datafile;
    global_datafile = NULL;
  }
}

model_parameters::model_parameters(int sz,int argc,char * argv[]) : 
 model_data(argc,argv) , function_minimizer(sz)
{
  initializationfunction();
  log_n_imm.allocate(1,size_n,1,30,1,"log_n_imm");
  log_n_mat.allocate(1,size_n,1,30,1,"log_n_mat");
  nat_m_dev.allocate(styr,endyr,-4,4,est_m_devs,"nat_m_dev");
  nat_m_mat_dev.allocate(styr,endyr,-4,4,est_m_mat_devs,"nat_m_mat_dev");
  q_dev.allocate(styr,endyr,-0.2,0.2,est_q_devs,"q_dev");
  q_mat_dev.allocate(styr,endyr,-0.2,0.2,est_q_mat_devs,"q_mat_dev");
  log_recruits.allocate(styr,endyr,1,30,1,"log_recruits");
  sigma_m.allocate(1,2,0.01,4,est_sigma_m,"sigma_m");
  log_m_mu.allocate(1,2,-5,3,est_log_m_mu,"log_m_mu");
  prop_rec.allocate(styr,endyr,0.01,1,2,"prop_rec");
  sigma_q.allocate(1,2,0.01,4,est_sigma_q,"sigma_q");
  imm_n_size_pred.allocate(styr,endyr,1,size_n,"imm_n_size_pred");
  #ifndef NO_AD_INITIALIZE
    imm_n_size_pred.initialize();
  #endif
  mat_n_size_pred.allocate(styr,endyr,1,size_n,"mat_n_size_pred");
  #ifndef NO_AD_INITIALIZE
    mat_n_size_pred.initialize();
  #endif
  nat_m.allocate(styr,endyr,1,size_n,"nat_m");
  #ifndef NO_AD_INITIALIZE
    nat_m.initialize();
  #endif
  nat_m_mat.allocate(styr,endyr,1,size_n,"nat_m_mat");
  #ifndef NO_AD_INITIALIZE
    nat_m_mat.initialize();
  #endif
  selectivity.allocate(styr,endyr,1,size_n,"selectivity");
  #ifndef NO_AD_INITIALIZE
    selectivity.initialize();
  #endif
  selectivity_mat.allocate(styr,endyr,1,size_n,"selectivity_mat");
  #ifndef NO_AD_INITIALIZE
    selectivity_mat.initialize();
  #endif
  temp_imm.allocate(1,size_n,"temp_imm");
  #ifndef NO_AD_INITIALIZE
    temp_imm.initialize();
  #endif
  temp_mat.allocate(1,size_n,"temp_mat");
  #ifndef NO_AD_INITIALIZE
    temp_mat.initialize();
  #endif
  trans_imm.allocate(1,size_n,"trans_imm");
  #ifndef NO_AD_INITIALIZE
    trans_imm.initialize();
  #endif
  sum_imm_numbers_obs.allocate(styr,endyr,"sum_imm_numbers_obs");
  #ifndef NO_AD_INITIALIZE
    sum_imm_numbers_obs.initialize();
  #endif
  sum_mat_numbers_obs.allocate(styr,endyr,"sum_mat_numbers_obs");
  #ifndef NO_AD_INITIALIZE
    sum_mat_numbers_obs.initialize();
  #endif
  imm_numbers_pred.allocate(styr,endyr,"imm_numbers_pred");
  #ifndef NO_AD_INITIALIZE
    imm_numbers_pred.initialize();
  #endif
  mat_numbers_pred.allocate(styr,endyr,"mat_numbers_pred");
  #ifndef NO_AD_INITIALIZE
    mat_numbers_pred.initialize();
  #endif
  imm_num_like.allocate("imm_num_like");
  #ifndef NO_AD_INITIALIZE
  imm_num_like.initialize();
  #endif
  mat_num_like.allocate("mat_num_like");
  #ifndef NO_AD_INITIALIZE
  mat_num_like.initialize();
  #endif
  imm_like.allocate("imm_like");
  #ifndef NO_AD_INITIALIZE
  imm_like.initialize();
  #endif
  mat_like.allocate("mat_like");
  #ifndef NO_AD_INITIALIZE
  mat_like.initialize();
  #endif
  nat_m_like.allocate("nat_m_like");
  #ifndef NO_AD_INITIALIZE
  nat_m_like.initialize();
  #endif
  nat_m_mat_like.allocate("nat_m_mat_like");
  #ifndef NO_AD_INITIALIZE
  nat_m_mat_like.initialize();
  #endif
  nat_m_mu_like.allocate("nat_m_mu_like");
  #ifndef NO_AD_INITIALIZE
  nat_m_mu_like.initialize();
  #endif
  nat_m_mat_mu_like.allocate("nat_m_mat_mu_like");
  #ifndef NO_AD_INITIALIZE
  nat_m_mat_mu_like.initialize();
  #endif
  smooth_q_like.allocate("smooth_q_like");
  #ifndef NO_AD_INITIALIZE
  smooth_q_like.initialize();
  #endif
  smooth_m_like.allocate("smooth_m_like");
  #ifndef NO_AD_INITIALIZE
  smooth_m_like.initialize();
  #endif
  q_like.allocate("q_like");
  #ifndef NO_AD_INITIALIZE
  q_like.initialize();
  #endif
  q_mat_like.allocate("q_mat_like");
  #ifndef NO_AD_INITIALIZE
  q_mat_like.initialize();
  #endif
  f.allocate("f");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
}

void model_parameters::userfunction(void)
{
  f =0.0;
 for(int size=1;size<=size_n;size++)
  {
   imm_n_size_pred(styr,size) = exp(log_n_imm(size));
   mat_n_size_pred(styr,size) = exp(log_n_mat(size));
   }
 for(int year=styr;year<=endyr;year++)
  for(int size=1;size<=size_n;size++)
  {
  nat_m(year,size) = log_m_mu(1) + nat_m_dev(year);
  nat_m_mat(year,size) = log_m_mu(2) + nat_m_dev(year);
  selectivity(year,size) = survey_sel(size) + q_dev(year);
  selectivity_mat(year,size) = survey_sel(size) + q_dev(year);
  //==option for estimating separate devs by maturity state
  if(est_q_mat_devs>0)
   selectivity_mat(year,size) = survey_sel(size) + q_mat_dev(year);  
  if(est_m_mat_devs>0)
   nat_m_mat(year,size) = log_m_mu(2) + nat_m_mat_dev(year);
  }
 for(int year=styr;year<endyr;year++)
  {
  for (int size=1;size<=size_n;size++) 
      {
		temp_imm(size) = imm_n_size_pred(year,size) * exp(-1*0.5*exp( nat_m(year,size)));
		temp_mat(size) = mat_n_size_pred(year,size) * exp(-1*0.5*exp( nat_m_mat(year,size)));
	   }	
	  // growth
	   trans_imm = size_trans * temp_imm;
       trans_imm(1) += exp(log_recruits(year))*prop_rec(year);
       trans_imm(2) += exp(log_recruits(year))*(1-prop_rec(year));
	  // maturity
	   for (int size=1;size<=size_n;size++) 
	   {
	  	 imm_n_size_pred(year+1,size) = (trans_imm(size) * (1-prop_term_molt(year,size))) * exp(-1*0.5*exp( nat_m(year,size)));	
	  	 mat_n_size_pred(year+1,size) = (trans_imm(size) * prop_term_molt(year,size) + temp_mat(size)) *  exp(-1*0.5*exp(nat_m_mat(year,size)));
       }
    }
   evaluate_the_objective_function();
}

void model_parameters::evaluate_the_objective_function(void)
{
  // make total numbers by maturity state from obs and preds
  imm_numbers_pred.initialize();
  mat_numbers_pred.initialize();
  sum_imm_numbers_obs.initialize();
  sum_mat_numbers_obs.initialize(); 
  for (int year=styr;year<=endyr;year++)
   for (int size=1;size<=size_n;size++)
   {
    imm_numbers_pred(year) += selectivity(year,size)*imm_n_size_pred(year,size);
    mat_numbers_pred(year) += selectivity_mat(year,size)*mat_n_size_pred(year,size);
	sum_imm_numbers_obs(year) += imm_n_size_obs(year,size);
	sum_mat_numbers_obs(year) += mat_n_size_obs(year,size);
   }  
  // likelihoods
  imm_num_like = 0;
  for (int year=styr;year<=endyr;year++)
   if (year!=2020)
    imm_num_like += square( log(imm_numbers_pred(year)) - log(imm_n_obs(year))) / (2.0 * square(sigma_numbers_imm(year)));
  mat_num_like = 0;
  for (int year=styr;year<=endyr;year++)
   if (year!=2020)
    mat_num_like += square( log(mat_numbers_pred(year)) - log(mat_n_obs(year))) / (2.0 * square(sigma_numbers_mat(year)));
  // immature numbers at size data
  imm_like = 0;
  for (int year=styr;year<=endyr;year++)
   for (int size=1;size<=size_n;size++)
    if (imm_n_size_obs(year,size) >0)
     imm_like += imm_eff_samp*(imm_n_size_obs(year,size)/sum_imm_numbers_obs(year)) * log( (selectivity(year,size)*imm_n_size_pred(year,size)/imm_numbers_pred(year)) / (imm_n_size_obs(year,size)/sum_imm_numbers_obs(year)));
  imm_like = -1*imm_like;
  // mature numbers at size data
  mat_like = 0;
  for (int year=styr;year<=endyr;year++)
   for (int size=1;size<=size_n;size++)
    if (mat_n_size_obs(year,size) >0)
     mat_like += mat_eff_samp*(mat_n_size_obs(year,size)/sum_mat_numbers_obs(year)) * log( (selectivity_mat(year,size)*mat_n_size_pred(year,size)/mat_numbers_pred(year)) / (mat_n_size_obs(year,size)/sum_mat_numbers_obs(year)));
  mat_like = -1*mat_like;
 //penalties on m 
  nat_m_mu_like =0;
   nat_m_mu_like += pow((log_m_mu(1)-log_mu_m_prior(1))/ (sqrt(2)*sqrt(sigma_m_mu(1))),2.0);
  nat_m_mat_mu_like =0;
   nat_m_mat_mu_like += pow((log_m_mu(2)-log_mu_m_prior(2))/ (sqrt(2)*sqrt(sigma_m_mu(2))),2.0); 
  if(est_m_devs>0)
  {
	  nat_m_like =0;
  for (int year=styr;year<=endyr;year++)
   nat_m_like += pow((nat_m(year,1)-log_m_mu(1))/ (sqrt(2)*sqrt(sigma_m(1))),2.0);
  nat_m_mat_like =0;
  for (int year=styr;year<=endyr;year++)
   nat_m_mat_like += pow((nat_m_mat(year,1)-log_m_mu(2))/ (sqrt(2)*sqrt(sigma_m(2))),2.0);
  }
  if(est_q_devs>0)
  {
  q_like =0;
  for (int year=styr;year<=endyr;year++)
   q_like += pow((selectivity(year,4)-survey_sel(4))/ (sqrt(2)*sqrt(sigma_q(1))),2.0);
  q_mat_like =0;
  for (int year=styr;year<=endyr;year++)
   q_mat_like += pow((selectivity_mat(year,4)-survey_sel(4))/ (sqrt(2)*sqrt(sigma_q(2))),2.0);
  }
  smooth_q_like = 0;
  smooth_q_like = smooth_q_weight* (norm2(first_difference(first_difference(q_dev))) +norm2(first_difference(first_difference(q_mat_dev)))) ;
  smooth_m_like = 0;
  smooth_m_like = smooth_m_weight* (norm2(first_difference(first_difference(nat_m_dev))) +norm2(first_difference(first_difference(nat_m_mat_dev)))) ;
  f = imm_num_like + mat_num_like + imm_like + mat_like + nat_m_like + nat_m_mat_like + nat_m_mu_like + nat_m_mat_mu_like + smooth_q_like + smooth_m_like + q_like + q_mat_like;
}

void model_parameters::report(const dvector& gradients)
{
 adstring ad_tmp=initial_params::get_reportfile_name();
  ofstream report((char*)(adprogram_name + ad_tmp));
  if (!report)
  {
    cerr << "error trying to open report file"  << adprogram_name << ".rep";
    return;
  }
  report<<"$likelihoods"<<endl;
  report<<f<<" "<<imm_num_like<<" "<<mat_num_like<<" "<<imm_like<<" "<<mat_like<<" "<<nat_m_like<<" "<<nat_m_mat_like<<" "<<nat_m_mu_like<<" "<<nat_m_mat_mu_like<<" "<<smooth_q_like<<" "<<smooth_m_like<<" "<<q_like<<" "<<q_mat_like<<endl;
  report<<"$log recruits"<<endl;
  report<<log_recruits<<endl;
  report <<"$natural mortality" << endl;
  for(int i=styr; i<=endyr; i++)
  {
    report << mfexp(nat_m(i))<<endl;
  }
  report <<"$recruits" << endl;
  report << mfexp(log_recruits)<<endl;
  report <<"$immature numbers at size" << endl;
  for(int i=styr; i<=endyr; i++)
  {
    report << (elem_prod(selectivity(i),imm_n_size_pred(i)))/imm_numbers_pred(i)<<endl;
  }
  report <<"$mature numbers at size" << endl;
  for(int i=styr; i<=endyr; i++)
  {
    report << (elem_prod(selectivity_mat(i),mat_n_size_pred(i)))/mat_numbers_pred(i)<<endl;
  }
  report<<"$imm_numbers_pred"<<endl;
  report<<imm_numbers_pred<<endl;
  report<<"$mat_numbers_pred"<<endl;
  report<<mat_numbers_pred<<endl;
  report<<"$imm_n_obs"<<endl;
  report<<imm_n_obs<<endl;
  report<<"$mat_n_obs"<<endl;
  report<<mat_n_obs<<endl;
  report<<"$survey_sel_input"<<endl;
  report<<survey_sel<<endl;
    report <<"$mature natural mortality" << endl;
  for(int i=styr; i<=endyr; i++)
  {
    report << mfexp(nat_m_mat(i))<<endl;
  }
  report <<"$survey selectivity" << endl;
  for(int i=styr; i<=endyr; i++)
  {
    report << selectivity(i)<<endl;
  }
  report <<"$mature survey selectivity" << endl;
  for(int i=styr; i<=endyr; i++)
  {
    report << selectivity_mat(i)<<endl;
  }
  report <<"$pred_imm_pop_num" << endl;
  for(int i=styr; i<=endyr; i++)
  {
    report << (imm_n_size_pred(i))<<endl;
  }
  report <<"$pred_mat_pop_num" << endl;
  for(int i=styr; i<=endyr; i++)
  {
    report << (mat_n_size_pred(i))<<endl;
  }
}

void model_parameters::set_runtime(void)
{
  dvector temp1("{10000}");
  maximum_function_evaluations.allocate(temp1.indexmin(),temp1.indexmax());
  maximum_function_evaluations=temp1;
  dvector temp("{1e-3}");
  convergence_criteria.allocate(temp.indexmin(),temp.indexmax());
  convergence_criteria=temp;
}

void model_parameters::preliminary_calculations(void){
#if defined(USE_ADPVM)

  admaster_slave_variable_interface(*this);

#endif
}

model_data::~model_data()
{}

model_parameters::~model_parameters()
{}

void model_parameters::final_calcs(void){}

#ifdef _BORLANDC_
  extern unsigned _stklen=10000U;
#endif


#ifdef __ZTC__
  extern unsigned int _stack=10000U;
#endif

  long int arrmblsize=0;

int main(int argc,char * argv[])
{
    ad_set_new_handler();
  ad_exit=&ad_boundf;
    gradient_structure::set_NO_DERIVATIVES();
#ifdef DEBUG
  #ifndef __SUNPRO_C
std::feclearexcept(FE_ALL_EXCEPT);
  #endif
  auto start = std::chrono::high_resolution_clock::now();
#endif
    gradient_structure::set_YES_SAVE_VARIABLES_VALUES();
    if (!arrmblsize) arrmblsize=15000000;
    model_parameters mp(arrmblsize,argc,argv);
    mp.iprint=10;
    mp.preliminary_calculations();
    mp.computations(argc,argv);
#ifdef DEBUG
  std::cout << endl << argv[0] << " elapsed time is " << std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start).count() << " microseconds." << endl;
  #ifndef __SUNPRO_C
bool failedtest = false;
if (std::fetestexcept(FE_DIVBYZERO))
  { cerr << "Error: Detected division by zero." << endl; failedtest = true; }
if (std::fetestexcept(FE_INVALID))
  { cerr << "Error: Detected invalid argument." << endl; failedtest = true; }
if (std::fetestexcept(FE_OVERFLOW))
  { cerr << "Error: Detected overflow." << endl; failedtest = true; }
if (std::fetestexcept(FE_UNDERFLOW))
  { cerr << "Error: Detected underflow." << endl; }
if (failedtest) { std::abort(); } 
  #endif
#endif
    return 0;
}

extern "C"  {
  void ad_boundf(int i)
  {
    /* so we can stop here */
    exit(i);
  }
}
