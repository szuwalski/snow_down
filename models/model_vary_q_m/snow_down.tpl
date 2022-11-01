DATA_SECTION
 //Bering Sea snow crab model
  
  init_int styr  																							
  init_int endyr 
  init_int dat_yr
  init_ivector years(1,dat_yr)
  init_int size_n
  init_vector sizes(1,size_n)
  init_vector imm_n_obs(styr,endyr)
  init_vector mat_n_obs(styr,endyr)
  init_matrix imm_n_size_obs(styr,endyr,1,size_n)
  init_matrix mat_n_size_obs(styr,endyr,1,size_n)
  init_matrix prop_term_molt(styr,endyr,1,size_n)
  init_matrix size_trans(1,size_n,1,size_n) 
  init_vector sigma_numbers_imm(styr,endyr)
  init_vector sigma_numbers_mat(styr,endyr)
  init_number mat_eff_samp
  init_number imm_eff_samp
  init_vector survey_sel(1,size_n)
  init_vector log_mu_m_prior(1,2)
  init_number est_m_devs
  init_number est_q_devs
  init_number est_m_mat_devs
  init_number est_q_mat_devs
  init_number est_sigma_m
  init_vector sigma_m_mu(1,2)
  init_number smooth_q_weight
  init_number smooth_m_weight
  init_number est_log_m_mu
  init_number est_sigma_q
  !!cout<<survey_sel<<endl;
  !!cout<<size_trans<<endl;
  !!cout<<sigma_m_mu<<endl;

PARAMETER_SECTION
  init_bounded_vector log_n_imm(1,size_n,1,30,1)
  init_bounded_vector log_n_mat(1,size_n,1,30,1)
  init_bounded_dev_vector nat_m_dev(styr,endyr,-4,4,est_m_devs)
  init_bounded_dev_vector nat_m_mat_dev(styr,endyr,-4,4,est_m_mat_devs)
  init_bounded_vector q_dev(styr,endyr,-0.2,0.2,est_q_devs)
  init_bounded_vector q_mat_dev(styr,endyr,-0.2,0.2,est_q_mat_devs)
  init_bounded_vector log_recruits(styr,endyr,1,30,1)
  init_bounded_vector sigma_m(1,2,0.01,4,est_sigma_m)
  init_bounded_vector log_m_mu(1,2,-5,3,est_log_m_mu)
  init_bounded_vector prop_rec(styr,endyr,0.01,1,2)
  init_bounded_vector sigma_q(1,2,0.01,4,est_sigma_q)
	
  matrix imm_n_size_pred(styr,endyr,1,size_n)
  matrix mat_n_size_pred(styr,endyr,1,size_n)
  matrix nat_m(styr,endyr,1,size_n)
  matrix nat_m_mat(styr,endyr,1,size_n)
  matrix selectivity(styr,endyr,1,size_n)
  matrix selectivity_mat(styr,endyr,1,size_n)
  
  vector temp_imm(1,size_n)
  vector temp_mat(1,size_n)
  vector trans_imm(1,size_n)

  vector sum_imm_numbers_obs(styr,endyr)
  vector sum_mat_numbers_obs(styr,endyr)
  vector imm_numbers_pred(styr,endyr)
  vector mat_numbers_pred(styr,endyr)
   
  number imm_num_like
  number mat_num_like
  number imm_like
  number mat_like
  number nat_m_like
  number nat_m_mat_like
  number nat_m_mu_like
  number nat_m_mat_mu_like
  number smooth_q_like
  number smooth_m_like
  number q_like
  number q_mat_like
    
  objective_function_value f
 
//==============================================================================
PROCEDURE_SECTION

// initial year numbers at size and selectivity
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
   
//==============================================================================
FUNCTION evaluate_the_objective_function
  
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
  
// ========================y==================================================   
REPORT_SECTION
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
  
RUNTIME_SECTION
//one number for each phase, if more phases then uses the last number
  maximum_function_evaluations 10000
  convergence_criteria 1e-3

