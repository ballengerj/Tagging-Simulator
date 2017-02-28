
GLOBALS_SECTION

  #include "admodel.h"  // make sure this file is in folder, used to create r compatible output
  #define EOUT(var) cout <<#var<<" "<<var<<endl

TOP_OF_MAIN_SECTION

  arrmblsize=500000000;  // memory allocation variables, should not need to be changed
  gradient_structure::set_MAX_NVAR_OFFSET(5000000);
  gradient_structure::set_NUM_DEPENDENT_VARIABLES(5000000);

DATA_SECTION

  ////////////////////////////////////////////////////////////////////
  // inputs for the model
  //////////////////////////////////////////////////////////////\
  init_int ncohorts
  init_int nages
  init_int max_age
  init_int nyrs
  init_int nstocks
  init_vector yrs(1,nyrs)
  init_vector ages(1,nages)
  init_vector tagging_month(1,nstocks) // ==1 tag jan 1st; used to define fractional  mortality in first year of release
  init_vector fish_start_month(1,nstocks) // start date (begin of month) used to define the fishing season length in comparison to the tag season (i.e., fraction of year that tags encounter fishing) 
  init_vector fish_end_month(1,nstocks) // end date (end of month) used to define the fishing season length in comparison to the tag season (i.e., fraction of year that tags encounter fishing) 
  init_ivector max_age_recap(1,nyrs)
  init_number phase_T_est_CNST
  init_number phase_T_est
  init_number phase_T_est_AGE
  init_number phase_T_est_AGE_no_age1
  init_number phase_T_res
  init_number phase_T_res_age
  init_number phase_T_res_no_age1
  init_number phase_report_rate
  init_number phase_F
  init_number phase_M_region
  init_number phase_M_age
  init_number phase_M_CNST
  init_number phase_M_year
  init_number report_rate_switch
  init_number initial_report_rate
  init_vector report_rate_fixed(1,nstocks)
  init_number report_rate_pen_switch  //==0 no prior, ==1 normal prior (use Report_rate_sigma, report_rate_ave), ==2 penalty (use report_rate_pen_low, report_rate_pen_hi)
  init_number report_rate_sigma
  init_number report_rate_ave
  init_number report_rate_pen_hi
  init_number report_rate_pen_low
  init_number report_rate_pen_mult
  init_number F_switch  //==1 F est for all ages/years/regions with a plus group at max_age, ==2 no plus group (estimate F for each year and resulting age of model)
  init_number initial_F
  init_3darray Fixed_F(1,nstocks,1,nyrs,1,nages)
  init_number F_pen_switch  //==0 no prior, ==1 normal prior (use F_sigma, F_ave), ==2 penalty (F_pen_hi)
  init_number F_sigma
  init_number F_ave
  init_number F_pen_hi
  init_number F_pen_low
  init_number F_pen_mult
  init_number M_switch  //==0 est M all ages all regions, ==1 est  for only 1 region
  init_number initial_M
  init_vector Fixed_M(1,nages)
  init_3darray Fixed_M_region(1,nstocks,1,nyrs,1,nages)
  init_number M_pen_switch  //==0 no prior, ==1 normal prior (use M_sigma, M_ave), ==2 penalty (use M_pen_low, M_pen_hi)
  init_number M_sigma
  init_number M_ave
  init_number M_pen_low     // low M values to avoid
  init_number M_pen_high   //high M values to avoid
  init_number M_pen_mult
  init_number movement_switch // ==-1, fix at SIM_T,==0 fix at 100% residency, ==1 est 1 Constant T, ==2 region T, ==3 est age and region, ==4 estimate residency, ==5 estimate residency by age
  init_number natal_homing //==0 markovian movement (random, no memory), ==1 natal homing (follows movement pattern of home stock at all times)
  init_number no_move_first_yr //==0 allow movement in year of release, ==1 do not allow movement in year of release
  init_matrix Move_Fract(1,nstocks,1,nstocks)  //matrix defining partitioning of movement for natal homing and symmetric movement scenarios, diagnols==1 and offdiagnols must sum to 1!
  init_number initial_T
  init_5darray Fixed_T(1,nstocks,1,ncohorts,1,nages,1,nyrs,1,nstocks)
  init_number move_pen_switch // ==0 no prior, ==1 normal prior (use T_ave and T_sigma), ==2 penalty (use T_low, T_hi)
  init_number T_sigma
  init_number T_ave
  init_number T_low            // value for T below which movement penalty becomes active (default=-10....effectively equivalent to T=0)
  init_number T_hi
  init_number Tpen_mult
  init_number tag_neff_switch // ==1 set tag_neff=ntags, else tag_neff==neff
  init_3darray neff(1,nstocks,1,ncohorts,1,nages)
  init_number TagCst              // constant for residuals to avoid log(0)
  init_number tag_loss_switch
  init_number handling_mortality_switch
  /////////////////////////////////////
  // Observed Data
  //////////////////////////////
  init_vector OBS_Tag_Loss(1,nages)
  init_vector OBS_Hand_Mort(1,nages)
  init_3darray ntags(1,nstocks,1,ncohorts,1,nages)  // tags at age by region
  init_5darray OBS_tag_rec(1,nstocks,1,ncohorts,1,nages,1,nyrs,1,nstocks)
  init_3darray OBS_total_tag_rec(1,nstocks,1,ncohorts,1,nages)
  init_3darray OBS_tag_not_rec(1,nstocks,1,ncohorts,1,nages)
  //////////////////////////////////////////////////////////
  // data manipulation vectors (initial calcs)
  //////////////////////////////////////////////////
  5darray OBS_tag_prop(1,nstocks,1,ncohorts,1,nages,1,nyrs,1,nstocks)
  3darray OBS_tag_prop_not_rec(1,nstocks,1,ncohorts,1,nages)
  3darray OBS_total_rec(1,nstocks,1,ncohorts,1,nages)
  matrix OBS_tag_rec_temp(1,nyrs,1,nstocks)
  //////////////////////////////////////////////////
  // Simulation Values for estimated parameters
  ////////////////////////////////////////////////////
  init_number SIM_mixing_switch
  init_5darray SIM_T(1,nstocks,1,ncohorts,1,nages,1,nyrs,1,nstocks)
  3darray SIM_T_region(1,nstocks,1,nages,1,nstocks)
  init_3darray SIM_M(1,nstocks,1,nyrs,1,nages)
  init_3darray SIM_F(1,nstocks,1,nyrs,1,nages)
  init_vector SIM_report_rate(1,nstocks)
  ////////////////////////////////////////////////////
  // integers for loops
  ////////////////////////////////////////
  init_number dummy
  init_number phase_dummy1 //(should only be positive when all other parameter phases are negative)
  int a
  int y
  int z
  int k
  int j
  int i
  int w
  int s
  ///////////////////////////////////////////////////////
  // check to make sure data input file read in correctly
  //////////////////////////////////////////////////////////
 !! cout << "input read" << endl;

 LOCAL_CALCS
  EOUT(dummy);
  if (dummy!=-1345.)
  {
    cout<<"problem .dat file"<<endl;
    exit(1);
  }
 END_CALCS
PARAMETER_SECTION
  init_number dummy1(phase_dummy1)  //for use when only want to simulate dynamics based on input parameters (i.e., no estimation)
  init_3darray ln_T_est_AGE(1,nstocks,1,nages,1,nstocks-1,phase_T_est_AGE)
  init_3darray ln_T_est_AGE_no_age1(1,nstocks,1,nages-1,1,nstocks-1,phase_T_est_AGE_no_age1)
  init_bounded_matrix ln_T_est(1,nstocks,1,nstocks-1,-120,10.0,phase_T_est)  //movement parameters
  init_bounded_number ln_T_est_CNST(-120,10.0,phase_T_est_CNST)  //for symmetric movement only 1 parameter to estimate
  init_bounded_vector ln_T_res(1,nstocks,-120,10.0,phase_T_res)  //for symmetric movement only 1 parameter to estimate
  init_bounded_matrix ln_T_res_age(1,nstocks,1,nages,-120,10.0,phase_T_res_age)  //for symmetric movement only 1 parameter to estimate
  init_bounded_matrix ln_T_res_no_age1(1,nstocks,1,nages-1,-120,10.0,phase_T_res_no_age1)  //for symmetric movement only 1 parameter to estimate

  matrix G(1,nstocks,1,nstocks)
  vector G_temp(1,nstocks)
  3darray G_age(1,nstocks,1,nages,1,nstocks)
  matrix G_temp_age(1,nages,1,nstocks)
  vector G_temp_sum_age(1,nages)
  matrix G_sum(1,nstocks,1,nages)
  5darray T(1,nstocks,1,ncohorts,1,nages,1,nyrs,1,nstocks)
  3darray T_region(1,nstocks,1,nages,1,nstocks)

  init_bounded_vector ln_report_rate(1,nstocks,-30,20,phase_report_rate)
  !! int ny=nyrs;
  !! int ns=nstocks;
  !! ivector age_recap=max_age_recap;
  init_3darray ln_F(1,ny,1,ns,1,age_recap,phase_F)
  3darray F(1,nstocks,1,nyrs,1,nages) 

  init_bounded_matrix ln_M_region(1,nstocks,1,nages,-10,2,phase_M_region)
  init_3darray ln_M_year(1,ny,1,ns,1,age_recap,phase_M_year) 
  init_bounded_vector ln_M_age(1,nages,-10,2,phase_M_age)  
  init_bounded_number ln_M_est(-10,2,phase_M_CNST)
  3darray M(1,nstocks,1,nyrs,1,nages) 
  3darray Total_Survival(1,nstocks,1,nyrs,1,nages)
  3darray Mort_Fish(1,nstocks,1,nyrs,1,nages)
  
  vector Hand_Mort(1,nages)
  vector Tag_Loss(1,nages)
  matrix report_rate(1,nyrs,1,nstocks)
  vector report_rate_CNST(1,nstocks)

  5darray tags_avail(1,nstocks,1,ncohorts,1,nages,1,nyrs,1,nstocks)
  5darray pred_rec(1,nstocks,1,ncohorts,1,nages,1,nyrs,1,nstocks)
  matrix total_recap_temp_stock(1,nyrs,1,nstocks)
  vector tags_avail_temp(1,nstocks)
  3darray total_rec(1,nstocks,1,ncohorts,1,nages)
  3darray not_rec(1,nstocks,1,ncohorts,1,nages)
  
  5darray pred_tag_prop(1,nstocks,1,ncohorts,1,nages,1,nyrs,1,nstocks)
  3darray pred_tag_prop_not_rec(1,nstocks,1,ncohorts,1,nages)

  3darray M_pen(1,nstocks,1,nyrs,1,nages)
  5darray T_pen(1,nstocks,1,ncohorts,1,nages,1,nyrs,1,nstocks)
  vector report_rate_pen(1,nstocks)
  3darray F_pen(1,nstocks,1,nyrs,1,nages)
  3darray tag_neff(1,nstocks,1,ncohorts,1,nages)

  5darray res_tags(1,nstocks,1,ncohorts,1,nages,1,nyrs,1,nstocks)
  3darray res_tags_not_rec(1,nstocks,1,ncohorts,1,nages)
  5darray rss_tags_temp(1,nstocks,1,ncohorts,1,nages,1,nyrs,1,nstocks)
  3darray rss_tags_temp_not_rec(1,nstocks,1,ncohorts,1,nages)
  matrix rss_tags_temp_age(1,nyrs,1,nstocks)
  3darray rss_tags_temp_stock(1,nstocks,1,ncohorts,1,nages)
  number rss_tags

 ///////////////////////////////////////////////////////////////////////////////
  // REMOVE WHEN NEXT ADMB VERSION RELEASED
  4darray T_reg1(1,ncohorts,1,nages,1,nyrs,1,nstocks)
  4darray T_reg2(1,ncohorts,1,nages,1,nyrs,1,nstocks)
  4darray T_reg3(1,ncohorts,1,nages,1,nyrs,1,nstocks)
  4darray T_reg4(1,ncohorts,1,nages,1,nyrs,1,nstocks)
  4darray tags_avail1(1,ncohorts,1,nages,1,nyrs,1,nstocks)
  4darray tags_avail2(1,ncohorts,1,nages,1,nyrs,1,nstocks)
  4darray tags_avail3(1,ncohorts,1,nages,1,nyrs,1,nstocks)
  4darray tags_avail4(1,ncohorts,1,nages,1,nyrs,1,nstocks)
  4darray pred_rec1(1,ncohorts,1,nages,1,nyrs,1,nstocks)
  4darray pred_rec2(1,ncohorts,1,nages,1,nyrs,1,nstocks)
  4darray pred_rec3(1,ncohorts,1,nages,1,nyrs,1,nstocks)
  4darray pred_rec4(1,ncohorts,1,nages,1,nyrs,1,nstocks)
  4darray pred_tag_prop1(1,ncohorts,1,nages,1,nyrs,1,nstocks)
  4darray pred_tag_prop2(1,ncohorts,1,nages,1,nyrs,1,nstocks)
  4darray pred_tag_prop3(1,ncohorts,1,nages,1,nyrs,1,nstocks)
  4darray pred_tag_prop4(1,ncohorts,1,nages,1,nyrs,1,nstocks)
  4darray res_tags1(1,ncohorts,1,nages,1,nyrs,1,nstocks)
  4darray res_tags2(1,ncohorts,1,nages,1,nyrs,1,nstocks)
  4darray res_tags3(1,ncohorts,1,nages,1,nyrs,1,nstocks)
  4darray res_tags4(1,ncohorts,1,nages,1,nyrs,1,nstocks)
  /////////////////////////////////////////////////////////////////////////////


  objective_function_value f

INITIALIZATION_SECTION  //set initial values
     
  ln_T_est initial_T;
  ln_T_est_CNST initial_T;
  ln_T_est_AGE initial_T;
  ln_T_est_AGE_no_age1 initial_T;
  ln_T_res initial_T;
  ln_T_res_age initial_T;
  ln_T_res_no_age1 initial_T;
  ln_report_rate initial_report_rate;
  ln_F initial_F;
  ln_M_region initial_M;
  ln_M_age initial_M;
  ln_M_est initial_M;
  ln_M_year initial_M;

PRELIMINARY_CALCS_SECTION

  for(int i=1;i<=nstocks;i++) //release stock
    {
    for(int x=1;x<=ncohorts;x++)  //release year
     {
      for (int a=1;a<=nages;a++) //release age  // setup 3d array to hold age (first index..within each age is a matrix of year x stock where now use nstock x nstock to hold all recap data)
        {
         OBS_tag_rec_temp=0;
         for(int y=1;y<=nyrs;y++)  //recap year
          {         
           for(int j=1;j<=nstocks;j++) //recap stock
            {
             OBS_tag_rec_temp(y,j)=OBS_tag_rec(i,x,a,y,j);
            if(ntags(i,x,a)>0)
              {
               OBS_tag_prop(i,x,a,y,j)=OBS_tag_rec(i,x,a,y,j)/ntags(i,x,a);
               OBS_tag_prop_not_rec(i,x,a)=OBS_tag_not_rec(i,x,a)/ntags(i,x,a);
              }
            if(ntags(i,x,a)==0)
              {                   
               OBS_tag_prop(i,x,a,y,j)=0;
               OBS_tag_prop_not_rec(i,x,a)=0;
              }
             } 
            }
           OBS_total_rec(i,x,a)=sum(OBS_tag_rec_temp);
           }
          }
         }  
PROCEDURE_SECTION

  get_movement();

  get_mortality();

  get_tag_loss();

  get_tag_recaptures();
 
  evaluate_the_objective_function();

FUNCTION get_movement
     for (int i=1;i<=nstocks;i++)
       {
       for(int x=1;x<=ncohorts;x++)  //release year
        {
         for(int a=1;a<=nages;a++)
          {
           for(int y=1;y<=nyrs;y++)  //recap year
            {  
             for (int j=1;j<=nstocks;j++)
             {
              T(i,x,a,y,j)=0;
             }
            }
           }
         }
        }
 //////////// Fix at Assumed Value ///////////////////////////////////////////////////////////////////
    if(movement_switch==-2)  //Fix at Fixed_T
     {
      for (int i=1;i<=nstocks;i++)
       {
       for(int x=1;x<=ncohorts;x++)  //release year
        {
         for(int a=1;a<=nages;a++)
          {
           for(int y=1;y<=nyrs;y++)  //recap year
            {  
             for (int j=1;j<=nstocks;j++)
             {
             if(no_move_first_yr==1 && y==x  && i==j)
             {
             T(i,x,a,y,j)=1;
             }
             if(no_move_first_yr==1 && y==x  && i!=j)
             {
             T(i,x,a,y,j)=0;
             }             
             if(no_move_first_yr==1 && y!=x)
             {
             T(i,x,a,y,j)=Fixed_T(i,x,a,y,j);
             }
             if(no_move_first_yr==0)
             {
             T(i,x,a,y,j)=Fixed_T(i,x,a,y,j);
             }
             }
             }
            }
           }
         }
        }

 //////////// Fix at True Value ///////////////////////////////////////////////////////////////////
    if(movement_switch==-1)  //Fix at SIM_T
     {
      for (int i=1;i<=nstocks;i++)
       {
       for(int x=1;x<=ncohorts;x++)  //release year
        {
         for(int a=1;a<=nages;a++)
          {
           for(int y=1;y<=nyrs;y++)  //recap year
            {  
             for (int j=1;j<=nstocks;j++)
             {
             if(no_move_first_yr==1 && y==x  && i==j)
             {
             T(i,x,a,y,j)=1;
             }
             if(no_move_first_yr==1 && y==x  && i!=j)
             {
             T(i,x,a,y,j)=0;
             }             
             if(no_move_first_yr==1 && y>x)
             {
             T(i,x,a,y,j)=SIM_T(i,x,a,y,j);
             }
             if(no_move_first_yr==0)
             {
             T(i,x,a,y,j)=SIM_T(i,x,a,y,j);
             }
             }
             }
            }
           }
         }
        }
 ////////////////////////////////////////////////////////////////////////////////////////////////////

 ///////////// Fix at 100% Residency ////////////////////////////////////////////////////////////////
    if(movement_switch==0)  //Fix at 100% residency
     {
      for (int i=1;i<=nstocks;i++)
       {
        for(int x=1;x<=ncohorts;x++)  //release year
         {       
         for(int a=1;a<=nages;a++)
          {
           for(int y=1;y<=nyrs;y++)  //recap year
            {  
             for (int j=1;j<=nstocks;j++)
             {
             if(i==j)
             {
             T(i,x,a,y,j)=1;
             }
             if(i!=j)
             {
             T(i,x,a,y,j)=0;
             }
            }
           }
          }
         }
        }
       }
 ///////////////////////////////////////////////////////////////////////////////////////////////////////

 ////////// Constant Symmetric Movement ///////////////////////////////////////////////////
 if(movement_switch==1)  
  {     
  for (int i=1;i<=nstocks;i++)
   {
   for(int x=1;x<=ncohorts;x++)  //release year
    {
    for(int a=1;a<=nages;a++)
    {
     for(int y=1;y<=nyrs;y++)  //recap year
      {  
       for (int j=1;j<=nstocks;j++)
        {
             if(no_move_first_yr==1 && y==x  && i==j)
             {
             T(i,x,a,y,j)=1;
             }
             if(no_move_first_yr==1 && y==x  && i!=j)
             {
             T(i,x,a,y,j)=0;
             }  
        if(no_move_first_yr==1 && y>x && i==j)
        {
        T(i,x,a,y,j)=(mfexp(ln_T_est_CNST)/(mfexp(ln_T_est_CNST)+1));
        }
        if(no_move_first_yr==1 && y>x && i!=j)
        {
        T(i,x,a,y,j)=(1-(mfexp(ln_T_est_CNST)/(mfexp(ln_T_est_CNST)+1)))*(1/(nstocks-1));
        }
        if(no_move_first_yr==0 && y>=x && i==j)
        {
        T(i,x,a,y,j)=(mfexp(ln_T_est_CNST)/(mfexp(ln_T_est_CNST)+1));
        }
        if(no_move_first_yr==0 && y>=x && i!=j)
        {
        T(i,x,a,y,j)=(1-(mfexp(ln_T_est_CNST)/(mfexp(ln_T_est_CNST)+1)))*(1/(nstocks-1));
        }
        }
      }
     }
    }
   }
  }
 //////////////////////////////////////////////////////////////////////////////////////////////////////////

 ////////// Estimate T for Every Region stocks*(stocks-1) T values Estimated ////////////////////////////
    if(movement_switch==2)  //estimate movement for each region
     {
     for (int j=1;j<=nstocks;j++)
      {
      for (int i=1;i<=nstocks;i++) 
       {
            if(j==i)
            {
            G(j,i)=1;
            }
            if(i>j)
            {
            G(j,i)=mfexp(ln_T_est(j,i-1));
            }
            if(j!=i && i<j)
            {
            G(j,i)=mfexp(ln_T_est(j,i));
            }
        }
       }    
        G_temp=rowsum(G);     
  for (int i=1;i<=nstocks;i++)
   {
   for(int x=1;x<=ncohorts;x++)  //release year
    {
    for(int a=1;a<=nages;a++)
    {
     for(int y=1;y<=nyrs;y++)  //recap year
      {  
       for (int j=1;j<=nstocks;j++)
        {
             if(no_move_first_yr==1 && y==x  && i==j)
             {
             T(i,x,a,y,j)=1;
             }
             if(no_move_first_yr==1 && y==x  && i!=j)
             {
             T(i,x,a,y,j)=0;
             }  
        if(no_move_first_yr==1 && y>x)
        {
        T(i,x,a,y,j)=G(i,j)/G_temp(i);
        }
        if(no_move_first_yr==0 && y>=x)
        {
        T(i,x,a,y,j)=G(i,j)/G_temp(i);
        }
        }
       }
      }
     }
    }
   }
 /////////////////////////////////////////////////////////////////////////////////////////////////////////////

 ////////////////// Estimate by Region and Age Age*stocks*(stocks-1) T Values Estimated //////////////////////
    if(movement_switch==3)  //estimate movement for each region and age
     {
    if(no_move_first_yr==0)
     {
     for (int j=1;j<=nstocks;j++)
      {
       G_temp_age=0;
       G_temp_sum_age=0;
       for(int a=1;a<=nages;a++)
       {
        for (int i=1;i<=nstocks;i++) 
         {
            if(j==i)
            {
            G_age(j,a,i)=1;
            }
            if(i>j)
            {
            G_age(j,a,i)=mfexp(ln_T_est_AGE(j,a,i-1)); 
            }
            if(j!=i && i<j)
            {
            G_age(j,a,i)=mfexp(ln_T_est_AGE(j,a,i));
            }
          G_temp_age(a,i)=G_age(j,a,i);           
        }
        G_temp_sum_age=rowsum(G_temp_age);
        G_sum(j,a)=G_temp_sum_age(a);
       }
      }
  for (int i=1;i<=nstocks;i++)
   {
   for(int x=1;x<=ncohorts;x++)  //release year
    {
    for(int a=1;a<=nages;a++)
    {
     for(int y=1;y<=nyrs;y++)  //recap year
      {  
       for (int j=1;j<=nstocks;j++)
        {
        if(y>=x)
        {
        T(i,x,a,y,j)=G_age(i,min((a+y-x),nages),j)/G_sum(i,min((a+y-x),nages));
        }
        }
       }
      }
     }
    }
   }
  }
 /////////////////////////////////////////////////////////////////////////////////////////////////////////////

 ////////////////// Estimate by Region and Age but no movement in first year of release (Age-1)*stocks*(stocks-1) T Values Estimated //////////////////////

    if(movement_switch==4)  //estimate movement for each region and age but no movement in first year of release
     {
      if(no_move_first_yr==1) //if don't move in the first year then there is no estimate of movement for the first age of release/recature because fish never move at that age (first movement is at youngest age of release+1)
       {
       for (int j=1;j<=nstocks;j++)
        {
         G_temp_age=0;
         G_temp_sum_age=0;
         for(int a=2;a<=nages;a++)
         {
          for (int i=1;i<=nstocks;i++) 
           {
            if(j==i)
            {
            G_age(j,a,i)=1;
            }
            if(i>j)
            {
            G_age(j,a,i)=mfexp(ln_T_est_AGE_no_age1(j,a-1,i-1)); 
            }
            if(j!=i && i<j)
            {
            G_age(j,a,i)=mfexp(ln_T_est_AGE_no_age1(j,a-1,i));
            }
          G_temp_age(a,i)=G_age(j,a,i);           
        }
        G_temp_sum_age=rowsum(G_temp_age);
        G_sum(j,a)=G_temp_sum_age(a);
       }
      }
  for (int i=1;i<=nstocks;i++)
   {
   for(int x=1;x<=ncohorts;x++)  //release year
    {
    for(int a=1;a<=nages;a++)
    {
     for(int y=1;y<=nyrs;y++)  //recap year
      {  
       for (int j=1;j<=nstocks;j++)
        {
           if(y==x  && i==j)
           {
           T(i,x,a,y,j)=1;
           }
           if(y==x  && i!=j)
           {
           T(i,x,a,y,j)=0;
           }
        if(y>x)
        {
        T(i,x,a,y,j)=G_age(i,min((a+y-x),nages),j)/G_sum(i,min((a+y-x),nages));  //need to adjust release age (a) for actual age (release age+year-year of release; a+y-x)
        }
       }
      }
     }
    }
   }
  }
 }
 /////////////////////////////////////////////////////////////////////////////////////////////////////////
 
 ////////// Age-Independent Residency///////////////////////////////////////////////////
  if(movement_switch==5)  //estimate age constant residency
  {       
  for (int i=1;i<=nstocks;i++)
   {
   for(int x=1;x<=ncohorts;x++)  //release year
    {
    for(int a=1;a<=nages;a++)
    {
     for(int y=1;y<=nyrs;y++)  //recap year
      {  
       for (int j=1;j<=nstocks;j++)
       {
             if(no_move_first_yr==1 && y==x  && i==j)
             {
             T(i,x,a,y,j)=1;
             }
             if(no_move_first_yr==1 && y==x  && i!=j)
             {
             T(i,x,a,y,j)=0;
             }  
        if(no_move_first_yr==1 && y>x && i==j)
        {
        T(i,x,a,y,j)=mfexp(ln_T_res(i))/(mfexp(ln_T_res(i)+1));
        }
        if(no_move_first_yr==1 && y>x && i!=j)
        {
        T(i,x,a,y,j)=(1-(mfexp(ln_T_res(i))/(mfexp(ln_T_res(i)+1))))*Move_Fract(i,j);
        }
        if(no_move_first_yr==0 && y>=x && i==j)
        {
        T(i,x,a,y,j)=mfexp(ln_T_res(i))/(mfexp(ln_T_res(i)+1));
        }
        if(no_move_first_yr==0 && y>=x && i!=j)
        {
        T(i,x,a,y,j)=(1-(mfexp(ln_T_res(i))/(mfexp(ln_T_res(i)+1))))*Move_Fract(i,j);
        }
        }
      }
     }
    }
   }
  }

 /////////////////////////////////////////////////////////////////////////////////////////////////////////
 
 ////////// Age-Dependent Residency ///////////////////////////////////////////////////
  if(movement_switch==6)  //estimate age dependent residency (movement dependent on age, release stock, and destination stock)
  {
  if(no_move_first_yr==0)
  {
  for (int i=1;i<=nstocks;i++)
   {
   for(int x=1;x<=ncohorts;x++)  //release year
    {
    for(int a=1;a<=nages;a++)
    {
     for(int y=1;y<=nyrs;y++)  //recap year
      {  
       for (int j=1;j<=nstocks;j++)
       {
        if(y>=x && i==j)
        {
        T(i,x,a,y,j)=mfexp(ln_T_res_age(i,min((a+y-x),nages)))/(mfexp(ln_T_res_age(i,min((a+y-x),nages))+1));
        }
        if(y>=x && i!=j)
        {
        T(i,x,a,y,j)=(1-(mfexp(ln_T_res_age(i,min((a+y-x),nages)))/(mfexp(ln_T_res_age(i,min((a+y-x),nages))+1))))*Move_Fract(i,j);
        }
        }
      }
     }
    }
   }
  }
 }
 /////////////////////////////////////////////////////////////////////////////////////////////////////////
 
 ////////// Age-Dependent Residency ///////////////////////////////////////////////////
  if(movement_switch==7) //estimate age dependent residency but no movement in first year (movement dependent on age, release stock, and destination stock)
  {
  if(no_move_first_yr==1)
  {
  for (int i=1;i<=nstocks;i++)
   {
   for(int x=1;x<=ncohorts;x++)  //release year
    {
    for(int a=1;a<=nages;a++)
    {
     for(int y=1;y<=nyrs;y++)  //recap year
      {  
       for (int j=1;j<=nstocks;j++)
       {
             if(y==x  && i==j)
             {
             T(i,x,a,y,j)=1;
             }
             if(y==x  && i!=j)
             {
             T(i,x,a,y,j)=0;
             }  
        if(y>x && i==j)
        {
        T(i,x,a,y,j)=mfexp(ln_T_res_no_age1(i,min((a+y-x),nages)-1))/(mfexp(ln_T_res_no_age1(i,min((a+y-x),nages)-1)+1));
        }
        if(y>x && i!=j)
        {
        T(i,x,a,y,j)=(1-(mfexp(ln_T_res_no_age1(i,min((a+y-x),nages)-1))/(mfexp(ln_T_res_no_age1(i,min((a+y-x),nages)-1)+1))))*Move_Fract(i,j);
        }
        }
      }
     }
    }
   }
  }
 }
FUNCTION get_mortality
  for(int i=1;i<=nstocks;i++) 
    {
     for(int y=1;y<=nyrs;y++) 
      {
         for (int a=1;a<=nages;a++)
          {
 /////////////////////// reporting rate //////////////////////////////////
          if(report_rate_switch==-2)
           {
            report_rate(y,i)=report_rate_fixed(i);
           }
          if(report_rate_switch==-1)
           {
            report_rate(y,i)=SIM_report_rate(i);
           }
          if(report_rate_switch==1)
           {
            report_rate(y,i)=(mfexp(ln_report_rate(i)))/(mfexp(ln_report_rate(i))+1);
           }
           report_rate_CNST(i)=report_rate(1,i);
 /////////////////////////////////////////////////////////////////////////////////////

 /////////////// Natural Mortality /////////////////////////////////////////////////////
          if(M_switch==-3)  //fixed assumed M by region, year, age
           {
            M(i,y,a)=Fixed_M_region(i,y,a);
           }
          if(M_switch==-2)  //fixed assumed F by age
           {
            M(i,y,a)=Fixed_M(a);
           }
          if(M_switch==-1)  //fixed True
           {
            M(i,y,a)=SIM_M(i,y,a);
           }
          if(M_switch==1) //est one M
           {
            M(i,y,a)=mfexp(ln_M_est);
           }
          if(M_switch==2) //est for each age, but const by region
           {
            M(i,y,a)=mfexp(ln_M_age(a));
           }           
          if(M_switch==3)  //est for each region and age
           {
            M(i,y,a)=mfexp(ln_M_region(i,a));
           }
        if(M_switch==4)  // no plus group
          {
           if(a<=max_age_recap(y))
               {
                M(i,y,a)=mfexp(ln_M_year(y,i,a)); 
               }
           if(a>max_age_recap(y))
               {
                M(i,y,a)=SIM_M(i,y,a); 
               }
             }    
 //////////////////////////////////////////////////////////////////////////////////////////////
 //////////////// Fishing Mortality //////////////////////////////////////////////////////////////
        if(F_switch==-2)  //Fix at assumed F
          {
          F(i,y,a)=Fixed_F(i,y,a);
          }
        if(F_switch==-1)  //Fix at true value
          {
          F(i,y,a)=SIM_F(i,y,a);
          }
        if(F_switch==1) //estimate
          {
           if(a<=max_age_recap(y))
               {
                F(i,y,a)=mfexp(ln_F(y,i,a)); 
               }
           if(a>max_age_recap(y))
               {
                F(i,y,a)=0; 
               }
             }
            }
           }
          }

 /////////// Total Mortality and Survival //////////////////////////////////////////////////////////////
  for(int i=1;i<=nstocks;i++) 
    {
     for(int y=1;y<=nyrs;y++) 
      {
         for (int a=1;a<=nages;a++)
          {
          Total_Survival(i,y,a)=mfexp(-(F(i,y,a)+M(i,y,a)));
          if(F(i,y,a)!=0 || M(i,y,a)!=0)
          {
          Mort_Fish(i,y,a)=F(i,y,a)*(1.-Total_Survival(i,y,a))/(F(i,y,a)+M(i,y,a)); //fraction of fish that die*fraction of deaths due to fishing
          }
          if(F(i,y,a)==0 && M(i,y,a)==0)
          {
          Mort_Fish(i,y,a)=0;
          }
          }
       }
    }
 /////////////////////////////////////////////////////////////////////////////////////////////

FUNCTION get_tag_loss

    for (int a=1;a<=nages;a++)
       {
       if(tag_loss_switch==1)  ///assume a yearly tag loss rate that is a function of true age
       {
        Tag_Loss(a)=OBS_Tag_Loss(a);
       }
       if(handling_mortality_switch==1)
       {       
       Hand_Mort(a)=OBS_Hand_Mort(a);  // function of age at release, decreases with age
       }
       if(tag_loss_switch==0)
       {
        Tag_Loss(a)=0;
        }                                           
       if(handling_mortality_switch==0)
       {       
       Hand_Mort(a)=0;  // function of age at release, decreases with age
       }
       }

FUNCTION get_tag_recaptures
  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////
  // In some calcs (mortality and movement) need to account for true age (age+time at large),
  // because multiple cohorts this means need to factor in release year
  // so true age becomes age+year-cohort(release year)
  // using subscript notation this equals= (a+y-x)
  /////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  // Plus group inherently accounted for by nages, ie, the nages mortality is used once fish reach nages
  // and the tags_avail and total_recaps are accounted for explicitly for every year and age
  // (there is no plus group per se just a constant mortality once fish reach nages)
  //////////////////////////////////////////////////////////////////////////////
  
  for(int i=1;i<=nstocks;i++) //release stock
    {
    for(int x=1;x<=ncohorts;x++)  //release year
     {
      for (int a=1;a<=nages;a++) //release age  // setup 3d array to hold age (first index..within each age is a matrix of year x stock where now use nstock x nstock to hold all recap data)
        {
        total_recap_temp_stock=0; 
         for(int y=1;y<=nyrs;y++)  //recap year
          {         
           for(int j=1;j<=nstocks;j++) //recap stock
            {

              if(y==x) //year of release for a cohort, account for time of tagging by discounting mortality 
              {
               if(fish_start_month(j)>tagging_month(i)) // tags undergo full year of fishing, but need to adjust M for mortality prior beginning of fishing season
                {
                 tags_avail(i,x,a,y,j)=ntags(i,x,a)*(1-Hand_Mort(a))*(1-Tag_Loss(a))*T(i,x,a,y,j)*mfexp(-(M(j,y,(a+y-x))*((fish_start_month(j)-tagging_month(i))/12)));  //tags released*tags don't move*tagged fish that don't die due to handling*tag fish that don't lose tags
                 pred_rec(i,x,a,y,j)=report_rate(y,j)*tags_avail(i,x,a,y,j)*F(j,y,(a+y-x))*(1.-mfexp(-(F(j,y,(a+y-x))+(M(j,y,(a+y-x))*((fish_end_month(j)-fish_start_month(j)+1)/12)))))/(F(j,y,(a+y-x))+(M(j,y,(a+y-x))*((fish_end_month(j)-fish_start_month(j)+1)/12)));  //recaps=tags available*fraction of fish that die*fraction of mortality due to fishing*tags inspected (reporting)                 
                }
               if(fish_start_month(j)<=tagging_month(i)) // tags only undergo partial/full (season<tag/season=tag) fishing, but no need to adjust M because fishing starts immediately at tag event time
                {                
                 tags_avail(i,x,a,y,j)=ntags(i,x,a)*(1-Hand_Mort(a))*(1-Tag_Loss(a))*T(i,x,a,y,j);  //tags released*tags don't move*tagged fish that don't die due to handling*tag fish that don't lose tags
                 pred_rec(i,x,a,y,j)=report_rate(y,j)*tags_avail(i,x,a,y,j)*(F(j,y,(a+y-x))*((1-(tagging_month(i)-fish_start_month(j)))/(fish_end_month(j)-fish_start_month(j)+1)))*(1.-mfexp(-((F(j,y,(a+y-x))*((1-(tagging_month(i)-fish_start_month(j)))/(fish_end_month(j)-fish_start_month(j)+1)))+(M(j,y,(a+y-x))*((fish_end_month(j)-tagging_month(i)+1)/12)))))/((F(j,y,(a+y-x))*((1-(tagging_month(i)-fish_start_month(j)))/(fish_end_month(j)-fish_start_month(j)+1)))+(M(j,y,(a+y-x))*((fish_end_month(j)-tagging_month(i)+1)/12)));  //recaps=tags available*fraction of fish that die*fraction of mortality due to fishing*tags inspected (reporting)                                  
                }
               total_recap_temp_stock(y,j)=pred_rec(i,x,a,y,j);
              }

            if(y==(x+1)) //year after tagging must discount tags avail for fraction of year underwent mortality in year of release  // must account for the maximum age so use min function to ensure that not exceeding the max age
             {
              tags_avail_temp=0;
              for(int n=1;n<=nstocks;n++)
              {               
               if(natal_homing==0)
               {
                if(fish_start_month(n)>tagging_month(i)) // tags undergo full year of fishing, but need to adjust M for mortality prior beginning of fishing season
                 {
                  tags_avail_temp(n)=tags_avail(i,x,a,y-1,n)*T(n,x,a,y,j)*(1-Tag_Loss(min((a+y-x),nages)))*mfexp(-(M(n,y-1,min(((a+y-x)-1),nages))*(1-(fish_end_month(n)/12))))*mfexp(-(F(n,y-1,min(((a+y-x)-1),nages))+(M(n,y-1,min(((a+y-x)-1),nages))*((fish_end_month(n)-fish_start_month(n)+1)/12)))); //tags_temp holds all tags moving into zone j, min function takes min of true age-1 (because mortality occurs in prev year) and nages allowed (plus group)
                 }
                if(fish_start_month(n)<=tagging_month(i)) // tags only undergo partial/full (season<tag/season=tag) fishing, but no need to adjust M because fishing starts immediately at tag event time
                 {
                  tags_avail_temp(n)=tags_avail(i,x,a,y-1,n)*T(n,x,a,y,j)*(1-Tag_Loss(min((a+y-x),nages)))*mfexp(-(M(n,y-1,min(((a+y-x)-1),nages))*(1-(fish_end_month(n)/12))))*mfexp(-((F(n,y-1,min(((a+y-x)-1),nages))*((1-(tagging_month(i)-fish_start_month(n)))/(fish_end_month(j)-fish_start_month(j)+1)))+(M(n,y-1,min(((a+y-x)-1),nages))*((fish_end_month(n)-tagging_month(i)+1)/12)))); //tags_temp holds all tags moving into zone j, min function takes min of true age-1 (because mortality occurs in prev year) and nages allowed (plus group)
                 }                                                                                                                                                              
               }
                //#####################################################################################################
                //  TRUE NATAL HOMING  T(n,x,a,y,j) becomes T(i,x,a,y,j) because need to maintain your natal origin
                //  movement values so T doesn't depend on current stock only origin stock and destination stock
                //########################################################################################################              
               if(natal_homing==1)
               {
                if(fish_start_month(n)>tagging_month(i)) // tags undergo full year of fishing, but need to adjust M for mortality prior beginning of fishing season
                 {
                  tags_avail_temp(n)=tags_avail(i,x,a,y-1,n)*T(i,x,a,y,j)*(1-Tag_Loss(min((a+y-x),nages)))*mfexp(-(M(n,y-1,min(((a+y-x)-1),nages))*(1-(fish_end_month(n)/12))))*mfexp(-(F(n,y-1,min(((a+y-x)-1),nages))+(M(n,y-1,min(((a+y-x)-1),nages))*((fish_end_month(n)-fish_start_month(n)+1)/12)))); //tags_temp holds all tags moving into zone j, min function takes min of true age-1 (because mortality occurs in prev year) and nages allowed (plus group)
                 }
                if(fish_start_month(n)<=tagging_month(i)) // tags only undergo partial/full (season<tag/season=tag) fishing, but no need to adjust M because fishing starts immediately at tag event time
                 {
                  tags_avail_temp(n)=tags_avail(i,x,a,y-1,n)*T(i,x,a,y,j)*(1-Tag_Loss(min((a+y-x),nages)))*mfexp(-(M(n,y-1,min(((a+y-x)-1),nages))*(1-(fish_end_month(n)/12))))*mfexp(-((F(n,y-1,min(((a+y-x)-1),nages))*((1-(tagging_month(i)-fish_start_month(n)))/(fish_end_month(j)-fish_start_month(j)+1)))+(M(n,y-1,min(((a+y-x)-1),nages))*((fish_end_month(n)-tagging_month(i)+1)/12)))); //tags_temp holds all tags moving into zone j, min function takes min of true age-1 (because mortality occurs in prev year) and nages allowed (plus group)
                 }               
               }               
              }
               tags_avail(i,x,a,y,j)=sum(tags_avail_temp);
               
               if(fish_start_month(j)>1) // tags undergo full year of fishing, but need to adjust M for mortality prior beginning of fishing season
                {
                 tags_avail(i,x,a,y,j)=tags_avail(i,x,a,y,j)*mfexp(-(M(j,y,min((a+y-x),nages))*((fish_start_month(j)-1)/12)));  //tags released*tags don't move*tagged fish that don't die due to handling*tag fish that don't lose tags
                }
               pred_rec(i,x,a,y,j)=report_rate(y,j)*tags_avail(i,x,a,y,j)*F(j,y,min((a+y-x),nages))*(1.-mfexp(-(F(j,y,min((a+y-x),nages))+(M(j,y,min((a+y-x),nages))*((fish_end_month(j)-fish_start_month(j)+1)/12)))))/(F(j,y,min((a+y-x),nages))+(M(j,y,min((a+y-x),nages))*((fish_end_month(j)-fish_start_month(j)+1)/12)));  //recaps=tags available*fraction of fish that die*fraction of mortality due to fishing*tags inspected (reporting)                 
               total_recap_temp_stock(y,j)=pred_rec(i,x,a,y,j);
             }

              if(y>x+1) //all other years assume full mortality  // must account for the maximum age so use min function to ensure that not exceeding the max age
              {
               tags_avail_temp=0;
                for(int n=1;n<=nstocks;n++)
                {               
                 if(natal_homing==0)
                 {
                  tags_avail_temp(n)=tags_avail(i,x,a,y-1,n)*T(n,x,a,y,j)*(1-Tag_Loss(min((a+y-x),nages)))*mfexp(-(M(n,y-1,min(((a+y-x)-1),nages))*(1-(fish_end_month(n)/12))))*mfexp(-(F(n,y-1,min(((a+y-x)-1),nages))+(M(n,y-1,min(((a+y-x)-1),nages))*((fish_end_month(n)-fish_start_month(n)+1)/12)))); //tags_temp holds all tags moving into zone j, min function takes min of true age-1 (because mortality occurs in prev year) and nages allowed (plus group)
                 }
                //#####################################################################################################
                //  TRUE NATAL HOMING  T(n,x,a,y,j) becomes T(i,x,a,y,j) because need to maintain your natal origin
                //  movement values so T doesn't depend on current stock only origin stock and destination stock
                //########################################################################################################              
                 if(natal_homing==1)
                 {
                  tags_avail_temp(n)=tags_avail(i,x,a,y-1,n)*T(i,x,a,y,j)*(1-Tag_Loss(min((a+y-x),nages)))*mfexp(-(M(n,y-1,min(((a+y-x)-1),nages))*(1-(fish_end_month(n)/12))))*mfexp(-(F(n,y-1,min(((a+y-x)-1),nages))+(M(n,y-1,min(((a+y-x)-1),nages))*((fish_end_month(n)-fish_start_month(n)+1)/12)))); //tags_temp holds all tags moving into zone j, min function takes min of true age-1 (because mortality occurs in prev year) and nages allowed (plus group)             
                 }               
                }

                 tags_avail(i,x,a,y,j)=sum(tags_avail_temp);
               
                if(fish_start_month(j)>0) // tags undergo full year of fishing, but need to adjust M for mortality prior beginning of fishing season
                  {
                   tags_avail(i,x,a,y,j)=tags_avail(i,x,a,y,j)*mfexp(-(M(j,y,min((a+y-x),nages))*((fish_start_month(j)-1)/12)));  //tags released*tags don't move*tagged fish that don't die due to handling*tag fish that don't lose tags
                  }
                 pred_rec(i,x,a,y,j)=report_rate(y,j)*tags_avail(i,x,a,y,j)*F(j,y,min((a+y-x),nages))*(1.-mfexp(-(F(j,y,min((a+y-x),nages))+(M(j,y,min((a+y-x),nages))*((fish_end_month(j)-fish_start_month(j)+1)/12)))))/(F(j,y,min((a+y-x),nages))+(M(j,y,min((a+y-x),nages))*((fish_end_month(j)-fish_start_month(j)+1)/12)));  //recaps=tags available*fraction of fish that die*fraction of mortality due to fishing*tags inspected (reporting)                 
                 total_recap_temp_stock(y,j)=pred_rec(i,x,a,y,j);
               }
             
             if(y<x) //can't have recaps before fish are released
              {
               tags_avail(i,x,a,y,j)=0;  //tags released*tags don't move*tagged fish that don't die due to handling*tag fish that don't lose tags
               pred_rec(i,x,a,y,j)=0;  //recaps=tags available*fraction of fish that die*fraction of mortality due to fishing*tags inspected (reporting)
               total_recap_temp_stock(y,j)=0;
              }              
             }
            }
             total_rec(i,x,a)=sum(total_recap_temp_stock);
             not_rec(i,x,a)=ntags(i,x,a)-total_rec(i,x,a);  //for ntags  at a given age all entries represent all tags released so can just use any of the entries (hence the i,x,a,1 subscripts)
           }
          }
         }
  for(int i=1;i<=nstocks;i++) //release stock
    {
    for(int x=1;x<=ncohorts;x++)  //release year
     {
      for (int a=1;a<=nages;a++) //release age  // setup 3d array to hold age (first index..within each age is a matrix of year x stock where now use nstock x nstock to hold all recap data)
        {
         for(int y=1;y<=nyrs;y++)  //recap year
          {         
           for(int j=1;j<=nstocks;j++) //recap stock
            {
            if(ntags(i,x,a)>0)
              {
               pred_tag_prop(i,x,a,y,j)=pred_rec(i,x,a,y,j)/ntags(i,x,a);
               pred_tag_prop_not_rec(i,x,a)=not_rec(i,x,a)/ntags(i,x,a);
              }
            if(ntags(i,x,a)==0)
              {                   
               pred_tag_prop(i,x,a,y,j)=0;
               pred_tag_prop_not_rec(i,x,a)=0;
              }
             } 
            }
           }
          }
         }

FUNCTION evaluate_the_objective_function
   f=0.0;

     for (int i=1;i<=nstocks;i++)
      {
       for(int x=1;x<=ncohorts;x++)  //release year
        {
         for(int a=1;a<=nages;a++)
          {
           for(int y=1;y<=nyrs;y++)  //recap year
            {  
             for (int j=1;j<=nstocks;j++)
              {
                T_pen(i,x,a,y,j)=0;
                M_pen(i,y,a)=0;
                F_pen(i,y,a)=0;
                report_rate_pen(j)=0;
               }
             }
           }
         }
       }
 //#################### Natural Mortality Priors   ########################################################################  # 
 if(M_pen_switch>0)
  {
   if(active(ln_M_region) || active(ln_M_age) || active(ln_M_est))
    {
     for(int i=1;i<=nstocks;i++) 
      {
       for(int y=1;y<=nyrs;y++) 
        {
         for (int a=1;a<=nages;a++)
          {
           if(M_pen_switch==1) //normal prior
            {
             M_pen(i,y,a)=log(square(M_sigma))+(square(M(i,y,a)-M_ave)/(2*square(M_sigma)));             
            }
           if(M_pen_switch==2) //penalty
            {
             if(M(i,y,a)>M_pen_high)
              {
               M_pen(i,y,a)=M_pen_mult*square(M(i,y,a)-M_pen_high);
              }
             if(M(i,y,a)<M_pen_low)
              {
               M_pen(i,y,a)=M_pen_mult*square(M(i,y,a)-M_pen_low);
              }
             if(M(i,y,a)>M_pen_low && M(i,y,a)<M_pen_high)
              {
               M_pen(i,y,a)=0;
              }
            }
          }
        }
      }
     f+=sum(M_pen);
   }
 }
 //##########################################################################################################################################

 //#################### Movement Priors   ########################################################################  # 
 if(move_pen_switch>0)
  {
   if(active(ln_T_est_AGE) || active(ln_T_est_AGE_no_age1) || active(ln_T_est) || active(ln_T_est_CNST) || active(ln_T_res) || active(ln_T_res_age))
    {
     for (int i=1;i<=nstocks;i++)
      {
       for(int x=1;x<=ncohorts;x++)  //release year
        {
         for(int a=1;a<=nages;a++)
          {
           for(int y=1;y<=nyrs;y++)  //recap year
            {  
             for (int j=1;j<=nstocks;j++)
              {
              if(y>=x)
              {
               if(movement_switch==1) //normal prior
                {
                 T_pen(i,x,a,y,j)=log(square(T_sigma))+(square(T(i,x,a,y,j)-T_ave)/(2*square(T_sigma)));
                }
               if(movement_switch==2) //penalty
                {
                 if(T(i,x,a,y,j)<T_low)
                  {
                   T_pen(i,x,a,y,j)=Tpen_mult*square(T_low-T(i,x,a,y,j));
                  }
                 if(T(i,x,a,y,j)>T_hi)
                  {
                   T_pen(i,x,a,y,j)=Tpen_mult*square(T_hi-T(i,x,a,y,j));
                  }
                 if(T(i,x,a,y,j)<T_hi && T(i,x,a,y,j)>T_low)
                  {
                   T_pen(i,x,a,y,j)=0;
                  }
                  }
                 if(no_move_first_yr==1 && y==x)
                  {
                   T_pen(i,x,a,y,j)=0;
                  }
                }
              }
            }
          }
        }
      }
    f+=sum(T_pen);
   }
  }
 //#######################################################################################################################################################
 //#################### Fishing Mortality Priors   ########################################################################  #       
  if(F_pen_switch>0) 
    {
     if(active(ln_F))
      {
       for(int i=1;i<=nstocks;i++) 
        {
         for(int y=1;y<=nyrs;y++) 
          {
           for (int a=1;a<=nages;a++)
            {
             if(F_pen_switch==1)  //==normal prior
              {
               F_pen(i,y,a)=log(square(F_sigma))+(square(F(i,y,a)-F_ave)/(2*square(F_sigma)));
              }
             if(F_pen_switch==2)  //==penalty
              {
               if(F(i,y,a)>F_pen_hi)
               {
                F_pen(i,y,a)=F_pen_mult*square(F_pen_hi-F(i,y,a));
               }
               if(F(i,y,a)<F_pen_low)
               {
                F_pen(i,y,a)=F_pen_mult*square(F_pen_low-F(i,y,a));
               }               
              }
             }
            }
           }
          f+=sum(F_pen);
         }
        }
 //#######################################################################################################################################

 //#################### Reporting Rate Priors   ########################################################################  # 
  if(report_rate_pen_switch>0) //penalize large and low values of ln_report_rate to avoid flat reponse surfaces near 0 and 1
   {
    if(active(ln_report_rate))
     {
      for (int j=1;j<=nstocks;j++) 
       {
        if(report_rate_pen_switch==1) //normal prior
         {
          report_rate_pen(j)=log(square(report_rate_sigma))+(square(report_rate(1,j)-report_rate_ave)/(2*square(report_rate_sigma)));
         }
        if(report_rate_pen_switch==2) //penalty
         {         
         if(report_rate(1,j)<report_rate_pen_low)
          {
           report_rate_pen(j)=report_rate_pen_mult*square(report_rate_pen_low-report_rate(1,j));
          }
          if(report_rate(1,j)>report_rate_pen_hi)
           {
            report_rate_pen(j)=report_rate_pen_mult*square(report_rate_pen_hi-report_rate(1,j));
           }
          if(report_rate(1,j)<report_rate_pen_hi && report_rate(1,j)>report_rate_pen_low)
           {
            report_rate_pen(j)=0;
           }
          }
         }
        f+=sum(report_rate_pen);
       }
      }
 //#################################################################################################################################################

 //####################Tag likelihood   ########################################################################  # 
  for(int i=1;i<=nstocks;i++) //release stock
    {
    for(int x=1;x<=ncohorts;x++)  //release year
     {
      for (int a=1;a<=nages;a++) //release age  // setup 3d array to hold age (first index..within each age is a matrix of year x stock where now use nstock x nstock to hold all recap data)
        {
         rss_tags_temp_age=0; 
         for(int y=1;y<=nyrs;y++)  //recap year
          {         
           for(int j=1;j<=nstocks;j++) //recap stock
            {        
             res_tags(i,x,a,y,j)=log((OBS_tag_prop(i,x,a,y,j)+TagCst))-log((pred_tag_prop(i,x,a,y,j)+TagCst)); //adjusted multinomial likelihood for tagging data
             res_tags_not_rec(i,x,a)=-log((pred_tag_prop_not_rec(i,x,a)+TagCst))+log((OBS_tag_prop_not_rec(i,x,a)+TagCst));
             rss_tags_temp(i,x,a,y,j)=OBS_tag_prop(i,x,a,y,j)*res_tags(i,x,a,y,j);
             rss_tags_temp_not_rec(i,x,a)=OBS_tag_prop_not_rec(i,x,a)*res_tags_not_rec(i,x,a);
             rss_tags_temp_age(y,j)=rss_tags_temp(i,x,a,y,j);

             if (tag_neff_switch==1)
              {
               tag_neff(i,x,a)=ntags(i,x,a);               
              }  
              if (tag_neff_switch==0)
              {
               tag_neff(i,x,a)=neff(i,x,a);         
              }
            }
           }
          rss_tags_temp_stock(i,x,a)=sum(rss_tags_temp_age);
         }
       }
     }

  rss_tags=sum(elem_prod(tag_neff,(rss_tags_temp_stock+rss_tags_temp_not_rec)));

  f+=rss_tags;

REPORT_SECTION

  report<<"$max_gradient"<<endl;
  report<<objective_function_value::gmax<<endl;

 //////// Put T into matrix for report (report always seems to ruin 3darray on printing and can't use to read into R same thing done for tags_avail and Pred_proportions ////////////////////
  for (int i=1;i<=nstocks;i++)
   {
   for(int x=1;x<=ncohorts;x++)  //release year
    {
    for(int a=1;a<=nages;a++)
    {
     for(int y=1;y<=nyrs;y++)  //recap year
      {  
       for (int j=1;j<=nstocks;j++)
        {
        if(y>x) //don't take first year where movement is 0 if no_move_first_yr==1
        {
         T_region(i,a,j)=T(i,x,a,y,j);
         SIM_T_region(i,a,j)=SIM_T(i,x,a,y,j);
        }
        }
       }
      }
     }
    }

 ////////////////////////////////////////////////////////////////
  report<<"$SIM_T_region"<<endl;
  report<<SIM_T_region<<endl;
  report<<"$T_region"<<endl;
  report<<T_region<<endl;
  report<<"$SIM_M"<<endl;
  report<<SIM_M<<endl;
  report<<"$M"<<endl;
  report<<M<<endl;
  report<<"$SIM_F"<<endl;
  report<<SIM_F<<endl;
  report<<"$F"<<endl;
  report<<F<<endl;
  report<<"$SIM_report_rate"<<endl;
  report<<SIM_report_rate<<endl;
  report<<"$report_rate_CNST"<<endl;
  report<<report_rate_CNST<<endl;
  
  report<<"$nages"<<endl;
  report<<nages<<endl;
  report<<"$ncohorts"<<endl;
  report<<ncohorts<<endl;
  report<<"$nyrs"<<endl;
  report<<nyrs<<endl;
  report<<"$nstocks"<<endl;
  report<<nstocks<<endl;
  report<<"$ages"<<endl;
  report<<ages<<endl;
  report<<"$years"<<endl;
  report<<yrs<<endl;
  report<<"$tagging_month"<<endl;
  report<<tagging_month<<endl;
  report<<"$max_age_recap"<<endl;
  report<<max_age_recap<<endl;
  
  report<<"$phase_T_est"<<endl;
  report<<phase_T_est<<endl;
  report<<"$phase_T_est_CNST"<<endl;
  report<<phase_T_est_CNST<<endl;
  report<<"$phase_T_est_AGE"<<endl;
  report<<phase_T_est_AGE<<endl;
  report<<"$phase_T_est_AGE_no_age1"<<endl;
  report<<phase_T_est_AGE_no_age1<<endl;
  report<<"$phase_T_res"<<endl;
  report<<phase_T_res<<endl;
  report<<"$phase_T_res_age"<<endl;
  report<<phase_T_res_age<<endl;
  report<<"$phase_T_res_no_age1"<<endl;
  report<<phase_T_res_no_age1<<endl;
  report<<"$phase_F"<<endl;
  report<<phase_F<<endl;
  report<<"$phase_report_rate"<<endl;
  report<<phase_report_rate<<endl;
  report<<"$phase_M_region"<<endl;
  report<<phase_M_region<<endl;
  report<<"$phase_M_age"<<endl;
  report<<phase_M_age<<endl;
  report<<"$phase_M_CNST"<<endl;
  report<<phase_M_CNST<<endl;
  report<<"$phase_M_year"<<endl;
  report<<phase_M_year<<endl;

  report<<"$M_switch"<<endl;
  report<<M_switch<<endl;
  report<<"$F_switch"<<endl;
  report<<F_switch<<endl;
  report<<"$report_rate_switch"<<endl;
  report<<report_rate_switch<<endl;
  report<<"$SIM_mixing_switch"<<endl;
  report<<SIM_mixing_switch<<endl;
  report<<"$movement_switch"<<endl;
  report<<movement_switch<<endl;
  report<<"$no_move_first_yr"<<endl;
  report<<no_move_first_yr<<endl;
  report<<"$natal_homing"<<endl;
  report<<natal_homing<<endl;
  report<<"$move_pen_switch"<<endl;
  report<<move_pen_switch<<endl;
  report<<"$tag_neff_switch"<<endl;
  report<<tag_neff_switch<<endl;
  report<<"$tag_loss_switch"<<endl;
  report<<tag_loss_switch<<endl;
  report<<"$handling_mortality_switch"<<endl;
  report<<handling_mortality_switch<<endl;
  
  report<<"$T_low"<<endl;
  report<<T_low<<endl;
  report<<"$T_hi"<<endl;
  report<<T_hi<<endl;
  report<<"$neff"<<endl;
  report<<neff<<endl;

  report<<"$ntags"<<endl;
  report<<ntags<<endl;

  report<<"$OBS_tag_rec"<<endl;
  report<<OBS_tag_rec<<endl;
  report<<"$OBS_tag_not_rec"<<endl;
  report<<OBS_tag_not_rec<<endl;
  
  report<<"$OBS_total_rec"<<endl;
  report<<OBS_total_rec<<endl;
  report<<"$OBS_tag_not_rec"<<endl;
  report<<OBS_tag_not_rec<<endl;
  report<<"$OBS_tag_prop"<<endl;
  report<<OBS_tag_prop<<endl;
  report<<"$OBS_tag_prop_not_rec"<<endl;
  report<<OBS_tag_prop_not_rec<<endl;

  report<<"$ln_M_est"<<endl;
  report<<ln_M_est<<endl;
  report<<"$ln_M_age"<<endl;
  report<<ln_M_age<<endl;
  report<<"$ln_M_region"<<endl;
  report<<ln_M_region<<endl;
  report<<"$M"<<endl;
  report<<M<<endl;

  report<<"$SIM_T1"<<endl;
  report<<SIM_T<<endl;
  report<<"$T_region"<<endl;
  report<<T_region<<endl;
  report<<"$ln_T_est"<<endl;
  report<<ln_T_est<<endl;
  report<<"$ln_T_est_CNST"<<endl;
  report<<ln_T_est_CNST<<endl;
  report<<"$ln_T_est_AGE"<<endl;
  report<<ln_T_est_AGE<<endl;
  report<<"$ln_T_est_AGE_no_age1"<<endl;
  report<<ln_T_est_AGE_no_age1<<endl;
  report<<"$ln_T_res"<<endl;
  report<<ln_T_res<<endl;
  report<<"$ln_T_res_age"<<endl;
  report<<ln_T_res_age<<endl;
  report<<"$report_rate"<<endl;
  report<<report_rate<<endl;
  report<<"$ln_report_rate"<<endl;
  report<<ln_report_rate<<endl;  
 // report<<"$ln_F"<<endl;
 // report<<ln_F<<endl;
  report<<"$F"<<endl;
  report<<F<<endl;
  report<<"$Total_Survival"<<endl;
  report<<Total_Survival<<endl;
  report<<"$Mort_Fish"<<endl;
  report<<Mort_Fish<<endl;
  report<<"$report_rate"<<endl;
  report<<report_rate<<endl;
  report<<"$Tag_Loss"<<endl;
  report<<Tag_Loss<<endl;
  report<<"$Hand_Mort"<<endl;
  report<<Hand_Mort<<endl;

  report<<"$total_rec"<<endl;
  report<<total_rec<<endl;
  report<<"$not_rec"<<endl;
  report<<not_rec<<endl;


  report<<"$F_pen"<<endl;
  report<<F_pen<<endl;
  report<<"$M_pen"<<endl;
  report<<M_pen<<endl;
  report<<"$report_rate_pen"<<endl;
  report<<report_rate_pen<<endl;
  
  report<<"$tag_neff"<<endl;
  report<<tag_neff<<endl;


  report<<"$res_tags_not_rec"<<endl;
  report<<res_tags_not_rec<<endl;
  report<<"$rss_tags"<<endl;
  report<<rss_tags<<endl;

  report<<"$Obj_Func"<<endl;
  report<<f<<endl;


 ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 // ADD WHEN NEXT VERSION OF ADMB RELEASED
 //////////////////////////////////////////////////////////////////////////////////////////////////////
 // report<<"$T"<<endl;
 // report<<T<<endl;
 // report<<"$tags_avail"<<endl;
 //  report<<tags_avail<<endl;
 //  report<<"$pred_rec"<<endl;
 //  report<<pred_rec<<endl;

  //report<<"$pred_tag_prop"<<endl;
  //report<<pred_tag_prop<<endl;
  //report<<"$pred_tag_prop_not_rec"<<endl;
  //report<<pred_tag_prop_not_rec<<endl;

  //report<<"$T_pen"<<endl;
  //report<<T_pen<<endl;
  //report<<"$res_tags"<<endl;
  //report<<res_tags<<endl;
  
 ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 // END Addition
 ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 
 ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 // REMOVE WHEN NEXT VERSION OF ADMB RELEASED
 //////////////////////////////////////////////////////////////////////////////////////////////////////
 if(nstocks>3)
 {
  for(int i=1;i<=nstocks;i++) //release stock
    {
    for(int x=1;x<=ncohorts;x++)  //release year
     {
      for (int a=1;a<=nages;a++) //release age  // setup 3d array to hold age (first index..within each age is a matrix of year x stock where now use nstock x nstock to hold all recap data)
        {
         for(int y=1;y<=nyrs;y++)  //recap year
          {         
           for(int j=1;j<=nstocks;j++) //recap stock
            {
            T_reg1(x,a,y,j)=T(1,x,a,y,j);
            T_reg2(x,a,y,j)=T(2,x,a,y,j);
            T_reg3(x,a,y,j)=T(3,x,a,y,j);
            T_reg4(x,a,y,j)=T(4,x,a,y,j);

            tags_avail1(x,a,y,j)=tags_avail(1,x,a,y,j);
            tags_avail2(x,a,y,j)=tags_avail(2,x,a,y,j);
            tags_avail3(x,a,y,j)=tags_avail(3,x,a,y,j);
            tags_avail4(x,a,y,j)=tags_avail(4,x,a,y,j);
            pred_rec1(x,a,y,j)=pred_rec(1,x,a,y,j);
            pred_rec2(x,a,y,j)=pred_rec(2,x,a,y,j);
            pred_rec3(x,a,y,j)=pred_rec(3,x,a,y,j);
            pred_rec4(x,a,y,j)=pred_rec(4,x,a,y,j);
            pred_tag_prop1(x,a,y,j)=pred_tag_prop(1,x,a,y,j);
            pred_tag_prop2(x,a,y,j)=pred_tag_prop(2,x,a,y,j);
            pred_tag_prop3(x,a,y,j)=pred_tag_prop(3,x,a,y,j);
            pred_tag_prop4(x,a,y,j)=pred_tag_prop(4,x,a,y,j);
            res_tags1(x,a,y,j)=res_tags(1,x,a,y,j);
            res_tags2(x,a,y,j)=res_tags(2,x,a,y,j);
            res_tags3(x,a,y,j)=res_tags(3,x,a,y,j);
            res_tags4(x,a,y,j)=res_tags(4,x,a,y,j);
            }
          }
        }
      }
    }

  report<<"$T_reg1"<<endl;
  report<<T_reg1<<endl;
  report<<"$T_reg2"<<endl;
  report<<T_reg2<<endl;
  report<<"$T_reg3"<<endl;
  report<<T_reg3<<endl;
  report<<"$T_reg4"<<endl;
  report<<T_reg4<<endl;

  report<<"$tags_avail1"<<endl;
  report<<tags_avail1<<endl;
  report<<"$tags_avail2"<<endl;
  report<<tags_avail2<<endl;
  report<<"$tags_avail3"<<endl;
  report<<tags_avail3<<endl;
  report<<"$tags_avail4"<<endl;
  report<<tags_avail4<<endl;
  report<<"$pred_rec1"<<endl;
  report<<pred_rec1<<endl;
  report<<"$pred_rec2"<<endl;
  report<<pred_rec2<<endl;
  report<<"$pred_rec3"<<endl;
  report<<pred_rec3<<endl;
  report<<"$pred_rec4"<<endl;
  report<<pred_rec4<<endl;
  report<<"$pred_tag_prop1"<<endl;
  report<<pred_tag_prop1<<endl;
  report<<"$pred_tag_prop2"<<endl;
  report<<pred_tag_prop2<<endl;
  report<<"$pred_tag_prop3"<<endl;
  report<<pred_tag_prop3<<endl;
  report<<"$pred_tag_prop4"<<endl;
  report<<pred_tag_prop4<<endl;
  report<<"$res_tags1"<<endl;
  report<<res_tags1<<endl;
  report<<"$res_tags2"<<endl;
  report<<res_tags2<<endl;
  report<<"$res_tags3"<<endl;
  report<<res_tags3<<endl;
  report<<"$res_tags4"<<endl;
  report<<res_tags4<<endl;
  }
 ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 // END REMOVE
 ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 
RUNTIME_SECTION
  convergence_criteria .001,.0001, 1.0e-4
  maximum_function_evaluations 100000
  
