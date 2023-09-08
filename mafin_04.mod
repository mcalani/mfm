%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       A Macro Financial Model for the Chilean Economy         %
%                       	 MAFIN                          	%
%                       BCCH GEE-GEFIN                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

// Note: must use Dynare 4.6.0 or above to run code. This enables comments after macro-processing commands (in the same line).

%---------------------------------------------------------------%
% 0. Housekeeping and options
%---------------------------------------------------------------%

close all;

@#define use_estimation = 1	  // = 1 to estimate the model, =0 to perform a stochastic simulation
@#define fin_fric = 1         // = 1 includes financial frictions and is the baseline model. =0 disables all fin. frictions, modifying model block and steady state calculation
@#define NH = 6               // Time-to-build periods before authorized housing projects increase housing stock. =0 implies there is no time to build <-> housing producers are symmetric to capital producers
@#define QLC = 1              // = 1 to include Quadratic Labor Adjustment Costs as in Lechthaler and Snower (2010).
@#define hist_decomp    = 0   // = 1 to performs historical decomposition of shocks
@#define ccyb_rule      = 0   // = 0 base (no ccyb); = [1 - 5] specify the ccyb rule (save_output==2)
@#define taylor_rule    = 0   // = 0 base ; = [1 - 4] specify the taylor rule        (save_output==3) 

%%% Outputs
@#define save_output    = 99   // =0 high/low phi; =1 FF/No-FF ; = 2 CCyB exercises; = 3 Taylor rule exercises ; = 4 CCyB-Bands exercises
@#define ss_values      = 1   // =1 export ss-values for comparative statics analysis in '..\Estatica_comparativa\SS_values.mat' file
@#define export_moments = 0   // =1 to export theoretical moments and save them in 'PLOTS' folder
@#define welfare_analysis = 1     // = 0 no; 1 = si (lento)

%---------------------------------------------------------------%
% 1. Definitions
%---------------------------------------------------------------%

@#include "variables.inc"
@#include "parameters.inc"

%---------------------------------------------------------------%
% 2. Calibration
%---------------------------------------------------------------%

	// Options used as parameters
FinFric	= @{fin_fric};

@#include "calibration.inc"
@#if save_output==4
    @#include "calibration_ex4.inc"
@#endif


%---------------------------------------------------------------%
% 3. Model
%---------------------------------------------------------------%

@#include "model.inc"

%---------------------------------------------------------------%
% 4. Steady State
%---------------------------------------------------------------%
//ss_params_script

resid;

steady;

//model_diagnostics;
//check;

%---------------------------------------------------------------%
% 5. Computation
%---------------------------------------------------------------%

	// Estimation or Stochastic Simulation
@#if use_estimation

    @#include "estimation.inc"
    //@#include "estimation_rpf.inc"

	@#if hist_decomp
		@#include "hist_decomp.inc"
        OUTPUT.M_ = M_;
        OUTPUT.oo_ = oo_;
        OUTPUT.oo_recursive_ = oo_recursive_;
        OUTPUT.options_ = options_;
        // save('PLOTS\Exercises\Ex3_Shocks_decom\output_shocks.mat','-struct','OUTPUT') 
	@#endif

@#else

	shocks;

	@#for va in vars_exo
		var u_@{va}	=	SIGMA_@{va}^2;
	@#endfor

    end;

	%stoch_simul(order=1, nograph, nomoments, nocorr, nodecomposition, nofunctions, ar=0, periods=0);
	%stoch_simul(order=1,  irf=25, nograph) gdp, ppi, R;
	stoch_simul(order=1, periods=0, irf=40, irf_shocks =(u_e_req,u_e_m), graph) gdp, ppi, l_f,l_h,q_l,phi_f,phi_h,ccyb,e_req,e_m;
@#endif

stoch_simul(order=1, nograph, nomoments, nocorr, nodecomposition, nofunctions, ar=0, periods=0);

@#if save_output==0
	stoch_simul(order=1, periods=0, irf=40,nograph);
    OUTPUT.bayestopt_ = bayestopt_ ;
    OUTPUT.dataset_ = dataset_ ;
    OUTPUT.dataset_info = dataset_info;
    OUTPUT.estim_params_ = estim_params_;
    OUTPUT.estimation_info = estimation_info;
    OUTPUT.M_ = M_;
    OUTPUT.oo_ = oo_;
    //OUTPUT.oo_recursive_ = oo_recursive_;
    OUTPUT.options_ = options_;
    @#if fin_fric==1
        //save('PLOTS\Exercises\Ex0_phi\output_phi08_ccyb2_taylor0.mat','-struct','OUTPUT') ;
        //save('PLOTS\Exercises\Ex0_phi\output_phi045.mat','-struct','OUTPUT') ;
        //save('PLOTS\Exercises\Ex0_phi\output_phi105.mat','-struct','OUTPUT') ;
        //save('PLOTS\Exercises\Ex0_phi\output_phi145.mat','-struct','OUTPUT') ;
    @#elseif fin_fric==0
        save('TEMP\output_TEMPsinFF.mat','-struct','OUTPUT') ;
    @#endif
    //run('PLOTS\doing_plots.m')

@#elseif save_output==1
	stoch_simul(order=1, periods=0, irf=40,nograph);
    OUTPUT.bayestopt_ = bayestopt_ ;
    OUTPUT.dataset_ = dataset_ ;
    OUTPUT.dataset_info = dataset_info;
    OUTPUT.estim_params_ = estim_params_;
    OUTPUT.estimation_info = estimation_info;
    OUTPUT.M_ = M_;
    OUTPUT.oo_ = oo_;
    //OUTPUT.oo_recursive_ = oo_recursive_;
    OUTPUT.options_ = options_;
    @#if fin_fric==1
        save('PLOTS\Ex1_FFyNoFF\output_FF.mat','-struct','OUTPUT')
    @#elseif fin_fric==0
        save('PLOTS\Ex1_FFyNoFF\output_sinFF.mat','-struct','OUTPUT')
    @#endif
    //run('PLOTS\doing_plots.m')

@#elseif save_output==2 
	stoch_simul(order=1, periods=0, irf=40,nograph);
    OUTPUT.bayestopt_      = bayestopt_ ;
    OUTPUT.dataset_        = dataset_ ;
    OUTPUT.dataset_info    = dataset_info;
    OUTPUT.estim_params_   = estim_params_;
    OUTPUT.estimation_info = estimation_info;
    OUTPUT.M_ = M_;
    OUTPUT.oo_ = oo_;
    //OUTPUT.oo_recursive_ = oo_recursive_;
    OUTPUT.options_ = options_;
    @#if     ccyb_rule==0
        save('PLOTS\Exercises\Ex2_CCyB_specs\output_ccyb0.mat','-struct','OUTPUT')

    @#elseif ccyb_rule==1
        save('PLOTS\Exercises\Ex2_CCyB_specs\output_ccyb1.mat','-struct','OUTPUT')

    @#elseif ccyb_rule==2
        save('PLOTS\Exercises\Ex2_CCyB_specs\output_ccyb2.mat','-struct','OUTPUT')

    @#elseif ccyb_rule==3
        save('PLOTS\Exercises\Ex2_CCyB_specs\output_ccyb3.mat','-struct','OUTPUT')

    @#elseif ccyb_rule==4
        save('PLOTS\Exercises\Ex2_CCyB_specs\output_ccyb4.mat','-struct','OUTPUT')

    @#elseif ccyb_rule==5
        save('PLOTS\Exercises\Ex2_CCyB_specs\output_ccyb5.mat','-struct','OUTPUT')
    @#endif
    
@#elseif save_output==3 
	stoch_simul(order=1, periods=0, irf=40,nograph);
    OUTPUT.bayestopt_      = bayestopt_     ;
    OUTPUT.dataset_        = dataset_       ;
    OUTPUT.dataset_info    = dataset_info   ;
    OUTPUT.estim_params_   = estim_params_  ;
    OUTPUT.estimation_info = estimation_info;
    OUTPUT.M_ = M_;
    OUTPUT.oo_ = oo_;
    //OUTPUT.oo_recursive_ = oo_recursive_;
    OUTPUT.options_ = options_;
    @#if     taylor_rule==0

        @#if ccyb_rule    == 0 // No ccyb
            save('PLOTS\Exercises\Ex3_TaylorRule_specs\output_taylor_noCCyB.mat','-struct','OUTPUT')
        @#elseif ccyb_rule    == 3 // ccyb (PD) 
            save('PLOTS\Exercises\Ex3_TaylorRule_specs\output_taylor_CCyB_PD.mat','-struct','OUTPUT')
        @#elseif ccyb_rule    == 5 // ccyb (PD) 
            save('PLOTS\Exercises\Ex3_TaylorRule_specs\output_taylor_CCyB_PDFL.mat','-struct','OUTPUT')
        @#else // ccyb (credit/gap)
            save('PLOTS\Exercises\Ex3_TaylorRule_specs\output_taylor0.mat','-struct','OUTPUT')
        @#endif

    @#elseif taylor_rule==1
        save('PLOTS\Exercises\Ex3_TaylorRule_specs\output_taylor1.mat','-struct','OUTPUT')

    @#elseif taylor_rule==2
        save('PLOTS\Exercises\Ex3_TaylorRule_specs\output_taylor2.mat','-struct','OUTPUT')

    @#elseif taylor_rule==3
        save('PLOTS\Exercises\Ex3_TaylorRule_specs\output_taylor3.mat','-struct','OUTPUT')

    @#elseif taylor_rule==4
        @#if ccyb_rule    == 0 
            save('PLOTS\Exercises\Ex3_TaylorRule_specs\output_pers_taylor_noCCyB.mat','-struct','OUTPUT')
        
        @#elseif ccyb_rule    == 1 // ccyb (credit)
            save('PLOTS\Exercises\Ex3_TaylorRule_specs\output_pers_taylor0.mat','-struct','OUTPUT')

        @#elseif ccyb_rule    == 2 // ccyb (credit/gap)
            save('PLOTS\Exercises\Ex3_TaylorRule_specs\output_pers_taylor0.mat','-struct','OUTPUT')

        @#elseif ccyb_rule    == 3 // ccyb (PD) 
            save('PLOTS\Exercises\Ex3_TaylorRule_specs\output_pers_taylor_CCyB_PD.mat','-struct','OUTPUT')

        @#elseif ccyb_rule    == 4 // ccyb (shock) 
            save('PLOTS\Exercises\Ex3_TaylorRule_specs\output_pers_taylor_CCyB_shock.mat','-struct','OUTPUT')

        @#elseif ccyb_rule    == 5 // ccyb (PD FL) 
            save('PLOTS\Exercises\Ex3_TaylorRule_specs\output_pers_taylor_CCyB_PDFL.mat','-struct','OUTPUT')
        @#endif

    @#elseif taylor_rule==5
        save('PLOTS\Exercises\Ex3_TaylorRule_specs\output_taylor5.mat','-struct','OUTPUT')
    
    @#endif

@#elseif save_output==4 
	stoch_simul(order=1, periods=0, irf=40,nograph);
    OUTPUT.bayestopt_      = bayestopt_     ;
    OUTPUT.dataset_        = dataset_       ;
    OUTPUT.dataset_info    = dataset_info   ;
    OUTPUT.estim_params_   = estim_params_  ;
    OUTPUT.estimation_info = estimation_info;
    OUTPUT.M_ = M_;
    OUTPUT.oo_ = oo_;
    //OUTPUT.oo_recursive_ = oo_recursive_;
    OUTPUT.options_ = options_;

    index_aux=num2str( index_ex4 );
    if index2_ex4 == 0.08
        save(['PLOTS\Exercises\Ex4_CCyB_Bands\output_CCyB_' index_aux([1 3:end]) '.mat'],'-struct','OUTPUT');
    elseif index2_ex4==0.105
        save(['PLOTS\Exercises\Ex4_CCyB_Bands\output_CCyB_phi105_' index_aux([1 3:end]) '.mat'],'-struct','OUTPUT');
    end
@#endif


@#if ss_values==1
    SS_values.param_names =  M_.param_names
    SS_values.params      =  M_.params
    //save('..\Estatica_comparativa\SS_values.mat','-struct','SS_values')
    save('SS_values.mat','-struct','SS_values')
@#endif


@#if export_moments==1 // not finished
    OUTPUT.oo_      = oo_                             ;
    MOMENTS.Var     = OUTPUT.oo_.var                  ;
    MOMENTS.Sigma   = diag(diag(MOMENTS.Var).^(1/2))  ;
    MOMENTS.Corr    = inv(MOMENTS.Sigma)*MOMENTS.Var*inv(MOMENTS.Sigma)  ;
    MOMENTS.Mean    = oo_.mean  ;
    MOMENTS.Sigma   = diag(MOMENTS.Sigma)  ;
@#endif