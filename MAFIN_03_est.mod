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

%%% Inputs
@#define use_estimation = 1     // =1 to estimate the model, =0 to perform a stochastic simulation
@#define fin_fric = 1           // =1 includes financial frictions and is the baseline model. =0 disables all fin. frictions, modifying model block and steady state calculation
@#define NH = 6                 // Time-to-build periods before authorized housing projects increase housing stock. =0 implies there is no time to build <-> housing producers are symmetric to capital producers
@#define QLC = 1                // =1 to include Quadratic Labor Adjustment Costs as in Lechthaler and Snower (2010).
@#define hist_decomp = 0        // = 1 to performs historical decomposition of shocks
@#define ccyb_rule      = 0     // = 0 base (no ccyb); = [1 - 5] specify the ccyb rule (save_output==2)
@#define taylor_rule    = 0     // = 0 base ; = [1 - 4] specify the taylor rule        (save_output==3) 

%%% Outputs 
@#define save_output    = 99    // =0 high/low phi; =1 FF/No-FF ; = 2 CCyB exercises; = 3 Taylor rule exercises
@#define ss_values      = 1     // =1 export ss-values for comparative statics analysis in '..\Estatica_comparativa\SS_values.mat' file
@#define export_moments = 0     // =1 to export theoretical moments and save them in 'PLOTS' folder
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

%---------------------------------------------------------------%
% 3. Model
%---------------------------------------------------------------%

@#include "model.inc"

%---------------------------------------------------------------%
% 4. Steady State
%---------------------------------------------------------------%

resid;          % compute residuals of ss equations using .m file 


%steady;         % compute ss numerically (or correct) 

%---------------------------------------------------------------%
% 5. Computation
%---------------------------------------------------------------%

check;          % compute policy and transition functions

	// Estimation or Stochastic Simulation
@#if use_estimation
    @#if fin_fric==1
        @#include "estimation.inc"
        //@#include "estimation_rpf.inc"
    @#elseif fin_fric==0
        @#include "estimation_noFF.inc"
    @#endif

	@#if hist_decomp
		@#include "hist_decomp.inc"
	@#endif
@#else
	shocks;
	@#for va in vars_exo
		var u_@{va}	=	SIGMA_@{va}^2;
	@#endfor
    end;
    
	%stoch_simul(order=1, nograph, nomoments, nocorr, nodecomposition, nofunctions, ar=0, periods=0);
	%stoch_simul(order=1,  irf=25, nograph) gdp, ppi, R;
	%stoch_simul(order=1, periods=0, irf=40, nograph);
@#endif

stoch_simul(order=1, nograph, nomoments, nocorr, nodecomposition, nofunctions, ar=0, periods=0);

@#if save_output==1
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
        save('PLOTS\output_FF.mat','-struct','OUTPUT')
    @#elseif fin_fric==0
        save('PLOTS\output_sinFF.mat','-struct','OUTPUT')
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
        save('PLOTS\Exercises\Ex1_CCyB_specs\output_ccyb0.mat','-struct','OUTPUT')

    @#elseif ccyb_rule==1
        save('PLOTS\Exercises\Ex1_CCyB_specs\output_ccyb1.mat','-struct','OUTPUT')

    @#elseif ccyb_rule==2
        save('PLOTS\Exercises\Ex1_CCyB_specs\output_ccyb2.mat','-struct','OUTPUT')

    @#elseif ccyb_rule==3
        save('PLOTS\Exercises\Ex1_CCyB_specs\output_ccyb3.mat','-struct','OUTPUT')

    @#elseif ccyb_rule==4
        save('PLOTS\Exercises\Ex1_CCyB_specs\output_ccyb4.mat','-struct','OUTPUT')

    @#elseif ccyb_rule==5
        save('PLOTS\Exercises\Ex1_CCyB_specs\output_ccyb5.mat','-struct','OUTPUT')
    @#endif
    //run('PLOTS\doing_ccyb_plots.m')
@#endif

@#if export_moments==1
    OUTPUT.oo_      = oo_ ;
    MOMENTS.Var     = OUTPUT.oo_.var  ;
    MOMENTS.Sigma   = diag(diag(MOMENTS.Var).^(1/2))  ;
    MOMENTS.Corr    = inv(MOMENTS.Sigma)*MOMENTS.Var*inv(MOMENTS.Sigma)  ;
    MOMENTS.Mean    = oo_.mean  ;
    MOMENTS.Sigma   = diag(MOMENTS.Sigma)  ;
@#endif

@#if ss_values==1
disp('Grabando SS_valuesDOTmat')
    SS_values.param_names =  M_.param_names
    SS_values.params      =  M_.params
    save('.\Estatica_comparativa\SS_values.mat','-struct','SS_values')
    save('SS_values.mat','-struct','SS_values')
@#endif

