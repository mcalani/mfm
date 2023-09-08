
% =========================================================================
% MAFIN Estimation
% =========================================================================
addpath('./DynarePrunning2016')
addpath('./_aux')

% -------------------------------------------------------------------------
	// Declare observables
varobs	gam_YR_obs 
        gam_Y_obs 
		gam_YCo_obs
		gam_C_obs        
		gam_IH_obs
		gam_IK_obs
		gam_G_obs
		gam_Ystar_obs
		gam_WN_obs 
				
		piZ_obs
		piM_obs
		pistar_obs
		piCostar_obs
		// piQH_obs
		gam_lf_obs
		gam_lh_obs
		
		R_obs
		Rstar_obs
		RLG_obs
		RD_obs
		RI_obs
		RL_obs
		xi_obs
		Roe_obs
		rer_obs
		stb_obs
		gam_N_obs
        ccyb_obs
		;

% -------------------------------------------------------------------------
	// Structural shocks
shocks;
@#for va in vars_exo
	var u_@{va}=SIGMA_@{va}^2;
@#endfor

% -------------------------------------------------------------------------  
	// Measurement errors
var gam_C_obs		;stderr	sqrt(1.062022624^2*0.1);
var	gam_G_obs		;stderr	sqrt(1.477989947^2*0.1);
var	gam_IH_obs		;stderr	sqrt(4.182871609^2*0.1);
var	gam_IK_obs		;stderr	sqrt(4.423631416^2*0.1);
var	gam_N_obs		;stderr	sqrt(1.339816319^2*0.1);
var	gam_WN_obs		;stderr	sqrt(0.455175549^2*0.1);
var	gam_YCo_obs		;stderr	sqrt(3.187357826^2*0.1);
var	gam_YR_obs		;stderr	sqrt(0.895501894^2*0.1);
var	gam_Y_obs		;stderr	sqrt(0.895501894^2*0.1);
var	gam_Ystar_obs	;stderr	sqrt(0.502379368^2*0.1);
var	piZ_obs			;stderr	sqrt(0.476537278^2*0.1);
var	piCostar_obs	;stderr	sqrt(13.47942165^2*0.1);
var	piM_obs			;stderr	sqrt(4.611897685^2*0.1);
//var	piQH_obs		;stderr	sqrt(1.379183027^2*0.1);
var	pistar_obs		;stderr	sqrt(2.753768226^2*0.1);
var	rer_obs			;stderr	sqrt(4.874002666^2*0.1);
var	RLG_obs			;stderr	sqrt(0.223641762^2*0.1);
var	R_obs			;stderr	sqrt(0.389344692^2*0.1);
var	RD_obs			;stderr	sqrt(0.389732641^2*0.1);
var	RI_obs			;stderr	sqrt(0.226138484^2*0.1);
var	RL_obs			;stderr	sqrt(0.552410384^2*0.1);
var	Rstar_obs		;stderr	sqrt(0.378975051^2*0.1);
var	stb_obs			;stderr	sqrt(4.844406331^2*0.1);
var	xi_obs			;stderr	sqrt(0.133787913^2*0.1);
var gam_lf_obs	    ;stderr	sqrt(1.796154931^2*0.1);
var gam_lh_obs	    ;stderr	sqrt(1.300762541^2*0.1);
var Roe_obs			;stderr	sqrt(0.654331862^2*0.1);
end;

% -------------------------------------------------------------------------
	// Priors
estimated_params;

% -------------------------------------------------------------------------------------------------------
%   PARAM_NAME, INIT_VAL,  L_BOUND,  U_BOUND,	PRIOR_SHAPE,    PRIOR_MEAN,	PRIOR_STANDARD_ERROR
% -------------------------------------------------------------------------------------------------------

	ALPHA_ppi,		,			1,		,		normal_pdf,		1.7,        0.1;
	ALPHA_R,		,			,		,		beta_pdf,		0.85,       0.025;
	ALPHA_W,		,			,		,		beta_pdf,		0.25,       0.075;
    ALPHA_y,		,			0,		,		normal_pdf,		0.125,      0.075;
    ETA,			0.99,		,		,		gamma_pdf,		1.0,        0.25; 
	ETA_CHAT,		0.99,       ,		,		gamma_pdf,		1.0,        0.25;
	ETA_ast,		,			,		,		gamma_pdf,		0.25,       0.075;
    GAMA_L,			0.05,			,		,		gamma_pdf,		0.1,       0.09;   
    GAMA_h,			,			,		,		gamma_pdf,		3.0,        0.25;   
	GAMA_k,			,			,		,		gamma_pdf,		3.0,        0.25; 	
	GAMA_n,			,			,		,		gamma_pdf,		3.0,        0.25;	   
    KAPPA_F,		,			,		,		beta_pdf,		0.5,        0.075;
	KAPPA_H,		,			,		,		beta_pdf,		0.5,        0.075;
    KAPPA_W,		,			,		,		beta_pdf,		0.85,       0.025;	
	PHI_ast,		0.99,			,		,	inv_gamma_pdf,	1,          inf;
	PHI_c,			,			,		,		beta_pdf,		0.85,       0.025;	
    PHI_HH,			,			,		,		beta_pdf,		0.85,       0.025;	 
	THETA_F,		,			,		,		beta_pdf,		0.5,        0.075;
	THETA_H,		,			,		,		beta_pdf,		0.5,        0.025;
	THETA_W,		,			,		,		beta_pdf,		0.5,        0.075;
	VARPHI,			,			,		,		gamma_pdf,		7.5,        1.5;
	MU_e,			,			,		,		beta_pdf,		0.3,        0.05;
    MU_f,			,			,		,		beta_pdf,		0.3,        0.05;
    MU_h,			,			,		,		beta_pdf,		0.3,        0.05;
    MU_i,			,			,		,		beta_pdf,		0.3,        0.05;
	ETA_ZETA_L,		,			,		,		gamma_pdf,		15/100,		15/100/5;

	RHO_a,			,			,		,		beta_pdf,		0.25,		0.075;
    RHO_e_m,		,			,		,		beta_pdf,		0.15,		0.075; 
	RHO_g,			,			,		,		beta_pdf,		0.75, 		0.075;
	RHO_p_co,		,			,		,		beta_pdf,		0.75, 		0.075;
	RHO_ppi_ast,	,			,		,		beta_pdf,		0.75, 		0.075;
	RHO_R_W,		,			,		,		beta_pdf,		0.75, 		0.075;
	RHO_sigma_ee,	,			,		,		beta_pdf,		0.75, 		0.075;
	RHO_sigma_ff,	,			,		,		beta_pdf,		0.75, 		0.075;
	RHO_sigma_hh,	,			,       ,		beta_pdf,		0.75, 		0.075;
	RHO_sigma_ii,	,			,		,		beta_pdf,		0.75, 		0.075;
	RHO_varrho,		,			,		,		beta_pdf,		0.75, 		0.075;
	RHO_xi_h,		,			,		,		beta_pdf,		0.75, 		0.075;
	RHO_xi_i,		,			,		,		beta_pdf,		0.75,		0.075;	
	RHO_xi_ih,      ,			,		,		beta_pdf,		0.75,		0.075;	
	RHO_xi_m,		,			,		,		beta_pdf,		0.75, 		0.075;
	RHO_xi_n,		,			,		,		beta_pdf,		0.75, 		0.075;
	RHO_xi_R,		,			,		,		beta_pdf,		0.75, 		0.075;
    RHO_xi_CHI_b,  	,     		,     	,  		beta_pdf, 		0.75, 		0.075;
    RHO_xi_CHI_e,  	,     		,     	,  		beta_pdf, 		0.75, 		0.075;
    RHO_xi_roe_r,  	,     		,     	,  		beta_pdf, 		0.75, 		0.075;
	RHO_y_ast,		,			,		,		beta_pdf,		0.85,		0.075;
	RHO_y_co,		,			,		,		beta_pdf,		0.85,		0.075;
	RHO_z,			,			,		,		beta_pdf,		0.85,		0.075;
	RHO_zetau,  	,     		,     	,  		beta_pdf, 		0.75, 		0.075;
	RHO_epsilon_L,  ,     		,     	,  		beta_pdf, 		0.75, 		0.05;

stderr u_a,			,			,		,		inv_gamma_pdf,	0.5,		inf;
stderr u_e_m,		,			,		,		inv_gamma_pdf,	0.5,		inf;
stderr u_g,			,			,		,		inv_gamma_pdf,	0.5,		inf;
stderr u_p_co,		,			,		,		inv_gamma_pdf,	0.5,		inf;
stderr u_ppi_ast,	,			,		,		inv_gamma_pdf,	0.5,		inf;
stderr u_R_W,		,			,		,		inv_gamma_pdf,	0.5,		inf;
stderr u_sigma_ee,	,			0,		100,		inv_gamma_pdf,	0.5,		inf;
stderr u_sigma_ff,	,			0,		100,		inv_gamma_pdf,	0.5,		inf;
stderr u_sigma_hh,	,			0,		100,		inv_gamma_pdf,	0.5,		inf;
stderr u_sigma_ii,	,			0,		100,		inv_gamma_pdf,	0.5,		inf;
stderr u_varrho,	,			,		,		inv_gamma_pdf,	0.5,		inf;
stderr u_xi_h,		,			,		,		inv_gamma_pdf,	0.5,		inf;
stderr u_xi_i,		,			,		,		inv_gamma_pdf,	0.5,		inf;
stderr u_xi_ih,		,			,		,		inv_gamma_pdf,	0.5,		inf;
stderr u_xi_m,		,			,		,		inv_gamma_pdf,	0.5,		inf;
stderr u_xi_n,		,			,		,		inv_gamma_pdf,	0.5,		inf;
stderr u_xi_R,		,			,		,		inv_gamma_pdf,	0.5,		inf;
stderr u_xi_CHI_b,	,     		,     	,  		inv_gamma_pdf, 	0.5, 		inf;
stderr u_xi_CHI_e,	,     		,     	,  		inv_gamma_pdf, 	0.5, 		inf;
stderr u_xi_roe_r,  ,     		,     	,  		inv_gamma_pdf, 	0.5, 		inf;
stderr u_y_ast,		,			,		,		inv_gamma_pdf,	0.5,		inf;
stderr u_y_co,		,			,		,		inv_gamma_pdf,	0.5,		inf;
stderr u_z,			,			,		,		inv_gamma_pdf,	0.5,		inf;
stderr u_zetau,		,     		,     	,  		inv_gamma_pdf, 	0.5, 		inf;
stderr u_epsilon_L,	,     		,     	10,  	inv_gamma_pdf, 	0.2, 		inf;

end;

% -------------------------------------------------------------------------
	// Estimation commands
estimation (datafile	=	Data_MAFIN_2023Q1_Pol, // Data_MAFIN_2023Q1_Pol,
			xls_sheet	=	DataDEMEAN,				
			plot_priors	=	1,
			prior_trunc	=	0,
			lyapunov	=	doubling,
			mode_compute=	0,
			mh_replic	=	0,
			order 		= 	1, 	
            //endogenous_prior,
			mode_file	=	mfm_ds01_v290423_2,
			//mode_check,
			kalman_algo	=	1,
            //fast_kalman_filter,
            //filter_covariance,             
			nograph,
            //mh_jscale = 0.125,
            //mh_nblocks = 1,
            //mh_drop = 0.5,
            filtered_vars ,
            smoother             
			);

%do_plot_priors_posteriors_laplace
%forecast(periods=11, conf_sig = 0.5);

//stoch_simul(order = 2, periods = 5000, pruning, noprint, nomoments, irf = 0);
//do_plot_simulation_stoch_steadystate

//stoch_simul(order=2, nograph, noprint, nomoments, periods=400, irf = 0);

//stoch_simul(order=2, nograph, noprint, nomoments,  irf = 0);
 