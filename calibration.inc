
% =========================================================================
% MAFIN Calibration
% =========================================================================

% -------------------------------------------------------------------------
	// CALIBRATED deep macro parameters
ALPHA		= 1-0.66;				//From XMAS: Labor share of 66% (Medina and Soto, 2007)
ALPHA_E		= 0.5;					//From XMAS: 50% expected inflation weight (XMAS assumption)
ALPHA_BSG 	= -0.4;     			//Short-term Bonds as Percentage of GDP
ALPHA_BLG 	= -4.5;     			//Long-term Bonds as Percentage of GDP
BETA_up		= 0.99997; // 0.99997;		//From XMAS: Annual real rate around 1.5%
BETA_i 		= 0.98;				// 
CHI			= 0.33;					//From XMAS: Codelco production / total commodity production, 2001-2015
DELTA_h		= 0.01;					//Assumption: same as capital depreciation
DELTA_k		= 0.01;					//From XMAS: Annual depreciation rate (Adolfson et al., 2008)
EPSILON_F	= 11;					//From XMAS: 10% ss price markup (Medina and Soto, 2007)
EPSILON_H	= 11;					//From XMAS: 10% ss price markup (Medina and Soto, 2007)
EPSILON_W	= 11;					//From XMAS: 10% ss price markup (Medina and Soto, 2007)
OMEGA		= 0.79;					//From XMAS: home bias in domestic demand, 2008-2016
N_H			= @{NH};				//Set in mod file. For now, =0 to shut down time-to-build
ppi_T		= 1.03^.25;				//From XMAS: Annual inflation target of 3%
RHO_VARPHI_H= 1+eps;				//=1 implies even expenditure across time to build of housing, <1 implies front-loading, and >1 back-loading. Warning: can't be exactly =1 due to division by 0
SIGMA		= 1;					//From XMAS: Log utility
V_PF		= 10^-5; 		 		//Strenght of wealth effect, =10^-5 means no wealth effect
ETA_CHAT	= 0.117528669045283; 					//CES Calibration
OMEGA_UP	= 0.7;					//Fraction of unrestricted patient households
UPSILON_H	= 1;					//Ratio of restricted patient housing to unrestricted patient housing
OMEGA_BL    = 0.822;    			//Ratio of long term bonds to short term bonds
//ETA_ZETA_L  = 0.144936296648825;				//Elasticity on transaction cost zeta_L by trading long-term governement bonds, from CCF
ETA_ZETA_L  = 0.0001;				//Elasticity on transaction cost zeta_L by trading long-term governement bonds, from CCF
EPSILON_TAU = 0.1;     	 			//Tax adjustment parameter
Duration_LOAN	= 10*4;				//Duration long-term housing debt, quarters
KAPPA_LOAN		= (1-1/Duration_LOAN);	//From Carlstrom et al., (2017)
Duration_BL		= 10*4;					//Duration long-term government bonds, quarters
KAPPA_B			= (1-1/Duration_BL);   	//From Carlstrom et al., (2017)
Duration_BB		= 5*4;					//Duration long-term bank bonds, quarters
KAPPA_BH		= (1-1/Duration_BB);    //From Carlstrom et al., (2017)
@#for j in 0:NH
	VARPHI_H_@{j}=RHO_VARPHI_H^(@{j})*(1-RHO_VARPHI_H)/(1-(RHO_VARPHI_H)^(N_H+1));	//Determined by choice of N_H and RHO_VARPHI_H
@#endfor
	// CALIBRATED financial parameters
CHI_b		= 0.04;					//From Clerc et al. (2015)
CHI_e		= 0.05;					//From Clerc et al. (2015)
GAMMA_bh	= 0.1;					//From Clerc et al. (2015)
GAMMA_d		= 0.1;					//From Clerc et al. (2015)
GAMMA_c     = 0.25;                 // (base=0.25) For CCyB process: map 10% of credit gap into 2.5% of ccyb
MU_e		= 0.118;					//From Mendoza
MU_f		= 0.118;					//
MU_h		= 0.118;					//
MU_i		= 0.118;					//
//phi_f		= 0.123;				//Chilean capital requirements of 8%, with 100% risk weight for corporate loans, +4.3% of excess capital in data
//phi_h		= 0.091;				//Chilean capital requirements of 8%, with 60% risk weight for housing loans, +4.3% of excess capital in data

ALPHA_bbval = 0.75;
ALPHA2_bbval = 0.5;

phi_h_base 	= 0.0975;					//112Chilean capital requirements of 8%/10.5%, with 100% risk weight for corporate loans, +4.3% of excess capital in data
phi_f_base 	= 0.0975;					//112Chilean capital requirements of 8%, with 60% risk weight for corporate loans, +4.3% of excess capital in data
ccybneutral = 0.005;

% -------------------------------------------------------------------------
	// Estimated deep macro parameters
ALPHA_ppi   = 1.914122213705633;
ALPHA_R     = 0.717899565126331;          
ALPHA_W     = 0.182626571836242;          
ALPHA_y     = 0.296923233254939;     
ETA         = 0.743630458321263;
ETA_ast     = 0.122187522459919;           
//GAMA_L      = 0.052267893189102;   
GAMA_L      = 0.88;        
GAMA_h      = 2.759096995495109;          
GAMA_k      = 2.982440179092768;           
GAMA_n      = 2.163525599003024;          
KAPPA_F     = 0.609483586048955;          
KAPPA_H     = 0.743582637035332;          
KAPPA_W     = 0.853615939887190;          
PHI_ast     = 0.318507419700392;            
PHI_c       = 0.813193654791216;           
PHI_HH      = 0.805880484097945;          
THETA_F     = 0.668798364544585;          
THETA_H     = 0.766206348812832;          
THETA_W     = 0.498980913334460;           
VARPHI      = 5.676599397226525;







% -------------------------------------------------------------------------
	// Steady state of macro shocks
a_ss		= (1+(0.033-0.0181))^.25;	//From XMAS: Annual real growth around 1.5% (per labor force)
bl_cb_ss    = 0;										//Normalization
//blg_ss    = 1;         								//Normalization
e_m_ss		= 1;										//Normalization
e_req_ss	= 1;										//Normalization
p_co_ss		= 1;										//Normalization
varrho_ss	= 1;										//Normalization
xi_h_ss		= 0.9; 										//Normalization
xi_i_ss		= 1;										//Normalization
xi_ih_ss	= 1;										//Normalization
xi_m_ss		= 1;										//Normalization
xi_CHI_b_ss	= 1;										//Normalization  
xi_CHI_e_ss = 1;										//Normalization	 
xi_roe_r_ss = 1;										//Normalization  
xi_R_ss		= 1.015^.25;								//From XMAS: EMBI Chile, 2001-2016 of 1.5% annual
z_ss		= 1;										//Normalization
zetau_ss  	= 1;										//Normalization
	// Calibrated SS macro moments
n_ss		= 0.3;					//Calibrated to XMAS: Hours worked (out of 1 unit of time)
r_h_k		= 0.65;					//Ratio of housing to capital stock, data average
s_co		= 0.12;					//Calibrated to XMAS: Share of commodity production in GDP
s_g			= 0.12;					//Calibrated to XMAS: Share of government expenditure in GDP
s_tb		= 0.05;					//Calibrated to XMAS: Share of trade balance in GDP
s_bast		= -0.14;				//Calibrated to XMAS: Net foreign asset position as share of GDPtive means debtor)
	// Normalizations of endogenous variables
ppi_S_ss		= 1;
p_H_ss			= 1;
	// Calibrated SS financial spreads
@#if fin_fric
  tpm_ss	= ppi_T*a_ss^SIGMA/BETA_up;
	R_D_ss	= (1.0028^.25-1) + tpm_ss;	//Quarterly spread over TPM from 2002-2019 average + R_ss
	R_L_ss	= (1.0448^.25-1) + tpm_ss;	//Quarterly spread over TPM from 2002-2019 average + R_ss
	Rnom_i_ss	= (1.0381^.25-1) + tpm_ss;  //Quarterly spread over TPM from 2002-2019 average + R_ss
	Rnom_BL_ss = (1.0175^.25-1) + tpm_ss;	//Long term rate, Quarterly spread over TPM from 2002-2019 average + R_ss
@#endif





% -------------------------------------------------------------------------
	// Autocorrelation of macro shocks
RHO_a		= 0.363656819945829;	//
RHO_bl_g    = 0.80;              	//Persistence long-term government bonds
RHO_bl_cb   = 0.80;					//Persistence of CB long-term government bonds purchases
RHO_e_m		= 0.108867133971157;	//Estimated without financial frictions
RHO_e_req	= 0.5; // 0.999;  0.65;	
//RHO_epsilon_L=0.823760348761645;      //Persistence shock on transaction costs
RHO_epsilon_L=0;      //Persistence shock on transaction costs
RHO_g		= 0.770377629939295;	//
RHO_p_co	= 0.880367938558405;	
RHO_ppi_ast	= 0.408731823582096;	//Estimated without financial frictions
RHO_R_W		= 0.885096185173148;	//Estimated without financial frictions
RHO_varrho	= 0.438322896543878;	//Estimated without financial frictions
RHO_xi_h	= 0.777020606233316;	//Estimated without financial frictions
RHO_xi_i	= 0.628527245439942;	//Estimated without financial frictions
RHO_xi_ih	= 0.607750982011759;	//Estimated without financial frictions
RHO_xi_m	= 0.768741282879974;	//Estimated without financial frictions
RHO_xi_n	= 0.768475072532623;	//Estimated without financial frictions
RHO_xi_R	= 0.795039102806919;	//Estimated without financial frictions
RHO_y_ast	= 0.902262857779588;	//Estimated without financial frictions
RHO_y_co	= 0.796226377303659;	//Estimated without financial frictions
RHO_z		= 0.893364545092234;	//Estimated without financial frictions
RHO_zetau   = 0.777020810710361;    				//XMAS



	// Autocorrelation of financial shocks
RHO_sigma_ee	= 0.957785569566603;				//Prior mean
RHO_sigma_ff	= 0.644431674401923;				//Prior mean
RHO_sigma_hh	= 0.769341784870793;				//Prior mean
RHO_sigma_ii	= 0.965302626572274;				//Prior mean
RHO_xi_CHI_b	= 0.819404180560220;				//Prior mean
RHO_xi_CHI_e	= 0.982210502564640;				//Prior mean
RHO_xi_roe_r	= 0.809663202046836;				//Prior mean

	// Variance of macro shocks
SIGMA_a			= 0.318732078641522;	//Estimated without financial frictions
SIGMA_bl_g		= 0;					//TURNED OFF FOR NOW
SIGMA_bl_cb		= 0;					//TURNED OFF FOR NOW
SIGMA_e_m		= 0.158440038352574;	//Estimated without financial frictions
SIGMA_e_req		= 0.01;					//to be estimated... (shock to credit gap)
//SIGMA_epsilon_L = 9.528607280128941;					//to be estimated...
SIGMA_epsilon_L = 0;					//to be estimated...
SIGMA_g			= 1.286765272375807;	//Estimated without financial frictions
SIGMA_p_co		= 1.299815465053385;		//Estimated without financial frictions
SIGMA_ppi_ast	= 2.302254121593250;		//Estimated without financial frictions
SIGMA_R_W		= 0.126945868676583;	//Estimated without financial frictions
SIGMA_varrho	= 7.747012306441073;		//Estimated without financial frictions
SIGMA_xi_h		= 0.246729668036235;	//Estimated without financial frictions
SIGMA_xi_i		= 7.689562975505806;		//Estimated without financial frictions
SIGMA_xi_ih		= 1.410157317152713;		//Estimated without financial frictions
SIGMA_xi_m		= 2.209355899685924;		//Estimated without financial frictions
SIGMA_xi_n		= 1.716805337778568;		//Estimated without financial frictions
SIGMA_xi_R		= 0.086680559975012;	//Estimated without financial frictions
SIGMA_y_ast		= 0.240231899798079;	//Estimated without financial frictions
SIGMA_y_co		= 2.717683281920666;		//Estimated without financial frictions
SIGMA_z			= 1.102240621247445;	//Estimated without financial frictions
SIGMA_zetau     = 0.517539249842122;     			//XMAS
	// Variance of financial shocks
SIGMA_sigma_ee	= 0.309551973922142;					//Shut down for now
SIGMA_sigma_ff	= 0.924929956837554;					//Shut down for now
SIGMA_sigma_hh	= 0.246053494773249;					//Shut down for now
SIGMA_sigma_ii	= 3.351607706507753;					//Shut down for now
SIGMA_xi_CHI_b	= 0.566435213645889;					//Shut down for now
SIGMA_xi_CHI_e	= 0.749539620619115;					//Shut down for now
SIGMA_xi_roe_r	= 0.749539620619115;					//Shut down for now








% -------------------------------------------------------------------------
	// Steady state of observables
@#for va in vars_obs
	@{va}_obs_ss = 0;
@#endfor




% -------------------------------------------------------------------------
	// No financial frictions
@#if !fin_fric
	BAR_l_h 	= 0.0421893; //0.0001;	 			//Subject to sensitivity.
	PHI_i 		= 0.0001;				//The lower, the less financial frictions. Subject to sensitivity.
	sigma_ee_ss	= 1;					//Irrelevant
	sigma_ff_ss	= 1;					//Irrelevant
	sigma_hh_ss	= 1;					//Irrelevant
	sigma_ii_ss	= 1;					//Irrelevant
	SIGMA_xi_CHI_b	= 1;				//Irrelevant
    SIGMA_xi_CHI_e	= 1;				//Irrelevant
	SIGMA_xi_roe_r	= 1;				//Irrelevant
	OMEGA_UP	= 1;					//Fraction of Unrestricted patient households
@#else
	BAR_l_h 	= 0.0421893; //0.0001;	 			//Subject to sensitivity.
	PHI_i 		= 0.0001;				//The lower, the less financial frictions. Subject to sensitivity.

	//IOTA		= a_ss^SIGMA*ppi_T*(1/BETA_i-1/BETA_up);
	IOTA		= a_ss^SIGMA*(1/BETA_i-ppi_T/BETA_up);
@#endif


