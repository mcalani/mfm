
% =========================================================================
% MAFIN MODEL
% =========================================================================

model;

% =========================================================================
%	Definicion reglas capital contracíclico 
% =========================================================================

@#if fin_fric
    @#if ccyb_rule==0       //No CCyB 
        buffer_guide = 0;
        ccyb = .76*(ccyb(-1))+ (e_req-1);  //9175=8, 958=16, 
        //ccyb = 0 ;  
    @#elseif ccyb_rule==1   //Credit gap 
        buffer_guide = (log(l_f + q_hat_l*l_h) - log(l_f_ss + q_l_ss*l_h_ss)) ;                   
        ccyb= GAMMA_c*(buffer_guide)*(e_req);

    @#elseif ccyb_rule==2   //Credit to GDP gap 
        buffer_guide = (log((l_f + q_l*l_h)/gdpn) - log((l_f_ss + q_l_ss*l_h_ss)/gdpn)) ;   
        ccyb= (1-ALPHA_R)*ccyb(-1)  + ALPHA_R*GAMMA_c*(buffer_guide) + GAMMA_c*log(e_req) *0.277665/0.0783217;

     @#elseif ccyb_rule==3   //Default 1per 
        //buffer_guide = PD_d-PD_d_ss + log(e_req);        
       buffer_guide = 3*log((1+PD_d(+1))/(1+PD_d_ss)) ;                                          
       ccyb= (1-ALPHA_R)*ccyb(-1)  + ALPHA_R*GAMMA_c*(buffer_guide) + GAMMA_c*log(e_req) *0.277665/0.173855;
 
     @#elseif ccyb_rule==4   //Shock 
        buffer_guide = 0 ;
		//ccyb = (1-1)*ccyb(-1) + 1*(GAMMA_c*(buffer_guide)) + (e_req(-2)-1) - (e_req(-10)-1) ; 
		//ccyb = (1-ALPHA_R)*ccyb(-1) + ALPHA_R*(GAMMA_c*(buffer_guide) + (e_req(-2)-1) - (e_req(-6)-1) ) ; 
		ccyb =  (e_req(-2)-1) - (e_req(-6)-1) ; 
        
     @#elseif ccyb_rule==5   //Default forward looking                                                                       
        buffer_guide = 10*((PD_d(+1)*(1-ALPHA_E)+PD_d(+4)*ALPHA_E) - PD_d_ss) ;
        ccyb= (1-ALPHA_R)*ccyb(-1)  + ALPHA_R*(buffer_guide) + GAMMA_c*log(e_req)*.277665/.0441636;
     @#endif

@#endif


@#if !fin_fric
buffer_guide=0;
ccyb = 0;
@#endif

% =========================================================================
%	Equilibrium Equations
% =========================================================================

	// zetal de capital total 

phi_f = 0.0433 + (phi_f_base+ccyb+ccybneutral) ;
phi_h = 0.0433 + 0.6*(phi_h_base+ccyb+ccybneutral);


% -------------------------------------------------------------------------
	// Unrestricted Patient households

c_hat_u 		= ((1-O_CHAT)^(1/ETA_CHAT)*(c_u-PHI_c*c_u(-1)/a)^((ETA_CHAT-1)/ETA_CHAT)
		+(O_CHAT)^(1/ETA_CHAT)*(xi_h*(h_u(-1)/a-PHI_HH*h_u(-2)/(a*a(-1))))^
		((ETA_CHAT-1)/ETA_CHAT))^(ETA_CHAT/(ETA_CHAT-1));												//(1)CES
la_u 			= c_hat_u^-SIGMA*((1-O_CHAT)*c_hat_u/(c_u-PHI_c*c_u(-1)/a))^(1/ETA_CHAT);				//(2)FOC consumption
varrho*la_u*q_h	= BETA_up*varrho(+1)*((c_hat_u(+1)*a(+1))^-SIGMA*xi_h(+1)*
		(O_CHAT*c_hat_u(+1)*a(+1)/(xi_h(+1)*(h_u-PHI_HH*h_u(-1)/a)))^(1/ETA_CHAT)
		+(1-DELTA_h)*la_u(+1)*a(+1)^-SIGMA*q_h(+1));													//(3)FOC housing
varrho*la_u		= BETA_up*R*varrho(+1)*la_u(+1)/ppi(+1)*a(+1)^-SIGMA;									//*****FOC government bonds local short term
varrho*la_u		= BETA_up*R_D_tilde(+1)/ppi(+1)*varrho(+1)*la_u(+1)*a(+1)^-SIGMA;						//*****FOC deposits
varrho*la_u 	= BETA_up*R_ast*varrho(+1)*la_u(+1)*ppi_S(+1)/ppi(+1)*a(+1)^-SIGMA;						//*****FOC bonds foreign
varrho*la_u*(1+zeta_L)*q_BL = BETA_up*varrho(+1)*la_u(+1)*a(+1)^-SIGMA*R_BL(+1)*q_BL(+1);				//*****FOC long-term government bonds
varrho*la_u*(1+zeta_L)*q_BB = BETA_up*varrho(+1)*la_u(+1)*a(+1)^-SIGMA*R_BB_tilde(+1)*q_BB(+1);			//*****FOC long-term banking bonds

% -------------------------------------------------------------------------
	// Restricted Patient households
@#if fin_fric
c_hat_r			= ((1-O_CHAT)^(1/ETA_CHAT)*(c_r-PHI_c*c_r(-1)/a)^((ETA_CHAT-1)/ETA_CHAT)
		+(O_CHAT)^(1/ETA_CHAT)*(xi_h*(h_r(-1)/a-PHI_HH*h_r(-2)/(a*a(-1))))^((
		ETA_CHAT-1)/ETA_CHAT))^(ETA_CHAT/(ETA_CHAT-1));													//*****CES
la_r			= c_hat_r^(-SIGMA)*((1-O_CHAT)*c_hat_r/(c_r-PHI_c*c_r(-1)/a))^(1/ETA_CHAT); 			//*****FOC consumption
varrho*la_r*q_h = BETA_rp*varrho(+1)*((c_hat_r(+1)*a(+1))^(-SIGMA)*xi_h(+1)*
		(O_CHAT*c_hat_r(+1)*a(+1)/(xi_h(+1)*(h_r-PHI_HH*h_r(-1)/a)))^(1/ETA_CHAT)
		+(1-DELTA_h)*la_r(+1)*a(+1)^(-SIGMA)*q_h(+1));													//*****FOC housing
varrho*la_r*q_BL = BETA_rp*varrho(+1)*la_r(+1)*R_BL(+1)*q_BL(+1)*a(+1)^-SIGMA;							//*****FOC long-term government bonds
q_BL*bl_r+c_r+q_h*h_r = q_BL*R_BL*bl_r(-1)/(a)+w*n_r+
		q_h*(1-DELTA_h)*h_r(-1)/(a);																	//*****Budget constraint restricted patient
@#endif

% -------------------------------------------------------------------------
	// Impatient households
R_h/ppi			= q_h*(1-DELTA_h)/q_h(-1);				 												//*****Return of "housing project"
c_hat_i			= ((1-O_CHAT)^(1/ETA_CHAT)*(c_i-PHI_c*c_i(-1)/a)^((ETA_CHAT-1)/ETA_CHAT)
					+(O_CHAT)^(1/ETA_CHAT)*(xi_h*(h_i(-1)/a-PHI_HH*h_i(-2)/(a*a(-1))))
					^((ETA_CHAT-1)/ETA_CHAT))^(ETA_CHAT/(ETA_CHAT-1));									//*****CES
la_i 			= c_hat_i^-SIGMA*((1-O_CHAT)*c_hat_i/(c_i-PHI_c*c_i(-1)/a))^(1/ETA_CHAT);				//*****FOC consumption
Rnom_i 	 		= R_i*ppi;

q_hat_l         =(1-KAPPA_LOAN*l_h(-1)/(l_h*a))*q_l+(KAPPA_LOAN*l_h(-1)/(l_h*a))*q_hat_l(-1);           //******New price of mortage --------

@#if fin_fric
	om_bar_i 	= R_hat_i*q_hat_l*l_h(-1)*ppi/(R_h*q_h(-1)*h_i(-1));                  					//*****Idiosyncratic shock threshold -------
	R_i 	 	= 1/q_l + KAPPA_LOAN;									   								//*****Mortgage Yield to Maturity ---------
    R_hat_i     = (1+KAPPA_LOAN*q_hat_l(-1))/q_hat_l;                                                   //***** Mortage efective-----------
	varrho*la_i*q_h = BETA_i*varrho(+1)* ( (c_hat_i(+1)*a(+1))^(-SIGMA)*xi_h(+1)
		*(O_CHAT*c_hat_i(+1)*a(+1)/(xi_h(+1)*(h_i-PHI_HH*h_i(-1)/a)))^(1/ETA_CHAT)
		+la_i(+1)*a(+1)^(-SIGMA)*Gam_i(+1)*R_h(+1)/ppi(+1)*q_h)
		+varrho*la_h*Gam_h(+1)*Gam_G_i(+1)*R_h(+1)*q_h;			                        				//*****FOC housing -----------
	BETA_i = varrho*la_h*ppi(+1)/(varrho(+1)*la_i(+1)*a(+1)^(-SIGMA))
		*Gam_h(+1)*Gam_G_der_i(+1)/Gam_der_i(+1);														//FOC leverage
	
		c_i+q_h*h_i-q_l*(l_h-KAPPA_LOAN*l_h(-1)/a)*(1-GAMA_L/2*(Dl-a_ss)^2)
	  -KAPPA_LOAN*l_h(-1)*q_hat_l(-1)/a= w*n/2 + Gam_i*R_h*q_h(-1)*h_i(-1)/(a*ppi);	                    //*****Budget constraint new 

		// sin costo ajuste 
	//varrho*q_l*(la_i-la_h*rho_h(+1)*phi_h)=BETA_i*(varrho(+1)*la_i(+1)									//FOC L_h 
	//*a(+1)^(-SIGMA)*KAPPA_LOAN*(q_l(+1)-q_hat_l)); 	
		
		// con costo ajuste 
	varrho*la_i*q_l*(1-GAMA_L/2*(Dl-a_ss)^2+GAMA_L*Dl*(Dl-a_ss))-varrho*la_h*rho_h(+1)*phi_h*q_l = 
	BETA_i*varrho(+1)*la_i(+1)*a(+1)^(-SIGMA)*(KAPPA_LOAN*q_l(+1)*(1-GAMA_L/2*(Dl(+1)-a_ss)^2)-
	q_l(+1)*Dl(+1)*GAMA_L*(Dl(+1)-a_ss)*(KAPPA_LOAN+Dl(+1))-KAPPA_LOAN*q_hat_l);						//FOC L_h 
	
	Dl = (l_h-KAPPA_LOAN*l_h(-1)/a)*a/(l_h(-1)-KAPPA_LOAN*l_h(-2)/a(-1));								//Def Dl 

	PD_i		= normcdf(aux_i);																		//*****Probability of default on mortgage loans
@#endif
Theta_i			= Xi_Tilde_i*(c_hat_i)^-SIGMA; 															//*****Preference Shifter, Patient
Xi_Tilde_i 		= Xi_Tilde_i(-1)^(1-V_PF)*(c_hat_i)^(SIGMA*V_PF); 										//*****Preference Shifter, Patient

% -------------------------------------------------------------------------
	// Entrepreneurs
q_k*k	= n_e+l_f;																						//(30)Balance sheet
R_e/ppi			= (r_k+(1-DELTA_k)*q_k)/q_k(-1);														//*****Return of capital project
@#if fin_fric
om_bar_e		= R_L(-1)*l_f(-1)/(R_e*q_k(-1)*k(-1));													//*****Idiosyncratic shock threshold
	c_e			= CHI_e*xi_CHI_e*psi_e;																	//*****Dividends
	n_e			= (1-CHI_e*xi_CHI_e)*psi_e;																//*****Retained profits (bequests or net worth)
	psi_e*a*ppi	=	Gam_e*R_e*q_k(-1)*k(-1); 															//*****Final wealth
	Gam_e(+1)	= la_e*(rho_f(+1)*phi_f/R_e(+1)- Gam_f(+1)*Gam_G_e(+1));								//*****FOC capital
	Gam_der_e(+1)	= la_e*Gam_f(+1)*Gam_G_der_e(+1);													//*****FOC leverage
	PD_e 		= normcdf(aux_e);																		//*****Probability of default on corporate loans
@#endif

% -------------------------------------------------------------------------
	// F Banks

d_f + e_f 		= l_f;																					//(40)Balance sheet
@#if fin_fric
om_bar_f 		= (1-phi_f)*R_D(-1)/R_f_tilde;															//*****Idiosyncratic shock threshold
	e_f			= phi_f*l_f;																			//*****Capital requirements
	rho_f 		= Gam_f*R_f_tilde/phi_f;																//*****Return on equity
	R_f_tilde	= Gam_G_e*R_e*q_k(-1)*k(-1)/l_f(-1);													//*****Return on assets
	PD_f 		= normcdf(aux_f);																		//*****Default probability of bank
@#endif

% -------------------------------------------------------------------------
	// H Banks
l_h*q_l			= q_BB*bb_tot + e_h;																	//*****Balance sheet------ Revisar
//l_h*q_l + (KAPPA_LOAN*l_h(-1)/a)*(q_hat_l(-1)-q_l)		= q_BB*bb_tot + e_h;																	//*****Balance sheet------ Revisar
@#if fin_fric
	om_bar_h	= (1-phi_h)*R_BB*q_BB/(R_h_tilde*q_BB(-1))*ppi;											//*****Idiosyncratic shock threshold
	e_h	 		= phi_h*q_l*l_h;
    //e_h	 		= phi_h*(l_h*q_l + (KAPPA_LOAN*l_h(-1)/a)*(q_hat_l(-1)-q_l)	);																		//*****Capital requirements
	rho_h		= Gam_h*R_h_tilde/phi_h;																//*****Return on equity
	R_h_tilde	= Gam_G_i*R_h*q_h(-1)*h_i(-1)/(q_l(-1)*l_h(-1));										//*****Return on assets
	PD_h		= normcdf(aux_h);																		//(50)Default probability of bank
@#endif

% -------------------------------------------------------------------------
	// Bankers and banking system
@#if fin_fric
	rho_tilde_h = (1-KAPPA_LOAN)*rho_h + KAPPA_LOAN*rho_tilde_h(+1);		   							//*****Expected rate
	rho_f(+1)	= xi_roe_r*rho_tilde_h(+1);																//*****Equal expected ROE
	c_b			= CHI_b*xi_CHI_b*psi_b;																	//*****Dividends
	n_b			= (1-CHI_b*xi_CHI_b)*psi_b;																//*****Retained profits (bequests or net worh)
	psi_b*a*ppi	= rho_f*e_f(-1)+rho_tilde_h*e_h(-1);													//*****Final wealth
	n_b			= e_f+e_h;																				//*****Equity distribution	
@#endif
PD_d 		= (q_BB(-1)*bb_tot(-1)*PD_h+d_f(-1)*PD_f)/(q_BB(-1)*bb_tot(-1)+d_f(-1));					//*****Probability of default on deposits

% -------------------------------------------------------------------------
	// Capital and housing producers
k	= (1-DELTA_k)*k(-1)/a+(1-GAMA_k/2*(i/i(-1)*a-a_ss)^2)*xi_i*i;										//*****Capital law of motion
1	= q_k*(1-GAMA_k/2*(i/i(-1)*a-a_ss)^2-GAMA_k*(i/i(-1)*a-a_ss)*i/i(-1)*a)*xi_i
		+BETA_up*varrho(+1)*la_p(+1)/(varrho*la_p)*a(+1)^(-SIGMA)*q_k(+1)
		*GAMA_k*(i(+1)/i*a(+1)-a_ss)*(i(+1)/i*a(+1))^2*xi_i(+1);										//*****FOC capital producer
h	= (1-DELTA_h)*h(-1)/a + (1-GAMA_h/2*(i_ah(-@{NH})/i_ah(-@{NH}-1)*a-a_ss)^2)
		*xi_ih(-@{NH})*i_ah(-@{NH})
		@#for i in 0:NH-1
		/a(@{-i})
		@#endfor
		;																								//(60)Housing law of motion
0	=
		@#for j in 0:NH
		+BETA_up^@{j}*varrho(@{j})*la_p(@{j})*VARPHI_H_@{j}
			@#for i in j+1:NH
			*a(@{i})^SIGMA
			@#endfor
		@#endfor
		-BETA_up^N_H*varrho(@{NH})*la_p(@{NH})*q_h(@{NH})
		*((1-GAMA_h/2*(i_ah/i_ah(-1)*a-a_ss)^2)-GAMA_h*(i_ah/i_ah(-1)*a-a_ss)*
		i_ah/i_ah(-1)*a)*xi_ih - BETA_up^(N_H+1)*varrho(@{NH}+1)*la_p(@{NH}+1)*
		q_h(@{NH}+1)*a(@{NH}+1)^(-SIGMA)
		*GAMA_h*(i_ah(+1)/i_ah*a(+1)-a_ss)*(i_ah(+1)/i_ah*a(+1))^2*xi_ih(+1);							//*****FOC housing producer
i_h	=
		@#for j in 0:NH
		+VARPHI_H_@{j}*i_ah(@{-j})
			@#for i in 0:j-1
			/a(@{-i})
			@#endfor
		@#endfor
		;																								//*****Housing investment aggregation

% -------------------------------------------------------------------------
	// Final goods
y_C			= (OMEGA^(1/ETA)*x_H^(1-1/ETA)+(1-OMEGA)^(1/ETA)*
					x_F^(1-1/ETA))^(ETA/(ETA-1));														//*****Final consumption good production function
x_F			= (1-OMEGA)*p_F^(-ETA)*y_C;																	//*****FOC foreign good
x_H			= OMEGA*p_H^(-ETA)*y_C;																		//*****FOC home good

% -------------------------------------------------------------------------
	//  Home goods
f_H			= (EPSILON_H-1)/EPSILON_H*p_tilde_H^(1-EPSILON_H)*y_H
					+BETA_up*THETA_H*varrho(+1)*la_p(+1)/(varrho*la_p)
					*a(+1)^(1-SIGMA)/ppi(+1)*(p_tilde_H*ppi_IH(+1)/
					p_tilde_H(+1))^(1-EPSILON_H)*ppi_H(+1)^EPSILON_H*f_H(+1);							//*****FOC (1) Calvo by price setter of variety j
f_H			= p_tilde_H^(-EPSILON_H)*mc_H*y_H+BETA_up*THETA_H*varrho(+1)
					*la_p(+1)/(varrho*la_p)*a(+1)^(1-SIGMA)/ppi(+1)*(p_tilde_H*ppi_IH(+1)
					/p_tilde_H(+1))^(-EPSILON_H)*ppi_H(+1)^(1+EPSILON_H)*f_H(+1);						//*****FOC (2) Calvo by price setter of variety j
1			= (1-THETA_H)*p_tilde_H^(1-EPSILON_H)+THETA_H*(ppi_IH/ppi_H)^(1-EPSILON_H);					//*****Price aggregation
ppi_IH		= ppi_H(-1)^KAPPA_H*ppi_T^(1-KAPPA_H);														//*****Price indexation
p_H*mc_H	= p_Z;																						//(70)Optimal marginal cost

% -------------------------------------------------------------------------
	// Wholesale Domestic goods
@#if QLC == 1
	mc_Z 	= 	1/(ALPHA^ALPHA*(1-ALPHA)^(1-ALPHA))*(r_k^ALPHA/z)*
				(w + GAMA_n*(n_tilde/n_tilde(-1) - 1)*(1/n_tilde(-1))*(y_Z*p_Z)-
				BETA_up*varrho(+1)*la_p(+1)/(varrho*la_p)*a(+1)^(1-SIGMA)*GAMA_n*
				(n_tilde(+1)/n_tilde - 1)*(n_tilde(+1)/n_tilde^2)
				*y_Z(+1)*p_Z(+1))^(1-ALPHA);															//*****Optimal marginal cost, with Labour Quadratic Adj Costs
	k(-1)	= n_tilde*a*ALPHA/((1-ALPHA)*r_k)*(w + GAMA_n*(n_tilde/n_tilde(-1) - 1)
				*(1/n_tilde(-1))*(y_Z*p_Z)-BETA_up*varrho(+1)*la_p(+1)/(varrho*la_p)
				*a(+1)^(1-SIGMA)*GAMA_n*(n_tilde(+1)/n_tilde - 1)
				*(n_tilde(+1)/n_tilde^2)*y_Z(+1)*p_Z(+1));												//*****Optimal capital/labor ratio, with Labour Quadratic Adj Costs
@#else
	mc_Z	= 1/(ALPHA^ALPHA*(1-ALPHA)^(1-ALPHA))*r_k^ALPHA*w^(1-ALPHA)/z;								//+++++Optimal marginal cost
	k(-1)	= n_tilde*ALPHA/(1-ALPHA)*w/r_k*a;															//+++++Optimal capital/labor ratio
@#endif
p_Z 		= mc_Z;																						//*****Optimal price

% -------------------------------------------------------------------------
	// Foreign goods
p_F*mc_F	= rer*xi_m;																					//(69)Marginal cost
f_F			= (EPSILON_F-1)/EPSILON_F*p_tilde_F^(1-EPSILON_F)*y_F
				+BETA_up*THETA_F*varrho(+1)*la_p(+1)/(varrho*la_p)*a(+1)^(1-SIGMA)
				/ppi(+1)*(p_tilde_F*ppi_IF(+1)/p_tilde_F(+1))^(1-EPSILON_F)
				*ppi_F(+1)^(EPSILON_F)*f_F(+1);															//*****FOC (1) Calvo by price setter of variety j
f_F			= p_tilde_F^(-EPSILON_F)*y_F*mc_F + BETA_up*THETA_F*varrho(+1)
				*la_p(+1)/(varrho*la_p)*a(+1)^(1-SIGMA)/ppi(+1)*(p_tilde_F*ppi_IF(+1)
				/p_tilde_F(+1))^-EPSILON_F*ppi_F(+1)^(1+EPSILON_F)*f_F(+1);								//*****FOC (2) Calvo by price setter of variety j
1			= (1-THETA_F)*p_tilde_F^(1-EPSILON_F)+THETA_F*(ppi_IF/ppi_F)^(1-EPSILON_F);					//*****Price aggregation
ppi_IF		= ppi_F(-1)^KAPPA_F*ppi_T^(1-KAPPA_F);														//*****Price indexation

% -------------------------------------------------------------------------
	// Wages
la_W		= (la_p+la_i)/2;																			//*****Union lagrangian multiplier
la_p		= OMEGA_UP*la_u + (1-OMEGA_UP)*la_r;														//*****Patient households' lagrangian multiplier
Theta		= ((OMEGA_UP*Theta_up+(1-OMEGA_UP)*Theta_rp) + Theta_i)/2; 									//*****Preference Shifter
mc_W		= Theta*xi_n*n_tilde^VARPHI/(la_W*w); 														//****Marginal cost, with Preference Shifter
f_W			= (EPSILON_W-1)/EPSILON_W*w_tilde^(1-EPSILON_W)*n_tilde
				+((OMEGA_UP*BETA_up+(1-OMEGA_UP)*BETA_rp)+BETA_i)/2*THETA_W
				*a(+1)^(-SIGMA)*varrho(+1)*la_W(+1)/(varrho*la_W)
				*ppi_W(+1)/ppi(+1)*(ppi_Wtilde(+1)/ppi_IW(+1))^(EPSILON_W-1)*f_W(+1);					//*****FOC (1) Calvo by wage setter of variety i
f_W			= w_tilde^(-EPSILON_W*(1+VARPHI))*mc_W*n_tilde
				+((OMEGA_UP*BETA_up+(1-OMEGA_UP)*BETA_rp)+BETA_i)/2*THETA_W
				*a(+1)^(-SIGMA)*varrho(+1)*la_W(+1)/(varrho*la_W)*ppi_W(+1)/ppi(+1)
				*(ppi_Wtilde(+1)/ppi_IW(+1))^(EPSILON_W*(1+VARPHI))*f_W(+1);							//*****FOC (2) Calvo by wage setter of variety i
1			= (1-THETA_W)*w_tilde^(1-EPSILON_W)+THETA_W*(ppi_IW/ppi_W)^(1-EPSILON_W);					//*****Wage aggregation
ppi_IW 		= a(-1)^ALPHA_W*a_ss^(1-ALPHA_W)*ppi(-1)^KAPPA_W*ppi_T^(1-KAPPA_W);							//*****Wage indexation
Theta_up 	= Xi_Tilde_up*(c_hat_u)^-SIGMA; 															//Preference Shifter, Unrestricted Patient
Xi_Tilde_up	= Xi_Tilde_up(-1)^(1-V_PF)*(c_hat_u)^(SIGMA*V_PF); 											//*****Preference Shifter, Unrestricted Patient

@#if fin_fric
	Theta_rp	= Xi_Tilde_rp*(c_hat_r)^-SIGMA; 														//*****Preference Shifter, Restricted Patient
	Xi_Tilde_rp	= Xi_Tilde_rp(-1)^(1-V_PF)*(c_hat_r)^(SIGMA*V_PF); 										//*****Preference Shifter, Restricted Patient
@#endif
	
% -------------------------------------------------------------------------
	// Monetary policy and rest of the world
@#if taylor_rule == 0
	R/R_ss		= (R(-1)/R_ss)^ALPHA_R*(((ppi*(1-ALPHA_E)+ppi(+4)*ALPHA_E)/ppi_T)^ALPHA_ppi
				*(gdp/gdp(-1))^ALPHA_y)^(1-ALPHA_R)*e_m;
@#elseif taylor_rule == 1
	R/R_ss		= (R(-1)/R_ss)^ALPHA_R*(((ppi*(1-ALPHA_E)+ppi(+4)*ALPHA_E)/ppi_T)^ALPHA_ppi //***** Dovish (incremented ALPHA_y)
				*(gdp/gdp(-1))^(3*ALPHA_y))^(1-ALPHA_R)*e_m;

@#elseif taylor_rule == 2
	R/R_ss		= (R(-1)/R_ss)^ALPHA_R*(((ppi*(1-ALPHA_E)+ppi(+4)*ALPHA_E)/ppi_T)^(3*ALPHA_ppi)  //***** Hawkish (incremented ALPHA_ppi)
				*(gdp/gdp(-1))^ALPHA_y)^(1-ALPHA_R)*e_m;

@#elseif taylor_rule == 3
	R/R_ss = (ppi/ppi_T)^ALPHA_ppi*e_m ;							     //***** Simple Rule

@#elseif taylor_rule == 4
	R/R_ss		= (R(-1)/R_ss)^ALPHA_R*(((ppi*(1-ALPHA_E)+ppi(+4)*ALPHA_E)/ppi_T)^ALPHA_ppi //***** Taylor rule (base)
				*(gdp/gdp(-1))^ALPHA_y)^(1-ALPHA_R)*e_m;
@#elseif taylor_rule == 5
	R/R_ss = (ppi/ppi_T)^ALPHA_ppi * ( (l_f + q_l*l_h)/(l_f_ss + q_l_ss*l_h_ss) )^(-0.36) ;  //***** Credit spread augmented
@#endif

rer/rer(-1)	= ppi_S*ppi_ast/ppi;																		//*****Real exchange rate
R_ast		= R_W*exp(-PHI_ast/100*(rer*b_ast_Tot/gdpn-s_bast))*xi_R*zetau;     						//*****Foreign interest rate and country premium observable a unobservable shocks
x_H_ast		= (p_H/rer)^-ETA_ast*y_ast;																	//*****Demand of home good by foreigners

% -------------------------------------------------------------------------
    // Goverment spending and taxes
tau + R(-1)*bs_g(-1)/(a*ppi) + R_BL*q_BL*bl_g(-1)/a + CHI*rer*p_co*y_co =
			g + q_BL*bl_g + bs_g + GAMMA_d*PD_f*R_D(-1)*d_f(-1)/(a*ppi)
			+GAMMA_bh*PD_h*R_BB*q_BB*bb_tot(-1)/a;	 													//(90)Government budget constraint
tau 		= ALPHA_TAU*gdpn + EPSILON_TAU*(bs_g_ss +
				q_BL_ss*bl_g_ss - bs_g - q_BL*bl_g);													//*****Definitions of Taxes

% -------------------------------------------------------------------------
	// Aggregation and market clearing
y_C			= c_i+c_p+i+i_h+g+nu;																		//*****Final goods market
c_p			= OMEGA_UP*c_u + (1-OMEGA_UP)*c_r; 		 													//*****Consumption patient households
nu*a*ppi 	= GAMMA_d*PD_f*R_D(-1)*d_f(-1)+ GAMMA_bh*PD_h*R_BB*q_BB*bb_tot(-1)
				+mu_G_e*R_e*q_k(-1)*k(-1)+mu_G_i*R_h*q_h(-1)*h_i(-1)
				+mu_G_h*R_h_tilde*l_h(-1)*q_l(-1)+mu_G_f*R_f_tilde*l_f(-1)
@#if QLC == 1
		+(GAMA_n/2)*(n_tilde/n_tilde(-1) - 1)^2*(y_Z*p_Z*a*ppi)											//*****Total default costs (CSV and deposit recovery)
@#endif
;
y_H			= x_H+x_H_ast;																				//*****Home goods market
y_F			= x_F;											 											//*****Foreign goods market
h			= h_i+h_p;																					//*****Housing market
h_p			= OMEGA_UP*h_u + (1-OMEGA_UP)*h_r;			 												//*****Housing patient household
bl_priv 	= OMEGA_UP*bl_u + (1-OMEGA_UP)*bl_r;				 										//*****Total private holdings of long-term government bonds
bs_priv		= OMEGA_UP*bs_u;						 													//(100)*****Total private holdings of short-term government bonds
bb_tot		= OMEGA_UP*bb_u;																			//*****Bank bonds total vs unrestricted holdings
b_ast_Tot	= OMEGA_UP*b_ast_u;																			//*****Foreign Currency Bonds
0        	= bl_priv + (bl_cb-1) + bl_g;					 											//*****Aggregate net holdings in zero net supply
0			= bs_priv + bs_g;							 												//*****Aggregate net holdings in zero net supply
d_f			= OMEGA_UP*d_u;																				//*****Deposits markets
@#if fin_fric
	//zeta_L 	= (((q_BL*bl_u+q_BB*bb_u)/(bs_u+d_u+rer*b_ast_u))^(ETA_ZETA_L)*epsilon_L);					//*****Transaction cost on long-term government bonds //yy_mcc
	zeta_L 	= ((bb_u + bl_u)/(bl_u_ss+bb_u_ss))^ETA_ZETA_L*epsilon_L;
	R_D_tilde 	= R_D(-1)*(1-GAMMA_d*PD_f);																//*****Definition R_D_tilde
	R_BB_tilde 	= R_BB*(1-GAMMA_bh*PD_h);																//*****Definition R_BB_tilde
@#endif
R_BL 		= 1/q_BL+KAPPA_B;																			//*****Long-term rate
R_BB 		= 1/q_BB+KAPPA_BH;                                                  						//*****Long-term banking rate
Rnom_BL    	= R_BL*ppi;
p_H/p_H(-1) = ppi_H/ppi;																				//*****Relative prices p^H
p_F/p_F(-1)	= ppi_F/ppi;																				//*****Relative prices p^F
ppi_W		= w/w(-1)*a*ppi;																			//*****Average wage inflation
ppi_Wtilde	= w_tilde/w_tilde(-1)*ppi_W;																//*****Reoptimizing wage inflation
y_H*Xim_h 	= x_Z;																						//*****Aggregation of home goods
y_Z 		= z*(k(-1)/a)^ALPHA*n_tilde^(1-ALPHA); 														//*****Aggregation of wholesale domestic goods
y_Z 		= x_Z;																						//*****Wholesale domestic goods market
Xim_h		= (1-THETA_H)*p_tilde_H^(-EPSILON_H)+
				THETA_H*(ppi_IH/ppi_H)^(-EPSILON_H)*Xim_h(-1);											//*****Price dispersion home good
imp			= y_F*Xi_F;																					//*****Imports
Xi_F		= (1-THETA_F)*p_tilde_F^(-EPSILON_F)+
				THETA_F*(ppi_IF/ppi_F)^(-EPSILON_F)*Xi_F(-1);											//*****Price dispersion foreign good
n			= n_tilde*Xi_W;																				//*****Aggregation of labor services
Xi_W		= (1-THETA_W)*w_tilde^(-EPSILON_W)+
				THETA_W*(ppi_IW/ppi_W)^(-EPSILON_W)*Xi_W(-1);											//*****Wage dispersion
n			= n_p + n_i;																				//*****Total Household labor supply
n_p			= n_i;																						//*****Equality labor patient & impatient
n_p			= OMEGA_UP*n_u + (1-OMEGA_UP)*n_r;															//*****Total Patient Household labor supply
gdp			= c_i+c_p+i+i_h+g+x_H_ast+y_co-imp;			 												//*****GDP
gdpn		= c_i+c_p+i+i_h+g+tb;																		//*****Nominal GDP
tb			= p_H*x_H_ast+rer*p_co*y_co-rer*xi_m*imp;													//*****Trade balance
rer*b_ast_Tot	= rer/a/ppi_ast*b_ast_Tot(-1)*R_ast(-1)
					+tb+rer*ren_ast-(1-CHI)*rer*p_co*y_co;												//(129)*****Balance of payments
@#if fin_fric
n_u			= n_r;																						//(130)*****Equality of labor unrestricted & restricted
@#endif
% =========================================================================
%	Equilibrium Equations - No Financial Frictions
% =========================================================================

@#if !fin_fric

% -------------------------------------------------------------------------
	// Unrestricted Patient households
zeta_L          = 0;                                													//(No)Transaction cost on long-term government bonds

% -------------------------------------------------------------------------
	// Restricted Patient households
c_hat_r			= 0;																					//CES
la_r			= 0;																					//FOC consumption
c_r				= 0;																					//FOC consumption
h_r				= 0;																					//FOC housing
bl_r			= 0;																					//FOC long-term government bonds
n_r 			= 0;																					//Labor of the Restricted Patient Household
Theta_rp		= 0; 																					//Preference Shifter, Restricted Patient
Xi_Tilde_rp		= 0; 																					//Preference Shifter, Restricted Patient

% -------------------------------------------------------------------------
	// Impatient households
om_bar_i		= 0;																					//Idiosyncratic shock threshold
R_i 	 		= 1/q_l + KAPPA_LOAN;									   								//Mortgage Yield to Maturity
varrho*la_i*q_h = BETA_i*varrho(+1)*((c_hat_i(+1)*a(+1))^(-SIGMA)*xi_h(+1)*
	(O_CHAT*c_hat_i(+1)*a(+1)/(xi_h(+1)*(h_i-PHI_HH*h_i(-1)/a)))^(1/ETA_CHAT)
	+(1-DELTA_h)*la_i(+1)*a(+1)^(-SIGMA)*q_h(+1));														//FOC housing
varrho*la_i*q_l	= BETA_i*varrho(+1)*la_i(+1)*a(+1)^(-SIGMA)*q_l(+1)*R_i(+1);							//FOC loans
c_i+q_h*h_i-q_l*l_h	= w*n_i+q_h*(1-DELTA_h)*h_i(-1)/a-R_i*q_l*l_h(-1)/(a);								//Budget constraint
PD_i			= 0;																					//Probability of default on mortgage loans (imposed)

% -------------------------------------------------------------------------
	// Entrepreneurs
om_bar_e		= 0;																					//Idiosyncratic shock threshold	
c_e				= psi_e;																				//Dividends
n_e				= 0;																					//No retained profits (imposed)
psi_e			= (R_e-R_L(-1))*l_f(-1)/(a*ppi);														//Final wealth
R_e(+1)			= R_L;																					//FOC capital
la_e			= 0;																					//FOC leverage (useless)
PD_e			= 0;																					//Probability of default on corporate loans (imposed)

% -------------------------------------------------------------------------
	// Bankers and banking system
e_f				= 0;																					//No equity for banks
e_h				= 0;																					//No equity for banks
c_b				= xi_CHI_b*psi_b;																		//Dividends
n_b				= 0;																					//No retained profits (imposed)
psi_b			= (R_i-R)*q_l*l_h(-1)/(a);																//Final wealth

% -------------------------------------------------------------------------
	// F Banks
om_bar_f 		= 0;																					//*****Idiosyncratic shock threshold	
R_L				= R;																					//Competitive banks
rho_f			= R;																					//Return on equity is actually undefined as there is no equity... =0 for simplicity
R_f_tilde		= R_L(-1);																				//Return on assets
PD_f			= 0;																					//Default probability of bank (imposed)

% -------------------------------------------------------------------------
	// H Banks
rho_h			= R;																					//Return on equity is actually undefined as there is no equity... =0 for simplicity
om_bar_h		= 0;																					//*****Idiosyncratic shock threshold
rho_tilde_h		= R;																					//Irrelevant
R_h_tilde		= R_i(-1);																				//Return on assets
R_i				= R+IOTA*exp(PHI_i*(q_l*l_h-q_l*BAR_l_h));												//Spread-charging banks
PD_h			= 0;																					//Default probability of bank (imposed)

% -------------------------------------------------------------------------
	// Aggregation + Others
R_D_tilde 		= R_D(-1);																				//*****Definition R_D_tilde
R_BB_tilde 		= R_BB;																					//*****Definition R_BB_tilde
@#endif


% =========================================================================
%	Other Equations
% =========================================================================

% -------------------------------------------------------------------------
	// Additional definitions
c				= c_i+c_p;																				//(131)*****Total consumption
c_p_i			= c_p/c_i;																				//(132)*****Consumption of patients over impatient HHs
i_agg			= i+i_h;																				//(133)*****Total investment
l				= l_h*q_hat_l+l_f;																		//(134)*****Total loans --------Revisar
spread_RBL_R	= R_BL*ppi-R;																			//Spread: long rate - risk-free rate
spread_RL_RD	= R_L-R_D;
spread_RL_R     = R_L-R;																				//Spread: corporate loans - deposits
spread_Ri_RD	= R_i*ppi-R_D;																			//Spread: housing loans - deposits
spread_Ri_RBL	= R_i-R_BL;																			//Spread: housing loans - deposits
spread_RD_R		= R_D-R;																				//Spread: deposits - risk-free rate
lev_e			= R_L*l_f/(q_k*k);																		//Leverage of entrepreneurs
nu_gdp			= nu/gdp;																				//Default costs as fraction of GDP
lab_income      = w*n;                  																//Labor Income
gdpR=(gdp-y_co);
mort_tot = q_l*l_h ;
mortVpar_tot = q_hat_l*l_h ;
ltot = mortVpar_tot+l_f ; 
brechacred = 100*((mortVpar_tot+l_f)/(gdpn+gdpn(-1)+gdpn(-2)+gdpn(-3)) - (mortVpar_tot_ss+l_f_ss)/(4*gdpn_ss));
lev_tot         = n_b / (l_f_ss + 0.6*mortVpar_tot_ss);
lev_f           = e_f / l_f_ss; 
lev_h           = e_h / mortVpar_tot_ss; 

% -------------------------------------------------------------------------
% Observables 
% -------------------------------------------------------------------------
	// Observables: dlog dmean 
gam_YR_obs		= 100*log((gdpR/(gdpR(-1)))*a/a_ss);
gam_Y_obs		= 100*log((gdp/(gdp(-1)))*a/a_ss);
gam_YCo_obs		= 100*log(y_co/y_co(-1)*a/a_ss);
gam_C_obs		= 100*log(c/c(-1)*a/a_ss);
gam_I_obs		= 100*log(i_agg/i_agg(-1)*a/a_ss);
gam_IK_obs		= 100*log(i/i(-1)*a/a_ss);
gam_IH_obs		= 100*log(i_h/i_h(-1)*a/a_ss);
gam_G_obs		= 100*log(g/g(-1)*a/a_ss);
gam_Ystar_obs	= 100*log(y_ast/y_ast(-1)*a/a_ss);
gam_WN_obs		= 100*log(w/w(-1)*a*ppi/a_ss/ppi_ss);
pi_obs			= 100*log(ppi/ppi_ss);
piZ_obs         = pi_obs;
piM_obs			= 100*log(ppi_S*ppi_ast*xi_m/xi_m(-1)/ppi_S_ss/ppi_ast_ss);
pistar_obs		= 100*log(ppi_ast/ppi_ast_ss);
piCostar_obs	= 100*log(p_co/p_co(-1)*ppi_ast/ppi_ast_ss);
piQH_obs		= 100*log(q_h/q_h(-1));
gam_lf_obs      = 100*log((l_f)/(l_f(-1))*a*ppi/a_ss/ppi_ss);  
gam_lh_obs      = 100*log((q_hat_l*l_h)/(q_hat_l(-1)*l_h(-1))*a*ppi/a_ss/ppi_ss);  
gam_ltot_obs	= 100*log((l_f+q_hat_l*l_h)/(l_f(-1)+q_hat_l(-1)*l_h(-1))*a*ppi/a_ss/ppi_ss);
ccyb_obs        = ccyb;
// Observables: tasas demean 
R_obs			= 100*log(R/R_ss);
Rstar_obs		= 100*log(R_W/R_W_ss);
RLG_obs			= 100*log(R_BL/R_BL_ss);
RD_obs			= 100*log(R_D/R_D_ss);
RI_obs			= 100*log(R_i/R_i_ss);
RL_obs			= 100*log(R_L/R_L_ss);
xi_obs			= 100*log(exp(-PHI_ast/100*(rer*b_ast_Tot/gdpn-s_bast))*(xi_R/xi_R_ss));

// Observables: demean 
@#if fin_fric
	Roe_aux			= 100*log(
							( (rho_f   *e_f(-1)+rho_tilde_h   *e_h(-1))/n_b(-1) ) /
							( (rho_f_ss*e_f_ss +rho_tilde_h_ss*e_h_ss) /n_b_ss  )
							);
	Roe_obs			=	(Roe_aux+Roe_aux(-1)+Roe_aux(-2)+Roe_aux(-3))/4;
@#else
	Roe_obs			= 0;
	Roe_aux			= 0;
@#endif

 //Observables: demean  
rer_obs			= 100*log(rer/rer_ss);
stb_obs			= 100*(tb/gdpn-tb_ss/gdpn_ss); 	 														//(166)
gam_N_obs		= 100*log(n/n(-1));
gam_ner_obs 	= 100*log(ppi_S/ppi_S_ss);
gam4_YR_obs     = gam_YR_obs + gam_YR_obs(-1) + gam_YR_obs(-2) + gam_YR_obs(-3);
piZ4_obs        = piZ_obs + piZ_obs(-1) +piZ_obs(-2) + piZ_obs(-3);

 
% -------------------------------------------------------------------------
	// BGG Auxiliaries
@#define agent = ["i", "e", "f", "h"]
@#for ag in agent
@#if fin_fric
	aux_@{ag}		= (log(om_bar_@{ag})+0.5*(sigma_@{ag}@{ag})^2)/sigma_@{ag}@{ag};					//Auxiliary: normally distributed transformation of omega_bar
	Gam_@{ag}		= 1-normcdf(aux_@{ag}-sigma_@{ag}@{ag})-om_bar_@{ag}*(1-normcdf(aux_@{ag}));		//Borrowers's share of the project (1-Gamma(omega))
	Gam_G_@{ag}		= (1-MU_@{ag})*normcdf(aux_@{ag}-sigma_@{ag}@{ag})+om_bar_@{ag}*(1-normcdf(aux_@{ag}));	//Lender's share of the project net of monitoring cost (Gamma(omega)-mu*G(omega))
	mu_G_@{ag}		= 1-(Gam_@{ag}+Gam_G_@{ag});														//Monitoring cost (mu*G(omega))
	Gam_der_@{ag}	= 1-(normcdf(aux_@{ag}));															//First derivative of Gamma(omega)
	Gam_G_der_@{ag}	= 1-(normcdf(aux_@{ag}))-MU_@{ag}*normpdf(aux_@{ag})/sigma_@{ag}@{ag};				//Derivative of (Gamma(omega)-mu*G(omega)) wrt omega
@#else
	aux_@{ag}		= 0; 																				//Undefined object
	Gam_@{ag}		= 0; 																				//Undefined object
	Gam_G_@{ag}		= 0; 																				//Undefined object
	mu_G_@{ag}		= 0; 																				//Undefined object
	Gam_der_@{ag}	= 0;																				//Undefined object
	Gam_G_der_@{ag}	= 0; 																				//Undefined object
@#endif
@#endfor

% -------------------------------------------------------------------------
	// Exogenous processes

@#if fin_fric
    @#define mediumshocks = ["a", "g", "ppi_ast", "xi_m", "y_ast", "z", "y_co", "epsilon_L", "sigma_ee"]
	@#define bigshocks =["xi_R", "R_W" , "e_m",  "zetau"] //
@#else
	@#define bigshocks =[ ]
@#endif


@#for va in vars_exo
    @#if va in mediumshocks
        log(@{va}/@{va}_ss)=RHO_@{va}*log((@{va}(-1))/@{va}_ss)+u_@{va}/100;
    @#elseif va in bigshocks 
        log(@{va}/@{va}_ss)=RHO_@{va}*log((@{va}(-1))/@{va}_ss)+u_@{va}/1000; 
    @#else 
        log(@{va}/@{va}_ss)=RHO_@{va}*log((@{va}(-1))/@{va}_ss)+u_@{va}/10; 
    @#endif
@#endfor


end;
