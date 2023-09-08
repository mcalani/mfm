%-------------------------------------------------------------------------%

function [ys,params,check] = MAFIN_03_est_steadystate(ys,exo,M_,options_)
%keyboard
%-------------------------------------------------------------------------%
% This function computes the SS for MAFIN
% Inputs: 
%   - ys        [vector] vector of initial values for the steady state of
%                   the endogenous variables
%   - exo       [vector] vector of values for the exogenous variables
%   - M_        [structure] Dynare model structure
%   - options_  [structure] Dynare options structure
%
% Output: 
%   - ys        [vector] vector of steady state values for the the endogenous variables
%   - params    [vector] vector of parameter values
%   - check     [scalar] set to 0 if steady state computation worked and to
%                    1 of not (allows to impose restrictions on parameters)
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
%	Intro
%-------------------------------------------------------------------------%

	% Read out parameters to access them with their name
%NumberOfParameters = M_.param_nbr;
% for ii = 1:NumberOfParameters
% 	paramname = M_.param_names{ii};
%     %string=[ paramname ' = M_.params(' int2str(ii) ');'];
%     string=sprintf('%s=M_.params(%d);',paramname,ii);
%     eval(string)
% 	%eval([ paramname ' = M_.params(' int2str(ii) ');'])
% end

  % Read out parameters to access them with their name
NumberOfParameters = M_.param_nbr;
for ii = 1:NumberOfParameters
	paramname = M_.param_names{ii};
	eval([ paramname ' = M_.params(' int2str(ii) ');'])
end
	% Read out shocks to access them with their name
NumberOfExogenous = M_.exo_nbr;
for ii = 1:NumberOfExogenous
	exoname = M_.exo_names{ii}(3:end);
	eval([ exoname ' = M_.params(loc(M_.param_names,''' exoname '_ss''));'])
end

% get_params

	% Initialize indicator
check = 0;

	% Calibrated SS values
p_H		= M_.params(loc(M_.param_names,'p_H_ss'));
ppi_S	= M_.params(loc(M_.param_names,'ppi_S_ss'));
n		= M_.params(loc(M_.param_names,'n_ss'));

if FinFric
  R_D		= M_.params(loc(M_.param_names,'R_D_ss'));  
  R_L		= M_.params(loc(M_.param_names,'R_L_ss'));
  Rnom_i	= M_.params(loc(M_.param_names,'Rnom_i_ss'));
  Rnom_BL	= M_.params(loc(M_.param_names,'Rnom_BL_ss'));
end


if FinFric   
%-------------------------------------------------------------------------%
%	Steady State Equations - Financial Frictions
%-------------------------------------------------------------------------%
%disp('STEADY STATE WITH FINANCIAL FRICTIONS')
buffer_guide 	= 0;
ccyb         	= 0.0;

phi_f       = (phi_f_base+ccyb) + 0.0433;
phi_h       = (phi_h_base+ccyb)*0.6 + 0.0433;


%% Bloque 1: sin soluciones numericas 
ppi = ppi_T;
R = ppi*a^SIGMA/BETA_up;
R_D_tilde = R;
R_ast	= R/ppi_S;
R_W = R_ast/xi_R;
ppi_ast = ppi/ppi_S;
ppi_H = ppi;
ppi_F = ppi;
ppi_IH = ppi;
ppi_IF = ppi;
ppi_W = a*ppi;
ppi_Wtilde = a*ppi;
ppi_IW = a*ppi;
p_tilde_H = 1;
p_tilde_F = 1;
w_tilde = 1;
mc_H = (EPSILON_H-1)/EPSILON_H;
mc_F = (EPSILON_F-1)/EPSILON_F;
mc_W = (EPSILON_W-1)/EPSILON_W;
Xim_h = 1;
Xi_F = 1;
Xi_W = 1;
q_k = 1/xi_i;
q_h = a^(N_H*SIGMA)*VARPHI_H_0/(BETA_up^N_H*xi_ih)*...
				(1-(BETA_up*RHO_VARPHI_H/a^SIGMA)^(N_H+1))/(1-BETA_up*RHO_VARPHI_H/a^SIGMA);
R_h = ppi*(1-DELTA_h);
n_tilde = n; 
rho_f = a*ppi/(1-CHI_b);
rho_h = rho_f;
rho_tilde_h = rho_h; 
PD_d = (1-R_D_tilde/R_D)/GAMMA_d;
PD_h = PD_d;
PD_f = PD_d;
R_BL = Rnom_BL/ppi; 
BETA_rp = a^SIGMA/R_BL; 
R_i = Rnom_i/ppi;     
R_hat_i = R_i;
q_hat_l =(R_hat_i-KAPPA_LOAN)^(-1);                
q_l = q_hat_l ; 
R_BB_tilde = R_BL; 
R_BB = R_BB_tilde/(1-GAMMA_bh*PD_h);
q_BB = 1/(R_BB-KAPPA_BH); 
q_BL = 1/(R_BL-KAPPA_B);
Dl  = a;

	% Bank numerical solutions
options		= optimoptions('fsolve','FunctionTolerance',1e-10,'OptimalityTolerance',1e-10,'MaxIterations',1000,'MaxFunctionEvaluations',1000,'Display','off');
x_guess		= [0.9,0.05];
[x,~]		= fsolve(@SS_num_banks,x_guess,options,phi_f,R_D,rho_f,PD_f);
om_bar_f	= x(1);
sigma_ff	= x(2);

x_guess		= [0.9052,0.0401];
[x,~]		= fsolve(@SS_num_banks_h,x_guess,options,phi_h,R_BB,rho_h,PD_h,ppi);
om_bar_h	= x(1);
sigma_hh	= x(2);

	% Bank BGG auxiliaries and ROA
aux_f = (log(om_bar_f)+0.5*sigma_ff^2)/sigma_ff;
Gam_f = 1-normcdf(aux_f-sigma_ff)-om_bar_f*(1-normcdf(aux_f));
Gam_G_f = (1-MU_f)*normcdf(aux_f-sigma_ff)+om_bar_f*(1-normcdf(aux_f));
Gam_der_f = 1-normcdf(aux_f);
Gam_G_der_f = 1-normcdf(aux_f)-MU_f*normpdf(aux_f)/sigma_ff;
mu_G_f = 1-Gam_f-Gam_G_f;
aux_h = (log(om_bar_h)+0.5*sigma_hh^2)/sigma_hh;
Gam_h = 1-normcdf(aux_h-sigma_hh)-om_bar_h*(1-normcdf(aux_h));
Gam_G_h = (1-MU_h)*normcdf(aux_h-sigma_hh)+om_bar_h*(1-normcdf(aux_h));
Gam_der_h = 1-normcdf(aux_h);
Gam_G_der_h = 1-normcdf(aux_h)-MU_h*normpdf(aux_h)/sigma_hh;
mu_G_h = 1-Gam_h-Gam_G_h;

R_f_tilde = rho_f*phi_f/Gam_f;
R_h_tilde = rho_h*phi_h/Gam_h;
  
	% Entrepreneur numerical solutions
options		= optimoptions('fsolve','FunctionTolerance',1e-10,'OptimalityTolerance',1e-10,'MaxIterations',10000,'MaxFunctionEvaluations',10000,'Display','off');
x_guess		= [0.1,0.9];
[x,ff]		= fsolve(@SS_num_e,x_guess,options,MU_e,R_f_tilde,R_L,CHI_e,a,ppi);
om_bar_e	= x(1);
sigma_ee	= x(2);
aux_e		= (log(om_bar_e)+0.5*sigma_ee^2)/sigma_ee;
Gam_e		= 1-normcdf(aux_e-sigma_ee)-om_bar_e*(1-normcdf(aux_e));
Gam_G_e		= (1-MU_e)*normcdf(aux_e-sigma_ee)+om_bar_e*(1-normcdf(aux_e));
Gam_der_e	= 1-normcdf(aux_e);
Gam_G_der_e	= 1-normcdf(aux_e)-MU_e*normpdf(aux_e)/sigma_ee;
mu_G_e		= 1-Gam_e-Gam_G_e;
PD_e		= normcdf(aux_e);

	% Impatient household numerical solutions
x_guess		= [0.8,.05];                                     
[x,ff]		= fsolve(@SS_num_i,x_guess,options,MU_i,a,SIGMA,ppi,R_h_tilde,BETA_i,R_i);
om_bar_i 	= exp(x(1))/(1+exp(x(1)));
sigma_ii 	= exp(x(2))/(1+exp(x(2)));
aux_i		= (log(om_bar_i)+0.5*sigma_ii^2)/sigma_ii;
Gam_i		= 1-normcdf(aux_i-sigma_ii)-om_bar_i*(1-normcdf(aux_i));
Gam_G_i		= (1-MU_i)*normcdf(aux_i-sigma_ii)+om_bar_i*(1-normcdf(aux_i));
Gam_der_i	= 1-normcdf(aux_i);
Gam_G_der_i	= 1-normcdf(aux_i)-MU_i*normpdf(aux_i)/sigma_ii;
mu_G_i		= 1-Gam_i-Gam_G_i;
PD_i		= normcdf(aux_i);

	% Second block of analytical solutions
R_e			= R_f_tilde*a*ppi/(a*ppi*Gam_G_e+Gam_e*(1-CHI_e)*R_f_tilde);
r_k			= q_k*(R_e/ppi-(1-DELTA_k));
p_Z			= p_H*mc_H;
mc_Z		= p_Z;
w			= (ALPHA^ALPHA*(1-ALPHA)^(1-ALPHA)*mc_Z*z/r_k^ALPHA)^(1/(1-ALPHA));
k			= ALPHA/(1-ALPHA)*n_tilde*w*a/r_k;
y_Z 		= z*(k/a)^ALPHA*n_tilde^(1-ALPHA);
x_Z 		= y_Z;
i 			= k*(1-(1-DELTA_k)/a)/xi_i;
y_H 		= x_Z/Xim_h;
psi_e		= Gam_e*R_e*q_k*k/(a*ppi);
n_e			= (1-CHI_e*xi_CHI_e)*psi_e;
c_e			= CHI_e*xi_CHI_e*psi_e;
la_e		= Gam_der_e/(Gam_f*Gam_G_der_e);
l_f			= q_k*k-n_e;
e_f			= phi_f*l_f;
d_f			= l_f-e_f;
d_u     	= d_f/OMEGA_UP;
h			= r_h_k*q_k*k/q_h;
i_ah		= h*a^N_H/xi_ih*(1-(1-DELTA_h)/a);
i_h			= i_ah*VARPHI_H_0*(1-(RHO_VARPHI_H/a)^(N_H+1))/(1-RHO_VARPHI_H/a);
p_F			= ((1-OMEGA*p_H^(1-ETA))/(1-OMEGA))^(1/(1-ETA));
rer			= mc_F*p_F/xi_m;


	% Numerical solution for l_h
options		= optimoptions('fsolve','FunctionTolerance',1e-10,'OptimalityTolerance',1e-10,'MaxIterations',1000,'MaxFunctionEvaluations',1000,'Display','off');
l_h_guess	= 0.03;
[l_h,~]		= fsolve(@SS_num_l_h,l_h_guess,options,R_i,om_bar_i,R_h,q_h,phi_h,d_f,a,ppi,GAMMA_d,PD_d,R_D,...
    mu_G_e,R_e,q_k,k,mu_G_i,mu_G_h,R_h_tilde,mu_G_f,R_f_tilde,...
    l_f,p_H,y_H,p_F,ETA,rer,xi_m,Xi_F,OMEGA,s_co,s_tb,h,w,n,...
    Gam_i,PHI_c,SIGMA,PHI_HH,BETA_i,Gam_G_i,BETA_up,DELTA_h,i,i_h,...
	s_g,xi_h,ETA_CHAT,OMEGA_UP,BETA_rp,q_l,q_BB,R_BB,GAMMA_bh,PD_h, ALPHA_BLG, ...
  ALPHA_BSG, q_BL,s_bast,OMEGA_BL,d_u,R_BL);

  % Third block of closed form solutions 
h_i = R_hat_i*q_hat_l*l_h*ppi/(om_bar_i*R_h*q_h);
e_h = phi_h*q_l*l_h;
n_b = e_f+e_h; 
psi_b = n_b/(1-CHI_b*xi_CHI_b);
c_b = CHI_b*xi_CHI_b*psi_b;
bb_tot = (1-phi_h)*q_l*l_h/q_BB;
nu = 1/(a*ppi)*(GAMMA_d*PD_d*R_D*d_f +GAMMA_bh*PD_h*R_BB*q_BB*bb_tot +mu_G_e*R_e*q_k*k...
				+mu_G_i*R_h*q_h*h_i + mu_G_h*R_h_tilde*l_h*q_l + mu_G_f*R_f_tilde*l_f) ;
gdpn		= (p_H*y_H+p_F^(-ETA)*(p_F-rer*xi_m*Xi_F)*(1-OMEGA)*nu-nu)/...
				(1-s_co-(1-s_tb)*p_F^(-ETA)*(p_F-rer*xi_m*Xi_F)*(1-OMEGA));
tb			= s_tb*gdpn;
g			= s_g*gdpn;
y_co		= s_co*gdpn/(p_co*rer);
b_ast_Tot = s_bast*gdpn/rer;
y_C			= gdpn+nu-tb;
x_F			= (1-OMEGA)*p_F^-ETA*y_C;
x_H			= OMEGA*p_H^-ETA*y_C;
x_H_ast		= y_H-x_H;
y_ast		= x_H_ast*(p_H/rer)^ETA_ast;
y_F			= x_F;
imp			= y_F*Xi_F;
h_p			= max(0,h-h_i);
c_i			= w*n/2+q_h*h_i*(Gam_i*R_h/(a*ppi)-1)+q_l*l_h;

O_CHAT		= (a^(-SIGMA*ETA_CHAT)*xi_h^(ETA_CHAT-1)* (a*c_i*(1-PHI_c/a)/(h_i*(1-PHI_HH/a))) ...
				* (1/BETA_i*(q_h-(Gam_G_i)*R_h*q_h/R_h_tilde)-(a^(-SIGMA)*Gam_i*R_h)/ppi*q_h)^(-ETA_CHAT) + 1)^-1;
           
c_hat_i		= ((1-O_CHAT)^(1/ETA_CHAT)*(c_i*(1-PHI_c/a))^((ETA_CHAT-1)/ETA_CHAT)...
				+(O_CHAT)^(1/ETA_CHAT)*(xi_h*h_i/a*(1-PHI_HH/a))^((ETA_CHAT-1)/ETA_CHAT))^(ETA_CHAT/(ETA_CHAT-1));
la_i		= c_hat_i^(-SIGMA)*((1-O_CHAT)*c_hat_i/(c_i*(1-PHI_c/a)))^(1/ETA_CHAT);
la_h        =la_i/(rho_tilde_h*phi_h); 	
    
  % Bond positions 
bl_g 		= ALPHA_BLG*gdpn/q_BL;
bs_g    	= ALPHA_BSG*gdpn;
bl_priv 	= - bl_g;
bs_priv 	= - bs_g;
bs_u 		= bs_priv/OMEGA_UP;
b_ast_u 	= b_ast_Tot/OMEGA_UP;
bb_u 		= bb_tot/OMEGA_UP;
bl_u 		= (OMEGA_BL*(bs_u + rer*b_ast_u + d_u ) - bb_u*q_BB)/q_BL;
bl_r 		= (bl_priv - OMEGA_UP*bl_u)/(1-OMEGA_UP); 
bl_cb 		= 1; 

  % Find h_r using eq 14 instead of upsilon 
corcho 		= a^(SIGMA*ETA_CHAT-1)*xi_h^(1-ETA_CHAT)*(q_h/BETA_rp-(1-DELTA_h)*a^(-SIGMA)*q_h)^ETA_CHAT*(1-O_CHAT)/O_CHAT*(1-PHI_HH/a)/(1-PHI_c/a);
h_r 		= (q_BL*bl_r*(R_BL/a-1)+w*n/2)/(q_h-q_h/a*(1-DELTA_h)+corcho);
%c_r			= ((a^SIGMA)*(xi_h*h_r*(1-PHI_HH/a)*(1-O_CHAT)/((1-PHI_c/a)*O_CHAT*a))^(1/ETA_CHAT)...
%				*(q_h/BETA_rp-(1-DELTA_h)*a^-SIGMA*q_h))^ETA_CHAT;     
c_r         = corcho*h_r;
c_hat_r		= ((1-O_CHAT)^(1/ETA_CHAT)*(c_r*(1-PHI_c/a))^((ETA_CHAT-1)/ETA_CHAT)...
				+(O_CHAT)^(1/ETA_CHAT)*(xi_h*h_r/a*(1-PHI_HH/a))^((ETA_CHAT-1)/ETA_CHAT))^(ETA_CHAT/(ETA_CHAT-1));
la_r		= c_hat_r^(-SIGMA)*((1-O_CHAT)*c_hat_r/(c_r*(1-PHI_c/a)))^(1/ETA_CHAT);
h_u  		= (h_p-(1-OMEGA_UP)*h_r)/OMEGA_UP;
c_u   		= ((a^SIGMA/BETA_up-(1-DELTA_h))*q_h/xi_h)^ETA_CHAT*(1-O_CHAT)/a/O_CHAT*h_u*xi_h*(1-PHI_HH/a)/(1-PHI_c/a);    
c_hat_u		= ((1-O_CHAT)^(1/ETA_CHAT)*(c_u*(1-PHI_c/a))^((ETA_CHAT-1)/ETA_CHAT)...
				+(O_CHAT)^(1/ETA_CHAT)*(xi_h*h_u/a*(1-PHI_HH/a))^((ETA_CHAT-1)/ETA_CHAT))^(ETA_CHAT/(ETA_CHAT-1));
la_u		= c_hat_u^-SIGMA*((1-O_CHAT)*c_hat_u/(c_u*(1-PHI_c/a)))^(1/ETA_CHAT);
la_p		= OMEGA_UP*la_u + (1-OMEGA_UP)*la_r;
c_p			= OMEGA_UP*c_u + (1-OMEGA_UP)*c_r;
c  			= c_p+c_i;
n_p			= n/2;
n_i			= n_p;
n_u			= n_p;
n_r			= n_u;
Xi_Tilde_up	= (c_hat_u)^(SIGMA);
Theta_up	= Xi_Tilde_up*(c_hat_u)^(-SIGMA);
Xi_Tilde_i	= (c_hat_i)^(SIGMA);
Theta_i		= Xi_Tilde_i*(c_hat_i)^(-SIGMA);	
Xi_Tilde_rp	= (c_hat_r)^(SIGMA);
Theta_rp	= Xi_Tilde_rp*(c_hat_r)^(-SIGMA);
Theta		= ((OMEGA_UP*Theta_up+(1-OMEGA_UP)*Theta_rp)+Theta_i)/2;
la_W		= (la_p+la_i)/2;
xi_n		= mc_W*la_W*w/(Theta*n_tilde^VARPHI);
gdp			= c_i+c_p+i+i_h+g+x_H_ast+y_co-imp;
ren_ast		= b_ast_Tot*(1-R_ast/a/ppi)-tb/rer+(1-CHI)*p_co*y_co;
epsilon_L   = BETA_up*R_BL*a^(-SIGMA)-1;
zeta_L      = epsilon_L; 
tau 		= g + GAMMA_d*PD_d*R_D*d_f/a/ppi+ GAMMA_bh*PD_h*R_BB*q_BB*bb_tot/a...
				-(R/(a*ppi)-1)*bs_g - (R_BL/(a)-1)*q_BL*bl_g - CHI*rer*p_co*y_co;
ALPHA_TAU   = tau/gdpn;
f_H			= p_tilde_H^-EPSILON_H*y_H*mc_H/(1-THETA_H*BETA_up*a^(1-SIGMA));
f_F			= p_tilde_F^-EPSILON_F*y_F*mc_F/(1-THETA_F*BETA_up*a^(1-SIGMA));
f_W			= w_tilde^(-EPSILON_W*(1+VARPHI))*mc_W*n_tilde/(1-((OMEGA_UP*BETA_up+(1-OMEGA_UP)*BETA_rp)+BETA_i)/2*THETA_W*a^(1-SIGMA));
mort_tot = q_l*l_h ;
mortVpar_tot =  q_hat_l*l_h;
ltot = q_l*l_h+l_f ; 

else 
  
%-------------------------------------------------------------------------%
%	Steady State Equations - No Financial Frictions
%-------------------------------------------------------------------------%
%disp('STEADY STATE WITHOUT FINANCIAL FRICTIONS')
buffer_guide = 0;
ccyb        = 0.0;
phi_f       = (0.08+ccyb) + 0.0433;
phi_h       = (0.08+ccyb)*0.6 + 0.0433;

%phi_f		= 0.123;				%Chilean capital requirements of 8%, with 100% risk weight for corporate loans, +4.3% of excess capital in data
%phi_h		= 0.091;				%Chilean capital requirements of 8%, with 60% risk weight for housing loans, +4.3% of excess capital in data

	% First block of non-financial variables
ppi			= ppi_T;
R			= ppi*a^SIGMA/BETA_up;
R_ast		= R/ppi_S;
R_W			= R_ast/xi_R;
ppi_ast		= ppi/ppi_S;
ppi_H		= ppi;
ppi_F		= ppi;
ppi_IH		= ppi;
ppi_IF		= ppi;
ppi_W		= a*ppi;
ppi_Wtilde	= a*ppi;
ppi_IW		= a*ppi;
p_tilde_H	= 1;
p_tilde_F	= 1;
w_tilde		= 1;
mc_H		= (EPSILON_H-1)/EPSILON_H;
mc_F		= (EPSILON_F-1)/EPSILON_F;
mc_W		= (EPSILON_W-1)/EPSILON_W;
Xim_h		= 1;
Xi_F		= 1;
Xi_W		= 1;
q_k			= 1/xi_i;
q_h			= a^(N_H*SIGMA)*VARPHI_H_0/(BETA_up^N_H*xi_ih)*...
				(1-(BETA_up*RHO_VARPHI_H/a^SIGMA)^(N_H+1))/(1-BETA_up*RHO_VARPHI_H/a^SIGMA);
R_h			= ppi*(1-DELTA_h);
n_tilde		= n;

	% BGG auxiliaries
aux_i		= 0;
Gam_i		= 0;
Gam_G_i		= 0;
Gam_der_i	= 0;
Gam_G_der_i	= 0;
mu_G_i		= 0;
PD_i		= 0;
aux_e		= 0;
Gam_e		= 0;
Gam_G_e		= 0;
Gam_der_e	= 0;
Gam_G_der_e	= 0;
mu_G_e		= 0;
PD_e		= 0;
aux_h		= 0;
Gam_h		= 0;
Gam_G_h		= 0;
Gam_der_h	= 0;
Gam_G_der_h	= 0;
mu_G_h		= 0;
PD_h		= 0;
aux_f		= 0;
Gam_f		= 0;
Gam_G_f		= 0;
Gam_der_f	= 0;
Gam_G_der_f	= 0;
mu_G_f		= 0;
PD_f		= 0;
PD_d		= 0;

	% Trivial financial and housing variables
%R_i			= ppi*a^SIGMA/BETA_i;
R_i			= a^SIGMA/BETA_i;
	
%q_l			= 1;
q_l			= 1/(R_i - KAPPA_LOAN) ;	
q_hat_l 	= q_l 
	
xi_CHI_b	= 1;
xi_CHI_e	= 1;
xi_roe_r	= 1;	
xi_h		= 1;
R_D_tilde	= R;
R_D			= R_D_tilde;
nu			= 0;
om_bar_i	= 0;

c_e			= 0;
n_e			= 0;
om_bar_e	= 0;
psi_e		= 0;
la_e		= 0;
e_f			= 0;
e_h			= 0;
n_b			= 0;
R_L			= R;
R_f_tilde	= R_L;

om_bar_f	= 0;


rho_f		= R;

R_e			= R_L;
R_BL	    = R/ppi;
Rnom_BL    	= R_BL*ppi; 
R_BB_tilde 	= R_BL;
R_BB		= R_BB_tilde;
q_BB		= 1/(R_BB - KAPPA_BH);
BETA_rp     = ppi*a^SIGMA/R_BL;
l_h			= BAR_l_h;

%IOTA		= a_ss^SIGMA*ppi_T*(1/BETA_i-1/BETA_up);
IOTA		= a_ss^SIGMA*(1/BETA_i-ppi_T/BETA_up);

%psi_b		= IOTA*l_h/(a*ppi);
psi_b		= IOTA*q_l*l_h/(a);

c_b			= xi_CHI_b*psi_b;
bb_tot		= (q_l*l_h-e_h)/q_BB;

%rho_h		= 0;
rho_h		= R;

R_h_tilde	= R_i;

%om_bar_h	= (1-phi_h)*(R_BB/R_h_tilde)*ppi;
om_bar_h	= 0;

rho_tilde_h	= R; 
Rnom_i 		= R_i*ppi;  


	% Second block of non-financial variables
r_k			= q_k*(R_e/ppi-(1-DELTA_k));
p_Z			= p_H*mc_H;
mc_Z		= p_Z;
w			= (ALPHA^ALPHA*(1-ALPHA)^(1-ALPHA)*mc_Z*z/r_k^ALPHA)^(1/(1-ALPHA));
k			= ALPHA/(1-ALPHA)*n_tilde*w/r_k*a;
y_Z			= z*(k/a)^ALPHA*n_tilde^(1-ALPHA);
x_Z			= y_Z;
i			= k*(1-(1-DELTA_k)/a)/xi_i;
y_H			= x_Z/Xim_h;
p_F			= ((1-OMEGA*p_H^(1-ETA))/(1-OMEGA))^(1/(1-ETA));
rer			= mc_F*p_F/xi_m;
gdpn		= (p_H*y_H-nu+nu*(1-OMEGA)*p_F^-ETA*(p_F-rer*xi_m*Xi_F))/(1-s_co-(p_F-rer*xi_m*Xi_F)*(1-OMEGA)*p_F^-ETA*(1-s_tb));
tb			= s_tb*gdpn;
g			= s_g*gdpn;
b_ast_Tot	= s_bast*gdpn/rer;
b_ast_u 	= b_ast_Tot/OMEGA_UP;
y_co		= s_co*gdpn/(p_co*rer);
y_C			= gdpn+nu-tb;
x_F			= (1-OMEGA)*p_F^-ETA*y_C;
x_H			= OMEGA*p_H^-ETA*y_C;
x_H_ast		= y_H-x_H;
y_ast		= x_H_ast*(p_H/rer)^ETA_ast;
y_F			= x_F;
imp			= y_F*Xi_F;
h			= r_h_k*q_k*k/q_h;
i_ah		= h*a^N_H/xi_ih*(1-(1-DELTA_h)/a);
i_h			= i_ah*VARPHI_H_0*(1-(RHO_VARPHI_H/a)^(N_H+1))/(1-RHO_VARPHI_H/a);
n_p			= n/2;
n_i			= n_p;
n_u			= n_p/OMEGA_UP;	
n_r        	= 0;

	% Numerical solution
x_guess		= 0.3;
options 	= optimoptions('fsolve','FunctionTolerance',1e-10,'OptimalityTolerance',1e-10,'MaxIterations',1000,'MaxFunctionEvaluations',1000,'Display','off');
[x,~] 		= fsolve(@SS_num_NF,x_guess,options,h,q_h,DELTA_h,a,w,n_i,l_h,R_i,ppi,y_C,i,i_h,g,BETA_up,BETA_i,SIGMA,OMEGA_UP,ETA_CHAT,q_l);
c_i			= x;

%h_i			= (w*n_i+l_h*(1-R_i/(a*ppi))-c_i)/(q_h*(1-(1-DELTA_h)/a));
h_i			= (w*n_i+q_l*l_h*(1-R_i/(a))-c_i)/(q_h*(1-(1-DELTA_h)/a));

c_p			= y_C-c_i-i-i_h-g;
h_p			= h-h_i;
h_u			= h_p/OMEGA_UP;
h_r			= 0;

	% Third block of non-financial variables
O_CHAT		= (a^(-SIGMA*ETA_CHAT)*(a*c_i*(1-PHI_c/a)/(xi_h*(h_i*(1-PHI_HH/a))))...
				*(q_h/BETA_i- (1-DELTA_h)*(a^-SIGMA)*q_h)^(-ETA_CHAT) + 1)^-1;
c_hat_i		= ((1-O_CHAT)^(1/ETA_CHAT)*(c_i*(1-PHI_c/a))^((ETA_CHAT-1)/ETA_CHAT)...
				+(O_CHAT)^(1/ETA_CHAT)*(xi_h*h_i/a*(1-PHI_HH/a))^((ETA_CHAT-1)/ETA_CHAT))^(ETA_CHAT/(ETA_CHAT-1));
la_i		= c_hat_i^-SIGMA*((1-O_CHAT)*c_hat_i/(c_i*(1-PHI_c/a)))^(1/ETA_CHAT);
c_u			= ((a^SIGMA)*(xi_h*h_u*(1-PHI_HH/a)*(1-O_CHAT)/((1-PHI_c/a)*O_CHAT*a))^(1/ETA_CHAT)...
				*(q_h/BETA_up-(1-DELTA_h)*a^-SIGMA*q_h))^ETA_CHAT;
c_hat_u		= ((1-O_CHAT)^(1/ETA_CHAT)*(c_u*(1-PHI_c/a))^((ETA_CHAT-1)/ETA_CHAT)...
				+(O_CHAT)^(1/ETA_CHAT)*(xi_h*h_u/a*(1-PHI_HH/a))^((ETA_CHAT-1)/ETA_CHAT))^(ETA_CHAT/(ETA_CHAT-1));
la_u		= c_hat_u^-SIGMA*((1-O_CHAT)*c_hat_u/(c_u*(1-PHI_c/a)))^(1/ETA_CHAT);
Xi_Tilde_up	= (c_hat_u)^(SIGMA);
Theta_up	= Xi_Tilde_up*(c_hat_u)^(-SIGMA);
Xi_Tilde_i	= (c_hat_i)^(SIGMA);
Theta_i		= Xi_Tilde_i*(c_hat_i)^(-SIGMA);	
c_r			= 0;
c_hat_r		= 0;
la_r		= 0;
Xi_Tilde_rp	= 0;
Theta_rp	= 0;
Theta		= ((OMEGA_UP*Theta_up+(1-OMEGA_UP)*Theta_rp)+Theta_i)/2;
la_p		= OMEGA_UP*la_u + (1-OMEGA_UP)*la_r;
la_W		= (la_p+la_i)/2;
xi_n		= mc_W*la_W*w/n_tilde^VARPHI;
c_p			= OMEGA_UP*c_u + (1-OMEGA_UP)*c_r;
gdp			= c_i+c_p+i+i_h+g+x_H_ast+y_co-imp;
ren_ast		= b_ast_Tot*(1-R_ast/a/ppi)-tb/rer+(1-CHI)*p_co*y_co;
f_H			= p_tilde_H^-EPSILON_H*y_H*mc_H/(1-THETA_H*BETA_up*a^(1-SIGMA));
f_F			= p_tilde_F^-EPSILON_F*y_F*mc_F/(1-THETA_F*BETA_up*a^(1-SIGMA));
f_W			= w_tilde^(-EPSILON_W*(1+VARPHI))*mc_W*n_tilde/(1-((OMEGA_UP*BETA_up+(1-OMEGA_UP)*BETA_rp)+BETA_i)/2*THETA_W*a^(1-SIGMA));
l_f			= q_k*k-n_e;
d_f			= l_f-e_f;
d_u     	= d_f/OMEGA_UP;
q_BL        = 1/(R_BL-KAPPA_B);
bl_g 		= ALPHA_BLG*gdpn/q_BL;
bs_g    	= ALPHA_BSG*gdpn;
bs_priv		= -bs_g;
bl_priv 	= -bl_g;
bs_u		= bs_priv/OMEGA_UP;	
bb_u 		= bb_tot/OMEGA_UP;	
bl_r		= 0;
bl_u        = (bl_priv - (1-OMEGA_UP)*bl_r)/OMEGA_UP;
zeta_L		= 0;
epsilon_L	= 1;
bl_cb 		= 1; 
tau 		= g + GAMMA_d*PD_d*R_D*d_f/a/ppi+ GAMMA_bh*PD_h*R_BB*q_BB*bb_tot/a...
				-(R/(a*ppi)-1)*bs_g - (R_BL/(a)-1)*q_BL*bl_g - CHI*rer*p_co*y_co;
ALPHA_TAU   = tau/gdpn;
mort_tot = q_l*l_h ;
mortVpar_tot =  q_hat_l*l_h;

end  


%-------------------------------------------------------------------------%
%	Steady State Equations - Definitions and Observables
%-------------------------------------------------------------------------%

	% Additional definitions
c			= c_i+c_p;
c_p_i		= c_p/c_i;
i_agg		= i+i_h;
l			= l_f+l_h*q_l;
spread_RD_R	= R_D-R;
spread_RL_RD= R_L-R_D;
spread_RL_R= R_L-R;
spread_Ri_RD= R_i*ppi-R_D;
spread_RBL_R= R_BL*ppi-R;

spread_Ri_RBL = R_i-R_BL; 

lev_e		= R_L*l_f/(q_k*k);
nu_gdp		= nu/gdp;
lab_income  = w*n;
gdpR=(gdp-y_co);
brechacred = 0 ; 
ltot = q_l*l_h+l_f ; 
lev_tot = n_b / (l_f + 0.6*mortVpar_tot);
lev_f = e_f / (l_f);
lev_h = e_h / (mortVpar_tot);

	% Observables
gam_C_obs		= 0;
gam_G_obs		= 0;
gam_I_obs		= 0;
gam_lh_obs		= 0;
gam_lf_obs		= 0;
gam_IH_obs		= 0;
gam_IK_obs		= 0;
gam_N_obs		= 0;
gam_WN_obs		= 0;
gam_YCo_obs		= 0;
gam_YR_obs		= 0;
gam_Y_obs		= 0;
gam_Ystar_obs	= 0;
gam_ner_obs 	= 0;
pi_obs			= 0;
piZ_obs			= 0;
piCostar_obs	= 0;
piM_obs			= 0;
piQH_obs		= 0;
pistar_obs		= 0;
rer_obs			= 0;
RLG_obs			= 0;
R_obs			= 0;
RD_obs			= 0;
RI_obs			= 0;
RL_obs			= 0;
Rstar_obs		= 0;
stb_obs			= 0;
xi_obs			= 0;
gam_ltot_obs	= 0;
Roe_obs			= 0;
Roe_aux         = Roe_obs;
gam4_YR_obs     = 0;
piZ4_obs        = 0;
ccyb_obs        =0; 
%-------------------------------------------------------------------------%
%	Outro
%-------------------------------------------------------------------------%

	% Set SS values as parameters
NumberOfEndogenous = M_.orig_endo_nbr; %auxiliary variables are set automatically
for ii = 1:NumberOfEndogenous
	varname = M_.endo_names{ii};
  % eval(sprintf('%s_ss=%s;',varname,varname))
	eval([ varname '_ss = ' varname ';'])
end

	% Update Parameters
params=zeros(NumberOfParameters,1);
for iter = 1:length(M_.params) %update parameters set in the file
   % params(iter)=eval(M_.param_names{iter});
    parname=M_.param_names{iter};
    params(iter)=eval(parname);
	%eval([ 'params(' num2str(iter) ') = ' M_.param_names{iter} ';' ]);
end

params=NaN(NumberOfParameters,1);
for iter = 1:length(M_.params) %update parameters set in the file
	eval([ 'params(' num2str(iter) ') = ' M_.param_names{iter} ';' ])
end

	% All Endogenous Variables
for ii = 1:NumberOfEndogenous
	varname = M_.endo_names{ii};
  ys(ii)=eval(varname);
	%eval(['ys(' int2str(ii) ') = ' varname ';'])
end

%save_params

end %mcc::: 

%-------------------------------------------------------------------------%
%	Numerical solution functions
%-------------------------------------------------------------------------%

function F = SS_num_banks(x,phi,R_D,rho,PD)

	% Guesses
om_bar	= x(1);
sigma	= x(2);

	% BGG Auxiliaries
aux	= (log(om_bar)+0.5*sigma^2)/sigma;
Gam	= 1-normcdf(aux-sigma)-om_bar*(1-normcdf(aux));

	% Target equations
F =[om_bar	- Gam*(1-phi)/phi*R_D/rho;
	PD		- normcdf(aux)];

end

%-------------------------------------------------------------------------%

function F = SS_num_banks_h(x,phi,R_BB,rho,PD,ppi)

	% Guesses
om_bar	= x(1);
sigma	= x(2);

	% BGG Auxiliaries
aux	= (log(om_bar)+0.5*sigma^2)/sigma;
Gam	= 1-normcdf(aux-sigma)-om_bar*(1-normcdf(aux));

	% Target equations
F =[om_bar*phi*rho	- Gam*(1-phi)*R_BB*ppi;
	PD		- normcdf(aux)];
end

%-------------------------------------------------------------------------%

function F = SS_num_e(x,MU,R_f_tilde,R_L,CHI_e,a,ppi)
	% Guesses
om_bar	= x(1);
sigma	= x(2);
	% BGG Auxiliaries
aux			= (log(om_bar)+0.5*sigma^2)/sigma;
Gam			= 1-normcdf(aux-sigma)-om_bar*(1-normcdf(aux));
Gam_G		= (1-MU)*normcdf(aux-sigma)+om_bar*(1-normcdf(aux));
Gam_der		= 1-normcdf(aux);
Gam_G_der	= 1-normcdf(aux)-MU*normpdf(aux)/sigma;
mu_G		= 1-Gam-Gam_G;
	% Target equations
F =[R_f_tilde			- Gam_G*R_L/om_bar;
	Gam_G_der/Gam_der	- (1-CHI_e)*R_f_tilde/a/ppi];
end


%-------------------------------------------------------------------------%
                        
function F = SS_num_i(x,MU_i,a,SIGMA,ppi,R_h_tilde,BETA_i,R_i)
  % Guesses
om_bar		= exp(x(1))/(1+exp(x(1)));
sigma	 	= exp(x(2))/(1+exp(x(2)));

  % BGG Auxiliaries
aux			= (log(om_bar)+0.5*sigma^2)/sigma;
Gam_G_i		= (1-MU_i)*normcdf(aux-sigma)+om_bar*(1-normcdf(aux));
Gam_der_i	= 1-normcdf(aux);
Gam_G_der_i	= 1-normcdf(aux)-MU_i*normpdf(aux)/sigma;
PD_i		= normcdf(aux);

  % Target equations
F =[R_h_tilde -	Gam_G_i*R_i*ppi/om_bar;
   ppi*a^(SIGMA)*Gam_G_der_i - Gam_der_i*BETA_i*R_h_tilde];
  %1000*(PD_i - 0.0112/4)];
end

function F = SS_num_i_eval(x,MU_i,a,SIGMA,ppi,R_h_tilde,BETA_i,R_i)
  % Guesses
om_bar=x(1);
sigma=x(2);
  % BGG Auxiliaries
aux			= (log(om_bar)+0.5*sigma^2)/sigma;
Gam_G_i		= (1-MU_i)*normcdf(aux-sigma)+om_bar*(1-normcdf(aux));
Gam_der_i	= 1-normcdf(aux);
Gam_G_der_i	= 1-normcdf(aux)-MU_i*normpdf(aux)/sigma;
PD_i		= normcdf(aux);

  % Target equations
F =[R_h_tilde -	Gam_G_i*R_i*ppi/om_bar;
  ppi*a^(-SIGMA) - Gam_der_i/Gam_G_der_i*BETA_i*R_h_tilde];
end

%-------------------------------------------------------------------------%

function F = SS_num_l_h(l_h,R_i,om_bar_i,R_h,q_h,phi_h,d_f,a,ppi,GAMMA_d,PD_d,R_D,...
    mu_G_e,R_e,q_k,k,mu_G_i,mu_G_h,R_h_tilde,mu_G_f,R_f_tilde,...
    l_f,p_H,y_H,p_F,ETA,rer,xi_m,Xi_F,OMEGA,s_co,s_tb,h,w,n,...
    Gam_i,PHI_c,SIGMA,PHI_HH,BETA_i,Gam_G_i,BETA_up,DELTA_h,i,i_h,...
	s_g,xi_h,ETA_CHAT,OMEGA_UP,BETA_rp,q_l,q_BB,R_BB,GAMMA_bh,PD_h,...
	ALPHA_BLG, ALPHA_BSG, q_BL,s_bast, OMEGA_BL,d_u,R_BL)
	
l_h 	= max(0.001,l_h);
	
	% Updating equations
h_i			= R_i*q_l*l_h*ppi/(om_bar_i*R_h*q_h);
bb_tot 		= (1-phi_h)*q_l*l_h/q_BB;
nu			= 1/(a*ppi)*(GAMMA_d*PD_d*R_D*d_f +GAMMA_bh*PD_h*R_BB*q_BB*bb_tot +mu_G_e*R_e*q_k*k...
			+mu_G_i*R_h*q_h*h_i + mu_G_h*R_h_tilde*l_h*q_l + mu_G_f*R_f_tilde*l_f) ;
gdpn		= (p_H*y_H+p_F^(-ETA)*(p_F-rer*xi_m*Xi_F)*(1-OMEGA)*nu-nu)/...
  (1-s_co-(1-s_tb)*p_F^(-ETA)*(p_F-rer*xi_m*Xi_F)*(1-OMEGA));
h_p			= max(0,h-h_i);
c_i			= w*n/2+q_h*h_i*(Gam_i*R_h/(a*ppi)-1)+q_l*l_h;
O_CHAT		= (xi_h^ETA_CHAT*a^(-SIGMA*ETA_CHAT)*(a*c_i*(1-PHI_c/a)/(xi_h*(h_i*(1-PHI_HH/a))))...
			*(1/BETA_i*(q_h-(Gam_G_i)*R_h*q_h/R_h_tilde)-(a^-SIGMA*Gam_i*R_h)/ppi*q_h)^(-ETA_CHAT) + 1)^-1;
  % Bond positions 
bl_g 		= ALPHA_BLG*gdpn/q_BL;
bs_g    	= ALPHA_BSG*gdpn;
bl_priv 	= - bl_g;
bs_priv 	= - bs_g;
bs_u 		= bs_priv/OMEGA_UP;
bb_u 		= bb_tot/OMEGA_UP;
b_ast_Tot 	= s_bast*gdpn/rer;
b_ast_u 	= b_ast_Tot/OMEGA_UP;
bl_u 		= (OMEGA_BL*(bs_u + rer*b_ast_u + d_u ) - bb_u*q_BB)/q_BL;
bl_r 		= (bl_priv - OMEGA_UP*bl_u)/(1-OMEGA_UP); 

  % Find h_r using eq 14 instead of upsilon 
corcho 		= a^(SIGMA*ETA_CHAT-1)*xi_h^(1-ETA_CHAT)*(q_h/BETA_rp-(1-DELTA_h)*a^(-SIGMA)*q_h)^ETA_CHAT*(1-O_CHAT)/O_CHAT*(1-PHI_HH/a)/(1-PHI_c/a);
h_r 		= (q_BL*bl_r*(R_BL/a-1)+w*n/2)/(q_h-q_h/a*(1-DELTA_h)+corcho);
%c_r			= ((a^SIGMA)*(xi_h*h_r*(1-PHI_HH/a)*(1-O_CHAT)/((1-PHI_c/a)*O_CHAT*a))^(1/ETA_CHAT)...
%				*(q_h/BETA_rp-(1-DELTA_h)*a^-SIGMA*q_h))^ETA_CHAT;     
c_r         = corcho*h_r;
h_u   		= (h_p-(1-OMEGA_UP)*h_r)/OMEGA_UP;
c_u   		= ((a^SIGMA/BETA_up-(1-DELTA_h))*q_h/xi_h)^ETA_CHAT*(1-O_CHAT)/a/O_CHAT*h_u*xi_h*(1-PHI_HH/a)/(1-PHI_c/a);    
c_p			= OMEGA_UP*c_u + (1-OMEGA_UP)*c_r;
	
	% Target equation
F	= gdpn - (c_p+c_i+i+i_h+s_g*gdpn+s_tb*gdpn);
end


%-------------------------------------------------------------------------%

function F = SS_num_NF(x,h,q_h,DELTA_h,a,w,n_i,l_h,R_i,ppi,y_C,i,i_h,g,BETA_up,...
					   BETA_i,SIGMA,OMEGA_UP,ETA_CHAT,q_l)

% Guess
c_i = x;

% Target equation
%F = (h*q_h*(1-(1-DELTA_h)/a)-w*n/2-(l_h)*(1-R_i/(a*ppi))+c_i)/(y_C-c_i-i-i_h-g)...
%    -((1-BETA_i*(1-DELTA_h)*a^-SIGMA)/(1-BETA_up*(1-DELTA_h)*a^-SIGMA))^ETA_CHAT*(BETA_up/BETA_i)^ETA_CHAT*((w*n/2+(l_h)*(1-R_i/(a*ppi))-c_i)/c_i);
	
F = (h*q_h*(1-(1-DELTA_h)/a)-w*n_i-(q_l*l_h)*(1-R_i/(a))+c_i)/(y_C-c_i-i-i_h-g)...
    -((1-BETA_i*(1-DELTA_h)*a^-SIGMA)/(1-BETA_up*(1-DELTA_h)*a^-SIGMA))^ETA_CHAT*(BETA_up/BETA_i)^ETA_CHAT*((w*n_i+(q_l*l_h)*(1-R_i/(a))-c_i)/c_i);
	
end


%-------------------------------------------------------------------------%

function F = SS_num_NF_NH(x,w,n,y_C,i,g,l_h,R_i,ppi,a,VARPHI,SIGMA)

	% Guess
c_p = x;

	% Target equation
F = (w*n/(y_C-c_p-i-g+l_h*(R_i/(ppi*a)-1))-1)^VARPHI...
	-((y_C-i-g)/c_p-1)^SIGMA;

end

%-------------------------------------------------------------------------%














