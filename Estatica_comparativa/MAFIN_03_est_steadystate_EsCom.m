%-------------------------------------------------------------------------%
%	Steady State Equations - Financial Frictions
%-------------------------------------------------------------------------%
%
%     clc
%     clear
%   close all
%   tic
%   run Parametros_y_Exogenas;


phi_f		=	capreq*1+0.043	  ;
phi_h		=	capreq*0.6+0.043 ;

ppi			= ppi_T;
R			= ppi*a^SIGMA/BETA_up;
R_D_tilde	= R;
ppi_S		= ppi/ppi_ast;
R_ast		= R/ppi_S; 
%R_W			= R_ast/xi_R;       %se toma como exogena *no entra en nada
%ppi_ast =ppi/ppi_S                 %exogena, se ecuentra ppi_S 


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
n_tilde = n; % n tiene que ser endogena fija n=0.3.
rho_f		= a*ppi/(1-CHI_b);
rho_h		= rho_f;
rho_tilde_h = rho_h;
%PD_d		= (1-R_D_tilde/R_D)/GAMMA_d; %Primero se encuentra R_D
%PD_h		= PD_d;
%PD_f		= PD_d;
%R_BL    	= Rnom_BL/ppi; se encuentra primero Rnom_BL
%BETA_rp     = a^SIGMA/R_BL; % parametro, se da vuelta y se encuentra R_BL 
R_BL		=a^SIGMA/BETA_rp  ;
Rnom_BL 	=R_BL/ppi;
R_BB_tilde 	= R_BL;

%R_BB		= R_BB_tilde/(1-GAMMA_bh*PD_h); %primero se encuentra PD_h
%R_i     	= Rnom_i/ppi; %primero se encuentra Rnom_i
%q_l     	=(R_i-KAPPA_LOAN)^(-1); 

q_BL    	= 1/(R_BL-KAPPA_B);
%q_BB		= 1/(R_BB-KAPPA_BH); % se encuentra más abajo


%Encontrar R_D (esto es solo de EstCom).
options		= optimoptions('fsolve','FunctionTolerance',1e-10,'OptimalityTolerance',1e-10,'MaxIterations',10000,'MaxFunctionEvaluations',10000,'Display','off');
x_guess=[1.5];
 [x,~]=fsolve(@SS_RD_PD,x_guess,options,phi_f,rho_f,sigma_ff,MU_f,GAMMA_d,R_D_tilde); 

R_D=x(1);

	% Bank numerical solutions
options		= optimoptions('fsolve','FunctionTolerance',1e-10,'OptimalityTolerance',1e-10,'MaxIterations',10000,'MaxFunctionEvaluations',10000,'Display','off');    
x_guess		= [0.9];
[x,~]		= fsolve(@SS_num_banks,x_guess,options,phi_f,R_D,rho_f,sigma_ff);
om_bar_f	= x(1);
%sigma_ff	= x(2); Es un parametro


	% Bank BGG auxiliaries and ROA
aux_f		= (log(om_bar_f)+0.5*sigma_ff^2)/sigma_ff;
Gam_f		= 1-normcdf(aux_f-sigma_ff)-om_bar_f*(1-normcdf(aux_f));
Gam_G_f		= (1-MU_f)*normcdf(aux_f-sigma_ff)+om_bar_f*(1-normcdf(aux_f));
Gam_der_f	= 1-normcdf(aux_f);
Gam_G_der_f	= 1-normcdf(aux_f)-MU_f*normpdf(aux_f)/sigma_ff;
mu_G_f		= 1-Gam_f-Gam_G_f;
PD_f		= normcdf(aux_f); %Solo EstCom
PD_d=PD_f; %Solo EstCom

options		= optimoptions('fsolve','FunctionTolerance',1e-10,'OptimalityTolerance',1e-10,'MaxIterations',10000,'MaxFunctionEvaluations',10000,'Display','off');
x_guess=[1.5];
 [x,~]=fsolve(@SS_R_BB,x_guess,options,phi_h,rho_h,ppi,MU_h,R_BB_tilde,GAMMA_bh,sigma_hh); 
R_BB=x(1);
q_BB		= 1/(R_BB-KAPPA_BH); % %Nuevo...verificar

x_guess		= [0.9052];
[x,~]		= fsolve(@SS_num_banks_h,x_guess,options,phi_h,R_BB,rho_h,ppi,sigma_hh);
om_bar_h	= x(1);
%sigma_hh	= x(2); % es un parametro

aux_h		= (log(om_bar_h)+0.5*sigma_hh^2)/sigma_hh;
Gam_h		= 1-normcdf(aux_h-sigma_hh)-om_bar_h*(1-normcdf(aux_h));
Gam_G_h		= (1-MU_h)*normcdf(aux_h-sigma_hh)+om_bar_h*(1-normcdf(aux_h));
Gam_der_h	= 1-normcdf(aux_h);
Gam_G_der_h	= 1-normcdf(aux_h)-MU_h*normpdf(aux_h)/sigma_hh;
mu_G_h		= 1-Gam_h-Gam_G_h;
PD_h		= normcdf(aux_h); %Solo EstCom



R_f_tilde	= rho_f*phi_f/Gam_f;
R_h_tilde	= rho_h*phi_h/Gam_h;

% Entrepreneur anumerical solutions
options		= optimoptions('fsolve','FunctionTolerance',1e-10,'OptimalityTolerance',1e-10,'MaxIterations',1000,'MaxFunctionEvaluations',1000,'Display','off');
x_guess		= [0.6];
[x,~]		= fsolve(@SS_num_e,x_guess,options,MU_e,R_f_tilde,CHI_e,a,ppi,sigma_ee);

om_bar_e	= x(1);
%sigma_ee	= x(2); % es un parametro
aux_e		= (log(om_bar_e)+0.5*sigma_ee^2)/sigma_ee;
Gam_e		= 1-normcdf(aux_e-sigma_ee)-om_bar_e*(1-normcdf(aux_e));
Gam_G_e		= (1-MU_e)*normcdf(aux_e-sigma_ee)+om_bar_e*(1-normcdf(aux_e));
Gam_der_e	= 1-normcdf(aux_e);
Gam_G_der_e	= 1-normcdf(aux_e)-MU_e*normpdf(aux_e)/sigma_ee;
mu_G_e		= 1-Gam_e-Gam_G_e;
PD_e		= normcdf(aux_e);

R_L         = R_f_tilde*om_bar_e/Gam_G_e ; %Nuevo...Revisar.


% Impatient household BGG auxiliaries and PD

options		= optimoptions('fsolve','FunctionTolerance',1e-10,'OptimalityTolerance',1e-10,'MaxIterations',1000,'MaxFunctionEvaluations',1000,'Display','off');
x_guess		= [0.945379];
[x,~]		= fsolve(@SS_num_i,x_guess,options,MU_i,R_h_tilde,BETA_i,a,SIGMA,ppi,sigma_ii);
om_bar_i 	= x(1);
%om_bar_i 	= exp(x(1))/(1+exp(x(1)));
%sigma_ii 	= exp(x(2))/(1+exp(x(2))); %es parametro

aux_i		= (log(om_bar_i)+0.5*sigma_ii^2)/sigma_ii;
Gam_i		= 1-normcdf(aux_i-sigma_ii)-om_bar_i*(1-normcdf(aux_i));
Gam_G_i		= (1-MU_i)*normcdf(aux_i-sigma_ii)+om_bar_i*(1-normcdf(aux_i));
Gam_der_i	= 1-normcdf(aux_i);
Gam_G_der_i	= 1-normcdf(aux_i)-MU_i*normpdf(aux_i)/sigma_ii;
mu_G_i		= 1-Gam_i-Gam_G_i;
PD_i		= normcdf(aux_i);

R_i =(R_h_tilde*om_bar_i)/(Gam_G_i*ppi); % nuevo! verificar!
q_l     	=(R_i-KAPPA_LOAN)^(-1); % nuevo! verificar!
Rnom_i=R_i*ppi;

R_e			= R_f_tilde*a*ppi/(a*ppi*Gam_G_e+Gam_e*(1-CHI_e)*R_f_tilde);
r_k			= q_k*(R_e/ppi-(1-DELTA_k));



% encontrar numericamente p_H
options		= optimoptions('fsolve','FunctionTolerance',1e-10,'OptimalityTolerance',1e-10,'MaxIterations',1000,'MaxFunctionEvaluations',1000,'Display','off');
x_guess=[1.5];
[x,~]=fsolve(@SS_phn,x_guess,options,om_bar_i,R_h,q_h,phi_h,a,ppi,GAMMA_d,PD_d,R_D,...
						mu_G_e,R_e,q_k,mu_G_i,mu_G_h,R_h_tilde,mu_G_f,R_f_tilde,ETA,...
						xi_m,Xi_F,OMEGA,s_co,s_tb,Gam_i,PHI_c,SIGMA,PHI_HH,BETA_i,Gam_G_i,DELTA_h,s_g,...
                        ALPHA,mc_H,z,r_k,DELTA_k,xi_i,Xim_h,Gam_e,CHI_e,Gam_der_e,Gam_f,Gam_G_der_e,om_bar_e,...
                        phi_f,q_l,r_h_k,N_H,xi_ih,VARPHI_H_0,RHO_VARPHI_H,mc_F,CHI_b,s_bast,y_ast,...
                        ETA_ast,xi_h,ETA_CHAT,GAMMA_bh,PD_h,R_BB,q_BB,OMEGA_UP,BETA_up,BETA_rp,n,n_tilde,...
                        xi_CHI_e,R_i,ALPHA_BLG,ALPHA_BSG,q_BL,OMEGA_BL,R_BL,xi_CHI_b); 
p_H=x(1);
                  


p_Z         =p_H*mc_H;
mc_Z        = p_Z;
w			= (ALPHA^ALPHA*(1-ALPHA)^(1-ALPHA)*mc_Z*z/r_k^ALPHA)^(1/(1-ALPHA));
k			= ALPHA/(1-ALPHA)*n_tilde*w/r_k*a;
y_Z         = z*(k/a)^ALPHA*n_tilde^(1-ALPHA);
x_Z         = y_Z;
i			= k*(1-(1-DELTA_k)/a)/xi_i;
y_H			= x_Z/Xim_h;
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
l_h_guess	= 0.04;
[l_h,~]		= fsolve(@SS_num_l_h,l_h_guess,options,R_i,om_bar_i,R_h,q_h,phi_h,d_f,a,ppi,GAMMA_d,PD_d,R_D,...
    mu_G_e,R_e,q_k,k,mu_G_i,mu_G_h,R_h_tilde,mu_G_f,R_f_tilde,...
    l_f,p_H,y_H,p_F,ETA,rer,xi_m,Xi_F,OMEGA,s_co,s_tb,h,w,n,...
    Gam_i,PHI_c,SIGMA,PHI_HH,BETA_i,Gam_G_i,BETA_up,DELTA_h,i,i_h,...
	s_g,xi_h,ETA_CHAT,OMEGA_UP,BETA_rp,q_l,q_BB,R_BB,GAMMA_bh,PD_h, ALPHA_BLG, ...
  ALPHA_BSG, q_BL,s_bast,OMEGA_BL,d_u,R_BL);	% Third block of analytical solutions

h_i  		= R_i*q_l*l_h*ppi/(om_bar_i*R_h*q_h);




bb_tot 		= (1-phi_h)*q_l*l_h/q_BB;
e_h			= phi_h*q_l*l_h;
n_b			= e_f+e_h;
psi_b		= n_b/(1-CHI_b*xi_CHI_b);
c_b			= CHI_b*xi_CHI_b*psi_b;

nu			= 1/(a*ppi)*(GAMMA_d*PD_d*R_D*d_f +GAMMA_bh*PD_h*R_BB*q_BB*bb_tot +mu_G_e*R_e*q_k*k...
				+mu_G_i*R_h*q_h*h_i + mu_G_h*R_h_tilde*l_h*q_l + mu_G_f*R_f_tilde*l_f) ;
gdpn		= (p_H*y_H+p_F^(-ETA)*(p_F-rer*xi_m*Xi_F)*(1-OMEGA)*nu-nu)/...
				(1-s_co-(1-s_tb)*p_F^(-ETA)*(p_F-rer*xi_m*Xi_F)*(1-OMEGA));
tb			= s_tb*gdpn;

%g			= s_g*gdpn;   %g es exogeno, se deja fijo y se encuentra s_g
%s_g         = g/gdpn;

%y_co		= s_co*gdpn/(p_co*rer);  %y_co es exogena, se encuentra s_co
%s_co        = (y_co*p_co*rer)/gdpn;

b_ast_Tot = s_bast*gdpn/rer;
y_C			= gdpn+nu-tb;
x_F			= (1-OMEGA)*p_F^-ETA*y_C;
x_H			= OMEGA*p_H^-ETA*y_C;
x_H_ast		= y_H-x_H;
%y_ast		= x_H_ast*(p_H/rer)^ETA_ast; % exogena. se usa  en la funcion de p_H

y_F			= x_F;
imp			= y_F*Xi_F;
h_p			= max(0,h-h_i);
c_i			= w*n/2+q_h*h_i*(Gam_i*R_h/(a*ppi)-1)+q_l*l_h;
O_CHAT		= (xi_h^ETA_CHAT*a^(-SIGMA*ETA_CHAT)*(a*c_i*(1-PHI_c/a)/(xi_h*(h_i*(1-PHI_HH/a))))...
				*(1/BETA_i*(q_h-(Gam_G_i)*R_h*q_h/R_h_tilde)-(a^-SIGMA*Gam_i*R_h)/ppi*q_h)^(-ETA_CHAT) + 1)^-1;
c_hat_i		= ((1-O_CHAT)^(1/ETA_CHAT)*(c_i*(1-PHI_c/a))^((ETA_CHAT-1)/ETA_CHAT)...
				+(O_CHAT)^(1/ETA_CHAT)*(xi_h*h_i/a*(1-PHI_HH/a))^((ETA_CHAT-1)/ETA_CHAT))^(ETA_CHAT/(ETA_CHAT-1));
la_i		= c_hat_i^-SIGMA*((1-O_CHAT)*c_hat_i/(c_i*(1-PHI_c/a)))^(1/ETA_CHAT);    

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
%bl_cb 		= 1;

  % Find h_r using eq 14 instead of upsilon
corcho 		= a^(SIGMA*ETA_CHAT-1)*xi_h^(1-ETA_CHAT)*(q_h/BETA_rp-(1-DELTA_h)*a^(-SIGMA)*q_h)^ETA_CHAT*(1-O_CHAT)/O_CHAT*(1-PHI_HH/a)/(1-PHI_c/a);
h_r 		= (q_BL*bl_r*(R_BL/a-1)+w*n/2)/(q_h-q_h/a*(1-DELTA_h)+corcho);
c_r		= ((a^SIGMA)*(xi_h*h_r*(1-PHI_HH/a)*(1-O_CHAT)/((1-PHI_c/a)*O_CHAT*a))^(1/ETA_CHAT)...
					*(q_h/BETA_rp-(1-DELTA_h)*a^-SIGMA*q_h))^ETA_CHAT;   
c_hat_r	= ((1-O_CHAT)^(1/ETA_CHAT)*(c_r*(1-PHI_c/a))^((ETA_CHAT-1)/ETA_CHAT)...
					+(O_CHAT)^(1/ETA_CHAT)*(xi_h*h_r/a*(1-PHI_HH/a))^((ETA_CHAT-1)/ETA_CHAT))^(ETA_CHAT/(ETA_CHAT-1));
la_r		= c_hat_r^(-SIGMA)*((1-O_CHAT)*c_hat_r/(c_r*(1-PHI_c/a)))^(1/ETA_CHAT);
h_u  		= (h_p-(1-OMEGA_UP)*h_r)/OMEGA_UP;
c_u   		= ((a^SIGMA/BETA_up-(1-DELTA_h))*q_h/xi_h)^ETA_CHAT*(1-O_CHAT)/a/O_CHAT*h_u*xi_h*(1-PHI_HH/a)/(1-PHI_c/a);
c_hat_u	= ((1-O_CHAT)^(1/ETA_CHAT)*(c_u*(1-PHI_c/a))^((ETA_CHAT-1)/ETA_CHAT)...
					+(O_CHAT)^(1/ETA_CHAT)*(xi_h*h_u/a*(1-PHI_HH/a))^((ETA_CHAT-1)/ETA_CHAT))^(ETA_CHAT/(ETA_CHAT-1));   %%%%NUEVO.0121
la_u		= c_hat_u^-SIGMA*((1-O_CHAT)*c_hat_u/(c_u*(1-PHI_c/a)))^(1/ETA_CHAT); 
la_p		= OMEGA_UP*la_u + (1-OMEGA_UP)*la_r;
c_p			= OMEGA_UP*c_u + (1-OMEGA_UP)*c_r; 
c			= c_i+c_p;
n_p			= n/2;
n_i			= n_p;
n_u			= n_p;
n_r			= n_u;
Xi_Tilde_up	= (c_hat_u)^(SIGMA);
Theta_up		= Xi_Tilde_up*(c_hat_u)^(-SIGMA);
Xi_Tilde_i	= (c_hat_i)^(SIGMA);
Theta_i		= Xi_Tilde_i*(c_hat_i)^(-SIGMA);
Xi_Tilde_rp	= (c_hat_r)^(SIGMA);
Theta_rp	= Xi_Tilde_rp*(c_hat_r)^(-SIGMA); 
Theta		= ((OMEGA_UP*Theta_up+(1-OMEGA_UP)*Theta_rp)+Theta_i)/2; %%%%NUEVO.0121
la_W		= (la_p+la_i)/2;
%xi_n		= mc_W*la_W*w/(Theta*n_tilde^VARPHI); %exogena
gdp			= c_i+c_p+i+i_h+g+x_H_ast+y_co-imp;
ren_ast		= b_ast_Tot*(1-R_ast/a/ppi)-tb/rer+(1-CHI)*p_co*y_co;
%epsilon_L	= (R_BL*BETA_up*a^(-SIGMA)-1)/((q_BL*bl_u+q_BB*bb_u)/(bs_u+d_u+rer*b_ast_u))^ETA_ZETA_L; % exogena
zeta_L		= epsilon_L*((q_BL*bl_u+q_BB*bb_u)/(bs_u+d_u+rer*b_ast_u))^ETA_ZETA_L;
tau 		= g + GAMMA_d*PD_d*R_D*d_f/a/ppi+ GAMMA_bh*PD_h*R_BB*q_BB*bb_tot/a...
				-(R/(a*ppi)-1)*bs_g - (R_BL/(a)-1)*q_BL*bl_g - CHI*rer*p_co*y_co;
%ALPHA_TAU   = tau/gdpn; %exogena

f_H			= p_tilde_H^-EPSILON_H*y_H*mc_H/(1-THETA_H*BETA_up*a^(1-SIGMA));
f_F			= p_tilde_F^-EPSILON_F*y_F*mc_F/(1-THETA_F*BETA_up*a^(1-SIGMA));
f_W			= w_tilde^(-EPSILON_W*(1+VARPHI))*mc_W*n_tilde/(1-((OMEGA_UP*BETA_up+(1-OMEGA_UP)*BETA_rp)+BETA_i)/2*THETA_W*a^(1-SIGMA));




%-------------------------------------------------------------------------%
%	Steady State Equations - Definitions and Observables
%-------------------------------------------------------------------------%

	% Additional definitions

c_p_i		= c_p/c_i;
i_agg		= i+i_h;
l			= l_f+l_h*q_l;
spread_RD_R	= R_D-R;
spread_RL_RD= R_L-R_D;
spread_Ri_RD= R_i*ppi-R_D;
spread_RBL_R= R_BL*ppi-R;
lev_e		= R_L*l_f/(q_k*k);
nu_gdp		= nu/gdp;
lab_income  =n*w;
gdpR=(gdp-y_co);

% disp('done!')
% toc







