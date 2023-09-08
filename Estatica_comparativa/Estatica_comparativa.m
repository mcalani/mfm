%-------------------------------------------------------------------------%
%	Steady State Equations - Financial Frictions
%-------------------------------------------------------------------------%

%% Comentar este bloque antes de correr el loop
% clc
% clear all
% format longG
% display_pars = 1;
% run Load_params;
% ccyb        = 0.00;

%% Bloque 1: sin soluciones numericas 
options		= optimoptions('fsolve','FunctionTolerance',1e-10,'OptimalityTolerance'...
    ,1e-30,'MaxIterations',1000,'MaxFunctionEvaluations',1000,'Display','off');
phi_f       = (phi_f_base+ccyb) + 0.0433;
phi_h       = (phi_h_base+ccyb)*0.6 + 0.0433;
ppi         = ppi_T;
R       	= ppi*a^SIGMA/BETA_up;
R_D_tilde	= R;
ppi_S		= ppi/ppi_ast; 
R_ast		= R/ppi_S;
R_W			= R_ast/xi_R;
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
Dl          = a;
q_h			= a^(N_H*SIGMA)*VARPHI_H_0/(BETA_up^N_H*xi_ih)*...
				(1-(BETA_up*RHO_VARPHI_H/a^SIGMA)^(N_H+1))/(1-BETA_up*RHO_VARPHI_H/a^SIGMA);
R_h			= ppi*(1-DELTA_h);
n_tilde		= n;
rho_f		= a*ppi/(1-CHI_b);
rho_h		= rho_f;
rho_tilde_h = rho_h;
R_BL        = a^SIGMA/BETA_rp;
Rnom_BL     = R_BL*ppi ;
R_BB_tilde 	= R_BL;
q_BL    	= 1/(R_BL-KAPPA_B);


%% Bloque 2: Num sol. om_bar_f & R_D
x_guess		= [0.9,0.05];
[x,~]		= fsolve(@EC_num_banks,x_guess,options,phi_f,sigma_ff,rho_f,GAMMA_d,R_D_tilde);
om_bar_f	= x(1);
R_D         = x(2);
PD_f = (1-R_D_tilde/R_D)/GAMMA_d; 

    % Bank BGG auxiliaries and ROA
aux_f		= (log(om_bar_f)+0.5*sigma_ff^2)/sigma_ff;
Gam_f		= 1-normcdf(aux_f-sigma_ff)-om_bar_f*(1-normcdf(aux_f));
Gam_G_f		= (1-MU_f)*normcdf(aux_f-sigma_ff)+om_bar_f*(1-normcdf(aux_f));
Gam_der_f	= 1-normcdf(aux_f);
Gam_G_der_f	= 1-normcdf(aux_f)-MU_f*normpdf(aux_f)/sigma_ff;
mu_G_f		= 1-Gam_f-Gam_G_f;
R_f_tilde	= rho_f*phi_f/Gam_f;

 
%% Bloque 2: Num sol. om_bar_h and R_BB 
% Bank numerical solutions

x_guess		= [0.904860304388942,1.00874615637069];
[x,~]		= fsolve(@EC_num_banks_H,x_guess,options,phi_h,rho_h,ppi,sigma_hh,GAMMA_bh,R_BB_tilde);
om_bar_h	= x(1);
R_BB	    = x(2);

aux_h		= (log(om_bar_h)+0.5*sigma_hh^2)/sigma_hh;
Gam_h		= 1-normcdf(aux_h-sigma_hh)-om_bar_h*(1-normcdf(aux_h));
Gam_G_h		= (1-MU_h)*normcdf(aux_h-sigma_hh)+om_bar_h*(1-normcdf(aux_h));
Gam_der_h	= 1-normcdf(aux_h);
Gam_G_der_h	= 1-normcdf(aux_h)-MU_h*normpdf(aux_h)/sigma_hh;
mu_G_h		= 1-Gam_h-Gam_G_h;

R_h_tilde	= rho_h*phi_h/Gam_h;
PD_h		= normcdf(aux_h);
q_BB		= 1/(R_BB-KAPPA_BH);


%% Entrepreneur anumerical solutions
x_guess		= 0.7;
[x,~]		= fsolve(@SS_num_e,x_guess,options,MU_e,R_f_tilde,CHI_e,a,ppi,sigma_ee);
om_bar_e	= x(1);  
aux_e		= (log(om_bar_e)+0.5*sigma_ee^2)/sigma_ee;
Gam_e		= 1-normcdf(aux_e-sigma_ee)-om_bar_e*(1-normcdf(aux_e));
Gam_G_e		= (1-MU_e)*normcdf(aux_e-sigma_ee)+om_bar_e*(1-normcdf(aux_e));
Gam_der_e	= 1-normcdf(aux_e);
Gam_G_der_e	= 1-normcdf(aux_e)-MU_e*normpdf(aux_e)/sigma_ee;
mu_G_e		= 1-Gam_e-Gam_G_e;
    
    % Related variables 
PD_e		= normcdf(aux_e);
R_L         = R_f_tilde*om_bar_e/Gam_G_e;
la_e        = Gam_der_e/(Gam_f*Gam_G_der_e);
R_e         = rho_f*phi_f/(Gam_e/la_e+Gam_f*Gam_G_e); 
r_k			= q_k*(R_e/ppi-(1-DELTA_k));

%% Impatient household BGG auxiliaries and PD
x_guess		= 0.8 ; % [0.960729732994068]; %This guess must be close to 1 !!
[x,fval]    = fsolve(@SS_num_i3,x_guess,options,MU_i,R_h_tilde,BETA_i,a,SIGMA,ppi,sigma_ii);
om_bar_i 	= x(1);     %om_bar_i    = 0.960729732994068;

aux_i		= (log(om_bar_i)+0.5*sigma_ii^2)/sigma_ii;
Gam_i		= 1-normcdf(aux_i-sigma_ii)-om_bar_i*(1-normcdf(aux_i));
Gam_G_i		= (1-MU_i)*normcdf(aux_i-sigma_ii)+om_bar_i*(1-normcdf(aux_i));
Gam_der_i	= 1-normcdf(aux_i);
Gam_G_der_i	= 1-normcdf(aux_i)-MU_i*normpdf(aux_i)/sigma_ii;
mu_G_i		= 1-Gam_i-Gam_G_i;
PD_i		= normcdf(aux_i);

R_hat_i     = (R_h_tilde*om_bar_i)/(Gam_G_i*ppi);
R_i         = R_hat_i ;
q_hat_l     =(R_hat_i-KAPPA_LOAN)^(-1);                
q_l     	=q_hat_l; 
Rnom_i      = R_i*ppi ; 

%% SE NORMALIZA A PH=1 COMO EN SS03

p_H = 1 ;
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
l_f			= q_k*k-n_e;
e_f			= phi_f*l_f;
d_f			= l_f-e_f;
d_u     	= d_f/OMEGA_UP;
p_F			= ((1-OMEGA*p_H^(1-ETA))/(1-OMEGA))^(1/(1-ETA));
rer			= mc_F*p_F/xi_m;

%% Numerical solution for l_h
x_guess	= 0.05; 
[x,~] = fsolve(@SS_num_l_h,x_guess,options,R_i,om_bar_i,R_h,q_h,a,ppi,BETA_i,Gam_G_i,...
    R_h_tilde,w,n,Gam_i,PHI_c,SIGMA,PHI_HH,xi_h,ETA_CHAT,q_l,O_CHAT) ;
l_h = x; 

%% Third block of closed form solutions  
h_i  		= R_i*q_l*l_h*ppi/(om_bar_i*R_h*q_h);
c_i			= w*n/2+q_h*h_i*(Gam_i*R_h/(a*ppi)-1)+q_l*l_h;
c_hat_i		= ((1-O_CHAT)^(1/ETA_CHAT)*(c_i*(1-PHI_c/a))^((ETA_CHAT-1)/ETA_CHAT)...
				+(O_CHAT)^(1/ETA_CHAT)*(xi_h*h_i/a*(1-PHI_HH/a))^((ETA_CHAT-1)/ETA_CHAT))^(ETA_CHAT/(ETA_CHAT-1));
la_i		= c_hat_i^(-SIGMA)*((1-O_CHAT)*c_hat_i/(c_i*(1-PHI_c/a)))^(1/ETA_CHAT);            
la_h        = la_i/(rho_h*phi_h); 
e_h 		= phi_h*q_l*l_h;
n_b			= e_f+e_h;
psi_b		= n_b/(1-CHI_b*xi_CHI_b);
c_b			= CHI_b*xi_CHI_b*psi_b;
bb_tot      = (q_l*l_h-e_h)/q_BB;
PD_d 		= (q_BB*bb_tot*PD_h+d_f*PD_f)/(q_BB*bb_tot+d_f);

x_H_ast     = y_ast/(p_H/rer)^ETA_ast;  
x_H		    = y_H-x_H_ast;  
y_C         = x_H/(OMEGA*p_H^-ETA);	 
x_F			= (1-OMEGA)*p_F^-ETA*y_C; 
y_F			= x_F; 
imp			= y_F*Xi_F; 
tb          =  p_H*x_H_ast+rer*p_co*y_co-rer*xi_m*imp;  
nu			= 1/(a*ppi)*(GAMMA_d*PD_f*R_D*d_f +GAMMA_bh*PD_h*R_BB*q_BB*bb_tot +mu_G_e*R_e*q_k*k...
				+mu_G_i*R_h*q_h*h_i + mu_G_h*R_h_tilde*l_h*q_l + mu_G_f*R_f_tilde*l_f) ;
gdpn        = y_C-nu+tb; 

% Ratios by definition 
s_g         = g/gdpn; 
s_co        = (y_co*p_co*rer)/gdpn; 
s_tb        = tb/gdpn; 
bl_g 		= ALPHA_BLG*gdpn/q_BL ;
bs_g    	= ALPHA_BSG*gdpn ;
bl_priv 	= - bl_g ;
bs_priv 	= - bs_g ;
bs_u 		= bs_priv/OMEGA_UP ;
bb_u 		= bb_tot/OMEGA_UP ;
n_p			= n/2;
n_i			= n_p;
n_u			= n_p;
n_r			= n_u; 


    %% Solve for rhk and bast simultaneously  
x_guess	= [0.4 , -0.02] ;
[x,~]		= fsolve(@SS_sbast_rhk,x_guess,options, ...
    q_k, k, q_h, a, N_H, xi_ih, DELTA_h, VARPHI_H_0, RHO_VARPHI_H, h_i, gdpn, rer, ...
    OMEGA_UP, OMEGA_BL, bl_priv, SIGMA, ETA_CHAT, BETA_rp, BETA_up, O_CHAT,R_BL,i,g,tb,la_i, bs_u, d_u, bb_u, ...
    q_BB, q_BL, xi_h, xi_n, c_i, PHI_HH, PHI_c, w, n, c_hat_i, mc_W, n_tilde, VARPHI);

r_h_k   = x(1) ; 
s_bast  = x(2); 
clear x 
   %Depend only on rhk
h			= r_h_k*q_k*k/q_h;
i_ah		= h*a^N_H/xi_ih*(1-(1-DELTA_h)/a);
i_h			= i_ah*VARPHI_H_0*(1-(RHO_VARPHI_H/a)^(N_H+1))/(1-RHO_VARPHI_H/a);
h_p			= max(0,h-h_i);

    %Depend on bast 
b_ast_Tot   = s_bast*gdpn/rer;
b_ast_u 	= b_ast_Tot/OMEGA_UP;
bl_u        = (OMEGA_BL*(bs_u + rer*b_ast_u + d_u ) - bb_u*q_BB)/q_BL; 
bl_r 		= (bl_priv - OMEGA_UP*bl_u)/(1-OMEGA_UP); 
corcho_r	= a^(SIGMA*ETA_CHAT-1)*xi_h^(1-ETA_CHAT)*(q_h/BETA_rp-(1-DELTA_h)*...
                a^(-SIGMA)*q_h)^ETA_CHAT*(1-O_CHAT)/O_CHAT*(1-PHI_HH/a)/(1-PHI_c/a);
h_r 		= (q_BL*bl_r*(R_BL/a-1)+w*n/2)/(q_h-q_h/a*(1-DELTA_h)+corcho_r);
c_r         = corcho_r*h_r;
c_hat_r		= ((1-O_CHAT)^(1/ETA_CHAT)*(c_r*(1-PHI_c/a))^((ETA_CHAT-1)/ETA_CHAT)...
                +(O_CHAT)^(1/ETA_CHAT)*(xi_h*h_r/a*(1-PHI_HH/a))^((ETA_CHAT-1)/ETA_CHAT))^(ETA_CHAT/(ETA_CHAT-1));
la_r		= c_hat_r^(-SIGMA)*((1-O_CHAT)*c_hat_r/(c_r*(1-PHI_c/a)))^(1/ETA_CHAT);

    %Depend on both 
h_u  		= (h_p-(1-OMEGA_UP)*h_r)/OMEGA_UP;
c_u   		= ((a^SIGMA/BETA_up-(1-DELTA_h))*q_h/xi_h)^ETA_CHAT*...
                (1-O_CHAT)/a/O_CHAT*h_u*xi_h*(1-PHI_HH/a)/(1-PHI_c/a);    
c_hat_u		= ((1-O_CHAT)^(1/ETA_CHAT)*(c_u*(1-PHI_c/a))^((ETA_CHAT-1)/ETA_CHAT)...
                +(O_CHAT)^(1/ETA_CHAT)*(xi_h*h_u/a*(1-PHI_HH/a))^((ETA_CHAT-1)/ETA_CHAT))^(ETA_CHAT/(ETA_CHAT-1));
la_u		= c_hat_u^-SIGMA*((1-O_CHAT)*c_hat_u/(c_u*(1-PHI_c/a)))^(1/ETA_CHAT);
la_p		= OMEGA_UP*la_u + (1-OMEGA_UP)*la_r;
c_p			= OMEGA_UP*c_u + (1-OMEGA_UP)*c_r;
la_W		= (la_p+la_i)/2;
Xi_Tilde_up	= (c_hat_u)^(SIGMA);
Theta_up	= Xi_Tilde_up*(c_hat_u)^(-SIGMA);
Xi_Tilde_i	= (c_hat_i)^(SIGMA);
Theta_i		= Xi_Tilde_i*(c_hat_i)^(-SIGMA);	
Xi_Tilde_rp	= (c_hat_r)^(SIGMA);
Theta_rp	= Xi_Tilde_rp*(c_hat_r)^(-SIGMA);
Theta		= ((OMEGA_UP*Theta_up+(1-OMEGA_UP)*Theta_rp)+Theta_i)/2;
ren_ast     = b_ast_Tot*(1-R_ast/a/ppi) + (1-CHI)*p_co*y_co-tb/rer;

gdp			= c_i+c_p+i+i_h+g+x_H_ast+y_co-imp;
ren_ast		= b_ast_Tot*(1-R_ast/a/ppi_ast)-tb/rer+(1-CHI)*p_co*y_co;
epsilon_L	= (R_BL*BETA_up*a^(-SIGMA)-1);
zeta_L		= epsilon_L;
tau 		= g + GAMMA_d*PD_f*R_D*d_f/a/ppi+ GAMMA_bh*PD_h*R_BB*q_BB*bb_tot/a...
				-(R/(a*ppi)-1)*bs_g - (R_BL/(a)-1)*q_BL*bl_g - CHI*rer*p_co*y_co;
ALPHA_TAU   = tau/gdpn;
f_H			= p_tilde_H^-EPSILON_H*y_H*mc_H/(1-THETA_H*BETA_up*a^(1-SIGMA));
f_F			= p_tilde_F^-EPSILON_F*y_F*mc_F/(1-THETA_F*BETA_up*a^(1-SIGMA));
f_W			= w_tilde^(-EPSILON_W*(1+VARPHI))*mc_W*n_tilde/(1-((OMEGA_UP*BETA_up+(1-OMEGA_UP)*BETA_rp)+BETA_i)/2*THETA_W*a^(1-SIGMA));

%%	% Additional definitions
c			= c_i+c_p;
c_p_i		= c_p/c_i;
i_agg		= i+i_h;
l			= l_f+l_h*q_l;
spread_RD_R	= R_D-R;
spread_RL_RD= R_L-R_D;
spread_RL_R = R_L-R;
spread_Ri_RD= R_i*ppi-R_D;
spread_RBL_R= R_BL*ppi-R;
spread_Ri_RBL = R_i - R_BL;
lev_e		= R_L*l_f/(q_k*k);
nu_gdp		= nu/gdp;
lab_income  = w*n;
gdpR        = (gdp-y_co);
mort_tot    = q_l*l_h ;
mortVpar_tot = q_hat_l*l_h;
brechacred = 0 ; 
ltot = q_l*l_h+l_f ; 
lev_tot = n_b / (l_f + 0.6*mortVpar_tot);
lev_f = e_f / (l_f);
lev_h = e_h / (mortVpar_tot);


%% Functions
function F = EC_num_banks(x,phi,sigma,rho,GAMMA_d, R_D_tilde)
    	% Guesses
    om_bar	= x(1);
    R_D	= x(2);
       % BGG Auxiliaries
    aux	= (log(om_bar)+0.5*sigma^2)/sigma;
    Gam	= 1-normcdf(aux-sigma)-om_bar*(1-normcdf(aux));
    PD	= normcdf(aux);
    	% Target equations
    F =[-rho + Gam*(1-phi)/phi*R_D/om_bar;
    -R_D + R_D_tilde / (1-GAMMA_d*PD)];
end

function F = EC_num_banks_H(x,phi,rho,ppi,sigma,GAMMA_bh,R_BB_tilde)
        % Guesses
    om_bar	= x(1);
    R_BB = x(2); 
        % BGG Auxiliaries
    aux	= (log(om_bar)+0.5*sigma^2)/sigma;
    Gam	= 1-normcdf(aux-sigma)-om_bar*(1-normcdf(aux));
    PD	= normcdf(aux);
        % Target equations
    F =[ -om_bar*phi*rho	+ Gam*(1-phi)*R_BB*ppi;
         R_BB*(1-GAMMA_bh*PD)-R_BB_tilde];
end

function F = SS_num_e(x,MU,R_f_tilde,CHI_e,a,ppi,sigma)
    	% Guesses
    om_bar	= x(1);
        % BGG Auxiliaries
    aux			= (log(om_bar)+0.5*sigma^2)/sigma;
    Gam_der		= 1-normcdf(aux);
    Gam_G_der	= 1-normcdf(aux)-MU*normpdf(aux)/sigma;
        % Target equations
     F = [ Gam_G_der/Gam_der	- (1-CHI_e)*R_f_tilde/a/ppi];
end

function F = SS_num_i3(x,MU,R_h_tilde,BETA_i,a,SIGMA,ppi,sigma)
        % Guesses
    om_bar		= x(1);
        % BGG Auxiliaries
    aux			= (log(om_bar)+0.5*sigma^2)/sigma;
    Gam_G		= (1-MU)*normcdf(aux-sigma)+om_bar*(1-normcdf(aux));
    Gam_der		= 1-normcdf(aux);
    Gam_G_der	= 1-normcdf(aux)-MU*normpdf(aux)/sigma;
        % Target equations
    F =[Gam_G_der/Gam_der - BETA_i*R_h_tilde/(ppi*a^(SIGMA))];
end


function F = SS_num_l_h(l_h,R_i,om_bar_i,R_h,q_h,a,ppi,BETA_i,Gam_G_i,...
   R_h_tilde,w,n,Gam_i,PHI_c,SIGMA,PHI_HH,xi_h,ETA_CHAT,q_l,O_CHAT)
        %Guess
    l_h 	= max(0.001,l_h);
        % Updating equations
    h_i  		= R_i*q_l*l_h*ppi/(om_bar_i*R_h*q_h);
    c_i=(((O_CHAT^-1)-1)*(1/BETA_i*(q_h-(Gam_G_i)*R_h*q_h/R_h_tilde)-(a^-SIGMA*Gam_i*R_h)/ppi*q_h)^(ETA_CHAT)...
        *(xi_h*(h_i*(1-PHI_HH/a))))/(xi_h^ETA_CHAT*a^(-SIGMA*ETA_CHAT)*(1-PHI_c/a)*a) ;
        % Target equation
    F	= l_h-(c_i-w*n/2-q_h*h_i*(Gam_i*R_h/(a*ppi)-1))/q_l;
end


function root = bisect(f,a,b,tol)
    if nargin < 4
        tol = 1.5e-8;
    end
        % Uso: root = bisect(f,a,b)
    if f(a)*f(b)>0 
        error('La funcion debe ser de signo distinto en los extremos');
    else
        s = sign(f(a));
        x = (a+b)/2;
        d = (b-a)/2;

        while d>tol
            d = d/2;
            if s == sign(f(x))
                x = x+d;
            else
                x = x-d;
            end
        end
        root = x;
    end
end



function F = SS_sbast_rhk(x, q_k, k, q_h, a, N_H, xi_ih, DELTA_h, VARPHI_H_0, RHO_VARPHI_H, h_i, gdpn, rer, ...
    OMEGA_UP, OMEGA_BL, bl_priv, SIGMA, ETA_CHAT, BETA_rp, BETA_up, O_CHAT,R_BL,i,g,tb,la_i, bs_u, d_u, bb_u, ...
    q_BB, q_BL, xi_h, xi_n, c_i, PHI_HH, PHI_c, w, n, c_hat_i, mc_W, n_tilde, VARPHI)
        %Guess 
    r_h_k   = x(1) ; 
    s_bast  = x(2); 
        %Depend only on rhk
    h			= r_h_k*q_k*k/q_h;
    i_ah		= h*a^N_H/xi_ih*(1-(1-DELTA_h)/a);
    i_h			= i_ah*VARPHI_H_0*(1-(RHO_VARPHI_H/a)^(N_H+1))/(1-RHO_VARPHI_H/a);
    h_p			= max(0,h-h_i);
        %Depend on bast 
    b_ast_Tot   = s_bast*gdpn/rer;
    b_ast_u 	= b_ast_Tot/OMEGA_UP;
    bl_u        = (OMEGA_BL*(bs_u + rer*b_ast_u + d_u ) - bb_u*q_BB)/q_BL; 
    bl_r 		= (bl_priv - OMEGA_UP*bl_u)/(1-OMEGA_UP); 
    corcho_r	= a^(SIGMA*ETA_CHAT-1)*xi_h^(1-ETA_CHAT)*(q_h/BETA_rp-(1-DELTA_h)*...
                    a^(-SIGMA)*q_h)^ETA_CHAT*(1-O_CHAT)/O_CHAT*(1-PHI_HH/a)/(1-PHI_c/a);
    h_r 		= (q_BL*bl_r*(R_BL/a-1)+w*n/2)/(q_h-q_h/a*(1-DELTA_h)+corcho_r);
    c_r         = corcho_r*h_r;
    c_hat_r		= ((1-O_CHAT)^(1/ETA_CHAT)*(c_r*(1-PHI_c/a))^((ETA_CHAT-1)/ETA_CHAT)...
                    +(O_CHAT)^(1/ETA_CHAT)*(xi_h*h_r/a*(1-PHI_HH/a))^((ETA_CHAT-1)/ETA_CHAT))^(ETA_CHAT/(ETA_CHAT-1));
    la_r		= c_hat_r^(-SIGMA)*((1-O_CHAT)*c_hat_r/(c_r*(1-PHI_c/a)))^(1/ETA_CHAT);
        %Depend on both 
    h_u  		= (h_p-(1-OMEGA_UP)*h_r)/OMEGA_UP;
    c_u   		= ((a^SIGMA/BETA_up-(1-DELTA_h))*q_h/xi_h)^ETA_CHAT*...
                    (1-O_CHAT)/a/O_CHAT*h_u*xi_h*(1-PHI_HH/a)/(1-PHI_c/a);    
    c_hat_u		= ((1-O_CHAT)^(1/ETA_CHAT)*(c_u*(1-PHI_c/a))^((ETA_CHAT-1)/ETA_CHAT)...
                    +(O_CHAT)^(1/ETA_CHAT)*(xi_h*h_u/a*(1-PHI_HH/a))^((ETA_CHAT-1)/ETA_CHAT))^(ETA_CHAT/(ETA_CHAT-1));
    la_u		= c_hat_u^-SIGMA*((1-O_CHAT)*c_hat_u/(c_u*(1-PHI_c/a)))^(1/ETA_CHAT);
    la_p		= OMEGA_UP*la_u + (1-OMEGA_UP)*la_r;
    c_p			= OMEGA_UP*c_u + (1-OMEGA_UP)*c_r;
    la_W		= (la_p+la_i)/2;

    Xi_Tilde_up	= (c_hat_u)^(SIGMA);
    Theta_up	= Xi_Tilde_up*(c_hat_u)^(-SIGMA);
    Xi_Tilde_i	= (c_hat_i)^(SIGMA);
    Theta_i		= Xi_Tilde_i*(c_hat_i)^(-SIGMA);	
    Xi_Tilde_rp	= (c_hat_r)^(SIGMA);
    Theta_rp	= Xi_Tilde_rp*(c_hat_r)^(-SIGMA);
    Theta		= ((OMEGA_UP*Theta_up+(1-OMEGA_UP)*Theta_rp)+Theta_i)/2;
        %Objective Function 
    F(1) = mc_W - Theta*xi_n*n_tilde^VARPHI/(la_W*w);
    F(2) = - gdpn + c_p + c_i + i + i_h +g +tb  ;
end 


