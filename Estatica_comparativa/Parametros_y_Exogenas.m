%-------------------------------------------------------------------------%
%	Parametros
%-------------------------------------------------------------------------%
%   clc
%   clear
% close all
% 
% %CCyB        = 0.025; 
% %phi_f		=	(0.08+CCyB)*1+0.043	  ;
% %phi_h		=	(0.08+CCyB)*0.6+0.043 ;

capreq=0.08;

ALPHA		=	1-0.66 ;
ALPHA_BSG   =  -0.4    ; %%%%NUEVO.0121
ALPHA_BLG   = -1.5     ; %%%%NUEVO.0121
ALPHA_TAU   = 1        ; %%%%NUEVO.0121 se toma como ==1
%ALPHA_E		=	0.5     ;
%ALPHA_ppi      =	1.7698	;
%ALPHA_R		=	0.7879	;
%ALPHA_W		=	0.1581	;
%ALPHA_y		=	0.1861	;

BETA_up     =   0.99997 ; %%%%NUEVO.0121
BETA_rp     =   0.995690002 ; %%%%NUEVO.0121
BETA_i		=	0.98	;
BETA_p		=	0.99997 ;

CHI         =	0.33	;
CHI_b		=	0.05	;
CHI_e		=	0.05	;

coupon		=	0.0450	;

DELTA_h		=	0.01	;
DELTA_k		=	0.01	;
EPSILON_F   =	11  	;
EPSILON_H   =	11  	;
EPSILON_W	=	11  	;

%ETA         =	1.0104	;  (cal)
ETA         = 1.783258765; %(est)

ETA_ast		=	0.148826201; %(est) % 0.2442 (cal);
%GAMA_h		=	3.0029	;
%GAMA_k		=	2.0641	;
%GAMA_n		=	1.0907 	;

GAMMA_d		=	0.1 	;
GAMMA_bh    =   0.1     ; %%%%NUEVO.0121

%KAPPA_F		=	0.2099	;
%KAPPA_H		=	0.1979	;
%KAPPA_W		=	0.1917	;
KAPPA_LOAN	=	0.975 	;
KAPPA_BH    =   0.8     ; %%%%NUEVO.0121

MU_e		=	0.3 	;
MU_f		=	0.3 	;
MU_h		=	0.3  ;
MU_i		=	0.3  ;

%PHI_ast		=	0.2171	;
PHI_c		=	0.339775949	; % (est) % 0.8560 (cal)	;
PHI_HH		=	0.889558938	; % (est) % 0  (cal)  ;



OMEGA		=	0.79	;

N_H     	=	6       ;
RHO_VARPHI_H=	1       ;

SIGMA		=	1   	;

THETA_F		=	0.0628	;
THETA_H		=	0.7937	;
THETA_W		=	0.7829	;

VARPHI		=	6.244102896 ; %(est) % 4.5474 (cal)	;
VARPHI_H_0	=	0.1429       ; % hay que calibrarlo!  Rho_varPHI_H^0*(1-RHO_VARPHI_H)/(1-RHO_VARPHI_H^N_H+1)=0?



ETA_CHAT	= 0.625871283 ; %(est) %0.7 (cal)  ; 
UPSILON_H   = 1         ;%%%%NUEVO.0121
OMEGA_UP    = 0.7       ;%%%%NUEVO.0121
KAPPA_B     = 0.8       ;%%%%NUEVO.0121
OMEGA_BL    = 0.822     ;%%%%NUEVO.0121 
ETA_ZETA_L  = 0.14013585 ; %(est) %0.9402 (cal)    ;%%%%NUEVO.0121 


%-------------------------------------------------------------------------%
%	Variables Exogenas
%-------------------------------------------------------------------------%
ppi_T		=	1.0074 ;
a		=	1.0037	; 
%e_m		=	1	;
g		=	0.1165	;
p_co	=	1	;
ppi_ast	=	1.0074	;
R_W		=	1.0074	;
sigma_ee=	0.206045 ;
sigma_ff=	0.0552	;
sigma_hh=	0.0401	;
sigma_ii=	0.0154693	;
%varrho	=	1	;
xi_h	=	1	;
xi_i	=	1	;
xi_ih	=	1	;
xi_m	=	1	;
xi_n	=	1.3753e+03	;
xi_R	=	1.0037	;
y_ast	=	0.1140	;
y_co	=	0.1281	;
%zetau	=	1	;
z		=	1	;
xi_CHI_e=   1   ; 
xi_CHI_b=   1   ;
xi_roe_r=   1   ;
epsilon_L=0.0044;   %%%%NUEVO.0121
blcb = 1    ;
n=0.3;


r_h_k		=	0.65	;
s_co		=	0.12	;
s_g         =	0.12	;
s_tb		=	0.05	;
s_bast		=	-0.14	;