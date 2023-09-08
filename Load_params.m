%-------------------------------------------------------------------------%
%	Parametros
%-------------------------------------------------------------------------%
% clc
% clear all

%%MC
% Exceptions={'phi_f_base','phi_h_base','ALPHA','ALPHA_BSG','ALPHA_BLG'...
%             ,'BETA_up','BETA_rp','BETA_i','CHI','CHI_b','CHI_e','DELTA_h','DELTA_k','EPSILON_F'...
%             ,'EPSILON_H','EPSILON_W','ETA','ETA_ast','ETA_CHAT','ETA_ZETA_L','GAMMA_d','GAMMA_bh'...
%             ,'KAPPA_LOAN','KAPPA_BH','KAPPA_B','MU_e','MU_f','MU_h','MU_i','N_H','OMEGA','OMEGA_UP'...
%             ,'OMEGA_BL','PHI_c','PHI_HH','RHO_VARPHI_H','SIGMA','THETA_F','THETA_H','THETA_W'...
%             ,'UPSILON_H','VARPHI','VARPHI_H_0'...
%             ,'a','bl_cb','epsilon_L','g','n','r_h_k','ppi_T'...
%             ,'p_co','ppi_ast','R_W','sigma_ee','sigma_ff','sigma_hh','sigma_ii','xi_CHI_e','xi_CHI_b','xi_h'...
%             ,'xi_i','xi_ih','xi_m','xi_n','xi_R','xi_roe_r','y_ast','y_co','z' ...
%             'O_CHAT','bl_g','bs_g','b_ast_Tot','ren_ast'}; 
         
Exogenous={'phi_f_base','phi_h_base','ccyb','ALPHA','ALPHA_BSG','ALPHA_BLG'...
            ,'BETA_up','BETA_rp','BETA_i','CHI','CHI_b','CHI_e','DELTA_h','DELTA_k','EPSILON_F'...
            ,'EPSILON_H','EPSILON_W','ETA','ETA_ast','ETA_CHAT','ETA_ZETA_L','GAMMA_d','GAMMA_bh'...
            ,'KAPPA_LOAN','KAPPA_BH','KAPPA_B','MU_e','MU_f','MU_h','MU_i','N_H','OMEGA','OMEGA_UP'...
            ,'OMEGA_BL','O_CHAT','PHI_c','PHI_HH','RHO_VARPHI_H','SIGMA','THETA_F','THETA_H','THETA_W'...
            ,'UPSILON_H','VARPHI','VARPHI_H_0'...
            ,'a','bl_cb','epsilon_L','g','n','r_h_k','ppi_T'...
            ,'p_co','ppi_ast','R_W','sigma_ee','sigma_ff','sigma_hh','sigma_ii','xi_CHI_e','xi_CHI_b','xi_h'...
            ,'xi_i','xi_ih','xi_m','xi_n','xi_R','xi_roe_r','y_ast','y_co','z'...
            ,'bl_g','bs_g','b_ast_Tot'};

SS_values=load("SS_values.mat");
%SS_values=load("SS_values_FFago10.mat");


if ~exist('display_pars')
    display_pars=0 ; 
end

if display_pars ==1
    disp('Fixed to its ss-value:')
end
for ii=1:length(SS_values.param_names)
    
    %drop '_ss'
    if contains(SS_values.param_names{ii},'_ss')
        SS_values.param_names{ii} = SS_values.param_names{ii}(1:end-3);
    end
        
    %evaluate to ss-value if it isn in the list
    if any(strcmp(Exogenous,SS_values.param_names{ii}))
        eval(append(SS_values.param_names{ii},'=SS_values.params(ii);'));
        if display_pars ==1
            disp(SS_values.param_names{ii})
        end
    else
    end
end




