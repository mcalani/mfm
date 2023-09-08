clear all
close all
clc
tic

disp('runing!')
vec=zeros(1,1);
vecex=zeros(1,1);
a_i=1;

%cabiar acam y name(linea 20)
acam={'capreq'};
%acam={'MU_e'};

for num_var = 0.045:0.005:0.15
%for num_var = 0.15:0.01:0.45

    clearvars -except a_i num_var vec acam name vecex


run Parametros_y_Exogenas;

eval ([acam{1}  '=num_var']);
name=[capreq];
%name=[MU_e];
run MAFIN_03_est_steadystate_EsCom;

vars2save={a,ALPHA,ALPHA_BLG,ALPHA_BSG,ALPHA_TAU,aux_e,aux_f,aux_h,aux_i,b_ast_Tot,b_ast_u,bb_tot,bb_u,BETA_i,BETA_p,BETA_rp,BETA_up,bl_cb,bl_g,bl_priv,bl_r,bl_u,blcb,bs_g,bs_priv,bs_u,c,c_b,c_e,c_hat_i,c_hat_r,c_hat_u,c_i,c_p,c_p_i,c_r,c_u,capreq,CHI,CHI_b,CHI_e,corcho,coupon,d_f,d_u,DELTA_h,DELTA_k,e_f,e_h,EPSILON_F,EPSILON_H,epsilon_L,EPSILON_W,ETA,ETA_ast,ETA_CHAT,ETA_ZETA_L,f_F,f_H,f_W,g,Gam_der_e,Gam_der_f,Gam_der_h,Gam_der_i,Gam_e,Gam_f,Gam_G_der_e,Gam_G_der_f,Gam_G_der_h,Gam_G_der_i,Gam_G_e,Gam_G_f,Gam_G_h,Gam_G_i,Gam_h,Gam_i,GAMMA_bh,GAMMA_d,gdp,gdpn,gdpR,h,h_i,h_p,h_r,h_u,i,i_agg,i_ah,i_h,imp,k,KAPPA_B,KAPPA_BH,KAPPA_LOAN,l,l_f,l_h,l_h_guess,la_e,la_i,la_p,la_r,la_u,la_W,lab_income,lev_e,mc_F,mc_H,mc_W,mc_Z,MU_e,MU_f,mu_G_e,mu_G_f,mu_G_h,mu_G_i,MU_h,MU_i,n,n_b,n_e,N_H,n_i,n_p,n_r,n_tilde,n_u,nu,nu_gdp,O_CHAT,om_bar_e,om_bar_f,om_bar_h,om_bar_i,OMEGA,OMEGA_BL,OMEGA_UP,p_co,p_F,p_H,p_tilde_F,p_tilde_H,p_Z,PD_d,PD_e,PD_f,PD_h,PD_i,PHI_c,phi_f,phi_h,PHI_HH,ppi,ppi_ast,ppi_F,ppi_H,ppi_IF,ppi_IH,ppi_IW,ppi_S,ppi_T,ppi_W,ppi_Wtilde,psi_b,psi_e,q_BB,q_BL,q_h,q_k,q_l,R,R_ast,R_BB,R_BB_tilde,R_BL,R_D,R_D_tilde,R_e,R_f_tilde,R_h,r_h_k,R_h_tilde,R_i,r_k,R_L,R_W,ren_ast,rer,rho_f,rho_h,rho_tilde_h,RHO_VARPHI_H,Rnom_BL,s_bast,s_co,s_g,s_tb,SIGMA,sigma_ee,sigma_ff,sigma_hh,sigma_ii,spread_RBL_R,spread_RD_R,spread_Ri_RD,spread_RL_RD,tau,tb,Theta,THETA_F,THETA_H,Theta_i,Theta_rp,Theta_up,THETA_W,UPSILON_H,VARPHI,VARPHI_H_0,w,w_tilde,x,x_F,x_guess,x_H,x_H_ast,x_Z,xi_CHI_b,xi_CHI_e,Xi_F,xi_h,xi_i,xi_ih,xi_m,xi_n,xi_R,xi_roe_r,Xi_Tilde_i,Xi_Tilde_rp,Xi_Tilde_up,Xi_W,Xim_h,y_ast,y_C,y_co,y_F,y_H,y_Z,z,zeta_L};
varsName={'a','ALPHA','ALPHA_BLG','ALPHA_BSG','ALPHA_TAU','aux_e','aux_f','aux_h','aux_i','b_ast_Tot','b_ast_u','bb_tot','bb_u','BETA_i','BETA_p','BETA_rp','BETA_up','bl_cb','bl_g','bl_priv','bl_r','bl_u','blcb','bs_g','bs_priv','bs_u','c','c_b','c_e','c_hat_i','c_hat_r','c_hat_u','c_i','c_p','c_p_i','c_r','c_u','capreq','CHI','CHI_b','CHI_e','corcho','coupon','d_f','d_u','DELTA_h','DELTA_k','e_f','e_h','EPSILON_F','EPSILON_H','epsilon_L','EPSILON_W','ETA','ETA_ast','ETA_CHAT','ETA_ZETA_L','f_F','f_H','f_W','g','Gam_der_e','Gam_der_f','Gam_der_h','Gam_der_i','Gam_e','Gam_f','Gam_G_der_e','Gam_G_der_f','Gam_G_der_h','Gam_G_der_i','Gam_G_e','Gam_G_f','Gam_G_h','Gam_G_i','Gam_h','Gam_i','GAMMA_bh','GAMMA_d','gdp','gdpn','gdpR','h','h_i','h_p','h_r','h_u','i','i_agg','i_ah','i_h','imp','k','KAPPA_B','KAPPA_BH','KAPPA_LOAN','l','l_f','l_h','l_h_guess','la_e','la_i','la_p','la_r','la_u','la_W','lab_income','lev_e','mc_F','mc_H','mc_W','mc_Z','MU_e','MU_f','mu_G_e','mu_G_f','mu_G_h','mu_G_i','MU_h','MU_i','n','n_b','n_e','N_H','n_i','n_p','n_r','n_tilde','n_u','nu','nu_gdp','O_CHAT','om_bar_e','om_bar_f','om_bar_h','om_bar_i','OMEGA','OMEGA_BL','OMEGA_UP','p_co','p_F','p_H','p_tilde_F','p_tilde_H','p_Z','PD_d','PD_e','PD_f','PD_h','PD_i','PHI_c','phi_f','phi_h','PHI_HH','ppi','ppi_ast','ppi_F','ppi_H','ppi_IF','ppi_IH','ppi_IW','ppi_S','ppi_T','ppi_W','ppi_Wtilde','psi_b','psi_e','q_BB','q_BL','q_h','q_k','q_l','R','R_ast','R_BB','R_BB_tilde','R_BL','R_D','R_D_tilde','R_e','R_f_tilde','R_h','r_h_k','R_h_tilde','R_i','r_k','R_L','R_W','ren_ast','rer','rho_f','rho_h','rho_tilde_h','RHO_VARPHI_H','Rnom_BL','s_bast','s_co','s_g','s_tb','SIGMA','sigma_ee','sigma_ff','sigma_hh','sigma_ii','spread_RBL_R','spread_RD_R','spread_Ri_RD','spread_RL_RD','tau','tb','Theta','THETA_F','THETA_H','Theta_i','Theta_rp','Theta_up','THETA_W','UPSILON_H','VARPHI','VARPHI_H_0','w','w_tilde','x','x_F','x_guess','x_H','x_H_ast','x_Z','xi_CHI_b','xi_CHI_e','Xi_F','xi_h','xi_i','xi_ih','xi_m','xi_n','xi_R','xi_roe_r','Xi_Tilde_i','Xi_Tilde_rp','Xi_Tilde_up','Xi_W','Xim_h','y_ast','y_C','y_co','y_F','y_H','y_Z','z','zeta_L'};

graphs={name,gdp,c,i_agg,k,h,R,ppi,l_f,l_h,R_L,R_i,R_D,PD_e*100,PD_i*100,PD_h*100,PD_f*100,n_e,n_b};

if strcmp(acam{1},'capreq')
    graphs(1)={phi_f*100};  
end


loop_length=length(graphs);
cols=ceil(loop_length/3);


for i=1:1:loop_length
x=cell2mat(graphs(i)); 
    vec(a_i,i)=x;
end

b_i=1;
for Z=1:numel(vars2save)

tt=cell2mat(vars2save(Z));    
vecex(a_i,b_i)=tt;

b_i=b_i+1;
end
a_i=a_i+1;
end


names={acam{1} 'GDP' 'Consumption' 'Investment' 'Capital' 'Housing' 'MPR' 'Inflation' 'Cor. Loans' 'Mort Loans'  '$R^L$' '$R^I$' '$R^D$' '$$PD^e$' '$PD^I$' '$PD^H$' '$PD^F$'  '$N^e$' '$N^b$'};

if strcmp(acam{1},'capreq')
    names{1}='$\phi_f$';
end    
    
for Z=1:numel(vars2save)
namesex(Z)=cellstr(varsName(Z));
end


clearvars -except vec names vecex namesex loop_length acam



disp('plottig!')
%%


close all


t = tiledlayout(4,6); 
set(gcf,'PaperOrientation','landscape')
for i = 2:1:loop_length
nexttile
plot(vec(:,1),vec(:,i),'b','LineWidth',2)
    title(names(i),'interpreter','latex','FontSize',8)
    xlabel(names(1),'interpreter','latex','FontSize',6)
    set(gca,'FontSize', 8,'FontName', 'Times');
end
disp('Exporting!')
figure=acam{1};
print(figure,'-dpdf','-bestfit');
filename='resultados.xlsx';
writecell(namesex,filename,'Sheet',figure,'Range','A1','WriteMode','overwritesheet')
writematrix(vecex,filename,'Sheet',figure,'Range','A2')

disp('done!')
toc