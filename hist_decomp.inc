%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%						MAFIN Historical Decomposition						%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

// Historical decomposition ------------------------------------------------
 /*	
	shock_groups(name = group1);
		'Non-stationary prod.'			= u_a;
		'Monetary policy'				= u_e_m;
		//'Government exp.'				= u_g;
		'Copper price'					= u_p_co;
		'Foreign inflation'				= u_ppi_ast;
		'Foreign interest rate'			= u_R_W;
		// 'Volatility of entrepreneurial capital idiosyncratic shock'	=	u_sigma_ee;
		// 'Volatility of corporate loans idiosyncratic shock'			=	u_sigma_ff;
		// 'Volatility of mortgage loans idiosyncratic shock'			=	u_sigma_hh;
		// 'Volatility of housing idiosyncratic shock'					=	u_sigma_ii;
		'Preferences'					= u_varrho;
		//'Housing'						= u_xi_h;
		//'Capital investment eff.'		= u_xi_i;
		//'Housing investment eff.'		= u_xi_ih;
		//'Imports relative price'		= u_xi_m;
		//'Leisure'						= u_xi_n;
		'Obs. Country premium'			= u_xi_R;
		'Foreign GDP'					= u_y_ast;
		'Copper prod.'					= u_y_co;
		'Unobs. Country premium'		= u_zetau;
		'Stationary prod.'				= u_z;
	end;
*/

	shock_groups(name = group1);
		'Productivity'                  = u_z,u_a,u_xi_ih,u_xi_i,u_xi_n;
        'Demand'                        = u_varrho, u_xi_h, u_g ; //
		'Monetary policy'				= u_e_m;
        'Financial: Risk'               = u_sigma_ee,u_sigma_ff,u_sigma_hh,u_sigma_ii,u_xi_roe_r,u_xi_CHI_b,u_xi_CHI_e,u_epsilon_L;
        'External'                      = u_ppi_ast,u_R_W,u_xi_R,u_y_ast,u_y_co,u_zetau;
		

	end;
	
	shock_decomposition(use_shock_groups = group1, nograph)        gam_YR_obs  piZ_obs R_obs gam_ltot_obs RL_obs;
	plot_shock_decomposition(use_shock_groups = group1, write_xls, colormap=parula) gam_YR_obs  piZ_obs R_obs gam_ltot_obs RL_obs;


	