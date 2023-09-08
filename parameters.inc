
% =========================================================================
% MAFIN Parameters Declaration
% =========================================================================

% -------------------------------------------------------------------------
	// List deep parameters
@#define pars =	   ["ALPHA"			,\\ //Capital share in home production
					"ALPHA_E"		,\\ //Expected inflation weight in Taylor rule
                    "GAMMA_c"       ,\\ // ccyb process
					"ALPHA_ppi"		,\\ //Monetary policy rate elasticity to inflation (Taylor rule)
					"ALPHA_R"		,\\ //Monetary policy rate persistence (Taylor rule)
					"ALPHA_TAU"     ,\\ //Percentage of GDP paid as lump sum taxes
                    "ALPHA_BSG"     ,\\ //Short-term Bonds as Percentage of GDP
					"ALPHA_BLG"     ,\\ //Long-term Bonds as Percentage of GDP
					"ALPHA_W"		,\\ //Weight of past productivity growth (vs SS) in wage indexation
					"ALPHA_y"		,\\ //Monetary policy rate elasticity to GDP growth (Taylor rule)
					"BETA_i"		,\\ //Discount factor of impatient HHs
					"BETA_up"		,\\ //Discount factor of unrestricted patient HHs
					"BETA_rp"		,\\ //Discount factor of restricted patient HHs
					"CHI"			,\\ //Share of commodity income that is collected by the government
					"CHI_b"			,\\ //Share of banker wealth that is transferred to HHs as dividends
					"CHI_e"			,\\ //Share of entrepreneur wealth that is transferred to HHs as dividends
					"DELTA_h"		,\\ //Housing depreciation rate
					"DELTA_k"		,\\ //Capital depreciation rate
					"Duration_LOAN"	,\\ //Duration long-term housing debt, quarters
					"Duration_BL"	,\\ //Duration long-term governement bonds, quarters
					"Duration_BB"	,\\ //Duration long-term housing debt, quarters
					"EPSILON_F"		,\\ //Elasticity of substitution between varieties of foreign good
					"EPSILON_H"		,\\ //Elasticity of substitution between varieties of home good
					"EPSILON_W"		,\\ //Elasticity of substitution between varieties of labor services
					"EPSILON_TAU"   ,\\ //Tax adjustment parameter
					"ETA"			,\\ //Elasticity of substitution between home and foreign composite goods
					"ETA_ast"		,\\ //Price-elasticity of foreign demand for home composite goods
					"ETA_CHAT"		,\\ //CES
                    "GAMA_L"		,\\ //Adjustment cost for new credit
					"GAMA_h"		,\\ //Adjustment cost for housing investment
					"GAMA_k"		,\\ //Adjustment cost for capital investment
					"GAMA_n"		,\\ //Adjustment cost for capital investment
					"GAMMA_bh"		,\\ //Defaulted deposit recovery cost, long-term bank bond
					"GAMMA_d"		,\\ //Defaulted deposit recovery cost
					"KAPPA_B"       ,\\ //Geometric decay long-term government bonds
					"KAPPA_BH"      ,\\ //Geometric decay long-term bank bonds
					"KAPPA_F"		,\\ //Weight of past inflation (vs SS) in foreign price indexation
					"KAPPA_H"		,\\ //Weight of past inflation (vs SS) in home price indexation
					"KAPPA_W"		,\\ //Weight of past inflation (vs SS) in wage indexation
					"KAPPA_LOAN"	,\\ //Geometric decay long term bonds
					"MU_e"			,\\ //Monitoring cost of corporate loan default
					"MU_f"			,\\ //Monitoring cost of F bank default
					"MU_h"			,\\ //Monitoring cost of H bank default
					"MU_i"			,\\ //Monitoring cost of housing loan default
					"O_CHAT"		,\\ //CES
					"OMEGA_BL"      ,\\ //Ratio of long term bonds to short term bonds
					"OMEGA_UP"      ,\\ //Fraction of Unrestricted patient households
					"UPSILON_H"     ,\\ //Ratio of Restricted patient housing to Unrestricted patient housing
					"PHI_ast"		,\\ //Country premium elasticity to foreign debt over GDP
					"PHI_c"			,\\ //Consumption habits parameter
					"PHI_HH"		,\\ //Housing habits parameter
					"ppi_T"			,\\ //Inflation target
					"OMEGA"			,\\ //Share of home goods in final goods production
					"N_H"			,\\ //Time-to-build periods in housing production
					"ren_ast"		,\\ //Net foreign rents
					"ETA_ZETA_L"    ,\\ //Elasticity on transaction costs by trading long-term governement bonds
					"RHO_VARPHI_H"	,\\ //Front(back)-loading of expenditure in housing investment
					"SIGMA"			,\\ //Inverse of the intertemporal elasticity of substitution of HHs
					"THETA_F"		,\\ //Frequency of price adjustments á-la-Calvo in foreign goods
					"THETA_H"		,\\ //Frequency of price adjustments á-la-Calvo in home goods
					"THETA_W"		,\\ //Frequency of price adjustments á-la-Calvo in nominal wages
					"V_PF"			,\\ //Strenght of wealth effect
                    "phi_h_base"    ,\\
                    "phi_f_base"    ,\\
                    "ccybneutral"    ,\\
					"VARPHI"		]	//Inverse elasticity of labor supply

% -------------------------------------------------------------------------
	// Steady state moments
@#define pars_SS = ["r_h_k"			,\\ //Ratio of housing to capital stock
					"s_bast"		,\\ //Net foreign asset position as share of GDP
					"s_co"			,\\ //Share of commodity production in GDP
					"s_g"			,\\ //Share of government expenditure in GDP
					"s_tb"			]	//Share of trade balance in GDP

% -------------------------------------------------------------------------
	// Parameters for no financial frictions
@#define pars_NF = ["BAR_l_h"		,\\ //Steady state level of l_h
					"IOTA"			,\\ //Steady state spread of housing loans over risk-free rate
					"PHI_i"			]	//Elasticity of impatient household's interest rate to their debt level

% -------------------------------------------------------------------------
	// Options as parameters
@#define pars_opt= ["FinFric"]			//Financial frictions

% -------------------------------------------------------------------------
	// Declare deep parameters
@#for pa in pars + pars_NF
	parameters @{pa};
@#endfor

% -------------------------------------------------------------------------
	// Deep parameters declared with macro-processing
@#for j in 0:NH
	parameters VARPHI_H_@{j};			//Fraction of projects authorized j periods ago, that is outlaid today.
@#endfor

% -------------------------------------------------------------------------
	// SS moments and options
@#for pa in pars_SS + pars_opt
	parameters @{pa};
@#endfor

% -------------------------------------------------------------------------
	// SS values
@#for va in vars_endo + vars_exo
	parameters @{va}_ss;
@#endfor

% -------------------------------------------------------------------------
	// Shock parameters
@#for va in vars_exo
	parameters SIGMA_@{va} RHO_@{va};
@#endfor
