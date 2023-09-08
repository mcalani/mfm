
% =========================================================================
% MAFIN Variables Declaration
% =========================================================================

% -------------------------------------------------------------------------
	// Main macro variables
@#define vars_macro =  ["b_ast_u" 		,\\ //Foreign bonds
						"b_ast_Tot"     ,\\ //Foreign bonds
						"bs_g"       	,\\ //Short term government assets
                        "bl_priv"     	,\\ //Total Private holding of long-term government bonds
						"bs_priv"     	,\\ //Total Private holding of short-term government bonds
						"c"     		,\\ //Consumption
						"gdp"   		,\\ //GDP real
                        "gdpR"   		,\\ //GDP real (no minero)
						"gdpn"  		,\\ //GDP nominal
						"i_agg" 		,\\ //Aggregate investment (capital + housing)
						"n"     		,\\ //Houshold labor supply
						"n_p"     		,\\ //Total Patient Household labor supply
						"n_i"     		,\\ //Impatient Household labor supply
						"n_u"     		,\\ //Unrestricted Patient Household labor supply
						"n_r"     		,\\ //Restricted Patient Household labor supply
						"ppi"   		,\\ //Inflation t-1 to t
						"ppi_S" 		,\\ //Nominal exchange rate inflation t-1 to t
						"R"    			,\\ //Local bonds interest rate
						"R_ast" 		,\\ //Foreign bonds interest rate
						"rer"   		,\\ //Real exchange rate
						"tau"       	,\\ //Lump sum taxes
						"Theta" 		,\\ //Endogenous preference shifter
						"Theta_up" 		,\\ //Endogenous preference shifter, Patient
						"Theta_rp" 		,\\ //Endogenous preference shifter, Restricted Patient
						"Theta_i" 		,\\ //Endogenous preference shifter, Impatient
						"tb"    		,\\ //Trade balance
						"w"     		,\\ //Wage
						"Xi_Tilde_up"	,\\ //Endogenous preference shifter, Patient
						"Xi_Tilde_rp"	,\\ //Endogenous preference shifter, Restricted Patient
						"Xi_Tilde_i"	,\\ //Endogenous preference shifter, Impatient
						"y_C"   		]	//Final goods

% -------------------------------------------------------------------------
	// Main financial variables
@#define vars_fin =    ["l"             ,\\ //Loans
						"n_b"           ,\\ //Banker net worth (bequest)
						"n_e"           ,\\ //Entrepreneur net worth (bequest)
						"nu"            ,\\ //Default costs (CSV monitoring costs and DIA deposit recovery costs)
						"PD_d"          ,\\ //Default probability of banks (on deposits)
						"PD_e"          ,\\ //Default probability of entrepreneurs
						"PD_i"          ,\\ //Default probability of impatient housholds
						"spread_RBL_R"  ,\\ //Spread between long rate and deposit rate
						"spread_RL_RD"  ,\\ //Spread between corporate loan rate and deposit rate
						"spread_RL_R"  ,\\  //Spread between corporate loan rate and deposit rate
						"spread_Ri_RD"  ,\\ //Spread between mortgage loan rate and deposit rate
                        "spread_Ri_RBL" , \\
						"spread_RD_R"   ]	//Spread between deposit loan rate and local bonds rate

% -------------------------------------------------------------------------
	// Patient household
@#define vars_p_HH =   ["bl_u"	,\\ //Unrestricted Patient household holding of long-term government bonds
						"bb_u"      	,\\ //Unrestricted Patient household holding of long-term bank bonds
						"bs_u"      	,\\ //Unrestricted Patient household holding of short-term government bonds
						"d_u"			,\\ //Deposits per capita by unrestricted households
						"c_p"   		,\\ //Consumption
						"c_hat_u"  		,\\ //CES Unrestricted (*****)
						"c_u"   		,\\ //Consumption unrestricted patient (*****)
						"h_p"   		,\\ //Houses
						"h_u"   		,\\ //Houses unrestricted patient (*****)
						"la_p"      	,\\	//Lagrange multiplier of budget constraint patient
						"la_u"     		,\\ //Lagrange multiplier of budget constraint unrestricted patient
						"q_BB"     		,\\ //Price of long term bank bond
                        "R_BB"     		,\\ //Long term bank bonds gross rate
                        "R_BB_tilde"	,\\ //Long term bank bonds gross return
						"Rnom_BL"        ,\\
						"q_BL"      	,\\ //Price of long term government bond rate (*****)
						"R_BL"      	,\\ //Long term government bonds rate (*****)
						"c_hat_r" 		,\\ //CES Restricted (*****)
						"c_r"					,\\ //consumption restricted patient (*****)
						"h_r"					,\\ //housing restricted patient (*****)
						"la_r"				,\\ //Lagrange multiplier of budget constraint restricted patient (*****)
						"bl_r"				,\\ //Restricted Patient household holding of long-term government bonds (*****)
						"zeta_L"    	]		//Transaction cost on long-term government bonds (*****)

% -------------------------------------------------------------------------
	// Impatient household
@#define vars_i_HH =   ["c_hat_i"   	,\\ //CES
						"c_i"   		,\\ //Consumption
						"h_i"   		,\\ //Houses
						"la_i"  		,\\ //Lagrange multiplier of budget constraint "lev_i" 		,\\ //Loan-to-value of mortgage loan
                        "la_h"          ,\\ //Lagrange multiplier FOC: l_h --------------------
						"Dl" , \\
						"R_h"   		,\\ //Gross nominal return per unit of housing bought
						"Rnom_i"        ,\\
                        "R_hat_i"       ,\\ // Mortage efective--------
						"R_i"				]	//Period return mortgage rate

% -------------------------------------------------------------------------
	// Entrepreneurs
@#define vars_e =      ["c_e"   		,\\ //Dividends
						"la_e"  		,\\ //Lagrange multiplier
						"lev_e" 		,\\ //Leverage
						"R_e"   		,\\ //Gross nominal return per unit of capital
						"r_k"   		,\\ //Rental rate of capital
						"R_L"   		,\\ //Gross nominal interest rate for corporate loan
						"psi_e" 		]	//Wealth

% -------------------------------------------------------------------------
	// Bankers & banks
@#define vars_b =      ["c_b"       	,\\ //Dividends
						"d_f"       	,\\ //Deposits in F banks
						"e_f"       	,\\ //Equity of F banks
						"bb_tot"		,\\ //banking bonds
						"e_h"       	,\\ //Equity of H banks
						"l_f"       	,\\ //Loans of F banks
						"l_h"       	,\\ //Loans of H banks
						"PD_f"      	,\\ //Default probability of F banks
						"PD_h"      	,\\ //Default probability of H banks
						"psi_b"     	,\\ //Banker wealth
						"R_D"       	,\\ //Gross nominal interest rate for deposits
						"R_D_tilde" 	,\\ //Gross nominal return of deposits
						"R_f_tilde" 	,\\ //Return on assets for F banks
						"R_h_tilde" 	,\\ //Return on assets for H banks 
						"Roe_aux" 		,\\ //Roe aux used for estimation 
						"rho_f"     	,\\ //Return on equity for F banks
						"rho_tilde_h"	,\\ //Return ave. on equity for H banks
						"rho_h"     	]	//Return on equity for H banks

% -------------------------------------------------------------------------
	// Capital and housing
@#define vars_kh =     ["h"				,\\ //Agregate housing
						"i"				,\\ //Investment in capital
						"i_ah"			,\\ //New authorized investment projects in housing
						"i_h"			,\\ //Investment in housing - current expenditure
						"k"				,\\ //Capital
						"q_h"			,\\ //Housing price
						"q_l"			,\\	//Loan price
						"q_hat_l"	    ,\\	//Loan price new --------
						"q_k"			,\\	//Capital price
                        "phi_f"         ,\\ //Capital requirement on corporate loans
                        "phi_h"         ,\\ //Capital requirement on mortgage loans
                        "buffer_guide"  ,\\ //ccyb buffer_guide ( CCYB(buffer_guide) )
                        "mort_tot"      ,\\ //total value of mortgage debt (Q_L*L_H)
						"mortVpar_tot"  ,\\ //total value of mortgage at par value (Q_hat_l + L_H)
                        "ltot"          ,\\ //credito total 
                        "ccyb"        ]     // Countercyclical capital buffer

% -------------------------------------------------------------------------
	// Production
@#define vars_prod =   ["f_F"			,\\ //Auxiliary for foreign goods of variety j
						"f_H"			,\\ //Auxiliary for home goods of variety j
						"imp"   		,\\ //Imports
						"mc_F"			,\\ //Real marginal cost in units of final foreign good
						"mc_H"			,\\ //Real marginal cost in units of final domestic good
						"mc_Z"			,\\ //Real marginal cost wholesale goods, in units of final domestic good
						"p_F"			,\\ //Price of foreign composite goods
						"p_H"			,\\ //Price of home composite goods
						"p_Z"			,\\ //Price of wholesale domestic goods (####)
						"p_tilde_F"		,\\ //Price that maximizes the current market value (foreign goods of variety)
						"p_tilde_H"		,\\ //Price that maximizes the current market value (home goods of variety)
						"ppi_F"			,\\ //Inflation t-1 to t foreign composite goods
						"ppi_H"			,\\ //Inflation t-1 to t home composite goods
						"ppi_IF"		,\\ //Price adjustment of non-optimizers, foreign goods
						"ppi_IH"		,\\ //Price adjustment of non-optimizers, home goods
						"x_F"			,\\ //Foreign composite goods (demand)
						"x_H"			,\\ //Home composite goods (home demand)
						"x_Z"			,\\ //Wholesale domestic goods (home demand) (####)
						"x_H_ast"		,\\ //Foreign demand for the home composite good
						"Xi_F"			,\\ //Price dispersion of foreign good
						"Xim_h"			,\\ //Price dispersion of home good
						"y_F"			,\\ //Foreign composite goods (supply)
						"y_H"			,\\ //Home composite goods (supply)
						"y_Z"			]	//Wholesale domestic goods (supply) (####)

% -------------------------------------------------------------------------
	// Wages
@#define vars_wages =  ["mc_W"			,\\	//Household marginal cost of providing an extra hour of work
						"n_tilde"		,\\	//Composite labor services
						"la_W"			,\\	//Average lagrange multiplier of budget constraints among households
						"f_W"				,\\	//Auxiliary for wages of variety i
						"w_tilde"		,\\	//Wages of reoptimizing unions relative to average wages
						"ppi_W"			,\\	//Nominal wage inflation
						"ppi_Wtilde"	,\\	//Nominal wage inflation of reoptimizing unions
						"ppi_IW"		,\\	//Nominal wage indexation of non-reoptimizing unions
						"Xi_W"			]	//Wage dispersion

% -------------------------------------------------------------------------
	// BGG auxiliaries
@#define vars_BGG =    ["aux"			,\\ //Auxiliary
						"Gam"			,\\ //Gamma: Share of loan returns taken by lender
						"Gam_der"		,\\ //Derivative of Gamma with respect to om_bar
						"Gam_G"			,\\ //Gamma minus mu_G: Share of loan returns taken by lender after paying for monitoring costs
						"Gam_G_der"		,\\ //Derivative of Gam_G with respect to om_bar
						"mu_G"			,\\ //mu times G: Total monitoring costs arising from default and costly state verification
						"om_bar"		]	//Default threshold for idiosyncratic return of loans (an omega below would default)

@#define agent = ["i", "e", "f", "h"	]	//Impatient households, Entrepreneurs, F banks and H banks respectively

% -------------------------------------------------------------------------
	// Other variables
@#define vars_other =  ["c_p_i"			,\\ //Consumption of patients relative to consumption of impatient HHs
						"brechacred"  	,\\ //brecha de credito
                        "lev_tot"       ,\\ //leverage total banking system 
                        "lev_f"         ,\\ //leverage f bank 
                        "lev_h"         ,\\ //leverage h bank 
						"lab_income"	,\\ //Labour income
						"nu_gdp"		]	//Default costs over GDP

% -------------------------------------------------------------------------
	// Observable variables
@#define vars_obs =    ["gam_C"			,\\ //Private consumption growth
						"gam_G"			,\\ //Government expenditure growth
						"gam_I"			,\\ //Capital investment growth
						"gam_IK"		,\\
						"gam_lh"		,\\ //Credit housing growth 
						"gam_lf"		,\\ //Credit comercial growth 
						"gam_ltot"  ,\\ //Total credit 
						"gam_IH"		,\\ //Housing investment growth
						"gam_N"			,\\ //Total hours growth
						"gam_WN"		,\\ //Nominal wage growth
						"gam_YCo"		,\\ //Commodity production growth
						"gam_YR"		,\\ //Non-commodity GDP growth
						"gam_Y"         ,\\ //GDP growth
                        "gam4_YR"		,\\ //Non-commodity GDP growth(annual)
						"gam_Ystar"		,\\ //Commodity production growth
						"pi"			,\\ //Inflation
						"piZ"			,\\ //Core Inflation
						"gam_ner"		,\\ //Core Inflation
                        "piZ4"			,\\ //Core Inflation(annual)
						"piCostar"		,\\ //Commodity price inflation
						"piQH"			,\\ //Housing price inflation
						"piM"			,\\ //Imports inflation
						"pistar"		,\\ //Foreign inflation
						"rer"			,\\ //Real exchange rate
						"Roe"			,\\ //ROE
						"RLG"			,\\ //Long-term rate
						"R"				,\\ //Monetary policy rate
						"RD"			,\\ //Bank funding rate
						"RI"			,\\ //Housing loans rate
						"RL"			,\\ //Corporate loans rate
						"Rstar"			,\\ //Foreign interest rate
						"stb"			,\\ //Trade balance share of GDP
                        "xi"            ,\\
                        "ccyb"          ]	//Country risk premium

% -------------------------------------------------------------------------

	// Exogenous variables
@#define vars_exo =    ["a"				,\\ //Non-stationary productivity
						"bl_g"       	,\\ //Long-term government bonds supply
						"bl_cb"      	,\\ //Long term government bonds purchases by the CB
						"e_m"			,\\ //Monetary policy shock
                        "e_req"         ,\\ //Capital requirement shock
						"epsilon_L" 	,\\ //Transaction costs shock (*****)
						"g"				,\\ //Government expenditure
						"p_co"			,\\ //Copper price
						"ppi_ast"		,\\ //Foreign inflation
						"R_W"			,\\ //Foreign interest rate
						"sigma_ee"		,\\ //Volatility of entrepreneurial capital idiosyncratic shock
						"sigma_ff"		,\\ //Volatility of corporate loans idiosyncratic shock
						"sigma_hh"		,\\ //Volatility of mortgage loans idiosyncratic shock
						"sigma_ii"		,\\ //Volatility of housing idiosyncratic shock
						"varrho"		,\\ //Intertemporal preference
						"xi_h"			,\\ //Housing preference
						"xi_i"			,\\ //Capital investment efficiency
						"xi_ih"			,\\ //Housing investment efficiency
						"xi_m"			,\\ //Imports relative price
						"xi_n"			,\\ //Leisure preference
						"xi_R"			,\\ //Country premium (on foreign interest rate) shock
						"xi_roe_r"  	,\\ //Expected returns ratio
						"xi_CHI_b"  	,\\ //Dividend policy banks
                        "xi_CHI_e"  	,\\ //Dividend policy entrepreneurs
						"y_ast"			,\\ //Foreign GDP
						"y_co"			,\\ //Copper production            
                        "zetau"			,\\ //Unobservable country premium shock
						"z"				]	//Stationary productivity

% -------------------------------------------------------------------------
	// Declare endogenous variables
@#define vars_endo = vars_macro + vars_fin + vars_p_HH + vars_i_HH + vars_e + vars_b + vars_kh + vars_prod + vars_wages + vars_other //Most endogenous variables, BGG auxiliaries and observables are added below
@#for va in vars_endo + vars_exo
	var @{va};
@#endfor

//var credit_gap;

% -------------------------------------------------------------------------
	// Declare auxiliary BGG variables
@#for va in vars_BGG
	@#for ag in agent
		var @{va}_@{ag};
		@#define vars_endo = vars_endo + [va + "_" + ag]
	@#endfor
@#endfor

% -------------------------------------------------------------------------
	// Declare observable variables as endogenous
@#for va in vars_obs
	var @{va}_obs;
	@#define vars_endo = vars_endo + [va + "_obs"]
@#endfor

% -------------------------------------------------------------------------
	// Declare shocks
@#for va in vars_exo
	varexo u_@{va};
@#endfor
