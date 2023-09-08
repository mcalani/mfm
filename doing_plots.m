% Plot and save impulse responses
clear;
clc;
close all;

ModelNames  = {'MAFIN_comp'};	
           
for q=1:length(ModelNames)
     
    DIR=['figs_' ModelNames{q}];
    LOC_final = [pwd '\' DIR];
    
                  
    model_names = { 
                    %'mafin_04_t05',          'h.l.=1';    % phi_base=0.097, AR_ccyb=0.5, peakin=1 
                    %'mafin_04_t07',          'h.l.=2';    % phi_base=0.097, AR_ccyb=0.7, peakin=1    
                    %'mafin_04_t09',          'h.l.=7';    % phi_base=0.097, AR_ccyb=0.9, peakin=1  
                    %'mafin_04_i09175',       'Capital base 8%';    % phi_base=0.097, AR_ccyb=0.9175, halflife = 8   
                    %'mafin_04_t09175',       '9.5%';    % phi_base=0.097, AR_ccyb=0.9175, halflife = 8 
                    'MAFIN_03_est_results',       '12.5%';    % phi_base=0.097, AR_ccyb=0.9175, halflife = 8 
                    %'mafin_04_phi112_halflife2',       'Vida media 2Q';    % phi_base=0.112, AR_ccyb=0.9175, halflife = 8 
                    %'mafin_04_phi112_halflife8',       'Vida media 8Q';    % phi_base=0.112, AR_ccyb=0.9175, halflife = 8 
                    %'mafin_04_phi112_halflife16',      'Vida media 16Q';    % phi_base=0.112, AR_ccyb=0.9175, halflife = 8 
                    
                    %'mafin_04_phi08_halflife8',       'Capital base 8%';    % phi_base=0.112, AR_ccyb=0.9175, halflife = 8 
                    %'mafin_04_phi97_halflife8',       '9.5%';    % phi_base=0.112, AR_ccyb=0.9175, halflife = 8 
                    %'mafin_04_phi125_halflife8',      '12.5%';    % phi_base=0.112, AR_ccyb=0.9175, halflife = 8 
                    
                    %'mafin_04_t096',         'h.l.=16';    % phi_base=0.097, AR_ccyb=0.9175, peakin=1  
                    %'mafin_04_t099',         'h.l.=68';   % phi_base=0.097, AR_ccyb=0.99, peakin=1  
                    
					%'MAFIN_noFF_reestim',     'MAFIN NoFF(reestim)';
                    %'mafin_04_compare0',          'market';    % phi_base=0.097, AR_ccyb=0.5, peakin=1 
                    %'mafin_04_compare1',          'book';    % phi_base=0.097, AR_ccyb=0.5, peakin=1 
                    %'mafin_04_compare2',          'book2';    % phi_base=0.097, AR_ccyb=0.5, peakin=1 
                    
                  };                 
    
    %                Name        Latex name        shock scale      long_name
    Impulse = 	{	%'z',            'z',                1,          'Transitory Technology';
                    %'a',            'a',				1,          'Permanent Technology';
                    %'xi_i',     	'\xi^{i}',      	1,         	'Capital Investment Efficiency';
                    %'xi_ih',		'\xi^{ih}',			1,         	'Housing Investment Efficiency';
                    %'varrho',		'\varrho',			1,      	'Preference';
                    %'xi_n',      	'\xi^{n},\kappa',   1,       	'Labor Disutility';		
                    %'sigma_ii',    	'\sigma^{ee}',   	1,         	'risk shock(households)';
                    %'sigma_ee',    	'\sigma^{ee}',   	1,         	'risk shock(entrepreneurs)';
                    %'sigma_ff',    	'\sigma^{ee}',   	1,         	'risk shock(banks firms)';
                    %'sigma_hh',    	'\sigma^{ee}',   	1,         	'risk shock(banks housing)';
                    'e_m',      	'e^{m}',			1000/139/0.8,        	'Monetary Policy';
                    %'e_req',      	'e^{RCC}',			10,        	'Financial Policy';
                    %'g',  			'g',                1,      	'Government Consumption';
                    %'y_ast',		'y^{*}',			1,       	'Foreign Output';
                    %'ppi_ast',		'\pi^{*}',          1,       	'Foreign Inflation';
                    %'R_W',			'R^{W}',			1,        	'Foreign Interest Rate';
                    %'xi_m',        	'\xi^{m}',          1,        	'Import Prices';
                    %'xi_R',			'\xi^{R}',          1,        	'Country Risk Premium';		
                    %'p_co',        	'p^{Co*}',          1,         	'Commodity Price';
                    %'y_co',			'y^{Co}',           1,         	'Commodity Output';
					%'bl_g',         	'BL^{G}',           100,          'Long-term government bonds supply';
					%'bl_cb',            'BL^{BC}',         1,          'Long term government bonds purchases by the CB';
					%'e_m',              'BL^{BC}',         1,          'Monetary policy shock';
					%'epsilon_L', 	 '\epsilon_L',         1,          'Transaction costs shock';
					%'xi_roe_r',  	 '\xi_{ROE}',         1,          'Expected returns ratio';
					%'xi_CHI_b',  	 '\xi_{\chi_{B}}',         1,          'Dividend policy banks';
                    %'xi_CHI_e',  	 'BL^{\chi_{B}}',         1,          'Dividend policy entrepreneurs';
			};
      
    
    %                Name           Latex Name              Cusum	Factor      Div_ss
    Response ={  	%'ccyb',         'CCyB\: (\%)',          0,		1,          0;
                                %'gam_YR_obs',	'GDP',              1,		1/100,      0;
                    %'ltot',         'Credito\: (\%)',       0,		1,          1;            
                    %'l_f',         'Credito F\: (\%)',       0,		1,          1;            
                    'gdp',          'PIB\: (\%)',           0,		1,          1;
                    %'gdpn',          'PIBN\: (\%)',       0,		1,          1;
                    %'c' ,            'Consumo\: (\%)',  0,		1,          1;
                    %%'gam_C_obs',     'C',               1,		1/100,		0;
                    'i_agg',        'I+I^{H}\: (\%)',       0,		1,          1;
                    'ppi',          '\pi',              0,		4,          0;
                    %'nu',           'perdida\: (\%)',       0,		1,          0;
                    %'h_i' ,         'h^I\: (\%)',           0,		1,          1;
                    %'h_p' ,         'h^P\: (\%)',           0,		1,          1;
                    %'ltot',          'ltot',		0,		1/100,		1;   
                    'i' ,            'I',               0,		1,          1;
                    %'i_h' ,          'I^{H}',           0,		1,          1;
                    
                    'R',			'R\: (\%)',             0,		1,          0;
                    'R_i',          'R^i\: (\%)',           0,		1,          0;
                    'R_L',          'R^L\: (\%)',           0,		1,          0;
                    'PD_e',         'PD^e\: (\%)',          0,		1,          0;
                    %'PD_i',           'PD^I\: (\%)',           0,		1,          0;
                    'n_e',           'Entr NW',        0,		1,          1;
                    

                    
                    'spread_RL_RD','spread RL RD',	0,		1,          0;
                    %'spread_RD_R',	'spread RD R',	0,		4,          0;
                    %%'R_f_tilde',		'\tilde{R}^F',	0,		1,          1;
                    %%'R_h_tilde',		'\tilde{R}^H',		0,		1,          1;
                    %'Gam_G_e',			'\Gamma_e-\mu G_e',           0,		1/100,          0;
                    %'Gam_G_i',			'\Gamma_e-\mu G_i',           0,		1/100,          0;
                    %'tb',			'TB',			0,		1,          1;
                    %'gam_ltot_obs','Loans',		1,		1/100,		0;
                    'l_f',          'Corp.Loans\: (\%)',       0,      1,          1;
                    %'mortVpar_tot', 'Mort.Loans\: (\%)',    0,      1,          1;
                    %'mort_tot', 'Mort.Loans',       0,      1,          1;
                    'q_k',           'q^K',			0,		1,          1;
                    'k',           'K',			0,		1,          1;
                    %'q_h',           'q^H\: (\%)',			0,		1,          1;
                    %'k',           'k',			0,		1,          1;
                    %'lev_tot',       'Nb \: (\% RWA_{SS})',	0,		1,          0;
                    'R_e',			'R^e',			0,		1,          0;
                    'R_h',			'R^h',			0,		1,          0;
                    %'r_k',         'R^{k}',     	0,		1,          0;
                    %'rer',			'rer',          0,		1,          1;
                    %'ppi_S',       'S',			1,		1,          0;
			};
    
    plot_options.n_col      = 4;
    plot_options.n_row      = 4;
    plot_options.LineStyle	= {'-','--','-.','--'};
    plot_options.LineColor  = {[0 0 0.9],[0.9 0 0],[0 0.8 0],[0 0 0]};
    plot_options.LineWidth  = {2, 1.5, 1.5, 1.5};
    plot_options.horizon    = 40;
    plot_options.grid       = 1;
    plot_options.latex      = 1;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot_options.model_names= model_names(:,1)';
legend_names            = model_names(:,2)';

shocks                  = Impulse(:,1)';
shocks_latex            = Impulse(:,2)';
shocks_scale            = Impulse(:,3)';
shocks_long_name        = Impulse(:,4)';

v_sel                   = Response(:,1)';
v_name                  = Response(:,2)';
plot_options.v_do_cusum = cell2mat(Response(:,3)');
plot_options.v_adj      = cell2mat(Response(:,4)');
plot_options.v_div_ss   = cell2mat(Response(:,5)');

plot_options.v_do_cusum = [plot_options.v_do_cusum 	0];  	% shock doesn't sum over
plot_options.v_div_ss   = [plot_options.v_div_ss 	1];    	% shock divided by its ss
plot_options.v_adj_orig = [plot_options.v_adj  		1]*100;	% no adjustment for shock

if exist(LOC_final, 'file')==7
    rmdir(LOC_final,'s')
end
mkdir(LOC_final)

for ii=1:length(plot_options.model_names)
    %disp(['loading data for ' plot_options.model_names{ii} ' (' num2str(ii) '/'  num2str(length(plot_options.model_names)) ')']);
    
    %name_mat=[plot_options.model_names{ii} '.mat'];
    %copyfile(['..\' name_mat],name_mat)
    
    load([plot_options.model_names{ii} '.mat'])
    eval(['data.M_.' plot_options.model_names{ii} '=M_;'])
    eval(['data.oo_.' plot_options.model_names{ii} '=oo_;'])
end


for j=1:size(shocks,2)
    %disp(['plotting IRFs for ' shocks{j} ' (' num2str(j) '/'  num2str(size(shocks,2)) ')']);
    plot_options.u_sel=shocks(j);
    if size(plot_options.model_names,2)>1
        for i=2:size(plot_options.model_names,2)
            plot_options.u_sel=[plot_options.u_sel shocks(j)];
        end
    end
    plot_options.v_sel  = [v_sel shocks(j)];
    plot_options.u_name = shocks_latex(j);
    plot_options.v_name = [ v_name shocks_latex(j)];
    plot_options.v_adj  = shocks_scale{j} * plot_options.v_adj_orig;
    
    plots_for_dynare(plot_options,data)
    
    h=legend(legend_names);
    set(h,'Position',[0.4 0 0.23 0.03],...
        'Orientation','Horizontal',...
        'Fontsize',12,...
        'Interpreter','Latex');
    legend('boxoff')
    
    
    fig_ratio=16/9;
    trim_hor=0.6;
    trim_top=0.05;
    % saveplot(['irf_' shocks{j}],fig_ratio,trim_hor,trim_top,LOC_final)
    
end

% for ii=1:length(plot_options.model_names)
%     name_mat=[plot_options.model_names{ii} '.mat'];
%     delete(name_mat)
% end


% Generate pdf file with all the graphs
latex_filename = [DIR '.tex'];

LATEX_LINE = ['\\documentclass{article}', '\n\n'];
LATEX_LINE = [LATEX_LINE, '\\usepackage[utf8]{inputenc}', '\n'];
LATEX_LINE = [LATEX_LINE, '\\usepackage{geometry}', '\n'];
LATEX_LINE = [LATEX_LINE, '\\geometry{verbose,tmargin=2cm,bmargin=2cm,lmargin=2cm,rmargin=2cm}', '\n'];
LATEX_LINE = [LATEX_LINE, '\\usepackage{amsmath}', '\n'];
LATEX_LINE = [LATEX_LINE, '\\usepackage{amssymb}', '\n'];
LATEX_LINE = [LATEX_LINE, '\\usepackage{graphicx}', '\n'];
LATEX_LINE = [LATEX_LINE, '\\makeatletter', '\n\n'];
LATEX_LINE = [LATEX_LINE, '\\begin{document}'];

LATEX_LINE = [LATEX_LINE,	'\n\n'];

for ii=1:length(shocks(1,:))
    NAME_FIG   = ['irf_' shocks{1,ii}];
    TITLE_FIG  = shocks_long_name{1,ii};
    
    LATEX_LINE = [LATEX_LINE,	'\\begin{figure}'];
    LATEX_LINE = [LATEX_LINE,	'\n',       '\\caption{' TITLE_FIG ' Shock}'];
    LATEX_LINE = [LATEX_LINE,	'\n',       '\\begin{center}'];
    LATEX_LINE = [LATEX_LINE,	'\n',       '\\includegraphics{' NAME_FIG '.pdf}'];
    LATEX_LINE = [LATEX_LINE,	'\n',       '\\end{center}'];
    LATEX_LINE = [LATEX_LINE,	'\n',       '\\end{figure}'];
    
    LATEX_LINE = [LATEX_LINE,	'\n\n'];
end

LATEX_LINE = [LATEX_LINE, '\\end{document}'];

if exist(latex_filename,'file')==2
    delete(latex_filename);
end
fid = fopen(latex_filename,'at');    %%% open the file
fprintf(fid,LATEX_LINE);             %%% append the extra code
fclose(fid);                         %%% close the file (the new tex file is generated)
movefile(latex_filename,LOC_final)

%addpath('C:\Program Files\MiKTeX 2.9-1\miktex\bin')
%[status,cmdout] = system(['pdflatex ' latex_filename]); %%% compile the tex file into a pdf
%rmpath('C:\Program Files\MiKTeX 2.9-1\miktex\bin')
disp(['Done for ' ModelNames{q}])
end