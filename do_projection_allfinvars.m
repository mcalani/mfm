clear
clc
close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% SETTING  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('gathering data')
% Import data from excel file
excel.input_file    = 'Data_MAFIN_2023Q1.xlsx'; % Excel file containing raw data
excel.input_sheet   = 'DataRawSA';            % Sheet of the excel file that contains raw data  
excel.input_range   = 'A1:AW105';              % Range of the raw data (include dates)

% Export data to excel file
excel.output_file   = excel.input_file;   % Excel file where output will be saved 

% Defining sample (typically: 2001.5 to "last quarter available")
sample.first    = 2000;   % Initial sample date 
sample.last     = 2025.75;   % Last sample date 

% Defining variables to transform (First: excel names, Second: Dynare names without "_obs")


vars	= {%	VarName     2nd Argument    Transformation  New Name        Detrend method
                'yrM'       ''              'dlog'          'gam_YR'        'demean'
                'yM'        ''              'dlog'          'gam_YCo'       'demean'
                'cp'        ''              'dlog'          'gam_C'         'demean'
                'i'         ''              'dlog'          'gam_I'         'demean'
                'ih'        ''              'dlog'          'gam_IH'        'demean'
                'iK'        ''              'dlog'          'gam_IK'        'demean'
                'cG'        ''              'dlog'          'gam_G'         'demean'
                'ystar'     ''              'dlog'          'gam_Ystar'     'demean'
                'wn'        ''              'dlog'          'gam_WN'        'demean'
                'p'         ''              'dlog'          'pi'            'demean'
                'pZ'         ''             'dlog'          'piZ'           'demean'
                'pM'        ''              'dlog'          'piM'           'demean'
                'pstar'     ''              'dlog'          'pistar'        'demean'
                'pCostar'   ''              'dlog'          'piCostar'      'demean'
                'qh'        ''              'dlog'          'piQH'          'demean'             
                'lf'        ''              'dlog'          'gam_lf'        'demean'
                'lh'        ''              'dlog'          'gam_lh'        'demean'
                'ltot'      ''              'dlog'          'gam_ltot'      'demean'
                'ytot'       ''             'dlog'          'gam_Y'         'demean'
           %    VarName     Base            Transformation  New Name        Detrend method
                'R'         '100'           'qrate'         'R'             'demean'
                'Rstar'     '100'           'qrate'         'Rstar'         'demean'
                'RLG'       '100'           'qrate'         'RLG'           'demean'
                'RD'        '100'           'qrate'         'RD'            'demean'
                'RL'        '100'           'qrate'         'RL'            'demean'
                'RI'        '100'           'qrate'         'RI'            'demean'
                'xi'        '10000'         'qrate'         'xi'            'demean'
                'div'       ''              ''              'div'           'demean' 
           %    VarName     Base            Transformation  New Name        Detrend method     
                'roe'       '100'           'qrate'         'Roe'           'demean'
           %    Numerator   Denominator     Transformation  New Name        Detrend method    
                'tb'        'yN'            'ratio'         'stb'           'demean'
           %    Numerator   Denominator     Transformation  New Name        Detrend method  
                'ntot_X_h'  'ft'            'dlogratio2'    'gam_N'         'demean'
           %    VarName     2nd Argument    Transformation  New Name        Detrend method          
                'rer'       ''              'loglevel'      'rer'           'demean'
                'tcn'       ''              'dlog'          'tcn'           'demean'
          };


					   				   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% Do not change anything from here onwards %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% COMPUTATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('making computations')
%%%%%%%%%%%%%%%% Creating Data Sample %%%%%%%%%%%%%%%%
% Import data from excel file
[~, ~, raw] = xlsread(excel.input_file,excel.input_sheet,excel.input_range);
raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {NaN};

% Assigning Variables
variable_values = reshape([raw{2:end,:}],size(raw(2:end,:)));
variable_names  = raw(1,:);
for j=1:size(variable_names,2)
	eval([variable_names{j} '_orig= variable_values(:,j);']);
end 

%%%%%%%%%%%%%%%% Transforming Data %%%%%%%%%%%%%%%%
for j=1:size(vars,1)
    if strcmp(vars{j,3},'dlog')
        eval([vars{j,4} '=100*log(' vars{j,1} '_orig(2:end)./' vars{j,1} '_orig(1:end-1));']);
        eval([vars{j,4} '=[nan; ' vars{j,4} '];']);   
    elseif strcmp(vars{j,3},'ydlog')
        eval([vars{j,4} '=100*log(('...
                                    vars{j,1} '_orig(8:end)  +'...
                                    vars{j,1} '_orig(7:end-1)+'...
                                    vars{j,1} '_orig(6:end-2)+'...
                                    vars{j,1} '_orig(5:end-3)'...
                                  ')./('...
                                    vars{j,1} '_orig(4:end-4)+'...
                                    vars{j,1} '_orig(3:end-5)+'...
                                    vars{j,1} '_orig(2:end-6)+'...
                                    vars{j,1} '_orig(1:end-7)'...
                                    '));'
                                    ]);
        eval([vars{j,4} '=[nan(7,1); ' vars{j,4} '];']);        
        for jj=1:numel(dates_orig)
            if dates_orig(jj)-floor(dates_orig(jj))~=0.75             
            eval([vars{j,4} '(jj)=nan;']);
            end
        end
    elseif strcmp(vars{j,3},'qrate')        
        eval([vars{j,4} '=100*0.25*log(1+' vars{j,1} '_orig/' vars{j,2} ');']);  
     elseif strcmp(vars{j,3},'rate')        
        eval([vars{j,4} '=100*log(1+' vars{j,1} '_orig/' vars{j,2} ');']);    
    elseif strcmp(vars{j,3},'ratio')       
        eval([vars{j,4} '=100*(' vars{j,1} '_orig./' vars{j,2} '_orig);']);   
    elseif strcmp(vars{j,3},'logratio')      
        eval([vars{j,4} '=100*log(' vars{j,1} '_orig./' vars{j,2} '_orig);']);
    elseif strcmp(vars{j,3},'dlogratio')      
        eval([vars{j,4} '=100*log(' vars{j,1} '_orig(2:end)./' vars{j,1} '_orig(1:end-1))-100*log(' vars{j,2} '_orig(2:end)./' vars{j,2} '_orig(1:end-1));']);
        eval([vars{j,4} '=[nan; ' vars{j,4} '];']);  
    elseif strcmp(vars{j,3},'dlogratio2')
        pos_X_=strfind(vars{j,1},'_X_');
        numerator_string1=vars{j,1}(1:pos_X_-1);
        numerator_string2=vars{j,1}(pos_X_+3);       
        eval([vars{j,4} '=' '100*log(' numerator_string1 '_orig(2:end)./' numerator_string1 '_orig(1:end-1))' '+100*log(' numerator_string2 '_orig(2:end)./' numerator_string2 '_orig(1:end-1))' '-100*log(' vars{j,2} '_orig(2:end)./' vars{j,2} '_orig(1:end-1));']);
        eval([vars{j,4} '=[nan; ' vars{j,4} '];']);  
    elseif strcmp(vars{j,3},'loglevel')       
        eval([vars{j,4} '=100*log(' vars{j,1} '_orig);']);  
    elseif strcmp(vars{j,3},'')       
        eval([vars{j,4}  '=' vars{j,1} '_orig;']);
    else
        error(['unknown transformation "' vars{j,3} '" for variable "' vars{j,1} '".']  ) 
    end
end

%%%%%%%%%%%%%%%% Computing Trend Variables %%%%%%%%%%%%%%%%
% Modifying sample
y0   = find(dates_orig(:,1)==sample.first);
yend = find(dates_orig(:,1)==sample.last);
dates = dates_orig(y0:yend);
new_names = vars(:,4)';
for j=1:numel(new_names)
    eval(['var_aux =' new_names{j} ';']);
	eval([new_names{j} '= var_aux(y0:yend);']);
end


for j=1:size(vars,1) %%aux1=orig_var %%aux2=missing data %%aux3=observed data %%aux4=detrended data
    var_aux1 = eval(vars{j,4});
    var_aux2 = isnan(var_aux1);
    var_aux3 = var_aux1(var_aux2==0);
    if strcmp(vars{j,5},'linear_trend')       
        if size(var_aux1,1) ~= size(var_aux3,1)
            var_aux4a = detrend(var_aux3,'linear'); var_aux4=var_aux1; idx=find(var_aux2==0); var_aux4(idx)=var_aux4a; 
        else
            var_aux4 = detrend(var_aux3,'linear');
            eval([vars{j,4} '_T = var_aux1-var_aux4;']);
        end
    elseif strcmp(vars{j,5},'demean')
        if size(var_aux1,1) ~= size(var_aux3,1)            
            var_aux4a = detrend(var_aux3,'constant'); var_aux4=var_aux1; idx=find(var_aux2==0); var_aux4(idx)=var_aux4a;                
        else
            var_aux4 = detrend(var_aux3,'constant');            
        end 
        eval([vars{j,4} '_T = var_aux1-var_aux4;']);
    elseif strcmp(vars{j,5},'')
            eval([vars{j,4} '_T = var_aux1*0;']);
    else
    error(['unknown detrending method "' vars{j,5} '" for variable "' vars{j,4} '".']  )    
    end
end


%%%%%%%%%%%%%%%% Saving results %%%%%%%%%%%%%%%%
disp('saving')
var_names  = new_names;
for j=1:numel(var_names)
var_names_T{j}=[var_names{j} '_T'];
var_names_obs{j}=[var_names{j} '_obs'];
end
var_values = NaN(size(dates,1),size(var_names,2));
var_mean   = var_values;

for j=1:size(var_names,2)
    eval(['var_values(:,j) =' var_names{j}   ';']);
    eval(['var_mean(:,j)   =' var_names_T{j} ';']);
    eval(['var_demean(:,j) =' var_names{j}   '-' var_names_T{j} ';']);
end

% Saving to excel file
% Writing dates in each sheet of the excel file
xlswrite(excel.output_file,{'dates'},'Data','A1')
xlswrite(excel.output_file,{'dates'},'DataMean','A1')
xlswrite(excel.output_file,{'dates'},'DataDemean','A1')
xlswrite(excel.output_file,dates,'Data','A2')
xlswrite(excel.output_file,dates,'DataMean','A2')
xlswrite(excel.output_file,dates,'DataDemean','A2')

% Writing variables names in each sheet of the excel file
xlswrite(excel.output_file,var_names,     'Data','B1')
xlswrite(excel.output_file,var_names_T,   'DataMean','B1')
xlswrite(excel.output_file,var_names_obs, 'DataDemean','B1')

% Writing variables vlues in each sheet of the excel file
xlswrite(excel.output_file,var_values,'Data','B2')
xlswrite(excel.output_file,var_mean,'DataMean','B2')
xlswrite(excel.output_file,var_demean,'DataDemean','B2')

clearvars -except var_names var_values var_mean var_demean 
%% Ejecución Dynare 

disp('Making Kalman Filter')
dynare MAFIN_03_est.mod                 % Calculo del filtro de Kalman 

disp('Saving Results')

global oo_ M_
%oo_.FilteredVariables.PD_d
Resulst = [oo_.FilteredVariables.l oo_.FilteredVariables.l_h oo_.FilteredVariables.l_f oo_.FilteredVariables.q_hat_l ...
    oo_.FilteredVariables.R_i oo_.FilteredVariables.R_L oo_.FilteredVariables.PD_f oo_.FilteredVariables.PD_d ...
    oo_.FilteredVariables.PD_h oo_.FilteredVariables.a oo_.FilteredVariables.ppi ];
for i=1:length(Resulst(:,1))-12
    Resulst(1,:)=[];
end
xlswrite('DPM_Proyection.xlsx',Resulst,'mfm_forecast','B2');

std_ltot = (oo_.forecast.HPDsup.l-oo_.forecast.HPDinf.l)/2;
std_R_i = ((oo_.forecast.HPDsup.R_i-oo_.forecast.HPDinf.R_i)/2)*100;
std_R_l = (oo_.forecast.HPDsup.R_L-oo_.forecast.HPDinf.R_L)/2;
std_PD_d = (oo_.forecast.HPDsup.PD_d-oo_.forecast.HPDinf.PD_d)/2;
std_PD_f = (oo_.forecast.HPDsup.PD_f-oo_.forecast.HPDinf.PD_f)/2;
std_PD_h = (oo_.forecast.HPDsup.PD_h-oo_.forecast.HPDinf.PD_h)/2;

res_std = [std_ltot std_R_i std_R_l std_PD_d std_PD_f std_PD_h];
xlswrite('DPM_Proyection.xlsx',res_std,'mfm_forecast','N3');

xlswrite('DPM_Proyection.xlsx',oo_.FilteredVariables.PD_d,'Figuras','AC3');

params_est = [M_.params];
%xlswrite('DPM_Proyection.xlsx',params_est,'Parameters','D2');

%save_results

