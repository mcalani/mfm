
%   CONSTRUCCION DE ESCENARIOS PARA MERCADO DE CREDITO SEGUN ESCENARIOS FORECAST IPOM 
%   *********************************************************************************

%   0. CLEAN HOUSE 
%   Asegurarse que 
%
%       i. Hoja excel que se llama del .mod sea la correcta y tenga data
%       desestacionalizada para todas las variables + proyeccion para
%       series enviadas por la GAM (usualmente pib inc. minería, ner,
%       inflación, tpm) 
%
%       ii. Se usa estimacion_rpf.inc en vez de estimacion.inc 




%   1. CORRER MODELO ESCENARIO BASE 
close all 
disp('Making Kalman Filter')
dynare MAFIN_03_est.mod                 
disp('Saving Results')

sample.first    = 2000;   % Initial sample date 
sample.last     = 2025.75;   % Last sample date 
timeline = linspace(sample.first,sample.last,(sample.last-sample.first)*4+1);

global oo_ M_

% Check prob default sistema bancario 
figure; 
plot(timeline,100* oo_.FilteredVariables.PD_d,'b-','LineWidth',2); 
title('Prob. default banca')
hold on; 
plot(timeline,100* oo_.SmoothedVariables.PD_d,'r-','LineWidth',1.5);
xline(2008.5,'--',{'Lehman Bro'})
xline(2019.5,'--',{'2019Q3'})
xline(2023.25,'-',{'Forecast'})
xlim([sample.first+6  sample.last]); 
ylim([0 2.5])
hold off

figure; 
plot(timeline,oo_.FilteredVariables.gam_Y_obs,'b-','LineWidth',2); 
title('Prob. default banca')
hold on; 
xline(2008.5,'--',{'Lehman Bro'})
xline(2019.5,'--',{'2019Q3'})
xline(2023.25,'-',{'Forecast'})
xlim([sample.first+6  sample.last]); 
ylim([-2 2.5])
hold off

figure; 
plot(timeline,100* oo_.FilteredVariables.PD_h,'b-','LineWidth',2); 
title('Prob. default banca')
hold on; 
xline(2008.5,'--',{'Lehman Bro'})
xline(2019.5,'--',{'2019Q3'})
xline(2023.25,'-',{'Forecast'})
xlim([sample.first+6  sample.last]); 
ylim([0 7])
hold off



figure; plot(timeline,oo_.FilteredVariables.R_i,timeline,oo_.SmoothedVariables.R_i,'g')
figure; plot(timeline,oo_.SmoothedVariables.R_i,timeline,(oo_.SmoothedVariables.RI_obs/100+1),'g')




%   2. 	IMPRIMIR RESULTADOS EN HOJA DPM_Proyection 
Resulst = [oo_.SmoothedVariables.l oo_.SmoothedVariables.l_h oo_.SmoothedVariables.l_f oo_.SmoothedVariables.q_hat_l ...
    oo_.SmoothedVariables.R_i oo_.SmoothedVariables.R_L oo_.SmoothedVariables.PD_f oo_.SmoothedVariables.PD_d ...
    oo_.SmoothedVariables.PD_h oo_.SmoothedVariables.a oo_.SmoothedVariables.ppi ];
for i=1:length(Resulst(:,1))-12
    Resulst(1,:)=[];
end
xlswrite('DPM_Proyection.xlsx',Resulst,'mfm_forecast','B2');

std_ltot = (oo_.forecast.HPDsup.l-oo_.forecast.HPDinf.l)./(oo_.forecast.HPDsup.l+oo_.forecast.HPDinf.l)*0.5*100;
std_R_i = ((oo_.forecast.HPDsup.R_i)-(oo_.forecast.HPDinf.R_i))*400;
std_R_L = ((oo_.forecast.HPDsup.R_L)-(oo_.forecast.HPDinf.R_L))*400;
std_PD_d = (oo_.forecast.HPDsup.PD_d-oo_.forecast.HPDinf.PD_d)/2*100;
std_PD_f = (oo_.forecast.HPDsup.PD_f-oo_.forecast.HPDinf.PD_f)/2*100;
std_PD_h = (oo_.forecast.HPDsup.PD_h-oo_.forecast.HPDinf.PD_h)/2*100;

    % imprimir desvest forecast 
res_std = [std_ltot std_R_i std_R_L std_PD_d std_PD_f std_PD_h];
xlswrite('DPM_Proyection.xlsx',res_std,'mfm_forecast','B34');

    % imprimir historia entera de inobservable PD_D 
xlswrite('DPM_Proyection.xlsx',oo_.FilteredVariables.PD_d,'Figuras','AC3');
xlswrite('DPM_Proyection.xlsx',oo_.SmoothedVariables.PD_d,'Figuras','AO3');




%   3. 	ESCENARIO DE RIESGO IPOM 
%   Asegurarse que en "estimation.rpf" se lee .XLS "Data_MAFIN_2023Q1_R1"
%   que debe tener cargado el escenario de riesgo en vez del central 

close all 
disp('Making Kalman Filter')
dynare MAFIN_03_est.mod             % (!) Changed to Data_MAFIN_2023Q1_R1 ?


Resulst = [oo_.SmoothedVariables.l oo_.SmoothedVariables.l_h oo_.SmoothedVariables.l_f oo_.SmoothedVariables.q_hat_l ...
    oo_.SmoothedVariables.R_i oo_.SmoothedVariables.R_L oo_.SmoothedVariables.PD_f oo_.SmoothedVariables.PD_d ...
    oo_.SmoothedVariables.PD_h oo_.SmoothedVariables.a oo_.SmoothedVariables.ppi ];
for i=1:length(Resulst(:,1))-12
    Resulst(1,:)=[];
end
xlswrite('DPM_Proyection.xlsx',Resulst,'mfm_forecast_R1','B2');

std_ltot = (oo_.forecast.HPDsup.l-oo_.forecast.HPDinf.l)./(oo_.forecast.HPDsup.l+oo_.forecast.HPDinf.l)*0.5*100;
std_R_i = ((oo_.forecast.HPDsup.R_i)-(oo_.forecast.HPDinf.R_i))*400;
std_R_L = ((oo_.forecast.HPDsup.R_L)-(oo_.forecast.HPDinf.R_L))*400;
std_PD_d = (oo_.forecast.HPDsup.PD_d-oo_.forecast.HPDinf.PD_d)/2*100;
std_PD_f = (oo_.forecast.HPDsup.PD_f-oo_.forecast.HPDinf.PD_f)/2*100;
std_PD_h = (oo_.forecast.HPDsup.PD_h-oo_.forecast.HPDinf.PD_h)/2*100;

    % imprimir desvest forecast 
res_std = [std_ltot std_R_i std_R_L std_PD_d std_PD_f std_PD_h];
xlswrite('DPM_Proyection.xlsx',res_std,'mfm_forecast_R1','B34');

    % imprimir historia entera de inobservable PD_D 
xlswrite('DPM_Proyection.xlsx',oo_.FilteredVariables.PD_d,'mfm_forecast_R1','AC3');
%params_est = [M_.params];
%xlswrite('DPM_Proyection.xlsx',params_est,'Parameters','D2');






%   4. 	Cambios en CCYB
%   Asegurarse que en "estimation.rpf" se lee .XLS "Data_MAFIN_2023Q1_Pol"
%   que debe tener cargado el escenario de ccyb 

%   i) Escenario 0
close all 
disp('Making Kalman Filter')
dynare mafin_04.mod             % (!) Changed to Data_MAFIN_2023Q1_Pol, changed ccyb_obs accordingly 

    %Recuperar series 
xlswrite('simCCyB.xlsx',oo_.SmoothedVariables.ltot ,'_ltot','C4');
xlswrite('simCCyB.xlsx',oo_.SmoothedVariables.ppi ,'_ppi','C4');
xlswrite('simCCyB.xlsx',oo_.SmoothedVariables.PD_d ,'_pdd','C4');
xlswrite('simCCyB.xlsx',oo_.SmoothedVariables.gdpR ,'_gdpNoM','C4');
xlswrite('simCCyB.xlsx',oo_.SmoothedVariables.gdp  ,'_gdp','C4');
xlswrite('simCCyB.xlsx',oo_.SmoothedVariables.R  ,'_tpm','C4');
xlswrite('simCCyB.xlsx',oo_.SmoothedVariables.ccyb  ,'_ccyb','C4');
xlswrite('simCCyB.xlsx',oo_.SmoothedVariables.i_agg  ,'_invA','C4');
xlswrite('simCCyB.xlsx',oo_.SmoothedVariables.i  ,'_invK','C4');
xlswrite('simCCyB.xlsx',oo_.SmoothedVariables.c  ,'_c','C4');
xlswrite('simCCyB.xlsx',oo_.SmoothedVariables.a  ,'_a','C4');
xlswrite('simCCyB.xlsx',oo_.SmoothedVariables.spread_RL_RD  ,'_spread_RL_RD','C4');
xlswrite('simCCyB.xlsx',oo_.SmoothedVariables.spread_Ri_RD  ,'_spread_Ri_RD','C4');
xlswrite('simCCyB.xlsx',oo_.SmoothedVariables.q_h  ,'_q_h','C4');
xlswrite('simCCyB.xlsx',oo_.SmoothedVariables.R_i  ,'_Ri','C4');
xlswrite('simCCyB.xlsx',oo_.SmoothedVariables.spread_Ri_RBL  ,'_spread_Ri_RBL','C4');


    % Recuperar series para ccyb positivo 50pb 
xlswrite('simCCyB.xlsx',oo_.SmoothedVariables.l_f ,'_l_f','D4');
xlswrite('simCCyB.xlsx',oo_.SmoothedVariables.ltot ,'_ltot','D4');
xlswrite('simCCyB.xlsx',oo_.SmoothedVariables.ppi ,'_ppi','D4');
xlswrite('simCCyB.xlsx',oo_.SmoothedVariables.PD_d ,'_pdd','D4');
xlswrite('simCCyB.xlsx',oo_.SmoothedVariables.gdpR ,'_gdpNoM','D4');
xlswrite('simCCyB.xlsx',oo_.SmoothedVariables.gdp  ,'_gdp','D4');
xlswrite('simCCyB.xlsx',oo_.SmoothedVariables.R  ,'_tpm','D4');
xlswrite('simCCyB.xlsx',oo_.SmoothedVariables.ccyb  ,'_ccyb','D4');
xlswrite('simCCyB.xlsx',oo_.SmoothedVariables.i_agg  ,'_invA','D4');
xlswrite('simCCyB.xlsx',oo_.SmoothedVariables.i  ,'_invK','D4');
xlswrite('simCCyB.xlsx',oo_.SmoothedVariables.c  ,'_c','D4');
xlswrite('simCCyB.xlsx',oo_.SmoothedVariables.a  ,'_a','D4');
xlswrite('simCCyB.xlsx',oo_.SmoothedVariables.spread_RL_RD  ,'_spread_RL_RD','D4');
xlswrite('simCCyB.xlsx',oo_.SmoothedVariables.spread_Ri_RD  ,'_spread_Ri_RD','D4');
xlswrite('simCCyB.xlsx',oo_.SmoothedVariables.q_h  ,'_q_h','D4');
xlswrite('simCCyB.xlsx',oo_.SmoothedVariables.R_i  ,'_Ri','D4');
xlswrite('simCCyB.xlsx',oo_.SmoothedVariables.spread_Ri_RBL  ,'_spread_Ri_RBL','D4');


    % Recuperar series para ccyb positivo 100pb 
xlswrite('simCCyB.xlsx',oo_.SmoothedVariables.ltot ,'_ltot','E4');
xlswrite('simCCyB.xlsx',oo_.SmoothedVariables.ppi ,'_ppi','E4');
xlswrite('simCCyB.xlsx',oo_.SmoothedVariables.PD_d ,'_pdd','E4');
xlswrite('simCCyB.xlsx',oo_.SmoothedVariables.gdpR ,'_gdpNoM','E4');
xlswrite('simCCyB.xlsx',oo_.SmoothedVariables.gdp  ,'_gdp','E4');
xlswrite('simCCyB.xlsx',oo_.SmoothedVariables.R  ,'_tpm','E4');
xlswrite('simCCyB.xlsx',oo_.SmoothedVariables.ccyb  ,'_ccyb','E4');
xlswrite('simCCyB.xlsx',oo_.SmoothedVariables.i_agg  ,'_invA','E4');
xlswrite('simCCyB.xlsx',oo_.SmoothedVariables.i  ,'_invK','E4');
xlswrite('simCCyB.xlsx',oo_.SmoothedVariables.c  ,'_c','E4');
xlswrite('simCCyB.xlsx',oo_.SmoothedVariables.a  ,'_a','E4');
xlswrite('simCCyB.xlsx',oo_.SmoothedVariables.spread_RL_RD  ,'_spread_RL_RD','E4');
xlswrite('simCCyB.xlsx',oo_.SmoothedVariables.spread_Ri_RD  ,'_spread_Ri_RD','E4');
xlswrite('simCCyB.xlsx',oo_.SmoothedVariables.q_h  ,'_q_h','E4');
xlswrite('simCCyB.xlsx',oo_.SmoothedVariables.R_i  ,'_Ri','E4');
xlswrite('simCCyB.xlsx',oo_.SmoothedVariables.spread_Ri_RBL  ,'_spread_Ri_RBL','E4');


    % Recuperar series para excenario central 
xlswrite('ECentral.xlsx',oo_.SmoothedVariables.ltot ,'_ltot','D4');
xlswrite('ECentral.xlsx',oo_.SmoothedVariables.ppi ,'_ppi','D4');
xlswrite('ECentral.xlsx',oo_.SmoothedVariables.PD_d ,'_pdd','D4');
xlswrite('ECentral.xlsx',oo_.SmoothedVariables.gdpR ,'_gdpNoM','D4');
xlswrite('ECentral.xlsx',oo_.SmoothedVariables.gdp  ,'_gdp','D4');
xlswrite('ECentral.xlsx',oo_.SmoothedVariables.R  ,'_tpm','D4');
xlswrite('ECentral.xlsx',oo_.SmoothedVariables.ccyb  ,'_ccyb','D4');
xlswrite('ECentral.xlsx',oo_.SmoothedVariables.i_agg  ,'_invA','D4');
xlswrite('ECentral.xlsx',oo_.SmoothedVariables.i  ,'_invK','D4');
xlswrite('ECentral.xlsx',oo_.SmoothedVariables.c  ,'_c','D4');
xlswrite('ECentral.xlsx',oo_.SmoothedVariables.spread_RL_RD  ,'_spread_RL_RD','D4');
xlswrite('ECentral.xlsx',oo_.SmoothedVariables.spread_Ri_RD  ,'_spread_Ri_RD','D4');
xlswrite('ECentral.xlsx',oo_.SmoothedVariables.q_h  ,'_q_h','D4');
xlswrite('ECentral.xlsx',oo_.SmoothedVariables.R_i  ,'_Ri','D4');
xlswrite('ECentral.xlsx',oo_.SmoothedVariables.spread_Ri_RBL  ,'_spread_Ri_RBL','D4');


%   i) Escenario 1 
close all 
disp('Making Kalman Filter')
dynare mafin_04.mod             % (!) Changed to Data_MAFIN_2023Q1_Pol, changed ccyb_obs accordingly 
xlswrite('ECentral.xlsx',oo_.SmoothedVariables.ltot ,'_ltot','D4');
xlswrite('ECentral.xlsx',oo_.SmoothedVariables.R_L ,'_rl','D4');
xlswrite('ECentral.xlsx',oo_.SmoothedVariables.R_i ,'_ri','D4');
xlswrite('ECentral.xlsx',oo_.SmoothedVariables.PD_d ,'_pdd','D4');
xlswrite('ECentral.xlsx',oo_.SmoothedVariables.PD_f ,'_pdf','D4');
xlswrite('ECentral.xlsx',oo_.SmoothedVariables.PD_h ,'_pdh','D4');
xlswrite('ECentral.xlsx',oo_.SmoothedVariables.gdpR ,'_gdpNoM','D4');
xlswrite('ECentral.xlsx',oo_.SmoothedVariables.gdp  ,'_gdp','D4');
xlswrite('ECentral.xlsx',oo_.SmoothedVariables.R  ,'_tpm','D4');




%   i) Escenario 1pm
close all 
disp('Making Kalman Filter')
dynare mafin_04.mod             % (!) Changed to Data_MAFIN_2023Q1_Pol, changed ccyb_obs accordingly 
xlswrite('efectoCCyB.xlsx',oo_.SmoothedVariables.ltot ,'_ltot','E4');
xlswrite('efectoCCyB.xlsx',oo_.SmoothedVariables.R_L ,'_rl','E4');
xlswrite('efectoCCyB.xlsx',oo_.SmoothedVariables.R_i ,'_ri','E4');
xlswrite('efectoCCyB.xlsx',oo_.SmoothedVariables.PD_d ,'_pdd','E4');
xlswrite('efectoCCyB.xlsx',oo_.SmoothedVariables.PD_f ,'_pdf','E4');
xlswrite('efectoCCyB.xlsx',oo_.SmoothedVariables.PD_h ,'_pdh','E4');
xlswrite('efectoCCyB.xlsx',oo_.SmoothedVariables.gdpR ,'_gdpNoM','E4');
xlswrite('efectoCCyB.xlsx',oo_.SmoothedVariables.gdp  ,'_gdp','E4');
xlswrite('efectoCCyB.xlsx',oo_.SmoothedVariables.R  ,'_tpm','E4');




%   ii) Escenario 2
close all 
disp('Making Kalman Filter')
dynare mafin_04.mod             % (!) Changed to Data_MAFIN_2023Q1_Pol, changed ccyb_obs accordingly 
xlswrite('efectoCCyB.xlsx',oo_.SmoothedVariables.ltot ,'_ltot','F4');
xlswrite('efectoCCyB.xlsx',oo_.SmoothedVariables.R_L ,'_rl','F4');
xlswrite('efectoCCyB.xlsx',oo_.SmoothedVariables.R_i ,'_ri','F4');
xlswrite('efectoCCyB.xlsx',oo_.SmoothedVariables.PD_d ,'_pdd','F4');
xlswrite('efectoCCyB.xlsx',oo_.SmoothedVariables.PD_f ,'_pdf','F4');
xlswrite('efectoCCyB.xlsx',oo_.SmoothedVariables.PD_h ,'_pdh','F4');
xlswrite('efectoCCyB.xlsx',oo_.SmoothedVariables.gdpR ,'_gdpNoM','F4');
xlswrite('efectoCCyB.xlsx',oo_.SmoothedVariables.gdp  ,'_gdp','F4');
xlswrite('efectoCCyB.xlsx',oo_.SmoothedVariables.R  ,'_tpm','F4');



%   ii) Escenario 2pm
close all 
disp('Making Kalman Filter')
dynare mafin_04.mod             % (!) Changed to Data_MAFIN_2023Q1_Pol, changed ccyb_obs accordingly 
xlswrite('efectoCCyB.xlsx',oo_.SmoothedVariables.ltot ,'_ltot','G4');
xlswrite('efectoCCyB.xlsx',oo_.SmoothedVariables.R_L ,'_rl','G4');
xlswrite('efectoCCyB.xlsx',oo_.SmoothedVariables.R_i ,'_ri','G4');
xlswrite('efectoCCyB.xlsx',oo_.SmoothedVariables.PD_d ,'_pdd','G4');
xlswrite('efectoCCyB.xlsx',oo_.SmoothedVariables.PD_f ,'_pdf','G4');
xlswrite('efectoCCyB.xlsx',oo_.SmoothedVariables.PD_h ,'_pdh','G4');
xlswrite('efectoCCyB.xlsx',oo_.SmoothedVariables.gdpR ,'_gdpNoM','G4');
xlswrite('efectoCCyB.xlsx',oo_.SmoothedVariables.gdp  ,'_gdp','G4');
xlswrite('efectoCCyB.xlsx',oo_.SmoothedVariables.R  ,'_tpm','G4');




%   iii) Escenario 3
close all 
disp('Making Kalman Filter')
dynare mafin_04.mod             % (!) Changed to Data_MAFIN_2023Q1_Pol, changed ccyb_obs accordingly 
xlswrite('efectoCCyB.xlsx',oo_.SmoothedVariables.ltot ,'_ltot','H4');
xlswrite('efectoCCyB.xlsx',oo_.SmoothedVariables.R_L ,'_rl','H4');
xlswrite('efectoCCyB.xlsx',oo_.SmoothedVariables.R_i ,'_ri','H4');
xlswrite('efectoCCyB.xlsx',oo_.SmoothedVariables.PD_d ,'_pdd','H4');
xlswrite('efectoCCyB.xlsx',oo_.SmoothedVariables.PD_f ,'_pdf','H4');
xlswrite('efectoCCyB.xlsx',oo_.SmoothedVariables.PD_h ,'_pdh','H4');
xlswrite('efectoCCyB.xlsx',oo_.SmoothedVariables.gdpR ,'_gdpNoM','H4');
xlswrite('efectoCCyB.xlsx',oo_.SmoothedVariables.gdp  ,'_gdp','H4');
xlswrite('efectoCCyB.xlsx',oo_.SmoothedVariables.R  ,'_tpm','H4');


%   iii) Escenario 3b
close all 
disp('Making Kalman Filter')
dynare mafin_04.mod             % (!) Changed to Data_MAFIN_2023Q1_Pol, changed ccyb_obs accordingly 
xlswrite('efectoCCyB.xlsx',oo_.SmoothedVariables.ltot ,'_ltot','I4');
xlswrite('efectoCCyB.xlsx',oo_.SmoothedVariables.R_L ,'_rl','I4');
xlswrite('efectoCCyB.xlsx',oo_.SmoothedVariables.R_i ,'_ri','I4');
xlswrite('efectoCCyB.xlsx',oo_.SmoothedVariables.PD_d ,'_pdd','I4');
xlswrite('efectoCCyB.xlsx',oo_.SmoothedVariables.PD_f ,'_pdf','I4');
xlswrite('efectoCCyB.xlsx',oo_.SmoothedVariables.PD_h ,'_pdh','I4');
xlswrite('efectoCCyB.xlsx',oo_.SmoothedVariables.gdpR ,'_gdpNoM','I4');
xlswrite('efectoCCyB.xlsx',oo_.SmoothedVariables.gdp  ,'_gdp','I4');
xlswrite('efectoCCyB.xlsx',oo_.SmoothedVariables.R  ,'_tpm','I4');




%   iii) Escenario 4
close all 
disp('Making Kalman Filter')
dynare mafin_04.mod             % (!) Changed to Data_MAFIN_2023Q1_Pol, changed ccyb_obs accordingly 
xlswrite('efectoCCyB.xlsx',oo_.SmoothedVariables.ltot ,'_ltot','J4');
xlswrite('efectoCCyB.xlsx',oo_.SmoothedVariables.R_L ,'_rl','J4');
xlswrite('efectoCCyB.xlsx',oo_.SmoothedVariables.R_i ,'_ri','J4');
xlswrite('efectoCCyB.xlsx',oo_.SmoothedVariables.PD_d ,'_pdd','J4');
xlswrite('efectoCCyB.xlsx',oo_.SmoothedVariables.PD_f ,'_pdf','J4');
xlswrite('efectoCCyB.xlsx',oo_.SmoothedVariables.PD_h ,'_pdh','J4');
xlswrite('efectoCCyB.xlsx',oo_.SmoothedVariables.gdpR ,'_gdpNoM','J4');
xlswrite('efectoCCyB.xlsx',oo_.SmoothedVariables.gdp  ,'_gdp','J4');
xlswrite('efectoCCyB.xlsx',oo_.SmoothedVariables.R  ,'_tpm','J4');



%   iii) Escenario 4b
close all 
disp('Making Kalman Filter')
dynare mafin_04.mod             % (!) Changed to Data_MAFIN_2023Q1_Pol, changed ccyb_obs accordingly 
xlswrite('efectoCCyB.xlsx',oo_.SmoothedVariables.ltot ,'_ltot','K4');
xlswrite('efectoCCyB.xlsx',oo_.SmoothedVariables.R_L ,'_rl','K4');
xlswrite('efectoCCyB.xlsx',oo_.SmoothedVariables.R_i ,'_ri','K4');
xlswrite('efectoCCyB.xlsx',oo_.SmoothedVariables.PD_d ,'_pdd','K4');
xlswrite('efectoCCyB.xlsx',oo_.SmoothedVariables.PD_f ,'_pdf','K4');
xlswrite('efectoCCyB.xlsx',oo_.SmoothedVariables.PD_h ,'_pdh','K4');
xlswrite('efectoCCyB.xlsx',oo_.SmoothedVariables.gdpR ,'_gdpNoM','K4');
xlswrite('efectoCCyB.xlsx',oo_.SmoothedVariables.gdp  ,'_gdp','K4');
xlswrite('efectoCCyB.xlsx',oo_.SmoothedVariables.R  ,'_tpm','K4');

