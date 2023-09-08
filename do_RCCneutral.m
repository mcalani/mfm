%   * CALCULO SIMPLE EVALUACION RCC NEUTRAL *
%   =========================================

%   SUPUESTO:       Definicion del nivel neutral es exogenamente dado 
%   METODOLOGIA:    Condicional en un shock financiero compuesto, se libera el maximo de RCC disponible y se cuantifica su aporte 
%
%
%   Autores:        JMoreno, MCalani 
%   Fecha revision: 31 Julio 2023


%% CLEAN HOUSE 
cd C:\Users\mcalani\Desktop\RPF
clear all
close all 

%% CORRER MODELO CON ESTADO ESTACIONARIO ARBITRARIO 
%       Select en calibración: i) capital base y ii) AR regla ccyb y e_req 
dynare mafin_04                                 


%%  CALCULATE IRF 
misVar = ["gdpn" "ppi" "R_L" "R" "PD_d" "i_agg" "ltot" "ccyb"]; 
%misShock = ["sigma_ee" "sigma_ff" "sigma_hh" "sigma_ii" "R_W" "xi_R" "e_req"];    %(!) e_req debe ser el último y los otros shocks deben ser "malos"                              
%misShock = ["sigma_ee" "sigma_ff" "sigma_hh" "sigma_ii" "xi_roe_r" "xi_CHI_b" "xi_CHI_e" "epsilon_L" "e_req"];
misShock = ["sigma_ee" "sigma_ff" "sigma_hh" "sigma_ii" "xi_roe_r"  "xi_CHI_e" "epsilon_L" "e_req"];

numV = length(misVar);
numS = length(misShock);
nPer = size(oo_.irfs.c_u_e_req,2);
vars_ss = oo_.steady_state;

%   Rescatar IRF variables y shocks definidos
misIRF = zeros(numV,nPer,numS); 
escIRF = misIRF; 

for iVar=1:numV
    tempVarSS = vars_ss(loc(oo_.var_list,eval("'"+misVar(iVar)+"'" )));
    for jVar = 1:numS
        misIRF(iVar,:,jVar) = eval("oo_.irfs."+misVar(iVar)+"_u_"+misShock(jVar));
        escIRF(iVar,:,jVar) = eval("oo_.irfs."+misVar(iVar)+"_u_"+misShock(jVar)) / tempVarSS ;        
    end
end

%%  CONSTRUIR SHOCK BENCHMARK    
    % Alternativa 1 escenario inventado 
rescale = 1.65;     %   Escala de shock 1.65 veces 1std - iff 

    % Alternativa 2 escenario 2008 
DefCris = min(oo_.SmoothedVariables.gam_lf_obs);                    %Mayor caida del crecimiento credito comercial / gam4_YR_obs si producto 
PerCris = find(oo_.SmoothedVariables.gam_lf_obs==DefCris); 

mediumshocks = ["a", "g", "ppi_ast", "xi_m", "y_ast", "z", "y_co", "epsilon_L", "sigma_ee"];
bigshocks = ["xi_R", "R_W" , "e_m",  "zetau"];
ShockCrisis = zeros(1,numS); 
scale = zeros(1,numS);

for iS=1:numS-1
    if find(bigshocks==misShock(iS))>0
        scale(iS)=1000;
    elseif find(mediumshocks==misShock(iS))>0
        scale(iS)=100;
    else
        scale(iS)=10;
    end
    tempShockSS = vars_ss(loc(oo_.var_list,eval("'"+misShock(iS)+"'" )));
    tempShock   = eval("oo_.SmoothedShocks.u_"+misShock(iS)+"("+PerCris+")");    % u_{exovar} en ese periodo 
    tempStd     = eval("oo_.posterior_mode.shocks_std.u_"+misShock(iS));          % SIGMA(u_{exovar}) estimado 
        %ShockCrisis(iS) = (tempShockSS-tempShock)*scale(iS)/tempStd;
        %ShockCrisis(iS) = (tempShock)*scale(iS)/tempStd;
    ShockCrisis(iS) = (tempShock)/tempStd;
end
ShockCrisis = [ShockCrisis 0];  % Para escalar que no haya shock de e_req      


%Escalar shocks a la crisis
for iVar=1:numV
    tempVarSS = vars_ss(loc(oo_.var_list,eval("'"+misVar(iVar)+"'" )));
    for jShock=1:numS
        shocks_aux(iVar,:,jShock)     = misIRF(iVar,:,jShock)*ShockCrisis(jShock);
        shocks_aux_esc(iVar,:,jShock) = shocks_aux(iVar,:,jShock)/tempVarSS;
    end
end


%%  AGREGAR IRFs 
    % Agregar shocks malos 
for iV=1:numV
    for iN=1:nPer  
        if misVar(iV) == "PD_d" | misVar(iV)== "gam4_YR_obs"
            %v_crisis(iV,iN) = 100 * rescale * sum(misIRF(iV,iN,:)) ;
            v_crisis(iV,iN) = 100 * sum(shocks_aux(iV,iN,:)); 
        elseif misVar(iV) == "ccyb"
            v_crisis(iV,iN) = 0;
        else
            %v_crisis(iV,iN) = 100 * rescale * sum(escIRF(iV,iN,:)) ;
            v_crisis(iV,iN) = 100 * sum(shocks_aux_esc(iV,iN,:)) ; 
        end
    end
end

    % IRF a liberar CCYB
%impulsoInit = oo_.irfs.ccyb_u_e_req(1);  % = 0.1%
impulsoMax = max(oo_.irfs.ccyb_u_e_req);
rescNeu  = ccybneutral/impulsoMax ;     % = 10 en el caso de liberar 1% 

v_libCCyB   = zeros(numV,nPer);
for iV=1:numV
    for iN=1:nPer
        if misVar(iV) == "PD_d" | misVar(iV)=="ccyb" | misVar(iV)== "gam4_YR_obs"
            v_libCCyB(iV,iN) = - 100 * rescNeu * misIRF(iV,iN,numS) ;
        else
            v_libCCyB(iV,iN) = - 100 * rescNeu * escIRF(iV,iN,numS); 
        end
    end
end

    % IRF en variables endogenas post liberación 
v_crisisCCyB = v_crisis + v_libCCyB ; 
    

%%  GRAFICAR
figurename = 'Crisis con Liberación de RCC';
numsubplot = 8;
nr = 2; 
nc = 4;
Figure= figure('Name',figurename);

for i=1:numV
    nam = misVar(i);
    subplot(nr,nc,i)
    hh = plot(v_crisis(i,1:20),'r','linewidth',2); hold on
    %plot(v_libCCyB(i,:),'c-*','linewidth',1.5); hold on 
    plot(v_crisisCCyB(i,1:20),'g--','linewidth',2); hold offp
    title(nam,'Interpreter','none')
end
legend('Crisis','Liberación','CrisisYLib')

%% Tareas
% Revisar el mort_tot
% Limpiar shocks de crisis anterior
% Definir los ratios
% Tablita con beneficios y costos (c y gdp)



