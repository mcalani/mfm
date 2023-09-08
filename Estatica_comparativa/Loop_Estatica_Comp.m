clear all
close all
clc

%% NOTA: Comentar primer bloque en archivo Estatica_Comp.m


%% Settings
% Choose range of values for CCyB
acam={'ccyb'}; rang = -0.03:0.001:0.05;

% 0.08+min(rang)

% Choose figure's x-axis
%xaxis = 'ccyb';
xaxis = 'phi_f';

% Choose variables to plot (Variable*scale, Name in figure)
selected_figures={   acam{1},    xaxis
                     'gdp',      'PIB'
                      'c',        'Consumo'
                      'PD_d*100', 'PD Bancos'
                      'l',        'Credito total'   
%                       'h',        'Inv. Inmobiliaria'                        
%                      'l_f',      'Cartera Com.'
%                      'q_l*l_h',  'Cartera Hip.'
                     
                     
                                           
%                       'i_agg',    'Inversion'
%                      'imp',      'Importaciones'
%                      'k',        'Capital'
%                      'R',        'TPM' % no se mueve 
%                       'ppi',      'Inflacion'% no se mueve 
%                      %'R_f_tilde', '$\tilde{R}^f$' %se mueve   
%                       'R_e',      'Retorno capital'                     
%                      %'R_h',      '$R^H$'% no se mueve 
                      %'R_BB',      '$R^{BB}$'% no se mueve 
%                      %'R_h',      '$R^H$'% no se mueve 
%                      %'R_BL',      '$R^{BL}$'% no se mueve 
%                      'R_L',      'Tasa Comercial'
%                      'R_i',      'Tasa hipotecario'
%                      'R_D',      'Tasa Depositos'
%                       'PD_e*100', 'PD Comercial'
%                       'PD_i*100', 'PD Hipotecario'
%                       "bl_priv"     	,'Bonos de Gob. LP'
%  					 "bs_priv"     	,'Bonos de Gob. CP'
% 					 %"spread_RBL_R"  ,'Dif. tasa larga y de depositos'% no se mueve 
%  					 "spread_RL_RD"  ,'Spread Com.'
%  					 "spread_Ri_RD"  ,'Spread Hip.'
% 					 "spread_RD_R"   ,'Spread Dep.'
%                       'PD_h*100', 'PD Bancos H'
%                       'PD_f*100', 'PD Bancos F'
%                      %'n_e',      '$N^e$' %se mueve
%                      %'n_b',      '$N^b$' %se mueve
%                      %'p_H' ,     'Precio bien compuesto' %se mueve
%                      %'p_F',      'p_F' %se mueve
%                      %'w',        'Salario' %se mueve
%                      'r_k',      'Tasa arriendo de K'
%                      %'q_k' ,     '$q_k$'% no se mueve 
%                      %'q_h',      '$q_h$'% no se mueve 
%                       'q_BB',       'Bonos bancarios LP'
                      %'q_BL',       'q_BL'               % no se mueve       
                     };

% Change layout
layout_x = 2; %3
layout_y = 2; %6

%%
%preallocation
rang_length = length(rang);
vec   = zeros(rang_length,24);
vecex = zeros(rang_length);

for jj = 1:rang_length
    num_var=rang(jj);

    clearvars -except jj layout_x layout_y num_var vec acam name vecex graphs rang xaxis selected_figures

    run Load_params;                %carga parámetros y exógenos
    eval([acam{1}  '=num_var']);    %evalúa ccyb en el loop
    run Estatica_comparativa;  %calcula SS dado el ccyb evaluado en el loop
    
    save TEMP -regexp ^(?!(layout_x|layout_y|selected_figures|x_guess|jj|num_var|vec|vecex|graphs|options|Exogenous|SS_valuesacam|rang|xaxis)$). %guarda valores obtenidos en objeto TEMP
    vars2save = load('TEMP');       %carga valores del workspace en vars2save   
    varsName  = string(fieldnames(vars2save)); %guarda los nombres de las variables
      
    eval(['name =' acam{1} ';']);

    selected_graphs = cell(1,size(selected_figures,1));
    for ii = 2:(size(selected_graphs,2))
        aux= eval([ selected_figures{ii,1} ]);
        selected_graphs(ii) = {aux};
    end

    selected_graphs(1) = {name};

    if isequal(xaxis,'phi_f')
        selected_graphs(1)={phi_f*100};
    end


    %Variables para figuras
    numb_graphs = length(selected_graphs)     ;
    for i=1:1:numb_graphs
        vec(jj,i)= cell2mat(selected_graphs(i));
    end
    
    %Resto de variables
    for k=1:numel(varsName)
        if( isnumeric(vars2save.(varsName{k})) )
            vecex(jj,k)= vars2save.(varsName{k});
        end
    end
end 



%% PLOTS
disp('plottig!')
names= {selected_figures{:,2}};

if isequal(xaxis,'phi_f')
    names{1}='Capital Base ($\phi$)';
end    

set(gcf,'PaperOrientation','landscape')

%discount voluntary capital? yes =1; no=0
dis_cap =0;
if dis_cap == 0
    disc = 0.0433*100;
elseif dis_cap == 0
    disc = 0;
end

for i = 2:numb_graphs
    subplot(layout_x,layout_y,i-1)
    plot(vec(:,1)-disc,vec(:,i),'b','LineWidth',2)
    xticks((floor(min(vec(:,1)-disc))):2:(floor(max(vec(:,1)-disc))))
    set(gca,'FontSize', 10,'FontName', 'Times','XLim', [min(vec(:,1)-disc), max(vec(:,1)-disc)]);
    title(names(i),'interpreter','latex','FontSize',15)
    xlabel(names(1),'interpreter','latex','FontSize',15)
end

fig_ratio=16/9;
trim_hor=0.6;
trim_top=0.05;

% set(gcf, 'Units', 'Inches', 'Position', [0, 0, 7.25, 9.125], ...
%     'PaperUnits', 'Inches', 'PaperSize', [7.25, 9.125])


% set(gcf,'Units','Inches');
% pos = get(gcf,'Position');
% set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])


set(gcf, 'Units', 'Inches', 'Position', [0, 0, 14.125, 7.25], ...
    'PaperUnits', 'Inches', 'PaperSize', [12.125,7.25,])
saveas(gcf, 'estatica_comp.svg')
saveas(gcf, 'estatica_comp.pdf')


%% Find optimal phi

pib_ind = find(contains(selected_figures(:,2),'PIB'));

max_gdp=max(vec(:,pib_ind)) ; 
ind_maxgdp = find(vec(:,2)==max_gdp) ;

phi_optimal = vec(ind_maxgdp,1)- 0.0433*100 ;

%% EXPORT SS-VALUES (not finished)
% 
% for Z=1:numel(vars2save)
%     namesex(Z)=cellstr(varsName(Z));
% end
% disp('Exporting!')
% figure=acam{1};
% print(figure,'-dpdf','-bestfit');
% filename='resultados.xlsx';
% 
% writecell(namesex,filename,'Sheet',figure,'Range','A1','WriteMode','overwritesheet')
% writematrix(vecex,filename,'Sheet',figure,'Range','A2')
% 
% %writetable(mytable,filename)


disp('done!')