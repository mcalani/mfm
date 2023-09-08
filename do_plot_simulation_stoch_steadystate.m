function plot_simulation_variables

global oo_

numfig = 8;
figurename = 'Simulation and Deterministic Steady State';
numsubplot = 25;
nr = 5;
nc = 5;
lr = 5;
lc = 5;

for plt = 1:numfig 
        hplt(plt) = figure('Name',figurename);
        for index=1:numsubplot
            i = (plt-1)*numsubplot + index;
            nam = oo_.var_list(i);
            subplot(nr,nc,index)
            hh = plot(oo_.endo_simul(i,:)); hold on
            yline(oo_.steady_state(i),'r-','linewidth',2); hold on
            yline(mean(oo_.endo_simul(i,:)),'g-','linewidth',2); hold on
            yline([oo_.steady_state(i)*1.5 oo_.steady_state(i)*0.5],'--')
            title(nam,'Interpreter','none')
        end
        %saveas(hplt(plt),['C:/Users/jmoreno/Desktop/MAFIN/model_v032023/Estimation/SimOrder2Adjust/Simulation' num2str(plt) '.fig']);
end
