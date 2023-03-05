function[nHSmean, imgHSvTemp]= plot_nHSvTemp(ising, simulation)

if simulation.numTrials > 1
    nHSmean = squeeze(mean(simulation.nHS))';
else
    nHSmean = mean(simulation.nHS, 2);
end

plt_title = 'High Spin Fraction vs Temperature';

figure
plot(ising.T_K', nHSmean, '.-c')
hold on
title({plt_title},'Interpreter', 'tex', 'Color', 'white');
xlabel("Temperature T (K)");
ylabel({'n_H_S'},'Interpreter','tex');
%legend(legArr,'Location','southeast')
axis([-inf inf 0 1.0]);
set(gca, 'Color', simulation.APSslideColor);
set(gca, 'XColor', [1, 1, 1]);
set(gca, 'YColor', [1, 1, 1]);
grid on;
hold off;
set(gcf, 'InvertHardcopy', 'off');

imgHSvTemp = strcat(simulation.dir_name,'/',simulation.dat_str,'_',...
    'delt',ising.bD_nom,'_J',ising.J_nom,'_nHSvsT');

saveas(gcf, strcat(imgHSvTemp,'.png'));
saveas(gcf, strcat(imgHSvTemp,'.fig'));

end