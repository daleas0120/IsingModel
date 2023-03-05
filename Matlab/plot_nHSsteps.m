function imgHSvTime = plot_nHSsteps(ising, simulation)
plt_title = 'Spin High Fraction vs Time';
figure
hold on

for idx = 1:length(ising.T)
    plot(1:simulation.dataPts, simulation.nHS(idx, :), '.-c');
end
set(gca, 'Color', simulation.APSslideColor);
set(gca, 'XColor', [1, 1, 1]);
set(gca, 'YColor', [1, 1, 1]);
grid on;
ylabel({'n_H_S'},'Interpreter','tex');
xlabel("Time (MCIMS Step)");
title({plt_title}, 'Color', 'white');
set(gcf, 'InvertHardcopy', 'off');

imgHSvTime = strcat(simulation.dir_name,'/',simulation.dat_str,'timeAvg_',...
    'delt',ising.bD_nom,'_J',ising.J_nom,'_nHSvsTime');

saveas(gcf, strcat(imgHSvTime, '.png'));
saveas(gcf, strcat(imgHSvTime, '.fig'));

end