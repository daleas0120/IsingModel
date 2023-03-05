function plot_nHS(ising, simulation)
plt_title = 'Spin High Fraction vs Time';
figure
hold on
for idx = 1:length(ising.T)
    plot(1:simulation.evo, simulation.nHS_evo(idx, :), '.-c');
end
set(gca, 'Color', simulation.APSslideColor);
set(gca, 'XColor', [1, 1, 1]);
set(gca, 'YColor', [1, 1, 1]);
grid on
ylabel({'n_H_S'},'Interpreter','tex');
xlabel("Time (MCIMS Step)");
title({plt_title}, 'Color', 'white');
set(gcf, 'InvertHardcopy', 'off');
saveas(gcf, strcat(simulation.dir_name,'/',simulation.dat_str,'timeEvo_',...
    'delt',ising.bD_nom,'_J',ising.J_nom,'_nHSvsT.png'));
saveas(gcf, strcat(simulation.dir_name,'/',simulation.dat_str,'timeEvo_',...
    'delt',ising.bD_nom,'_J',ising.J_nom,'_nHSvsT.fig'));

end