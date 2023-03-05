function [nHSmean, imgHSvMag, trendline_fit, gof]= plot_nHSvMag(ising, simulation)

if simulation.numTrials > 1
    nHSmean = squeeze(mean(simulation.nHS))';
else
    nHSmean = mean(simulation.nHS, 2);
end

%%
num_pts = size(nHSmean, 1);
field = ising.h_field(1, 1:num_pts);
if num_pts > 1
    
    [trendline_fit, gof] = polyfit(field', nHSmean, 1);
    x = linspace(min(field), max(field), 100);
    y = polyval(trendline_fit, x);
end

%%

plt_title = 'High Spin Fraction vs Magnetic Field';

%figure
plot(field', nHSmean, '.-c')
hold on
if num_pts > 1
    plot(x, y)
end

title({plt_title},'Interpreter', 'tex', 'Color', 'white');
xlabel("Magnetic Field H (T)");
ylabel({'n_H_S'},'Interpreter','tex');
%legend(legArr,'Location','southeast')
axis([-inf inf -inf inf]);
set(gca, 'Color', simulation.APSslideColor);
set(gca, 'XColor', [1, 1, 1]);
set(gca, 'YColor', [1, 1, 1]);
grid on;
hold off;
set(gcf, 'InvertHardcopy', 'off');

imgHSvMag = strcat(simulation.dir_name,'/',simulation.dat_str,'_',...
    'delt',ising.bD_nom,'_J',ising.J_nom,'_nHSvsH');

saveas(gcf, strcat(imgHSvMag,'.png'));
saveas(gcf, strcat(imgHSvMag,'.fig'));

end