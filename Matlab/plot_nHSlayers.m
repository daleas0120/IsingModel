function plot_nHSlayers(nHS_by_layer1, nHS_by_layer2, T1, T2)

nHS = [nHS_by_layer1, nHS_by_layer2];
[m, n] = size(nHS);
nHS(m, :) = 0;

T_min = T1(1);
T_max = T1(n/2);

T_up = linspace(T_min, T_max, 5);
T_down = linspace(T_up(4), T_min, 4);

T_tics = linspace(1, m, 9);
T_labels = string([T_up, T_down]);

T = string([T1, T2]);

[X, Y] = meshgrid(1:n, m:-1:1);

figure;
colordef white
set(gcf, 'color', [1 1 1])

h = waterfall(X, Y, nHS);
set(h, 'FaceColor', 'flat')
set(h, 'FaceAlpha', 0.3)
set(h, 'EdgeColor', 'k')
set(h, 'FaceVertexCData', rand(m,3))
xlabel('Temperature T(K)')
xticks(1:n)
xticklabels(T)
ylabel({'Layer #','(last layer is boundary condition)'})
zlabel('nHS')
title('High Spin Fraction by Layer as Temperature Varies')

