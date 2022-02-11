function plot_nHSlayers_hyst(nHS_by_layer1, nHS_by_layer2, T1, T2)

[m1, n1] = size(nHS_by_layer1);
[m2, n2] = size(nHS_by_layer2);

nHS_by_layer1(m1, :) = 0;
nHS_by_layer2(m2, :) = 0;

T = string(T1);

[X1, Y1] = meshgrid(1:n1, m1:-1:1);
[X2, Y2] = meshgrid(n2:-1:1, m2:-1:1);

nHS_by_layer2(1:2, :) = 0;

% Create figure
figure1 = figure('Color',[1 1 1]);

% Create axes
axes1 = axes('Parent',figure1);

set(gcf, 'color', [1 1 1])

h1 = waterfall([X1, X2], [Y1, Y2], [nHS_by_layer1, nHS_by_layer2]);
hold on

set(h1, 'FaceColor', 'flat')
set(h1, 'FaceAlpha', 0.3)
set(h1, 'EdgeColor', 'k')
set(h1, 'FaceVertexCData', rand(m1,3))
view(axes1,[0 0]);
xlabel('Temperature T(K)')
xticks(1:n1)
xticklabels(T)
ylabel({'Layer #','(last layer is boundary condition)'})
zlabel('nHS')
title('High Spin Fraction by Layer as Temperature Varies')

