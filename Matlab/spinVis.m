%{
spinVis.m
Ashley Dale
Visualizes a spin lattice
%}
function leg = spinVis(spinD)
[m, n, p] = size(spinD);

%Set colors
spins = spinD;

%threshold spin lattice

%find where the spins=1, set these spins equal to 10
mask_HS = (spins == 1);
spins(mask_HS) = 10;
color10 = [1,0,1];

%find where  0.75<spins<1, set these spins equal to 9
mask_HSpt75 = ((spins < 1) & (spins >0.75));
spins(mask_HSpt75) = 9;
color9 = [245/255, 106/255, 255/255];

%find where  0.5<spins<=0.75, set these spins equal to 8
mask_HSpt5 = ((spins <= 0.75) & (spins >0.5));
spins(mask_HSpt5) = 8;
color8 = [238/255, 151/255, 255/255];

%find where  0.25<spins<=0.5, set these spins equal to 7
mask_HSpt25 = ((spins <= 0.5) & (spins >0.25));
spins(mask_HSpt25) = 7;
color7 = [238/255, 151/255, 255/255];

%find where 0<spins<=0.25, set these spins equal to 6
mask_HSpt0 = ((spins <= 0.25) & (spins >0));
spins(mask_HSpt0) = 6;
color6 = [235/255, 186/255, 255/255];

%find where the spins are 0, set these spins equal to 5
mask_border = (spins == 0);
spins(mask_border) = 5;
color5 = [0,1,1];

%find where -0.25<spins<0
mask_LSpt0 = ((spins >= -0.25) & (spins < 0));
spins(mask_LSpt0) = 4;
color4 = [213/255, 232/255, 255/255]; %#d5e8ff

%find where -0.5<spins<-0.25
mask_LS0pt25 = ((spins >= -0.5) & (spins < -0.25));
spins(mask_LS0pt25) = 3;
color3 = [195/255, 253/255, 255/255]; %#c3fdff

%find where -0.75<spins<-0.5
mask_LS0pt5 = ((spins >= -0.75) & (spins < -0.5));
spins(mask_LS0pt5) = 2;
color2 = [0, 255/255, 203/255]; %#00ffcb

%find where -1<spins<-0.75
mask_LS0pt75 = ((spins > -1) & (spins < -0.75));
spins(mask_LS0pt75) = 1;
color1 = [0, 255/255, 137/255]; %#00ff89

%find where the spins are -1, set these spins equal to 4
mask_LS = (spins == -1);
spins(mask_LS) = 0;
color0 = [0,1,0]; %#00ff00

%create a set of x y z vectors that describe the space

x = 1:m;
y = 1:n;
z = 1:p;

[cx, cy, cz] = ndgrid(x, y, z);

X = cx(:);
Y = cy(:);
Z = cz(:);

C = spins(:);

%separate out the different values to correctly generate legend
%HS
hs_c = (C == 10);
hs_x = X(hs_c);
hs_y = Y(hs_c);
hs_z = Z(hs_c);

hs9_c = (C==9);
hs9_x = X(hs9_c);
hs9_y = Y(hs9_c);
hs9_z = Z(hs9_c);

hs8_c = (C==8);
hs8_x = X(hs8_c);
hs8_y = Y(hs8_c);
hs8_z = Z(hs8_c);

hs7_c = (C==7);
hs7_x = X(hs7_c);
hs7_y = Y(hs7_c);
hs7_z = Z(hs7_c);

hs6_c = (C==6);
hs6_x = X(hs6_c);
hs6_y = Y(hs6_c);
hs6_z = Z(hs6_c);

%BC
bc_c = (C == 5);
bc_x = X(bc_c);
bc_y = Y(bc_c);
bc_z = Z(bc_c);

%LS

ls4_c = (C == 4);
ls4_x = X(ls4_c);
ls4_y = Y(ls4_c);
ls4_z = Z(ls4_c);

ls3_c = (C == 3);
ls3_x = X(ls3_c);
ls3_y = Y(ls3_c);
ls3_z = Z(ls3_c);

ls2_c = (C == 2);
ls2_x = X(ls2_c);
ls2_y = Y(ls2_c);
ls2_z = Z(ls2_c);

ls1_c = (C == 1);
ls1_x = X(ls1_c);
ls1_y = Y(ls1_c);
ls1_z = Z(ls1_c);

ls_c = (C == 0);
ls_x = X(ls_c);
ls_y = Y(ls_c);
ls_z = Z(ls_c);

%adjust marker size based on number of points:
if (m*n*p) > 50000
    markerSize = 15;
else
    markerSize = 30;

% plot the items
scatter3(hs_x, hs_y, hs_z, markerSize, color10)
hold on
scatter3(hs9_x, hs9_y, hs9_z, markerSize, color9)
scatter3(hs8_x, hs8_y, hs8_z, markerSize, color8)
scatter3(hs7_x, hs7_y, hs7_z, markerSize, color7)
scatter3(hs6_x, hs6_y, hs6_z, markerSize, color6)

scatter3(bc_x, bc_y, bc_z, markerSize, color5)

scatter3(ls4_x, ls4_y, ls4_z, markerSize, color4)
scatter3(ls3_x, ls3_y, ls3_z, markerSize, color3)
scatter3(ls2_x, ls2_y, ls2_z, markerSize, color2)
scatter3(ls1_x, ls1_y, ls1_z, markerSize, color1)
scatter3(ls_x, ls_y, ls_z, markerSize, color0)

grid on

set(gcf, 'Position',[100, 100, 1250, 750])
hold off
pause(0.5)

%title({strcat(num2str(m), 'x', num2str(n), 'x', num2str(p), ' Lattice')},...
%    'Interpreter', 'tex')

leg = {'HS','','','','', '0','','','','','LS'};
end