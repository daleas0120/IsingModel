%{
spinVis.m
Ashley Dale
Visualizes a spin lattice
%}
function spinVis(spinD)
[m, n, p] = size(spinD);

%Set colors
spins = spinD;

mask_border = (spins == 0);
spins(mask_border) = spins(mask_border) + 6;

mask_HS = (spins == 1);
spins(mask_HS) = spins(mask_HS) + 4;

mask_LS = (spins == -1);
spins(mask_LS) = spins(mask_LS) + 3;


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
hs_c = (C == 5);
hs_x = X(hs_c);
hs_y = Y(hs_c);
hs_z = Z(hs_c);

%LS
ls_c = (C == 2);
ls_x = X(ls_c);
ls_y = Y(ls_c);
ls_z = Z(ls_c);

%BC
bc_c = (C == 6);
bc_x = X(bc_c);
bc_y = Y(bc_c);
bc_z = Z(bc_c);

%adjust marker size based on number of points:
if (m*n*p) > 50000
    markerSize = 15;
else
    markerSize = 2;
end

%scatter3(X, Y, Z, 2, C)
scatter3(hs_x, hs_y, hs_z, markerSize, [1,0,1])
hold on
scatter3(ls_x, ls_y, ls_z, markerSize, [0,1,0])
scatter3(bc_x, bc_y, bc_z, markerSize, [0,1,1])

%{
for i = 1:m
    for j = 1:n
        for k = 1:p
            if spinD(i,j,k) == 1
                argC = 'r';
                argS = '^';
            elseif spinD(i,j,k) == (-1)
                argC = 'k';
                argS = 'v';
            else
                argC = 'b';
                argS = '*';
            end
            plot3(i,j,k,argS,'Color', argC, 'MarkerFaceColor',argC)
            hold on
        end
    end
end
%}
grid on
%axis equal
legend({'HS', 'LS', '0'})
set(gcf, 'Position',[100, 100, 1250, 750])
hold off
pause(1)
title(strcat(num2str(m), 'x', num2str(n), 'x', num2str(p), ' Lattice'))
end