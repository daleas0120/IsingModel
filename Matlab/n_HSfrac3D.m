function nHS = n_HSfrac3D(spins)
%{
nHS.m
Ashley Dale

Calculates the high-spin fraction for a 2D matrix of spins with the border
spins locked in high spin
%}

[N, ~, D] = size(spins);

sum_Si = sum(spins(2:N-1, 2:N-1, 2:D-1), 'all');
edge = (N-2);
corners = 8;
top = 4*edge + (N-2)^2; 
bottom = 4*edge + (N-2)^2; 
sides = (D-2)*4*(N-1);

mean = ((corners + top + bottom + sides + sum_Si))/(N*N*D);

nHS = (1+mean)/2;

end