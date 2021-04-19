function nHS = n_HSfrac3D(spins)
%{
nHS.m
Ashley Dale
Calculates the high-spin fraction for a 3D matrix of spins
%}

[N, M, D] = size(spins);

S = spins(2:N-1, 2:M-1, 4:D-1);

sum_Si = sum(S, 'all');

mean = sum_Si/((N-2)*(M-2)*(D-4));

nHS = (1+mean)/2;

end