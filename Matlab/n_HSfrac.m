function nHS = n_HSfrac(N, spins)
%{
nHS.m
Ashley Dale

Calculates the high-spin fraction for a 2D matrix of spins with the border
spins locked in high spin
%}

sum_Si = sum(spins(2:N-1, 2:N-1), 'all');

mean = ((4*(N - 1)) + sum_Si)/(N*N);

nHS = (1+mean)/2;

end