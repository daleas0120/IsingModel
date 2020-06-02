function nHS = n_HSfrac(spins)
%{
nHS.m
Ashley Dale

Calculates the high-spin fraction for a 2D matrix of spins with the border
spins locked in high spin
%}

N = max(size(spins));

%sum_Si = sum(spins(2:N-1, 2:N-1), 'all');
sum_Si = sum(spins, 'all');

%mean = ((4*(N - 1)) + sum_Si)/(N*N);

mean = sum_Si/((N-2)^2);

nHS = (1+mean)/2;

end