function H = magnetism(spins)
% return the sum of the lattice spins

[m, n, p] = size(spins);

H = sum(spins, 'all')/(m*n*p);

end