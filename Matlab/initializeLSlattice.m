function [spins] = initializeLSlattice(N, M, b)
%{

%initializeLattice.m
%Ashley Dale
%creates a randomly initialized NxN lattice
N: lattice dimension (NxN lattice)
b: boundary condition 
lock: [1, -1]
p: probability with which spin is locked

%}

spins = -1.*ones((N-2), (M-2));
spins = padarray(spins,[1 1], b, 'both');


end