function [spinD, listLS] = initializeLattice3D_ones(N,D,b,pLS,pHS)
%{
%initializeLattice3D.m
%Ashley Dale
%creates a randomly initialized NxN lattice

N: number of edge spins
D: number of stacked layers
b: [1, 0 -1]; lock edge spins into HS, open boundary, or LS
lock: determines if spin is to be locked in state lock
pLS: probability spin is locked in LS state
pHS: probability spin is locked in HS state
%}

%spinD = b.*ones(N);
spinD = -1*ones(N);

pt = 1;
listLS(pt, :) = [0 0 0];

for kdx = 1:(D-2)
    
    spins = ones(N-2); %decide how many ones there are
    spins = padarray(spins,[1 1],b,'both');
    spinD = cat(3,spinD, spins);
end

spinD = cat(3, spinD, b.*ones(N));

end

