
function [spins, locked] = initializeLattice(N, M, b, pLS, pHS)
%{

%initializeLattice.m
%Ashley Dale
%creates a randomly initialized NxN lattice
N: lattice dimension (NxN lattice)
b: boundary condition 
lock: [1, -1]
p: probability with which spin is locked

%}
spins = rand((N - 2), (M-2)); %decide how many ones there are
pt = 1;
locked(pt, :) = [0 0];

for idx = 1:N-2
    for jdx = 1:M-2
        if (spins(idx, jdx) > 0.5)
            spins(idx, jdx) = 1;
            if rand < 2*pHS
                locked(pt,:) = [(idx + 1) (jdx + 1)];
                pt = pt+1;
            end
        else
            spins(idx, jdx) = (-1);
            %add spin location to list of locked spins with probability p
            if rand < 2*pLS
                locked(pt, :) = [(idx + 1)  (jdx + 1)];
                pt = pt+1;
            end
        end
    end
end
%spins = padarray(spins,[2 2],0,'both');

spins = padarray(spins,[1 1], b, 'both');

end
