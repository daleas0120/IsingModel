
function [spins, locked] = initializeLattice_noPad(N, M, pLS, pHS)
%{
%initializeLattice.m
%Ashley Dale
%creates a randomly initialized NxN lattice
N: lattice dimension (NxN lattice)
lock: [1, -1]
p: probability with which spin is locked
%}

spins = rand(N, M); %decide how many ones there are
pt = 1;
locked(pt, :) = [0 0];

for idx = 1:N
    for jdx = 1:M
        if (spins(idx, jdx) > 0.5)
            spins(idx, jdx) = 1;
            if rand < 2*pHS
                locked(pt,:) = [idx jdx];
                pt = pt+1;
            end
        else
            spins(idx, jdx) = (-1);
            %add spin location to list of locked spins with probability p
            if rand < 2*pLS
                locked(pt, :) = [idx  jdx];
                pt = pt+1;
            end
        end
    end
end

end
