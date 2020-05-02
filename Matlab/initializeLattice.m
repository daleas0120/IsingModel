%initializeLattice.m
%Ashley Dale
%creates a randomly initialized NxN lattice

function spins = initializeLattice(N)
spins = rand(N - 2); %decide how many ones there are
for idx = 1:N-2
    for jdx = 1:N-2
        if (spins(idx, jdx) > 0.5)
            spins(idx, jdx) = 1;
        else
            spins(idx, jdx) = (-1);
        end
    end
end
%spins = padarray(spins,[2 2],0,'both');
spins = padarray(spins,[1 1], 0, 'both');
end