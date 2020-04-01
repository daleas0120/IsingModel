%initializeLattice.m
%Ashley Dale
%creates a randomly initialized NxN lattice

function spins = initializeLattice(N)
spins = rand(N); %decide how many ones there are
for idx = 1:N
    for jdx = 1:N
        if (spins(idx, jdx) > 0.5)
            spins(idx, jdx) = 1;
        else
            spins(idx, jdx) = (-1);
        end
    end
end
spins = padarray(spins,[2 2],0,'both');
end