function spinD = initializeLattice3D(N,D)
%{
%initializeLattice3D.m
%Ashley Dale
%creates a randomly initialized NxN lattice

%}

spinD = ones(N);
for kdx = 1:(D-2)
    spins = rand(N-2); %decide how many ones there are
    for idx = 1:(N-2)
        for jdx = 1:(N-2)
            if (spins(idx, jdx) > 0.5)
                spins(idx, jdx) = 1;
            else
                spins(idx, jdx) = (-1);
            end
        end
    end
    spins = padarray(spins,[1 1],1,'both');
    spinD = cat(3,spinD, spins);
end
    spinD = cat(3, spinD, ones(N));
end