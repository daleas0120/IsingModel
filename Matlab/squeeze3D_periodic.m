function S = squeeze3D_periodic(spins)

[m, n, p] = size(spins);
% reduced lattice size = ((size minus padding)/3) + add padding
m0 = (m/3);
n0 = (n/3);
p0 = ((p-2)/3)+2;

%create new lattice to hold results
S = zeros(m0, n0, p0);
% counters to iterate through new lattice positions
i = 1;
j = 1;
k = 2;
for idx = 2:3:(m-1)
    for jdx = 2:3:(n-1)
        for kdx = 2:3:(p-1)
            
            % value of new spin location, normalized by max value;
            sigma = (sumNNN3D(spins, idx, jdx, kdx, [1 1 1], true) + spins(idx, jdx, kdx))/27;
            
            S(i, j, k) = sigma;
            k = k+1;
        end
        j = j+1;
        k = 1;
    end
    i = i+1;
    j = 1;
end

end
