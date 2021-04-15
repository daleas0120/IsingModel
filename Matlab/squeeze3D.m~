function S = squeeze3D(spins)

[m, n, p] = size(spins);

m0 = ((m-2)/3)+2;
n0 = ((n-2)/3)+2;
p0 = ((p-2)/3)+2;

S = zeros(m0, n0, p0);
i = 2;
j = 2;
k = 2;
for idx = 3:3:(m - 2)
    for jdx = 3:3:(n-2)
        for kdx = 3:3:(p-2)
            
            sigma = (sumNNN3D(spins, idx, jdx, kdx, [1 1 1]) + spins(idx, jdx, kdx))/27;
            
            S(i, j, k) = sigma;
            k = k+1;
        end
        j = j+1;
        k = 2;
    end
    i = i+1;
    j = 2;
end

end
