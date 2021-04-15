function C_NN = crossCorrelation(spins)

[M, N, P] = size(spins);
for k = 2:3
    for i = 2:M-1
        for j = 2:N-1
            
            C_NN(i-1, j-1, k-1) = spins(i, j, k)*sumNN(spins, i, j, k);
            
        end
    end
end

end