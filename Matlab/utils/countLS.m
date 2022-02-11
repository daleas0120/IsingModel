function listLS = countLS(spins)

N = max(size(spins));

[t1, t2] = find(spins(2:N-1, 2:N-1) == (-1));

listLS = [t1 t2];

end