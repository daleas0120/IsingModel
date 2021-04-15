function S = sumNNN(spins, i, j)
%{

%}

[m n] = size(spins);

sNN_wt = 1;
sNNN_wt = 0.75;
sNNNN_wt = 0.25;

sNN = spins(i, j+1)+...
    spins(i+1, j)+...
    spins(i-1, j)+...
    spins(i, j-1);

sNNN = spins(i+1, j+1) +...
    spins(i-1, j+1)+...
    spins(i+1, j-1)+...
    spins(i-1, j-1);

if i>2 && i<m-2 && j>2 && j<n-2
    sNNNN = spins(i-2, j)+...
        spins(i+2, j)+...
        spins(i, j-2)+...
        spins(i, j+2);
else
    sNNNN = 0;
end

S = sNN_wt*sNN + sNNN_wt*sNNN + sNNNN_wt*sNNNN;

end