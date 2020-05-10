
function Snn = nearestN3D(spins)
%{
%nearestN3D.m
%Ashley Dale
%Sums over all nearest neighbors in 3D
%}
[N, ~, D] = size(spins);
Snn = 0;
for i = 2:N - 1
    for j = 2:N - 1
        for k = 2:D - 1
            Snn = Snn + ...
                spins(i, j, k)*(...
                spins(i-1, j, k) +...
                spins(i+1, j, k) +...
                spins(i, j-1, k) +...
                spins(i, j+1, k) +...
                spins(i, j, k-1) +...
                spins(i, j, k+1));
        end
    end
end
%because we visit each spin 6 times, need to divide total interaction by 6
Snn = Snn/6;
end