
N = 12;

J = 1;

k_b = 8.617333262*10^-5;%eV/K

T = [0.01:0.5:5];

beta = J./(k_b.*T);

T_c = 2.27;

steps = 1000;

[spins, ~] = initializeLattice(N, 0, 0, 0);

figure;
spinVis(spins)
pause(0.5)
close all

for n = 1:length(beta)
    for t = 1:steps
        
        for x = 2:N-1
            for y = 2:N-1
                
                delta_spin = -1*spins(x,y) - spins(x,y);
                
                Snn = spins(x+1, y) + ...
                    spins(x-1, y) + ...
                    spins(x, y+1) + ...
                    spins(x, y-1);
                
                spins(x, y) = -1*spins(x,y);
                
                dE = delta_spin*(-1*J*Snn);
                
                p = exp(-1*dE/beta(n));
                r = rand;
                
                if dE < 0 || p >= r
                    continue
                else
                    spins(x,y) = -1*spins(x,y);
                end
                
            end
        end
        
        % At the end of each step
        nHS(t, n) = n_HSfrac(spins);
    end
end

figure;

plot(T, nHS)