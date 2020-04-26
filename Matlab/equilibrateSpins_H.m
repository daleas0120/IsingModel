function [spins, E, B, nHS] = equilibrateSpins_H(...
    time, spins, k, T, mu, H, J, big_delta, ln_g, ...
    frameRate, dir_name, saveIntResults)
%{
%equilabrateSpins_H.m
%Ashley Dale

%Cools a matrix of spins to a given temperature, and at various times saves
%an image of the spin matrix to a file

%time: an integer that determines how long the system cools
%N: the square root of the number of spins
%k: 1/temperature
%mu: atomic magnetic moment
%H: external magnetic field
%J: spin exchange coupling constant
%frameRate: determines how frequently intermediate spin matrices are saved
%to an image
%spins_last: previous spin value
%dir_name: where to save results
%saveIntResults: boolean to control writing of frame samples

%}
    
set(0,'DefaultTextInterpreter','none')
k_b = 8.617333262*10^-5;%eV/K

E = zeros(time, 1);
B = zeros(time, 1);
nHS = zeros(time, 1);
N = max(size(spins));

for idx = 1:time% how many times to let the system evolve
    for row = 2:N-1
        for col = 2:N-1
            i = row;
            j = col;
            
            if N > 2      
                %pick spin and flip right away
                %spins(i,j) = -1*spins(i,j);
              
                S = -1*spins(i,j);
                delta_sig = S - spins(i,j);
                
                spins(i, j) = -1*spins(i,j);
                
                sum_nn = (spins((i-1),j) + spins((i+1),j) +...
                    spins(i,(j-1)) + spins(i,(j+1)));
                
                %then do change in energy with correct sign
                %dE = 2*spins(i,j) * (J*sum_nn + H*mu);
                dE = delta_sig*(-1*J*sum_nn + (big_delta/2 - T*ln_g/2));
                
                %E1 = spins(i,j)*J*sum_nn - (big_delta/2 - k_b*T*ln_g/2)*spins(i,j);
                %E2 = (-1*spins(i,j))*J*sum_nn - (big_delta/2 - k_b*T*ln_g/2)*(-1*spins(i,j));
                %dE = E2 - E1;
                
                p = exp(-1*dE/T);
                
                if dE < 0 || p > rand()
                    continue
                else
                    spins(i,j) = -1*spins(i,j);
                end
                
                %then check with random number; if state is acceptable,
                %keep and move on; if state is not acceptable then flip
                %sign back and try a different spin
                
            end
        end
    end
    
    Snn = nearestN(spins);
    sum_Si = sum(spins(2:N-1, 2:N-1), 'all');
    
    E(idx, 1) = -J*Snn;
    B(idx, 1) = mu*sum_Si;
    nHS(idx, 1) = n_HSfrac(N, spins);
    
    if mod(idx, frameRate) == 0
        pltTitle = strcat(num2str(N),'spins','_',num2str(T),'K_', num2str(idx));
        imagesc(spins)
        title(pltTitle)
        colorbar
        axis square;
        pause(0.05);
        if saveIntResults
            frame_name = strcat(dir_name,'\frames\',pltTitle,".png");
            saveas(gcf, frame_name)
        end
    end
    
end

E = mean(E);
B = (mean(B)./(N*N));
nHS = mean(nHS);

end