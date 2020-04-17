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

function [spins, E, B] = equilibrateSpins_H(...
    time, N, spins, k, T, mu, H, J, big_delta, ln_g, ...
    frameRate, spins_last, dir_name)

set(0,'DefaultTextInterpreter','none')
k_b = 8.617333262*10^-5;%eV/K
E = zeros(time, 1);
B = zeros(time, 1);

for idx = 1:time% how many times to let the system evolve
    for row = 1:N
        for col = 1:N
            %i = row;
            %j = col;
            
            i = randi([2, N-1]);
            j = randi([2, N-1]);
            
            %pick spin and flip right away
            spins(i,j) = -1*spins(i,j);
            
            sum_nn = (spins((i-1),j) + spins((i+1),j) +...
                spins(i,(j-1)) + spins(i,(j+1)));
            
            %then do change in energy with correct sign
            %dE = 2*spins(i,j) * (J*sum_nn + H*mu);
            dE = spins(i,j)*(2*J*sum_nn - (big_delta - k_b*T*ln_g));
            
            boltzConst = exp(dE*k);
            p = exp(-1*dE*k);
            
            if dE < 0
                continue
            elseif p < rand()
                spins(i,j) = -1*spins(i,j);
            end
            %then check with random number; if state is acceptable,
            %keep and move on; if state is not acceptable then flip
            %sign back and try a different spin
            
            %r = rand(); %r is between 0 and 1
            %if boltzConst < r
            %    spins(i,j) = -1*spins(i,j);
            %end
        end
    end
    
    Snn = nearestN(spins);
    sum_Si = sum(spins(2:N-1, 2:N-1), 'all');
    
    E(idx, 1) = -J*Snn;
    B(idx, 1) = mu*sum_Si;
    
    if mod(idx, frameRate) == 0
        frame_temp = num2str(J/(k*k_b));
        frame_name = num2str(idx);
        s = num2str(N);
        frame_name = strcat(dir_name,'\frames\',s,'spins','_',...
            frame_temp,'K_', frame_name, ".png");
        comp = spins;
        imagesc(comp)
        title(frame_name)
        colorbar
        axis square;
        saveas(gcf, frame_name)
    end
    
end
E = mean(E);
B = (mean(B)./(N*N));
s_ave = tanh(k*mu*H + k*4*J);
end