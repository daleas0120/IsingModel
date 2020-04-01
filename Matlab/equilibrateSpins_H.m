%equilabrateSpins_H.m
%Ashley Dale

%Cools a matrix of spins to a given temperature, and at various times saves
%an image of the spin matrix to a file

%time: an integer that determines how long the system cools
%N: the square root of the number of spins
%k: 1/temperature
%mu: atomic magnetic moment
%h:
%H: external magnetic field
%J: spin exchange coupling constant
%frameRate: determines how frequently intermediate spin matrices are saved
%to an image
%spins_last: previous spin value
%dir_name: where to save results

function [spins, E, B] = equilibrateSpins_H(...
    time, N, spins, k, mu, h, H, J, frameRate, spins_last, dir_name)



E = zeros(time, 1);
B = zeros(time, 1);

for idx = 1:time% how many times to let the system evolve
    for row = 3:N+2
        for col = 3:N+2
            i = row;
            j = col;
            
            % = randi([3, N + 2]);
            %j = randi([3, N + 2]);
            
            %pick spin and flip right away
            spins(i,j) = -1*spins(i,j);
            
            sum_nn = (spins((i-1),j) + spins((i+1),j) +...
                spins(i,(j-1)) + spins(i,(j+1)));
            
            %then do change in energy with correct sign
            dE = 2*spins(i,j) * (J*sum_nn + mu*h);
            
            boltzConst = exp(dE*k);
            
            %then check with random number; if state is acceptable,
            %keep and move on; if state is not acceptable then flip
            %sign back and try a different spin
            
            r = rand(); %r is between 0 and 1
            if boltzConst < r
                % if boltzConst = inf then flip
                spins(i,j) = -1*spins(i,j);
            end
        end
    end
    
    Snn = nearestN(spins);
    sum_Si = sum(spins, 'all');
    
    E(idx, 1) = -J*Snn;
    B(idx, 1) = mu*H*sum_Si;
    
    if mod(idx, frameRate) == 0
        frame_temp = num2str(k);
        frame_name = num2str(idx);
        frame_name = strcat(dir_name,'/frames/', frame_temp,'_', frame_name, ".png");
        comp = spins;
        imagesc(comp)
        title(frame_name)
        axis square;
        saveas(gcf, frame_name)
    end
    
end
end