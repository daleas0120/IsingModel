function [spins, E, nHS] = equilibrateSpins_3D(...
    time, spins, ~, T, ~, ~, J,omega, big_delta, ln_g, listLS, ...
    frameRate, dir_name, saveIntResults)
%{
equilabrateSpins_H.m
Ashley Dale

Cools a matrix of spins to a given temperature, and at various times saves
an image of the spin matrix to a file

time: an integer that determines how long the system cools
N: the square root of the number of spins
k: 1/temperature
mu: atomic magnetic moment
H: external magnetic field
J: spin exchange coupling constant
ln_g:
listLS:
    
frameRate: determines how frequently intermediate spin matrices are saved
to an image
spins_last: previous spin value
dir_name: where to save results
saveIntResults: boolean to control writing of frame samples

    %}
    
    set(0,'DefaultTextInterpreter','none')
    %omega = 0.001;
    E = zeros(time, 1);
    nHS = zeros(time, 1);
    [N, ~, D] = size(spins);
    
    %% some optimization
    
    longRange = omega*(big_delta - T*ln_g)/2;
    
    
    
    %%
    
    
    for idx = 1:time% how many times to let the system evolve
        if N > 2
            for xAxis = 2:N-1
                for yAxis = 2:N-1
                    for zAxis = 2:D-1
                        i = xAxis;
                        j = yAxis;
                        k = zAxis;
                        
                        if ismember([i j k], listLS, 'rows')
                            continue
                        else                            
                            delta_sig = -1*spins(i, j, k) - spins(i, j, k);
                            
                            %pick spin and flip right away
                            spins(i, j, k) = -1*spins(i, j, k);
                            
                            sum_nn = sumNN3D(spins, i, j, k);
                            
                            %avg_spin = (sum(spins,'all'))/(N*M*D);
                            
                            %then do change in energy with correct sign
                            %dE = 2*spins(i,j) * (J*sum_nn + H*mu);
                            %dE = delta_sig*(-1*J*sum_nn +...
                            %    (big_delta/2 - T*ln_g/2));
                            
                            dE = delta_sig*(-1*J*sum_nn + longRange);
                            
                            p = exp(-1*dE/T);
                            
                            if dE < 0 || p > rand()
                                continue
                            else
                                spins(i,j,k) = -1*spins(i,j,k);
                            end
                            
                            %then check with random number; if state is acceptable,
                            %keep and move on; if state is not acceptable then flip
                            %sign back and try a different spin
                        end
                    end
                end
            end
        end
        %%
        %Snn = nearestN3D(spins);
        
        %E(idx, 1) = -J*Snn;
        nHS(idx, 1) = n_HSfrac3D(spins);
        
        %%{
        if (mod(idx, frameRate) == 0) && saveIntResults
            pltTitle = strcat(num2str(N),'spins','_T=',num2str(T),'_', num2str(idx));
            spinVis(spins)
            title(pltTitle)
            %axis equal;
            pause(0.05);
            if saveIntResults
                frame_name = strcat(dir_name,'\frames\',pltTitle,".png");
                saveas(gcf, frame_name)
            end
        
        end
        %%}
    end
    
    %E = mean(E);
    E=0;
    nHS = mean(nHS);
    
end