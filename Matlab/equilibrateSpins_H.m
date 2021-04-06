function [spins, E, nHS] = equilibrateSpins_H(...
    time, spins, ~, T, mu, ~, J, big_delta, ln_g, G, listLS, ...
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
%ln_g: log of ratio of degeneracy HS to degeneracy LS
%G:
%listLS:
    
%frameRate: determines how frequently intermediate spin matrices are saved
%to an image
%spins_last: previous spin value
%dir_name: where to save results
%saveIntResults: boolean to control writing of frame samples

    %}
    
    set(0,'DefaultTextInterpreter','none')
    %k_b = 8.617333262*10^-5;%eV/K
    
    E = zeros(time, 1);
    %B = zeros(time, 1);
    nHS = zeros(time, 1);
    [M, N] = (size(spins));
    
    %% some optimization
    longRange = (big_delta/2 - T*ln_g/2);
    
    f = waitbar(0, 'Starting');
    
    
    for idx = 1:time% how many times to let the system evolve
        for row = 2:M-1
            for col = 2:N-1
                i = row;
                j = col;
                
                if N > 2
                    
                    
                    if ismember([i j], listLS, 'rows')
                        continue
                    else
                        %tmp1 = countLS(spins(2:N-1, 2:N-1));
                        %tmp2 = ismember(listLS, tmp1, 'rows');
                        %if max(size(listLS)) ~= max(size(tmp2))
                        %    print("ERROR")
                        %end
                        
                        spinsLast = spins;
                        %tmp3 = countLS(spins(2:N-1, 2:N-1));
                        
                        delta_sig = -1*spins(i,j) - spins(i,j);
                        
                        %pick spin and flip right away
                        spins(i, j) = -1*spins(i,j);
                        
                        %sum_nn = (spins((i-1),j) + spins((i+1),j) +...
                        %    spins(i,(j-1)) + spins(i,(j+1)));
                        
                        sum_nn = sumNNN(spins, i, j);
                        
                        %spin_avg = mean(mean(spins));
                        spin_avg = 1;
                        
                        
                        %then do change in energy with correct sign
                        %dE = 2*spins(i,j) * (J*sum_nn + H*mu);
                        %dE = delta_sig*(-1*J*sum_nn + (big_delta/2 - T*ln_g/2 - G*spin_avg));
                        dE = delta_sig*(-1*J*sum_nn + (big_delta/2 - T*ln_g/2));
                        
                        %E1 = spins(i,j)*J*sum_nn - (big_delta/2 - k_b*T*ln_g/2)*spins(i,j);
                        %E2 = (-1*spins(i,j))*J*sum_nn - (big_delta/2 - k_b*T*ln_g/2)*(-1*spins(i,j));
                        %dE = E2 - E1;
                        
                        p = exp(-1*dE/T);
                        r = rand;
                        
                        if dE < 0 || p >= r
                            continue
                            
                        else
                            spins(i,j) = -1*spins(i,j);
                        end
                        
                        %then check with random number; if state is acceptable,
                        %keep and move on; if state is not acceptable then flip
                        %sign back and try a different spin
                        
                        %tmp4 = countLS(spins(2:N-1, 2:N-1));
                        
                        %if max(size(listLS)) > max(size(tmp4))
                        %    tmp5 = setdiff(listLS, tmp4, 'rows');
                        %    bp = 5;
                        %end
                        
                    end
                    
                    
                end
            end
        end
        
        %%
        %Snn = nearestN(spins);
        %sum_Si = sum(spins(2:N-1, 2:N-1), 'all');
        
        %E(idx, 1) = -J*Snn;
        %B(idx, 1) = mu*sum_Si;
        
        %Collapse lattice
        
        reducedSpins = reduceLattice(spins, 9);
        
        nHS(idx, 1) = n_HSfrac(reducedSpins);
        
        %plot(nHS)
        %pause(0.05);
        
        
        
        %{
    
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
    
        %}
        waitbar(idx/time, f, sprintf('Progress: %d %%', floor(idx/time*100)));
        pause(0.1);
    end
    close(f)
    E = 0;
    %E = mean(E);
    %B = (mean(B)./(N*N));
    %nHS = mean(nHS);
    
end

