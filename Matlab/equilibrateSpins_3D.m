function [spins, E, nHS] = equilibrateSpins_3D(...
    time, spins, ~, T, omega, mu, h_field, ...
    weights, J, big_delta, ln_g, listLS, ...
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
    f = waitbar(0,'1','Name','Annealing ...',...
        'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
    
    setappdata(f, 'canceling', 0);
    
    set(0,'DefaultTextInterpreter','none')
    %omega = 0.001;
    E = zeros(time, 1);
    nHS = zeros(time, 1);
    H = zeros(time, 1);
    [N, ~, D] = size(spins);
    
    %% some optimization
    
    %longRange = omega*(big_delta/2 - T*ln_g/2 - mu*h_field);
    longRange = omega*(-1*mu*h_field);
    periodic = true;
    %%
    
    for idx = 1:time% how many times to let the system evolve
        waitbar(idx/time, f)
        
        if N > 2
            for xAxis = 2:N-1
                for yAxis = 2:N-1
                    for zAxis = 2:D-1
                        if getappdata(f, 'canceling')
                            break
                        end
                        
                        i = xAxis;
                        j = yAxis;
                        k = zAxis;
                        
                        tmp = find(listLS(:, 1) == i & listLS(:, 2) == j & listLS(:, 3)==k, 1);
                        bool = isempty(tmp);
                        
                        if ~bool
                            continue
                        else
                            
                            delta_sig = -1*spins(i, j, k) - spins(i, j, k);
                            
                            spins(i, j, k) = -1*spins(i, j, k);
                            
                            sum_nn = sumNNN3D(spins, i, j, k, weights, periodic);
                            
                            dE = delta_sig*(-1*J*sum_nn + longRange);
                            
                            p = exp(-1*dE/T);
                            
                            if dE < 0 || p > rand()
                                continue
                            else
                                spins(i,j,k) = -1*spins(i,j,k);
                            end
                            
                        end
                    end
                end
            end
        end
        %%
        
        nHS(idx, 1) = n_HSfrac3D(spins);
        H(idx, 1) = magnetism(spins(2:N-1, 2:N-1, 2:D-1));
        
        %%{
        if ((mod(idx, frameRate) == 0) && saveIntResults)
            
            pltTitle = strcat(num2str(N),'spins','_T=',num2str(T),'_', num2str(idx));
            [~] = spinVis(spins);
            set(gca,'xticklabel',[]);
            set(gca,'yticklabel',[]);
            set(gca,'zticklabel',[]);
            set(gca,'xtick',[]);
            set(gca,'ytick',[]);
            set(gca,'ztick',[]);
            set(gca, 'Color', [34/255, 42/255, 53/255]);
            set(gcf, 'InvertHardcopy', 'off');
            pause(0.05);
            if saveIntResults == 1
                frame_name = strcat(dir_name,'/frames/',pltTitle);
                saveas(gcf, strcat(frame_name,".png"));
                saveas(gcf, strcat(frame_name,".fig"));
            end
            close
            
            figure
            squeezeSpins = squeeze3D(spins);
            [~] = spinVis(squeezeSpins);
            set(gca,'xticklabel',[]);
            set(gca,'yticklabel',[]);
            set(gca,'zticklabel',[]);
            set(gca,'xtick',[]);
            set(gca,'ytick',[]);
            set(gca,'ztick',[]);
            set(gca, 'Color', [34/255, 42/255, 53/255]);
            set(gcf, 'InvertHardcopy', 'off');
            pause(0.05);
            if saveIntResults == 1
                frame_name = strcat(dir_name,'/frames/squeeze',pltTitle);
                saveas(gcf, strcat(frame_name,".png"));
                saveas(gcf, strcat(frame_name,".fig"));
            end
            close
            
        end
        %%}
    end
    E = 0;
    F = findall(0, 'type', 'figure', 'tag', 'TMWWaitbar');
    delete(F);
end