function [spins, E, nHS] = equilibrateSpins_fast(...
    time, spins, T, weights, J, big_delta, ln_g, listLS, ...
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
    f = waitbar(0,'1','Name','equilibrateSpins_H',...
        'CreateCancelBtn','setappdata(gcbf,''canceling'',Annealing)');
    
    setappdata(f, 'canceling', 0);
    
    set(0,'DefaultTextInterpreter','none');
    
    %k_b = 8.617333262*10^-5;%eV/K
    
    E = zeros(time, 1);
    %B = zeros(time, 1);
    nHS = zeros(time, 1);
    [~, M] = size(spins);
    
    %% some optimization
    longRange = (big_delta/2 - T*ln_g/2);
    periodic = true;
    
    KERNEL_WTS = ...
        [0.0000 0.0000 0.0000 0.3333 0.0000 0.0000 0.0000;
        0.0000 0.3536 0.4472 0.5000 0.4472 0.3536 0.0000;
        0.0000 0.4472 0.7071 1.0000 0.7071 0.4472 0.0000;
        0.3333 0.5000 1.0000 0.0000 1.0000 0.5000 0.3333;
        0.0000 0.4472 0.7071 1.0000 0.7071 0.4472 0.0000;
        0.0000 0.3536 0.4472 0.5000 0.4472 0.3536 0.0000;
        0.0000 0.0000 0.0000 0.3333 0.0000 0.0000 0.0000];

        % N6 Kernel
%     KERNEL = ...
%         [0 0 0 1 0 0 0;
%         0 0 1 1 1 0 0;
%         0 1 1 1 1 1 0;
%         1 1 1 1 1 1 1;
%         0 1 1 1 1 1 0;
%         0 0 1 1 1 0 0;
%         0 0 0 1 0 0 0];
% 

        KERNEL = ...
            [ 0 0 0 0 0 0 0;
            0 0 0 1 0 0 0;
            0 0 1 1 1 0 0;
            0 1 1 0 1 1 0;
            0 0 1 1 1 0 0;
            0 0 0 1 0 0 0;
            0 0 0 0 0 0 0];


    N = KERNEL.*KERNEL_WTS;
    
    J = J.*ones(M, M);
    
%     spinFig = figure;
%     nHSFig = figure;
%     
    for idx = 1:time% how many times to let the system evolve

        if getappdata(f, 'canceling')
            break
        end
        
        % apply periodic boundary conditions by circularly padding spin
        % array
        
        spinsOG = spins;
        
        % sum over all nearest neighbors
        
        spins_x = [spins spins spins;
            spins spins spins;
            spins spins spins];
        sNN = conv2(spins_x, N, 'same');
        sNN = sNN(M+1:2*M, M+1:2*M);
        
        deltaE = 2.*J.*(spins.*sNN);
        
        p = exp(-deltaE./T);
        r = rand(M);
        q = rand(M);
        A = ((r < p).*(q < 0.1).*(-2) + 1);
        
        spins = spins.*A;
        
%         figure(spinFig)
%         imagesc(spins)
%         axis square
%         xticklabels([])
%         yticklabels([])
%         caxis([-1 1])
        

     
        %%

        %E(idx, 1) = -J*Snn;
        %B(idx, 1) = mu*sum_Si;
        
        nHS(idx, 1) = n_HSfrac(spins);
%         figure(nHSFig)
%         plot(nHS)
        
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
                waitbar(idx/time, f)
        
    end
    delete(f)
    close all
    E = 0;
    %E = mean(E);
    %B = (mean(B)./(N*N));
    %nHS = mean(nHS);
    
end