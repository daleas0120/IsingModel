for temp = 1:length(k2)
    %copy spins for later comparison
    spins_last = spins;
    
    %let state reach equilibrium
    X = sprintf('Cooling %d x %d spins to temp %f ....',N, N, T_inv2(temp));
    disp(X)
    [spins, ~, ~] = equilibrateSpins_periodic(...
        evo, spins, k2(temp), T2(temp), weights, H2, J2,...
        big_delta1, ln_g2, G2, locked,...
        frameRate, dir_name, saveIntResults);
    %take data
    fprintf("Taking Data\n")
    [spins, E(p, temp, numSpins), n_HS2(:, temp)] = ...
        equilibrateSpins_periodic(...
        dataPts, spins, k2(temp), T2(temp), weights, H2, J2,...
        big_delta1, ln_g2, G2, locked,...
        frameRate, dir_name, saveIntResults);
    
    close;
    toc
    
    figure;
        subplot(2,2,1)
        imagesc(spins);
        axis square
        title('Binary Lattice')
        caxis([-1 1])
        
        nHS = n_HSfrac(spins)
        
        rSpins = reduceLattice_periodic(spins, 9);
        subplot(2,2,2)
        imagesc(rSpins)
        axis square
        colorbar
        title('Averaged Lattice')
        caxis([-1 1])
        
        subplot(2,2,3)
        histogram(rSpins)
        title('Site Value Distribution')
        
        subplot(2,2,4)
        plot(n_HS2)
        title('nHS vs Time')
        ylim([0 1])
        grid on
        grid minor
        
        saveas(gcf, strcat(trial_dir, '\', 'grid_', num2str(p), 'down.png'));
        
        imgName = strcat(trial_dir,'\',dat_str0,'_J',J_nom1, 'K_T', num2str(T1(temp)),'K_',...
            num2str(weights(1)), '_', num2str(weights(2)), '_', num2str(weights(3)),...
            '_S',num2str(S1));
        
        saveSpinsImg(spins, strcat(imgName,'_', num2str(p), '_OG.png'))
        saveSpinsImg(rSpins, strcat(imgName,'_', num2str(p),  '_squeeze9.png'))
        saveSpinsImg(imresize(rSpins, [1043 1025]), strcat(imgName,'_', num2str(p), '_squeeze9LG.png'))
        
        saveSpinImg(rSpins, strcat(imgName,'_', num2str(p),  'v1.png'));
        saveSpinImg(imresize(rSpins, [1043 1025]), strcat(imgName,'_', num2str(p),'v2.png'));
end