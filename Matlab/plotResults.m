%% PLOTTING
legArr = {strcat("T Inc, J=",J_nom1, "K, ln(g)=",num2str(ln_g1)),...
    strcat("T Dec, J=", J_nom2, "K, ln(g)=",num2str(ln_g2))};
%legArr = makeLegend(L);
set(0,'DefaultTextInterpreter','none')

if numTrials > 1
    
    meanE = squeeze(mean(E));
    mean_nHS = squeeze(mean(n_HS))';
    
    %%
    close all
    %{
    figure
    plot(T, meanE', "*-")
    hold on
    title("Energy vs Temperature")
    xlabel("Temperature T (K)")
    ylabel("Energy")
    hold off
    saveas(gcf, strcat(dir_name,'\',dat_str,'_',num2str(N),'netEvsT','.png'))
    %}
    plt_title = strcat('\rm ',' \Delta=',bD_nom1,'K');
    named = strcat(trial_dir,'\',dat_str0,'_',...
        'nHSvsT','_Jinc',J_nom1,'K_Jdec',J_nom2,'K_D',bD_nom1,'K_lnginc',num2str(ln_g1),'_lngdec',num2str(ln_g2));
    
    figure
    plot(T_inv1, mean_nHS,'b*-')
    hold on
    title(plt_title, 'Interpreter', 'tex')
    xlabel("Temperature T (K)")
    ylabel("n_H_S", 'Interpreter','tex')
    axis([-inf inf 0 1.01])
    legend(legArr,'Location','southeast')
    axis([-inf inf 0 1.0])
    hold off
    
    saveas(gcf, strcat(named, '.png'))
    T_out = [T_inv1 T_inv2];
    nHS_out = [n_HS1 n_HS2];
    
    writematrix([T_out nHS_out],strcat(named, '.txt'));
else
    n_HS1 = squeeze(n_HS1);
    n_HS2 = squeeze(n_HS2);
    
    nom = strcat(trial_dir,'\',dat_str0,'_',...
        'nHSvsT','_J',J_nom1,'K_D',bD_nom1,'K_pLS',num2str(pLS1),...
        '_pHS',num2str(pHS1),'_L',num2str(L));
    
    named = strcat(trial_dir,'\',dat_str0,'_',...
        'nHSvsT','_Jinc',J_nom1,'K_Jdec',J_nom2,'K_D',bD_nom1,'K_lnginc',num2str(ln_g1),'_lngdec',num2str(ln_g2));
    
    
    %figure
    %plot(T, E)
    %title("E")
    plt_title = strcat('\rm \Delta=',bD_nom1,'K \n',...
        num2str(L), 'x', num2str(L),'x', num2str(D), ' Lattice');
    figure
    plot(T_inv1, n_HS1,'r.-')
    hold on
    plot(T_inv2, n_HS2,'b.-')
    grid on
    title(plt_title, 'interpreter','tex')
    xlabel("Temperature T (K)")
    ylabel("n_H_S",'Interpreter','tex')
    axis([-inf inf 0 1.01])
    legend(legArr,'Location','southeast')
    hold off
    
    
    saveas(gcf, strcat(named, 'a.png'))
    T_out = [T_inv1 T_inv2];
    nHS_out = [n_HS1 n_HS2];
    
    writematrix([T_out; nHS_out]',strcat(named, '.txt'));
    
    plot_nHSlayers(nHS_by_layer1, nHS_by_layer2, T_inv1, T_inv2)
    saveas(gcf, strcat(named, 'BYLAYER.png'))
    plot_nHSlayers_hyst(nHS_by_layer1, nHS_by_layer2, T_inv1, T_inv2);
    saveas(gcf, strcat(named, 'BYLAYER_hyst.png'))
    
end

function legArr = makeLegend(L,D)
%returns a legend given an array of lattice sizes
legArr = cell(max(size(L)),1);
for idx = 1:max(size(L))
    legArr{idx} = strcat(num2str(L(idx)),'x',num2str(L(idx)),'x',num2str(D),' spins');
end
end
