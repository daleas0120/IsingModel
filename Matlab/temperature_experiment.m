function [rootName, sim] = temperature_experiment(p, numSpins, spins, ising, lattice, sim)
k = ising.k;
N = lattice.N;
D = lattice.D;
T_K = ising.T_K;
evo = sim.evo;
omega = ising.omega;
mu = ising.mu;
h_field = ising.h_field;
weights = lattice.weights;
J = ising.J;
big_delta = ising.big_delta;
ln_g = ising.ln_g;
listLS = lattice.listLS;
T = ising.T;
for temp = 1:length(k)
    
    %copy spins for later comparison
    %spins_last = spins;
    
    %let state reach equilibrium
    X = sprintf('Cooling %d x %d x %d spins to temp %f ....',...
        N, N, D, T_K(temp));
    fprintf(X)
    
    tic
    [spins, ~, sim.nHS_evo(temp, :)] = equilibrateSpins_3D(...
        evo, spins, k(temp), T(temp), omega, mu, h_field, weights, J, ...
        big_delta, ln_g, listLS, ...
        sim.frameRate, sim.dir_name, sim.saveIntResults);
    
    
    %take data
    fprintf("Taking Data\n")
    [spins, sim.E(p, temp, numSpins), sim.nHS(temp, :)] = ...
        equilibrateSpins_3D(...
        sim.dataPts, spins, k(temp), T(temp), omega, mu, h_field, weights, J, ...
        big_delta, ln_g, listLS, ...
        sim.frameRate, sim.dir_name, 'false');
    
    close;
    %%{
    rootName = strcat(sim.dat_str, sim.p_name{p}, num2str(N),...
        'spins_k_', num2str(T_K(temp)), 'K');
    if sim.saveIntResults == 1
        
        file_name = strcat(dir_name,'/txt/',rootName,'.txt');
        png_name = strcat(dir_name,'/png/',rootName, '.png');
        fig_name = strcat(dir_name,'/fig/',rootName, '.fig');
        %save spin matrix to text file
        writematrix(spins,file_name);
        
        %save final spin
        f = figure;
        spinVis(spins);
        title({strcat(num2str(T_K(temp)), 'K')},...
            'Interpreter', 'tex',...
            'Color', 'white');
        set(gca,'xticklabel',[]);
        set(gca,'yticklabel',[]);
        set(gca,'zticklabel',[]);
        set(gca,'xtick',[]);
        set(gca,'ytick',[]);
        set(gca,'ztick',[]);
        set(gca, 'Color', APSslideColor);
        set(gcf, 'Color', APSslideColor);
        set(gcf, 'InvertHardcopy', 'off');
        saveas(gcf, png_name);
        saveas(gcf, fig_name);
        close
    end
    %}
    toc
end
end
