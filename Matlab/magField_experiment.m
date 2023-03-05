function [rootName, sim] = magField_experiment(p, numSpins, spins, ising, lattice, sim)
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

nHS_fig = figure

for field = 1:length(h_field)
    
    %copy spins for later comparison
    %spins_last = spins;
    
    %let state reach equilibrium
    X = sprintf('Cooling %d x %d x %d spins to temp %f @ field %f ....',...
        N, N, D, T_K, h_field(field));
    fprintf(X)
    
    tic
    [spins, ~, sim.nHS_evo(field, :)] = equilibrateSpins_3D(...
        evo, spins, k, T, omega, mu, h_field(field), weights, J, ...
        big_delta, ln_g, listLS, ...
        sim.frameRate, sim.dir_name, sim.saveIntResults);
    
    
    %take data
    fprintf("Taking Data\n")
    [spins, sim.E(p, field, numSpins), sim.nHS(field, :)] = ...
        equilibrateSpins_3D(...
        sim.dataPts, spins, k, T, omega, mu, h_field(field), weights, J, ...
        big_delta, ln_g, listLS, ...
        sim.frameRate, sim.dir_name, 0);
    
    set(0,'CurrentFigure',nHS_fig)
    plot_nHSvMag(ising, sim);
    
    %%{
    rootName = strcat(sim.dat_str, sim.p_name{p}, num2str(N),...
        'spins_h_', num2str(h_field(field)), 'B');
    if sim.saveIntResults == 1
        
        file_name = strcat(sim.dir_name,'/txt/',rootName,'.txt');
        png_name = strcat(sim.dir_name,'/png/',rootName, '.png');
        fig_name = strcat(sim.dir_name,'/fig/',rootName, '.fig');
        %save spin matrix to text file
        writematrix(spins,file_name);
        
        %save final spin
%         f = figure;
%         spinVis(spins);
%         title({strcat(num2str(T_K(field)), 'K')},...
%             'Interpreter', 'tex',...
%             'Color', 'white');
%         set(gca,'xticklabel',[]);
%         set(gca,'yticklabel',[]);
%         set(gca,'zticklabel',[]);
%         set(gca,'xtick',[]);
%         set(gca,'ytick',[]);
%         set(gca,'ztick',[]);
%         set(gca, 'Color', APSslideColor);
%         set(gcf, 'Color', APSslideColor);
%         set(gcf, 'InvertHardcopy', 'off');
%         saveas(gcf, png_name);
%         saveas(gcf, fig_name);
%         close
    end
    %}
    toc
end
end
