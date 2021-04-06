%%
tic
%clear;
bd = 1700;
L = [82];

k_b = 8.617333262*10^-5;%eV/K
mu = 1; %atomic magnetic moment

%% Simulation Parameters
J1 = 1;%
T1 = [100:10:400];%K
big_delta1 = bd;%K
%ln_g1 = 44.7/8.31; %ratio of degeneracy HS to LS
ln_g1 = 81.9/8.31;
G1 = 0;%K
H1 = 0; %external magnetic field

pLS1 = 0; %percentage of interior spins locked in LS
pHS1 = 0; %percentage of interior spins locked in HS
boundCond1 = (0); %boundary condition

bD_nom1 = num2str(big_delta1);
J_nom1 = num2str(J1);

J_ev1 = J1*k_b; %coupling constant/exchange energy in eV
T_ev1 = T1.*k_b;
bD_ev1 = big_delta1*k_b;
G_ev1 = G1*k_b;

k1 = J_ev1./(k_b.*T1); % dimensionless inverse temperature

% DIMENSIONLESS UNITS

big_delta1 = (k_b*big_delta1)/abs(J_ev1);
T1 = (k_b.*T1)./abs(J_ev1);
J1 = J_ev1/abs(J_ev1);
G1 = G_ev1/abs(J_ev1);
T_inv1 = (abs(J_ev1).*T1)./k_b;

%%
evo = 50e1; %number of MC steps to let the system burn in; this is discarded
dataPts = 50e1; %number of MC steps to evaluate the system
numTrials = 1; %number of times to repeat the experiment
frameRate = 10; % provides a modulus to save snapshot of system

% naming system for the files and folders holding data from repeated trials
p_name = {'a_', 'b_', 'c_', 'd_', 'e_', 'f_', 'g_', 'h_', 'i_', 'j_', 'k_',...
    'l_', 'm_', 'n_', 'o_', 'p_', 'q_', 'r_', 's_', 't_', 'u_', 'v_', 'w_', ...
    'x_', 'y_', 'z_', 'A_', 'B_', 'C_', 'D_'};

% save intermediate results:
saveIntResults = false;

%Energy output variables
E = zeros(1, length(k1));
Snn = zeros(1, length(k1));

%Magnetism output variables
B = zeros(1, length(k1));

%Spin fraction output variables
n_HS1 = zeros(1, length(k1));


%%
for p = 1:numTrials
    
    if numTrials>1 || ~saveIntResults
        % save all trials in a single directory at highest level
        t = datetime('now');
        t.Format = "yyMMdd";
        tryName = num2str(numTrials);
        dat_str0 = string(t);
        trial_dir = strcat(dat_str0,'_',tryName,'trialRuns');
        mkdir(trial_dir)
    else
        % no group directory required
        trial_dir = '..\..';
    end
    
    %%
    
    for numSpins = 1:length(L)% square root of number of spins
        
        N = L(numSpins);
        
        %results folder for this particular data run
        if saveIntResults
            t = datetime('now');
            t.Format = "yyMMdd";
            dat_str = string(t);
            dir_name = strcat(trial_dir,'\',dat_str,p_name{p},'_',num2str(N),'spins');
            mkdir(dir_name)
            mkdir(dir_name,'frames')
        else
            dir_name = "";
        end
        %initialize 2D lattice
        [spins, listLS] = initializeLattice(N, boundCond1, pLS1, pHS1); %randomly initializes 2D lattice
        
        origSpins = spins;
        
        % View initial lattice
        %{
        figure;
        imagesc(spins)
        axis square
        title("Initial Lattice")
        %}
        figure;
        %%
        for temp = 1:length(k1)
            %copy spins for later comparison
            spins_last = spins;
            
            %let state reach equilibrium
            X = sprintf('Cooling %d x %d spins to temp %f ....',N, N, T_inv1(temp));
            disp(X)
            [spins, ~, ~] = equilibrateSpins_H(...
                evo, spins, k1(temp), T1(temp), mu, H1, J1,...
                big_delta1, ln_g1, G1, listLS,...
                frameRate, dir_name, saveIntResults);
            
            %take data
            fprintf("Taking Data\n")
            [spins, E(p, temp, numSpins), n_HS1(p, temp, numSpins)] = ...
                equilibrateSpins_H(...
                dataPts, spins, k1(temp), T1(temp), mu, H1, J1, ...
                big_delta1, ln_g1, G1, listLS, ...
                frameRate, dir_name, saveIntResults);
            
            close;
            toc
        end
    end
end

%% PLOTTING

legArr = {strcat("T Inc, J=",J_nom1, "K, ln(g)=",num2str(ln_g1))};
%legArr = makeLegend(L);
set(0,'DefaultTextInterpreter','tex')

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
        'nHSvsT','_Jinc',J_nom1,'K_D',bD_nom1,'K_lnginc',num2str(ln_g1));
    
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
    T_out = [T_inv1];
    nHS_out = [n_HS1];
    
    writematrix([T_out nHS_out],strcat(named, '.txt'));
else
    n_HS1 = squeeze(n_HS1);
    
    
    nom = strcat(trial_dir,'\',dat_str0,'_',...
        'nHSvsT','_J',J_nom1,'K_D',bD_nom1,'K_pLS',num2str(pLS1),...
        '_pHS',num2str(pHS1),'_L',num2str(L));
    
    named = strcat(trial_dir,'\',dat_str0,'_',...
        'nHSvsT','_Jinc',J_nom1,'K_D',bD_nom1,'K_lnginc',num2str(ln_g1));
    
    
    %figure
    %plot(T, E)
    %title("E")
    plt_title = strcat('\rm \Delta=',bD_nom1,'K');
    figure
    plot(T_inv1, n_HS1,'r.-')
    
    grid on
    title(plt_title, 'interpreter','tex')
    xlabel("Temperature T (K)")
    ylabel("n_H_S",'Interpreter','tex')
    axis([-inf inf 0 1.01])
    legend(legArr,'Location','southeast')
    hold off
    
    
    saveas(gcf, strcat(named, '.png'))
    T_out = [T_inv1];
    nHS_out = [n_HS1];
    
    writematrix([T_out; nHS_out]',strcat(named, '.txt'));
    save(strcat(named, '.mat'));
end



function legArr = makeLegend(L)
%returns a legend given an array of lattice sizes
legArr = cell(max(size(L)),1);
for idx = 1:max(size(L))
    legArr{idx} = strcat(num2str(L(idx)),'x',num2str(L(idx)),' spins');
end
end
