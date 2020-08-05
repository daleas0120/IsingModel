% Monte Carlo Ising Model
% Ashley Dale
% Calls the following matlab files: initializeLattice.m,
% equilibrateSpins_H.m

%%
tic
clear;

k_b = 8.617333262*10^-5;%eV/K 
mu = 1; %atomic magnetic moment

%% SIMULATION PARAMETERS

evo = 50e0; %number of MC steps to let the system burn in; this is discarded
dataPts = 50e0; %number of MC steps to evaluate the system
numTrials = 1; %number of times to repeat the experiment
frameRate = 10; % provides a modulus to save snapshot of system
saveIntResults = false;% save intermediate results:

%% LATTICE PARAMETERS
boundCond = (0); %boundary condition
%L = [4, 7, 10, 40];
L = [37];
D = 37;

%% MOLECULE PARAMETERS
bd = 1325;

%% Way UP (LS to HS)
J1 = 10;%
%T1 = 2200;
T1 = [0 80:20:320 2200];%K
big_delta1 = bd;%K
%ln_g1 = 44.7/8.31; %ratio of degeneracy HS to LS
ln_g1 = 67.5/8.31;
G1 = 0;%K
H1 = 0; %external magnetic field

pLS1 = 0; %percentage of interior spins locked in LS
pHS1 = 0; %percentage of interior spins locked in HS
boundCond1 = (0); %boundary condition

%% WAY DOWN (HS to LS)
J2 = 10;%
%T2 = 0;
T2 = [2200 320:-20:80 0];%K
big_delta2 = bd;%K
%ln_g2 = 47.4/8.31;
ln_g2 = 75/8.31; %ratio of degeneracy HS to LS
G2 = 0;%K
H2 = 0; %external magnetic field

pLS2 = 0; %percentage of interior spins locked in LS
pHS2 = 0; %percentage of interior spins locked in HS
boundCond2 = (0); %boundary condition

%% CONVERSION TO DIMENSIONLESS UNITS (Formatting for program)

k_b = 8.617333262*10^-5;%eV/K
mu = 1; %atomic magnetic moment

bD_nom1 = num2str(big_delta1);
J_nom1 = num2str(J1);

J_ev1 = J1*k_b; %coupling constant/exchange energy in eV
T_ev1 = T1.*k_b;
bD_ev1 = big_delta1*k_b;
G_ev1 = G1*k_b;

k1 = J_ev1./(k_b.*T1); % dimensionless inverse temperature

big_delta1 = (k_b*big_delta1)/abs(J_ev1);
T1 = (k_b.*T1)./abs(J_ev1);
J1 = J_ev1/abs(J_ev1);
G1 = G_ev1/abs(J_ev1);
T_inv1 = (abs(J_ev1).*T1)./k_b;

bD_nom2 = num2str(big_delta2);
J_nom2 = num2str(J2);

J_ev2 = J2*k_b; %coupling constant/exchange energy in eV
T_ev2 = T2.*k_b;
bD_ev2 = big_delta2*k_b;
G_ev2 = G2*k_b;

k2 = J_ev2./(k_b.*T2); % dimensionless inverse temperature

big_delta2 = (k_b*big_delta2)/abs(J_ev2);
T2 = (k_b.*T2)./abs(J_ev2);
J2 = J_ev2/abs(J_ev2);
G2 = G_ev2/abs(J_ev2);
T_inv2 = (abs(J_ev2).*T2)./k_b;

% naming system for the files and folders holding data from repeated trials
p_name = {'a_', 'b_', 'c_', 'd_', 'e_', 'f_', 'g_', 'h_', 'i_', 'j_', 'k_',...
    'l_', 'm_', 'n_', 'o_', 'p_', 'q_', 'r_', 's_', 't_', 'u_', 'v_', 'w_', ...
    'x_', 'y_', 'z_', 'A_', 'B_', 'C_', 'D_'};
%Energy output variables
E = zeros(1, length(k1));
Snn = zeros(1, length(k1));

%Magnetism output variables
B = zeros(1, length(k1));
%%
%Spin fraction output variables
n_HS1 = zeros(1, length(k1));
n_HS2 = zeros(1, length(k2));

nHS_by_layer1 = zeros(D, length(k1));
nHS_by_layer2 = zeros(D, length(k2));


%%
for p = 1:numTrials
    
    %%
    if numTrials>1 || ~saveIntResults
        % save all trials in a single directory at highest level
        t = datetime('now');
        t.Format = "yyMMdd";
        tryName = num2str(numTrials);
        dat_str0 = string(t);
        trial_dir = strcat('..\..\',dat_str0,'_',tryName,'_3DtrialRuns');
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
            dir_name = strcat(trial_dir,'\',dat_str,p_name{p},'_',num2str(N),'spins3D');
            mkdir(dir_name)
            mkdir(dir_name,'frames')
        else
            dir_name = "";
        end
        
        %initialize 3D lattice
        
        %[spins, listLS] = initializeLattice3D_ones(...
        %    N, D, boundCond, pLS1, pHS1);
        
        %[spins, listLS] = initializeLattice3D(...
        %    N, D, boundCond, pLS1, pHS1); %randomly initializes 3D lattice
        
        [spins, listLS] = initializeLattice3D_pin(...
            N, D, boundCond, pLS1, pHS1, 1); 
        
        origSpins = spins;
        
        % View initial lattice
        %%{
        figure
        spinVis(spins)
        axis equal
        pause(1)
        close
        %}
        figure;
        %%
        for temp = 1:length(k1)
            %copy spins for later comparison
            spins_last = spins;
            
            %let state reach equilibrium
            X = sprintf('Cooling %d x %d x %d spins to temp %f ....',...
                N, N, D, T_inv1(temp));
            disp(X)
            
            [spins, ~, ~] = equilibrateSpins_3D(...
                evo, spins, k1(temp), T1(temp), mu, H1, J1,...
                big_delta1, ln_g1, listLS,...
                frameRate, dir_name, saveIntResults);
            
            %take data
            fprintf("Taking Data\n")
            [spins, E(p, temp, numSpins), n_HS1(p, temp, numSpins)] = ...
                equilibrateSpins_3D(...
                dataPts, spins, k1(temp), T1(temp), mu, H1, J1, ...
                big_delta1, ln_g1, listLS, ...
                frameRate, dir_name, saveIntResults);
            
            nHS_by_layer1(:,temp) = nHS_layer(spins);
            
            %close;
            %spinVis(spins)
            %pause(1)
            %axis equal
            toc
        end
        %%
        %%{
        for temp = 1:length(k2)
            %copy spins for later comparison
            spins_last = spins;
            
            %let state reach equilibrium
            X = sprintf('Cooling %d x %d x %d spins to temp %f ....',...
                N, N, D, T_inv2(temp));
            disp(X)
            
            [spins, ~, ~] = equilibrateSpins_3D(...
                evo, spins, k2(temp), T2(temp), mu, H2, J2, big_delta2, ln_g2, listLS, ...
                frameRate, dir_name, saveIntResults);
            
            %take data
            fprintf("Taking Data\n")
            [spins, E(p, temp, numSpins), n_HS2(p, temp, numSpins)] = ...
                equilibrateSpins_3D(...
                dataPts, spins, k2(temp), T2(temp), mu, H2, J2, ...
                big_delta2, ln_g2, listLS, ...
                frameRate, dir_name, saveIntResults);
            
            nHS_by_layer2(:, temp) = nHS_layer(spins);
            
            %spinVis(spins)
            %axis equal
            %pause(1)
            
            toc
        end
        %}
        
    end
end
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
    plt_title = strcat('\rm \Delta=',bD_nom1,'K');
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
    
    
   % saveas(gcf, strcat(named, '.png'))
    T_out = [T_inv1 T_inv2];
    nHS_out = [n_HS1 n_HS2];
    
    writematrix([T_out; nHS_out]',strcat(named, '.txt'));
    
    plot_nHSlayers(nHS_by_layer1, nHS_by_layer2, T_inv1, T_inv2)
    %saveas(gcf, strcat(named, 'BYLAYER.png'))
    
end

function legArr = makeLegend(L,D)
%returns a legend given an array of lattice sizes
legArr = cell(max(size(L)),1);
for idx = 1:max(size(L))
    legArr{idx} = strcat(num2str(L(idx)),'x',num2str(L(idx)),'x',num2str(D),' spins');
end
end
