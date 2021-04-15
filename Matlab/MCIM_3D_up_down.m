% Monte Carlo Ising Model
% Ashley Dale
% Calls the following matlab files: initializeLattice.m,
% equilibrateSpins_H.m

%%
tic
clear;

k_b = 8.617333262*10^-5;%eV/K 
mu = 1; %atomic magnetic moment
T_c = 4.515;

%% SIMULATION PARAMETERS

evo = 200e0; %number of MC steps to let the system burn in; this is discarded
dataPts = 200e0; %number of MC steps to evaluate the system
numTrials = 1; %number of times to repeat the experiment
frameRate = 1000; % provides a modulus to save snapsh ot of system

saveIntResults = false;% save intermediate results:

%% LATTICE PARAMETERS
boundCond = (0); %boundary condition
%L = [4, 7, 10, 40];
L = [37];
D = 37;
omega = 0.1;

%% MOLECULE PARAMETERS
bd = 1550;

%% Way UP (LS to HS)

%J1 = 20;%
%T1 = 2200;
%T1 = [100:2:200];%K

J1 = 100;%
%T1 = 2200;
T1 = [250:10:750];%K
big_delta1 = bd;%K
%ln_g1 = 44.7/8.31; %ratio of degeneracy HS to LS
%ln_g1 = 67.5/8.31;
%ln_g1 = 48.2/8.31;
ln_g1 = 81.9/8.31;

G1 = 0;%K
H1 = 0; %external magnetic field

pLS1 = 0; %percentage of interior spins locked in LS
pHS1 = 0; %percentage of interior spins locked in HS
boundCond1 = (0); %boundary condition

%% WAY DOWN (HS to LS)

J2 = 20;%
%T2 = 0;
T2 = [200:-2:100];%K
big_delta2 = bd;%K
%ln_g2 = 47.4/8.31;
%ln_g2 = 49.8/8.31; %ratio of degeneracy HS to LS
%ln_g2 = 5.8002;

ln_g2 = ln_g1;
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
H1 = zeros(1, length(k1));
H2 = zeros(1, length(k2));
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
        t.Format = "yyMMddHHmm";
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
        
        [spins, listLS] = initializeLattice3D(...
            N, D, boundCond, pLS1, pHS1); %randomly initializes 3D lattice
        
        %[spins, listLS] = initializeLattice3D_pin(...
        %    N, D, boundCond, pLS1, pHS1, 2); 

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
                omega, big_delta1, ln_g1, listLS,...
                frameRate, dir_name, saveIntResults);
            
            %take data
            fprintf("Taking Data\n")
            [spins, E(p, temp, numSpins), n_HS1(p, temp, numSpins)] = ...
                equilibrateSpins_3D(...
                dataPts, spins, k1(temp), T1(temp), mu, H1, J1, ...
                omega, big_delta1, ln_g1, listLS, ...
                frameRate, dir_name, saveIntResults);
            
            nHS_by_layer1(:,temp) = nHS_layer(spins);
            H1(1, temp) = magnetism(spins(2:(L-1), 2:(L-1), 2:(D-1)));
            
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
                evo, spins, k2(temp), T2(temp), mu, H2, J2,...
                omega, big_delta2, ln_g2, listLS, ...
                frameRate, dir_name, saveIntResults);
            
            %take data
            fprintf("Taking Data\n")
            [spins, E(p, temp, numSpins), n_HS2(p, temp, numSpins)] = ...
                equilibrateSpins_3D(...
                dataPts, spins, k2(temp), T2(temp), mu, H2, J2, ...
                omega, big_delta2, ln_g2, listLS, ...
                frameRate, dir_name, saveIntResults);
            
            nHS_by_layer2(:, temp) = nHS_layer(spins);
            
            H2(1, temp) = magnetism(spins(2:(L-1), 2:(L-1), 2:(D-1)));
            
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
    plt_title = strcat('\rm \Delta=',bD_nom1,'K \n',...
        num2str(L), 'x', num2str(L),'x', num2str(D), ' Lattice');
    figure

    plot(T_inv1, n_HS1,'r.-', 'LineWidth', 1, 'MarkerSize', 10)
    hold on
    plot(T_inv2, n_HS2,'b.-', 'LineWidth', 1, 'MarkerSize', 10)

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
%% Save workspace variables
mat_name = strcat(nom, '.mat');
save(mat_name);

%%
function legArr = makeLegend(L,D)
%returns a legend given an array of lattice sizes
legArr = cell(max(size(L)),1);
for idx = 1:max(size(L))
    legArr{idx} = strcat(num2str(L(idx)),'x',num2str(L(idx)),'x',num2str(D),' spins');
end
end
