% Monte Carlo Ising Model
% Ashley Dale
% Calls the following matlab files: initializeLattice.m,
% equilibrateSpins_H.m

%%
tic
clear;

%% LATTICE PARAMETERS
probLock = 0.575; %percentage of spins locked
lock = (-1); %Locked in LS or HS
boundCond = (-1); %boundary condition

%L = [4, 7, 10, 40];
L = [10];
D = 3;

%% MOLECULE PARAMETERS

%% Way UP (LS to HS)
J1 = 170;%
T1 = [300:5:400];%K
big_delta1 = bd;%K
%ln_g1 = 44.7/8.31; %ratio of degeneracy HS to LS
ln_g1 = 42.75/8.31;
G1 = 0;%K
H1 = 0; %external magnetic field

pLS1 = 0; %percentage of interior spins locked in LS
pHS1 = 0; %percentage of interior spins locked in HS
boundCond1 = (0); %boundary condition

%% WAY DOWN (HS to LS)
J2 = 125;%
T2 = [400:-5:300];%K
big_delta2 = bd;%K
%ln_g2 = 47.4/8.31;
ln_g2 = 49.8/8.31; %ratio of degeneracy HS to LS
G2 = 0;%K
H2 = 0; %external magnetic field

pLS2 = 0; %percentage of interior spins locked in LS
pHS2 = 0; %percentage of interior spins locked in HS
boundCond2 = (0); %boundary condition

%% SIMULATION PARAMETERS

evo = 50e1; %number of MC steps to let the system burn in; this is discarded
dataPts = 50e1; %number of MC steps to evaluate the system
numTrials = 1; %number of times to repeat the experiment
frameRate = 10; % provides a modulus to save snapshot of system

k_b = 8.617333262*10^-5;%eV/K
mu = 1; %atomic magnetic moment

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

% save intermediate results:
saveIntResults = false;

%Energy output variables
E = zeros(1, length(k));
Snn = zeros(1, length(k));

%Magnetism output variables
B = zeros(1, length(k));

%Spin fraction output variables
n_HS = zeros(1, length(k));




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
        [spins, listLS] = initializeLattice3D(...
            N, D, boundCond, lock, probLock); %randomly initializes 3D lattice
        
        origSpins = spins;
        
        % View initial lattice
        %{
        figure
        spinVis(spins)
        axis equal
        pause(1)
        close
        %}
        figure;
        
        for temp = 1:length(k)
            %copy spins for later comparison
            spins_last = spins;
            
            %let state reach equilibrium
            X = sprintf('Cooling %d x %d x %d spins to temp %f ....',...
                N, N, D, T_inv(temp));
            disp(X)
         
            [spins, ~, ~] = equilibrateSpins_3D(...
                evo, spins, k(temp), T(temp), mu, H, J, big_delta, ln_g, listLS, ...
                frameRate, dir_name, saveIntResults);
            
            %take data
            fprintf("Taking Data\n")
            [spins, E(p, temp, numSpins), n_HS(p, temp, numSpins)] = ...
                equilibrateSpins_3D(...
                dataPts, spins, k(temp), T(temp), mu, H, J, ...
                big_delta, ln_g, listLS, ...
                frameRate, dir_name, saveIntResults);
            
            close;
            %{
            if saveIntResults
                
                file_name = strcat(dir_name,'\',dat_str, p_name{p}, num2str(N),...
                    'spins_k_', num2str(T(temp)), 'K.txt');
                image_name = strcat(dir_name,'\',dat_str, p_name{p},num2str(N),...
                    'spins_k_', num2str(T(temp)), 'K.png');
                
                %save spin matrix to text file
                writematrix(spins,file_name);
                
                %save final spin
                figure;
                spinVis(spins)
                axis square;
                saveas(gcf, image_name);
                close
            end
            %}
            toc
        end
    end
end

%% PLOTTING

legArr = makeLegend(L, D);
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
    
    plt_title = strcat('\rm ','J=',J_nom, 'K and ',' \Delta=',bD_nom,'K');
    
    figure
    plot(T_inv, mean_nHS,'*-k')
    hold on
    title(plt_title, 'Interpreter', 'tex')
    xlabel("Temperature T (K)")
    ylabel("n_H_S", 'Interpreter','tex')
    legend(legArr,'Location','southeast')
    axis([-inf inf 0 1.0])
    hold off
    saveas(gcf, strcat(trial_dir,'\',dat_str0,'_',...
        'delt',bD_nom,'_J',J_nom,'_nHSvsT','.png'))
    
    writematrix([T_inv' mean_nHS], strcat(trial_dir,'\',dat_str0,'_',...
        'delt',bD_nom,'_J',J_nom,'_nHSvsT','.txt'));
else
    
    %figure
    %plot(T, E)
    %title("E")
    
    figure
    plt_title = strcat('\rm ','J=',J_nom, 'K and ',' \Delta=',bD_nom,'K');
    
    plot(T_inv', n_HS(1,:), '*-k')
    hold on
    title(plt_title,'Interpreter', 'tex')
    xlabel("Temperature T (K)")
    ylabel("n_H_S"','Interpreter','tex')
    legend(legArr,'Location','southeast')
    axis([-inf inf 0 1.0])
    hold off
    saveas(gcf, strcat(trial_dir,'\',dat_str0,'_',...
        'delt',bD_nom,'_J',J_nom,'_nHSvsT','.png'))
    
    writematrix([T_inv' n_HS(1,:)'], strcat(trial_dir,'\',dat_str0,'_',...
        'delt',bD_nom,'_J',J_nom,'_nHSvsT','.txt'));
end

function legArr = makeLegend(L,D)
%returns a legend given an array of lattice sizes
legArr = cell(max(size(L)),1);
for idx = 1:max(size(L))
    legArr{idx} = strcat(num2str(L(idx)),'x',num2str(L(idx)),'x',num2str(D),' spins');
end
end
