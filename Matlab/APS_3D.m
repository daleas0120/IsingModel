% Monte Carlo Ising Model
% Ashley Dale
% Calls the following matlab files: initializeLattice.m,
% equilibrateSpins_H.m

%%
tic
clear;

%% lattice definitions
lattice = struct;
%lattice.weights = [1 0.5 0.3333];
lattice.weights = [1 0 0];
%L = [164];
lattice.L = [38];
%D = 29;
lattice.D = 14;
lattice.numPinnedLayers=0;
lattice.pinningVal = 1;
lattice.pLS = 0; %percentage of interior spins locked in LS
lattice.pHS = 0; %percentage of interior spins locked in HS
lattice.boundCond = (0); %boundary condition

ising = struct;
ising.k_b = 8.617333262*10^-5;%[eV/K]
ising.mu_B = 5.7883818060*10^-5; %[eV/T]
ising.mu_atom = ising.mu_B*4.89; %[eV/T]atomic magnetic moment for Fe2+
ising.J_K = 40;%
%T_K = [100:10:250 255:5:300 310:10:400];%K
ising.T_K = 273;
ising.big_delta_K = 2450;%K
ising.ln_g = 81.9/8.31;
ising.G = 0;%K
ising.h_field = [0:5:120]; % [T] applied magnetic field
ising.omega = 1;
ising.bD_nom = num2str(ising.big_delta_K);
ising.J_nom = num2str(ising.J_K);

ising.J_ev = ising.J_K*ising.k_b; %coupling constant/exchange energy in eV
ising.T_ev = ising.T_K.*ising.k_b; % [K]*[eV/K] = [eV]
ising.bD_ev = ising.big_delta_K*ising.k_b;
ising.k = ising.J_ev./(ising.k_b.*ising.T_K); % dimensionless inverse temperature
ising.H = 0; %external magnetic field

%% DIMENSIONLESS UNITS

ising.big_delta = (ising.k_b*ising.big_delta_K)/ising.J_ev;
ising.T = (ising.k_b.*ising.T_K)./ising.J_ev;
ising.J = ising.J_ev/ising.J_ev;
ising.G = (ising.k_b*ising.G)/ising.J_ev;
ising.T_inv = (ising.J_ev.*ising.T_K)./ising.k_b;
ising.mu = ising.mu_atom/ising.J_ev;
%%
simulation = struct;
simulation.evo = 5e1; %number of MC steps to let the system burn in; this is discarded
simulation.dataPts = 5e1; %number of MC steps to evaluate the system
simulation.frameRate = 100; % provides a modulus to save snapshot of system
simulation.numTrials = 1; %number of times to repeat the experiment

% naming system for the files and folders holding data from repeated trials
simulation.p_name = {'a_', 'b_', 'c_', 'd_', 'e_', 'f_', 'g_', 'h_', 'i_', 'j_', 'k_',...
    'l_', 'm_', 'n_', 'o_', 'p_', 'q_', 'r_', 's_', 't_', 'u_', 'v_', 'w_', ...
    'x_', 'y_', 'z_', 'A_', 'B_', 'C_', 'D_'};

% save intermediate results:
simulation.saveIntResults = 1;

%Energy output variables
simulation.E = zeros(1, length(ising.T_K));
simulation.Snn = zeros(1, length(ising.T_K));

%Magnetism output variables
simulation.B = zeros(1, length(ising.T_K));

%Spin fraction output variables
simulation.nHS = zeros(length(ising.T_K), simulation.dataPts);
simulation.nHS_evo = zeros(length(ising.T_K), simulation.evo);

simulation.APSslideColor = [34/255, 42/255, 53/255];
set(0,'DefaultFigureColor',simulation.APSslideColor)

%%
for p = 1:simulation.numTrials
    t = datetime('now');
    t.Format = "yyMMdd";
    dat_str0 = string(t);
    
    if simulation.numTrials>1 || simulation.saveIntResults
        % save all trials in a single directory at highest level
        tryName = num2str(simulation.numTrials);
        trial_dir = strcat(dat_str0,'_',tryName,'trialRuns');
        mkdir(trial_dir)
    else
        % no group directory required
        trial_dir = '..';
    end
    %%
    for numSpins = 1:length(lattice.L)% square root of number of spins
        
        lattice.N = lattice.L(numSpins);
        
        %results folder for this particular data run
        t = datetime('now');
        t.Format = "yyMMddhhmm";
        simulation.dat_str = string(t);
        if simulation.saveIntResults
            simulation.dir_name = strcat(trial_dir,'/',simulation.dat_str,simulation.p_name{p},'_',num2str(lattice.N),'spins3D');
            mkdir(simulation.dir_name)
            mkdir(simulation.dir_name,'frames')
            mkdir(simulation.dir_name,'png')
            mkdir(simulation.dir_name,'fig')
            mkdir(simulation.dir_name,'txt')
        else
            simulation.dir_name = "";
        end
        
        %% initialize 3D lattice
        [lattice.spinsOG, lattice.listLS] = initializeLattice3D_pin(lattice);
        
        % View initial lattice
        %%{
        figure
        spinVis(lattice.spinsOG);
        set(gca,'xticklabel',[]);
        set(gca,'yticklabel',[]);
        set(gca,'zticklabel',[]);
        set(gca,'xtick',[]);
        set(gca,'ytick',[]);
        set(gca,'ztick',[]);
        set(gca, 'Color', simulation.APSslideColor);
        pause(3);
        %close;
        %}
        %%
        
        figure;

        spins = lattice.spinsOG;
        %[rootName, simulation] = temperature_experiment(p, numSpins, spins, ising, lattice, simulation);
        [rootName, simulation] = magField_experiment(p, numSpins, spins, ising, lattice, simulation);

    end
end

save(strcat(simulation.dir_name,'/', rootName))
close all
%% PLOTTING

%legArr = makeLegend(lattice.L, lattice.D);
set(0,'DefaultTextInterpreter','none')

%Plot nHS vs Temperature
%[nHSmean, imgHSvTemp] = plot_nHSvTemp(ising, simulation);

% Plot nHS vs Mag Field
[nHSmean, imgHSvMag, trendline_fit, gof]= plot_nHSvMag(ising, simulation);


%Plot nHS vs steps
imgHSvTime = plot_nHSsteps(ising, simulation);

% Plot nHS evolution
plot_nHS(ising, simulation);

% Write Results
%writematrix([ising.T_K' nHSmean], strcat(imgHSvTemp,'.txt'));

% writematrix([(1:simulation.dataPts)' simulation.nHS(idx,:)'], strcat(imgHSvTime,'.txt'));
% 
% writematrix([(1:simulation.dataPts)' simulation.nHS(idx,:)'], ...
%     strcat(simulation.dir_name,'/',simulation.dat_str,'timeEvo_',...
%     'delt',ising.bD_nom,'_J',ising.J_nom,'_nHSvsTime','.txt'));


function legArr = makeLegend(L,D)
%returns a legend given an array of lattice sizes
legArr = cell(max(size(L)),1);
for idx = 1:max(size(L))
    legArr{idx} = strcat(num2str(L(idx)),'x',num2str(L(idx)),'x',num2str(D),' spins');
end
end

