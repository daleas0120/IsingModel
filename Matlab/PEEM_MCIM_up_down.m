% Monte Carlo Ising Model
% Ashley Dale
% Calls the following matlab files: initializeLattice.m,
% equilibrateSpins_H.m

addpath('utils/')

%%
tic
%clear;
bd = 3001;
%L = [209, 311];
%L = [1002, 1002];
L = [300, 300];
wts = [1, 0.7071, 0.5, 0, 0, 0];
%wts = [1, 0.5, 0.25, 0, 0, 0];
%wts = [1, 0, 0, 0, 0, 0];

%weights = [1, 0.7071, 0.5, 0.4472];
%wts = [1, 0.7071, 0.5, 0.4472, 0.3536, 0.3333]; %linear
%wts = [1, 0.5, 0.25, 0.125, 0.0625, 0.0313]; %quadratic
%weights = [1, 0, 0, 0];
%weights = [1 1 1 1];
%wts = wts./sum(wts);

%% Way UP (LS to HS)
J1 = 48.5;%
%T1 = 297;%K
T1 = [285:1:300];
%T1 = [200:10:400];
big_delta1 = bd;%K
%ln_g1 = 44.7/8.31; %ratio of degeneracy HS to LS
S1 = 83.9;
ln_g1 = S1/8.31;
G1 = 0;%K
H1 = 0; %external magnetic field

pLS1 = 0; %percentage of interior spins locked in LS
pHS1 = 0; %percentage of interior spins locked in HS
boundCond1 = (0); %boundary condition

%% WAY DOWN (HS to LS)
% J2 = J1;%
% T2 = [400:-10:100];%K
% big_delta2 = bd;%K
% %ln_g2 = 47.4/8.31;
% ln_g2 = ln_g1; %ratio of degeneracy HS to LS
% G2 = 0;%K
% H2 = 0; %external magnetic field
%
% pLS2 = 0; %percentage of interior spins locked in LS
% pHS2 = 0; %percentage of interior spins locked in HS
% boundCond2 = (0); %boundary condition

%%
evo = 250; %number of MC steps to let the system burn in; this is discarded
dataPts = 250; %number of MC steps to evaluate the system
numTrials = 1; %number of times to repeat the experiment
frameRate = 310001; % provides a modulus to save snapshot of system

%% DIMENSIONLESS UNITS (Formatting for program)

mu = 1; %atomic magnetic moment
k_b = 8.61733326200000e-05; %eV/K
bD_nom1 = num2str(big_delta1);
J_nom1 = num2str(J1);

J_ev1 = J1*k_b; %coupling constant/exchange energy in eV
T_ev1 = T1.*k_b;
bD_ev1 = big_delta1*k_b;
G_ev1 = G1*k_b;

k1 = J_ev1./(k_b.*T1); % dimensionless inverse temperature

% These are the values that actually get passed in
big_delta1 = (k_b*big_delta1)/abs(J_ev1);
T1 = (k_b.*T1)./abs(J_ev1);
J1 = J_ev1/abs(J_ev1);
G1 = G_ev1/abs(J_ev1);
T_inv1 = (abs(J_ev1).*T1)./k_b;

bD_Tnorm = bd./T1;
J_Tnorm = J1./T1;

%%

% bD_nom2 = num2str(big_delta2);
% J_nom2 = num2str(J2);
%
% J_ev2 = J2*k_b; %coupling constant/exchange energy in eV
% T_ev2 = T2.*k_b;
% bD_ev2 = big_delta2*k_b;
% G_ev2 = G2*k_b;
%
% k2 = J_ev2./(k_b.*T2); % dimensionless inverse temperature
%
% big_delta2 = (k_b*big_delta2)/abs(J_ev2);
% T2 = (k_b.*T2)./abs(J_ev2);
% J2 = J_ev2/abs(J_ev2);
% G2 = G_ev2/abs(J_ev2);
% T_inv2 = (abs(J_ev2).*T2)./k_b;
%%

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
n_HS1 = zeros(dataPts, length(k1));
%n_HS2 = zeros(dataPts, length(k2));

%%
for p = 1:numTrials
    
    if numTrials>1 || ~saveIntResults
        % save all trials in a single directory at highest level
        t = datetime('now');
        t.Format = "yyMMdd";
        tryName = num2str(numTrials);
        dat_str0 = string(t);
        trial_dir = strcat('..\..\',dat_str0,'_',tryName,'trialRuns');
        mkdir(trial_dir)
    else
        % no group directory required
        trial_dir = '';
    end
    
    %%
    
    for numSpins = 1:(length(L)-1) % square root of number of spins
        
        N = L(1);
        M = L(2);
        
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
        [spins, locked] = initializeLattice_noPad(N, M, pLS1, pHS1); %randomly initializes 2D lattice
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
            [spins, ~, ~] = equilibrateSpins_periodic(...
                evo, spins, T1(temp), wts, J1,...
                big_delta1, ln_g1, locked,...
                frameRate, dir_name, saveIntResults);
            
            %take data
            fprintf("Taking Data\n")
            
            [spins, E(p, temp, numSpins), n_HS1(:, temp)] = ...
                equilibrateSpins_periodic(...
                dataPts, spins, T1(temp), wts, J1, ...
                big_delta1, ln_g1, locked, ...
                frameRate, dir_name, saveIntResults);
            
            close;
            toc
            %%
            rSpins = reduceLattice_periodic(spins, 9);
            
            %             imgName = strcat(dat_str0,'_J',J_nom1, 'K_T', num2str(T_inv1(temp)),'K_',...
            %                 num2str(wts(1)), '_', num2str(wts(2)), '_', num2str(wts(3)),'_',...
            %                 num2str(wts(4)), '_', num2str(wts(5)), '_', num2str(wts(6)),...
            %                 '_S',num2str(S1));
            %
            %             saveSpinsImg(spins, strcat(imgName,'_', num2str(p), '_OGup.png'))
            %             saveSpinsImg(rSpins, strcat(imgName,'_', num2str(p),  '_squeeze9up.png'))
            %             saveSpinsImg(imresize(rSpins, [1043 1025]), strcat(imgName,'_', num2str(p), '_squeeze9LGup.png'))
            %
            %             saveSpinImg(rSpins, strcat(imgName,'_', num2str(p),  'v1up.png'));
            %             saveSpinImg(imresize(rSpins, [1043 1025]), strcat(imgName,'_', num2str(p),'v2up.png'));
            %             %%
            %             figure;
            %             subplot(2,2,1)
            %             imagesc(spins);
            %             axis square
            %             title('Binary Lattice')
            %             caxis([-1 1])
            %
            %             subplot(2,2,2)
            %             imagesc(rSpins)
            %             axis square
            %             colorbar
            %             title('Averaged Lattice')
            %             caxis([-1 1])
            %
            %             subplot(2,2,3)
            %             histogram(rSpins)
            %             title('Site Value Distribution')
            %
            %             subplot(2,2,4)
            %             plot(n_HS1)
            %             title('nHS vs Time')
            %             ylim([0 1])
            %             grid on
            %             grid minor
            %
            %             saveas(gcf, strcat(imgName, 'grid_', num2str(p),'_',num2str(T_inv1(temp)), 'up.png'));
            
        end
        
        %%
        %MCIM_2D_coolDown()
        
        
    end
end

%% PLOTTING

%legArr = {strcat("T Inc, J=",J_nom1, "K, ln(g)=",num2str(ln_g1)),...
%    strcat("T Dec, J=", J_nom2, "K, ln(g)=",num2str(ln_g2))};
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
   % named = strcat(trial_dir,'\',dat_str0,'_',...
   %     'nHSvsT','_Jinc',J_nom1,'K_Jdec',J_nom2,'K_D',bD_nom1,'K_lnginc',num2str(ln_g1),'_lngdec',num2str(ln_g2));
    
    named = strcat(dat_str0,'_',...
        'nHSvsT','_Jinc',J_nom1,'K_Jdec',J_nom2,'K_D',bD_nom1,'K_lnginc',...
        num2str(ln_g1),'_lngdec',num2str(ln_g2));
    
    figure
    plot(T_inv1, mean_nHS,'b*-')
    hold on
    title(plt_title, 'Interpreter', 'tex')
    xlabel("Temperature T (K)")
    ylabel("n_H_S", 'Interpreter','tex')
    axis([-inf inf 0 1.01])
    %legend(legArr,'Location','southeast')
    axis([-inf inf 0 1.0])
    hold off
    
    saveas(gcf, strcat(named, '.png'))
    T_out = [T_inv1 T_inv2];
    nHS_out = [n_HS1 n_HS2];
    
    writematrix([T_out nHS_out],strcat(named, '.txt'));
else
    n_HS1 = squeeze(n_HS1);
    %n_HS2 = squeeze(n_HS2);
    
    n_HS1 = mean(n_HS1, 1);
    %n_HS2 = mean(n_HS2, 1);
    
    nom = strcat(trial_dir,'\',dat_str0,'_',...
        'nHSvsT','_J',J_nom1,'K_D',bD_nom1,'K_pLS',num2str(pLS1),...
        '_pHS',num2str(pHS1),'_L',num2str(L));
    
    %named = strcat(trial_dir,'\',dat_str0,'_',...
    %    'nHSvsT','_Jinc',J_nom1,'K_','D',bD_nom1,'K_lnginc',num2str(ln_g1));
    named = strcat('nHSvsT','_Jinc',J_nom1,'K_','D',bD_nom1,'K_lnginc',num2str(ln_g1));
    
    %figure
    %plot(T, E)
    %title("E")
    plt_title = strcat('\rm \Delta=',bD_nom1,'K');
    figure
    plot(T_inv1, n_HS1,'r.-')
    hold on
    %plot(T_inv2, n_HS2,'b.-')
    grid on
    title(plt_title, 'interpreter','tex')
    xlabel("Temperature T (K)")
    ylabel("n_H_S",'Interpreter','tex')
    axis([-inf inf 0 1.01])
    %legend(legArr,'Location','southeast')
    hold off
    
    saveas(gcf, strcat(named, '.png'))
    %T_out = [T_inv1 T_inv2];
    %nHS_out = [n_HS1 n_HS2];
    
    writematrix([T_inv1; n_HS1]',strcat(named, '.txt'));
    
end



function legArr = makeLegend(L)
%returns a legend given an array of lattice sizes
legArr = cell(max(size(L)),1);
for idx = 1:max(size(L))
    legArr{idx} = strcat(num2str(L(idx)),'x',num2str(L(idx)),' spins');
end
end
