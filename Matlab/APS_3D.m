% Monte Carlo Ising Model
% Ashley Dale
% Calls the following matlab files: initializeLattice.m,
% equilibrateSpins_H.m

%%
tic
clear;
bd = 2450;
k_b = 8.617333262*10^-5;%eV/K
mu = 1; %atomic magnetic moment
weights = [1 0.5 0.3333];

J_K = 40;%
%T_K = [100:10:250 252:2:300 310:10:400];%K
T_K = 249;
big_delta_K = bd;%K
ln_g = 81.9/8.31;
G = 0;%K

pLS = 0; %percentage of interior spins locked in LS
pHS = 0; %percentage of interior spins locked in HS
boundCond = (0); %boundary condition
omega = 1;


bD_nom = num2str(big_delta_K);
J_nom = num2str(J_K);

J_ev = J_K*k_b; %coupling constant/exchange energy in eV
T_ev = T_K.*k_b;
bD_ev = big_delta_K*k_b;

k = J_ev./(k_b.*T_K); % dimensionless inverse temperature

H = 0; %external magnetic field

%% DIMENSIONLESS UNITS

big_delta = (k_b*big_delta_K)/J_ev;
T = (k_b.*T_K)./J_ev;
J = J_ev/J_ev;
G = (k_b*G)/J_ev;

T_inv = (J_ev.*T_K)./k_b;

%%
evo = 2e2; %number of MC steps to let the system burn in; this is discarded
dataPts = 2.5e1; %number of MC steps to evaluate the system
frameRate = 1; % provides a modulus to save snapshot of system
numTrials = 1; %number of times to repeat the experiment

% naming system for the files and folders holding data from repeated trials
p_name = {'a_', 'b_', 'c_', 'd_', 'e_', 'f_', 'g_', 'h_', 'i_', 'j_', 'k_',...
    'l_', 'm_', 'n_', 'o_', 'p_', 'q_', 'r_', 's_', 't_', 'u_', 'v_', 'w_', ...
    'x_', 'y_', 'z_', 'A_', 'B_', 'C_', 'D_'};

% save intermediate results:
saveIntResults = true;

%Energy output variables
E = zeros(1, length(T_K));
Snn = zeros(1, length(T_K));

%Magnetism output variables
B = zeros(1, length(T_K));

%Spin fraction output variables
nHS = zeros(length(T_K), dataPts);
nHS_evo = zeros(length(T_K), evo);

L = [164];
D = 29;
numPinnedLayers=0;
pinningVal = 1;

APSslideColor = [34/255, 42/255, 53/255];
%%
set(0,'DefaultFigureColor',[34/255, 42/255, 53/255])
for p = 1:numTrials
    t = datetime('now');
    t.Format = "yyMMdd";
    dat_str0 = string(t);
    
    if numTrials>1 || saveIntResults
        % save all trials in a single directory at highest level
        tryName = num2str(numTrials);
        trial_dir = strcat(dat_str0,'_',tryName,'trialRuns');
        mkdir(trial_dir)
    else
        % no group directory required
        trial_dir = '..';
    end
    %%
    for numSpins = 1:length(L)% square root of number of spins
        
        N = L(numSpins);
        
        %results folder for this particular data run
        if saveIntResults
            t = datetime('now');
            t.Format = "yyMMddhhmm";
            dat_str = string(t);
            dir_name = strcat(trial_dir,'/',dat_str,p_name{p},'_',num2str(N),'spins3D');
            mkdir(dir_name)
            mkdir(dir_name,'frames')
            mkdir(dir_name,'png')
            mkdir(dir_name,'fig')
            mkdir(dir_name,'txt')
        else
            dir_name = "";
        end
        
        %% initialize 3D lattice
        [spinsOG, listLS] = initializeLattice3D_pin(...
            N, D, boundCond, pLS, pHS, numPinnedLayers, pinningVal);
        
        % View initial lattice
        %%{
        figure
        spinVis(spinsOG);
        set(gca,'xticklabel',[])
        set(gca,'yticklabel',[])
        set(gca,'zticklabel',[])
        set(gca,'xtick',[])
        set(gca,'ytick',[])
        set(gca,'ztick',[])
        set(gca, 'Color', APSslideColor)
        pause(3)
        close
        %}
        
        figure;
        %%
        spins = spinsOG;
        for temp = 1:length(k)
            %copy spins for later comparison
            %spins_last = spins;
            
            %let state reach equilibrium
            X = sprintf('Cooling %d x %d x %d spins to temp %f ....',...
                N, N, D, T_K(temp));
            disp(X)
            
            tic
            [spins, ~, nHS_evo(temp, :)] = equilibrateSpins_3D(...
                evo, spins, k(temp), T(temp), omega, weights, J, ...
                big_delta, ln_g, listLS, ...
                frameRate, dir_name, saveIntResults);
            
            
            %take data
            fprintf("Taking Data\n")
            [spins, E(p, temp, numSpins), nHS(temp, :)] = ...
                equilibrateSpins_3D(...
                dataPts, spins, k(temp), T(temp), omega, weights, J, ...
                big_delta, ln_g, listLS, ...
                frameRate, dir_name, 'false');
            
            close;
            %%{
            rootName = strcat(dat_str, p_name{p}, num2str(N),...
                'spins_k_', num2str(T_K(temp)), 'K');
            if saveIntResults
                
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
                set(gca,'xticklabel',[])
                set(gca,'yticklabel',[])
                set(gca,'zticklabel',[])
                set(gca,'xtick',[])
                set(gca,'ytick',[])
                set(gca,'ztick',[])
                set(gca, 'Color', APSslideColor)
                set(gcf, 'Color', APSslideColor)
                set(gcf, 'InvertHardcopy', 'off')
                saveas(gcf, png_name);
                saveas(gcf, fig_name);
                close
            end
            %}
            toc
        end
    end
end
%%

save(strcat(dir_name,'/', rootName))

%% PLOTTING

legArr = makeLegend(L, D);
set(0,'DefaultTextInterpreter','none')

if numTrials > 1
    
    meanE = squeeze(mean(E));
    mean_nHS = squeeze(mean(nHS))';
    
    
    %%
    close all
    
    plt_title = 'High Spin Fraction vs Temperature';
    figure
    plot(T_inv, mean_nHS,'*-k')
    hold on
    title({plt_title}, 'Interpreter', 'tex')
    xlabel("Temperature T (K)")
    ylabel("n_H_S", 'Interpreter','tex')
    legend(legArr,'Location','southeast')
    axis([-inf inf 0 1.0])
    hold off
    saveas(gcf, strcat(trial_dir,'/',dat_str0,'_',...
        'delt',bD_nom,'_J',J_nom,'_nHSvsT','.png'))
    
    writematrix([T_inv' mean_nHS], strcat(trial_dir,'/',dat_str0,'_',...
        'delt',bD_nom,'_J',J_nom,'_nHSvsT','.txt'));
else
    nHSmean = mean(nHS, 2);
    
    plt_title = 'High Spin Fraction vs Temperature';
    
    figure
    plot(T_K', nHSmean, '.-c')
    hold on
    title({plt_title},'Interpreter', 'tex', 'Color', 'white')
    xlabel("Temperature T (K)")
    ylabel({'n_H_S'},'Interpreter','tex')
    %legend(legArr,'Location','southeast')
    axis([-inf inf 0 1.0])
    set(gca, 'Color', APSslideColor)
    set(gca, 'XColor', [1, 1, 1])
    set(gca, 'YColor', [1, 1, 1])
    grid on
    hold off
    set(gcf, 'InvertHardcopy', 'off')
    
    imgHSvTemp = strcat(dir_name,'/',dat_str,'_',...
        'delt',bD_nom,'_J',J_nom,'_nHSvsT');
    
    saveas(gcf, strcat(imgHSvTemp,'.png'))
    saveas(gcf, strcat(imgHSvTemp,'.fig'))
    
    writematrix([T_K(:)' nHSmean(:)], strcat(imgHSvTemp,'.txt'));
    
    %%
    %Plot nHS vs steps
    plt_title = 'Spin High Fraction vs Time';
    figure
    hold on
    
    for idx = 1:length(T)
        plot(1:dataPts, nHS(idx, :), '.-c')
    end
    set(gca, 'Color', APSslideColor)
    set(gca, 'XColor', [1, 1, 1])
    set(gca, 'YColor', [1, 1, 1])
    grid on
    ylabel({'n_H_S'},'Interpreter','tex')
    xlabel("Time (MCIMS Step)")
    title({plt_title}, 'Color', 'white')
    set(gcf, 'InvertHardcopy', 'off')
    
    imgHSvTime = strcat(dir_name,'/',dat_str,'timeAvg_',...
        'delt',bD_nom,'_J',J_nom,'_nHSvsTime');
    
    saveas(gcf,strcat(imgHSvTime, '.png'));
    saveas(gcf, strcat(imgHSvTime, '.fig'));
    
    writematrix([(1:dataPts)' nHS(idx,:)'], strcat(imgHSvTime,'.txt'));
    
    %%
    % Plot nHS evolution
    plt_title = 'Spin High Fraction vs Time';
    figure
    hold on
    for idx = 1:length(T)
        plot(1:evo, nHS_evo(idx, :), '.-c')
    end
    set(gca, 'Color', APSslideColor)
    set(gca, 'XColor', [1, 1, 1])
    set(gca, 'YColor', [1, 1, 1])
    grid on
    ylabel({'n_H_S'},'Interpreter','tex')
    xlabel("Time (MCIMS Step)")
    title({plt_title}, 'Color', 'white')
    set(gcf, 'InvertHardcopy', 'off')
    saveas(gcf, strcat(dir_name,'/',dat_str,'timeEvo_',...
        'delt',bD_nom,'_J',J_nom,'_nHSvsT.png'))
    saveas(gcf, strcat(dir_name,'/',dat_str,'timeEvo_',...
        'delt',bD_nom,'_J',J_nom,'_nHSvsT.fig'))
    
    writematrix([(1:dataPts)' nHS(idx,:)'], strcat(dir_name,'/',dat_str,'timeEvo_',...
        'delt',bD_nom,'_J',J_nom,'_nHSvsTime','.txt'));
    
end

function legArr = makeLegend(L,D)
%returns a legend given an array of lattice sizes
legArr = cell(max(size(L)),1);
for idx = 1:max(size(L))
    legArr{idx} = strcat(num2str(L(idx)),'x',num2str(L(idx)),'x',num2str(D),' spins');
end
end
