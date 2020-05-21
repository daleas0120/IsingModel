% Monte Carlo Ising Model
% Ashley Dale
% Calls the following matlab files: initializeLattice.m,
% equilibrateSpins_H.m

%%
tic
%clear;

L = [202];

k_b = 8.617333262*10^-5;%eV/K
mu = 1; %atomic magnetic moment

J = -75;%
T = [100:10:400 400:-10:100];%K
big_delta = 2300;%K
ln_g = 81.9/8.31; %ratio of degeneracy HS to LS
G = 0;%K
H = 0; %external magnetic field

%%
pLS = 0.3; %percentage of interior spins locked in LS
pHS = 0.12; %percentage of interior spins locked in HS
boundCond = (0); %boundary condition

%%
evo = 50e1; %number of MC steps to let the system burn in; this is discarded
dataPts = 50e1; %number of MC steps to evaluate the system
numTrials = 1; %number of times to repeat the experiment
frameRate = 10; % provides a modulus to save snapshot of system


%%
bD_nom = num2str(big_delta);
J_nom = num2str(J);

J_ev = J*k_b; %coupling constant/exchange energy in eV
T_ev = T.*k_b;
bD_ev = big_delta*k_b;
G_ev = G*k_b;

k = J_ev./(k_b.*T); % dimensionless inverse temperature



%% DIMENSIONLESS UNITS

big_delta = (k_b*big_delta)/abs(J_ev);
T = (k_b.*T)./abs(J_ev);
J = J_ev/abs(J_ev);
G = G_ev/abs(J_ev);
T_inv = (abs(J_ev).*T)./k_b;


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
        [spins, listLS] = initializeLattice(N, boundCond, pLS, pHS); %randomly initializes 2D lattice
        
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
        for temp = 1:length(k)
            %copy spins for later comparison
            spins_last = spins;
            
            %let state reach equilibrium
            X = sprintf('Cooling %d x %d spins to temp %f ....',N, N, T_inv(temp));
            disp(X)
            [spins, ~, ~] = equilibrateSpins_H(...
                evo, spins, k(temp), T(temp), mu, H, J,...
                big_delta, ln_g, G, listLS,...
                frameRate, dir_name, saveIntResults);
            
            %take data
            fprintf("Taking Data\n")
            [spins, E(p, temp, numSpins), n_HS(p, temp, numSpins)] = ...
                equilibrateSpins_H(...
                dataPts, spins, k(temp), T(temp), mu, H, J, ...
                big_delta, ln_g, G, listLS, ...
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
                imagesc(spins)
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

legArr = makeLegend(L);
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
    img_nom = strcat(trial_dir,'\',dat_str0,'_',...
        'nHSvsT','_J',J_nom,'K_D',bD_nom,'K_pLS',num2str(pLS),'_pHS',num2str(pHS),'.txt');
    figure
    plot(T_inv, mean_nHS,'b*-')
    hold on
    title(plt_title, 'Interpreter', 'tex')
    xlabel("Temperature T (K)")
    ylabel("n_H_S", 'Interpreter','tex')
    axis([-inf inf 0 1.01])
    legend(legArr,'Location','southeast')
    axis([-inf inf 0 1.0])
    hold off
    
    saveas(gcf, ...
        strcat(trial_dir,'\',dat_str0,'_','nHSvsT','_J',J_nom,'K_D',bD_nom,'K.png'))
    
    writematrix([T_inv mean_nHS],img_nom)
else
    n_HS = squeeze(n_HS);
    nom = strcat(trial_dir,'\',dat_str0,'_',...
        'nHSvsT','_J',J_nom,'K_D',bD_nom,'K_pLS',num2str(pLS),...
        '_pHS',num2str(pHS),'_L',num2str(L));
    %figure
    %plot(T, E)
    %title("E")
    plt_title = strcat('\rm ','J=',J_nom, 'K and ',' \Delta=',bD_nom,'K');
    figure
    plot(T_inv, n_HS,'.-')
    hold on
    grid on
    title(plt_title, 'interpreter','tex')
    xlabel("Temperature T (K)")
    ylabel("n_H_S",'Interpreter','tex')
    axis([-inf inf 0 1.01])
    legend(legArr,'Location','southeast')
    hold off
    saveas(gcf,...
        strcat(nom,'.png'))
    writematrix([T_inv' n_HS'],...
        strcat(nom,'.txt') )
    
end



function legArr = makeLegend(L)
%returns a legend given an array of lattice sizes
legArr = cell(max(size(L)),1);
for idx = 1:max(size(L))
    legArr{idx} = strcat(num2str(L(idx)),'x',num2str(L(idx)),' spins');
end
end
