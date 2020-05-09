% Monte Carlo Ising Model
% Ashley Dale
% Calls the following matlab files: initializeLattice.m,
% equilibrateSpins_H.m

%%
tic
clear;

k_b = 8.617333262*10^-5;%eV/K
mu = 1; %atomic magnetic moment

%%{
%%PRB 84 Constants
J = 10;%K
T = [100:10:400];%K
big_delta = 1125;%K
ln_g = 6; %ratio of degeneracy HS to LS
G=0;
%}

%%2013 Chiruta Constants
%{
J = 200;%K
T = [50:5:360 360:-5:50];
big_delta = 2247;%K
ln_g = 7.216;
G = 0;
%}

bD_nom = num2str(big_delta);
J_nom = num2str(J);

J_ev = J*k_b; %coupling constant/exchange energy in eV
T_ev = T.*k_b;
bD_ev = big_delta*k_b;

k = J_ev./(k_b.*T); % dimensionless inverse temperature

H = 0; %external magnetic field

%% DIMENSIONLESS UNITS

big_delta = (k_b*big_delta)/J_ev;
T = (k_b.*T)./J_ev;
J = J_ev/J_ev;
G = (k_b*G)/J_ev;

T_inv = (J_ev.*T)./k_b;

%%
evo = 5e2; %number of MC steps to let the system burn in; this is discarded
dataPts = 1e2; %number of MC steps to evaluate the system
frameRate = 250e7+1; % provides a modulus to save snapshot of system
numTrials = 2; %number of times to repeat the experiment

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

%L = [4, 7, 10, 40];
L = [20];
D = 10;

%L = [5];
%D = 5;

for p = 1:numTrials
    
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
        
        %initialize 2D lattice
        [spins, listLS] = initializeLattice3D(N,D,(-1),(-1), 0.3); %randomly initializes 3D lattice
        
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
            X = sprintf('Cooling %d spins to temp %f ....',N, T_inv(temp));
            disp(X)
            [spins, ~, ~, ~] = equilibrateSpins_3D(...
                evo, spins, k(temp), T(temp), mu, H, J, big_delta, ln_g, listLS, ...
                frameRate, dir_name, saveIntResults);
            
            %take data
            fprintf("Taking Data\n")
            [spins, E(p, temp, numSpins), ~, n_HS(p, temp, numSpins)] = ...
                equilibrateSpins_3D(...
                dataPts, spins, k(temp), T(temp), mu, H, J, ...
                big_delta, ln_g, listLS, ...
                frameRate, dir_name, saveIntResults);
            
            close;
            
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
        end
        toc
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
    figure
        plt_title = strcat('\rm ','J=',J_nom, 'K and ',' \Delta=',bD_nom,'K');
    
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
else
    
    %figure
    %plot(T, E)
    %title("E")
    
    figure
    plot(T, n_HS)
    hold on
    title("n_H_S")
    xlabel("Temperature T (K)")
    ylabel("n_H_S")
    legend(legArr,'Location','southeast')
    hold off
    saveas(gcf, strcat(trial_dir,'\',dat_str0,'_',...
        'delt',bD_nom,'_J',J_nom,'_nHSvsT','.png'))
end

function legArr = makeLegend(L,D)
%returns a legend given an array of lattice sizes

for idx = 1:max(size(L))
    legArr{idx} = strcat(num2str(L(idx)),'x',num2str(L(idx)),'x',num2str(D),' spins');
end
end
