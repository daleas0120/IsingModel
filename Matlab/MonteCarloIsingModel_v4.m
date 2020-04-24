% Monte Carlo Ising Model
% Ashley Dale
% Calls the following matlab files: initializeLattice.m,
% equilibrateSpins_H.m

%%
tic
clear;

k_b = 8.617333262*10^-5;%eV/K
J_ev = 160*k_b; %coupling constant/exchange energy in eV -- should also be unitless

mu = 1; %atomic magnetic moment

%T = J./(k_b.*k);

%T = [1:25:501 501:-25:1]; %degrees Kelvin
T = [1:5:101 101:-5:1];
k = J_ev./(k_b.*T); % dimensionless temperature

H = 0; %external magnetic field
g = 6; %
ln_g = log(g); %degeneracy
%big_delta = J_ev/(1300*k_b); %ERROR: difference in energy between HS and LS
%big_delta = 1300*k_b;
big_delta = 8.5;

beta = 1./(k_b.*T);
%J = big_delta/8.125;
J = big_delta/8;

evo = 5e3;
dataPts = 5e3;
frameRate = 5001;
numTrials = 3;
p_name = {'a_', 'b_', 'c_', 'd_', 'e_', 'f_', 'g_', 'h_', 'i_', 'j_', 'k_',...
    'l_', 'm_', 'n_', 'o_', 'p_', 'q_', 'r_', 's_', 't_', 'u_', 'v_', 'w_', ...
    'x_', 'y_', 'z_', 'A_', 'B_', 'C_', 'D_'};

%Energy output variables
E = zeros(1, length(k));

Snn = zeros(1, length(k));

%Magnetism output variables
B = zeros(1, length(k));

%L = [4, 7, 10, 40];
L = [4,7,10,40,100];

for p = 1:numTrials
    
    for numSpins = 1:length(L)% square root of number of spins
        
        N = L(numSpins);
        
        %results folder
        t = datetime('now');
        t.Format = "yyMMddHHmmss";
        dat_str = string(t);
        dir_name = strcat('..\..\',dat_str,'_',num2str(N),'spins');
        mkdir(dir_name)
        mkdir(dir_name,'frames')
        
        
        %initialize 2D lattice
        %spins = initializeLattice(N);
        
        spins = ones(N);
        spins(2:N-1, 2:N-1) = -1*spins(2:N-1, 2:N-1)
        meanS = aveS(N, spins);
        tmp = (1+meanS)/2;
        
        figure;
        imagesc(spins)
        axis square
        title("Initial Lattice")
        
        figure;
        for temp = 1:length(k)
            %create figure to view spins
            %figure;
            
            temp_name = num2str(T(temp));
            file_name = strcat(dir_name,'\',dat_str, p_name{p}, num2str(N),...
                'spins_k_', temp_name, 'K.txt')
            image_name = strcat(dir_name,'\',dat_str, p_name{p},num2str(N),...
                'spins_k_', temp_name, 'K.png');
            m = 1;
            
            %copy spins for later comparison
            spins_last = spins;
            
            %let state reach equilibrium
            fprintf("Cooling...\n")
            [spins, ~, ~] = equilibrateSpins_H(...
                evo, N, spins, k(temp), T(temp), mu, H, J, big_delta, ln_g, ...
                frameRate, spins_last, dir_name);
            %close all;
            
            %take data
            fprintf("Taking Data\n")
            [spins, E(p, temp, numSpins), B(p, temp, numSpins)] = equilibrateSpins_H(...
                dataPts, N, spins, k(temp), T(temp), mu, H, J, big_delta, ln_g, ...
                frameRate, spins_last, dir_name);
            close;
            
            %save spin matrix to text file
            %writematrix(spins,file_name);
            
            %show new spins matrix
            figure;
            imagesc(spins)
            axis square;
            %saveas(gcf, image_name)
            
            %add temp to legend for later plots
            %plt_legend{temp} = num2str(k(temp));
            
            meanS = aveS(N, spins);
            
            n_HS(p, temp, numSpins) = (1+meanS)/2;
            
        end
        toc
    end
end
%%
if numTrials > 1
    
    meanE = squeeze(mean(E));
    mean_nHS = squeeze(mean(n_HS))';
    
    %%
    close all
    
    figure
    plot(T, meanE', "*-")
    hold on
    title("Energy vs Temperature")
    xlabel("Temperature T (K)")
    ylabel("Energy")
    hold off
    saveas(gcf, strcat(dir_name,'\',dat_str,'_',num2str(N),'netEvsT','.png'))
    
    figure
    plot(T, mean_nHS,'*-')
    hold on
    title("Calculated thermal dependence of the HS fraction")
    xlabel("Temperature T (K)")
    ylabel("n_H_S")
    hold off
    saveas(gcf, strcat(dir_name,'\',dat_str,'_',num2str(N),'nHSvsT','.png'))
else
    
    figure
    plot(T, E)
    title("E")
    
    figure
    plot(T, n_HS)
    title("n_H_S")
end


%%
if length(L) > 1
    
    figure
    plot(T, meanE,"*-")
    hold on
    title('Averaged Energy vs Temperature')
    xlabel("Temperature T (K)")
    ylabel("Energy")
    legend('4 spins','7 spins', '10 spins', '40 spins')
    hold off
    saveas(gcf, strcat(dir_name,'\',dat_str,'_',num2str(N),'avgEvsT','.png'))
    
    figure
    plot(T, mean_nHS,"*-")
    hold on
    title('Averaged n_H_S vs Temperature')
    xlabel("Temperature T (K)")
    ylabel("n_H_S")
    legend('4 spins','7 spins', '10 spins', '40 spins')
    hold off
    saveas(gcf, strcat(dir_name,'\',dat_str,'_',num2str(N),'avgnHSvsT','.png'))
    
end
%%
function mean = aveS(N, spins)
sum_Si = sum(spins(2:N-1, 2:N-1), 'all');
mean = ((4*(N - 1)) + sum_Si)/(N*N);
end
