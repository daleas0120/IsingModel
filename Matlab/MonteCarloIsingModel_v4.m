% Monte Carlo Ising Model
% Ashley Dale
% Calls the following matlab files: initializeLattice.m,
% equilibrateSpins_H.m

%%
tic
clear;

k = 0.35:0.05:2;
k = [k 3.5];
%k = 0.44;

k_b = 8.617333262*10^-6;%eV/K
J = 0.001; %coupling constant/exchange energy in eV
mu = 1; %atomic magnetic moment

T = J./(k_b.*k);
H = 0; %external magnetic field
ln_g = 6; %degeneracy
big_delta = 1300;
beta = 1./(k_b.*T);

evo = 1e5;
dataPts = 1.2e5;
frameRate = 1.2e5 + 1;
numTrials = 20;
p_name = {'a_', 'b_', 'c_', 'd_', 'e_', 'f_', 'g_', 'h_', 'i_', 'j_', 'k_',...
    'l_', 'm_', 'n_', 'o_', 'p_', 'q_', 'r_', 's_', 't_', 'u_', 'v_', 'w_', ...
    'x_', 'y_', 'z_', 'A_', 'B_', 'C_', 'D_'};

%Energy output variables
E = zeros(1, length(k));

Snn = zeros(1, length(k));

%Magnetism output variables
B = zeros(1, length(k));

L = [4 7 10 40];

for numSpins = 1:length(L)% square root of number of spins
    
    N = L(numSpins);
    
    %results folder
    dat_str = '200408';
    dir_name = strcat('..\..\' , dat_str,'_',num2str(N),'spins');
    mkdir(dir_name)
    mkdir(dir_name,'frames')
    
    for p = 1:numTrials
        %initialize 2D lattice
        spins = initializeLattice(N);
        figure;
        imagesc(spins)
        axis square
        title("Initial Lattice")
        
        for temp = 1:length(k)
            %create figure to view spins
            %figure;
            
            temp_name = num2str(k(temp));
            file_name = strcat(dir_name,'/',dat_str, p_name{p}, num2str(N),...
                'spins_k_', temp_name, '.txt')
            image_name = strcat(dir_name,'/',dat_str, p_name{p},num2str(N),...
                'spins_k_', temp_name, '.png');
            m = 1;
            
            %copy spins for later comparison
            spins_last = spins;
            
            %let state reach equilibrium
            [spins, ~, ~] = equilibrateSpins_H(...
                evo, N, spins, k(temp), mu, H, J, frameRate, spins_last, dir_name);
            close all;
            
            %take data
            [spins, E(p, temp), B(p, temp)] = equilibrateSpins_H(...
                dataPts, N, spins, k(temp), mu, H, J, frameRate, spins_last, dir_name);
            close all;
            
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
            
            n_HS(p, temp) = (1+meanS)/2;
            
        end
        writematrix([T' n_HS'], strcat(dir_name,'/',dat_str,p_name{p},num2str(N),...
            'nHSvsT','_','.txt'));
        writematrix([T' E'], strcat(dir_name,'/',dat_str,p_name{p},num2str(N),...
            'EvsT','_','.txt'));
        writematrix([T' B'], strcat(dir_name,'/',dat_str,p_name{p},num2str(N),...
            'BvsT','_','.txt'));
    end
    toc
    %%
    meanB(numSpins, :) = abs(mean(B));
    meanE(numSpins, :) = mean(E);
    mean_nHS(numSpins, :) = mean(n_HS);
    
    %%
    close all
    
    figure
    plot(T, (meanB(numSpins, :)),"*-")
    hold on
    title("Magnetism vs Temperature")
    xlabel("Temperature T (K)")
    ylabel("Net Magnetism")
    hold off
    saveas(gcf, strcat(dir_name,'/',dat_str,num2str(N),'netMagvsT','.png'))
    
    figure
    plot(T, meanE(numSpins, :), "*-")
    hold on
    title("Energy vs Temperature")
    xlabel("Temperature T (K)")
    ylabel("Energy")
    hold off
    saveas(gcf, strcat(dir_name,'/',dat_str,num2str(N),'netEvsT','.png'))
    
    figure
    plot(T, mean_nHS(numSpins, :),'*-')
    hold on
    title("Calculated thermal dependence of the HS fraction")
    xlabel("Temperature T (K)")
    ylabel("n_H_S")
    hold off
    saveas(gcf, strcat(dir_name,'/',dat_str,num2str(N),'nHSvsT','.png'))
end
%%
figure
plot(T, meanB,"*-")
hold on
title('Averaged Magnetism vs Temperature')
xlabel("Temperature T (K)")
ylabel("Magnetism")
legend('4 spins','7 spins', '10 spins', '40 spins')
hold off
saveas(gcf, strcat(dir_name,'/',dat_str,num2str(N),'avgMagvsT','.png'))

figure
plot(T, meanE,"*-")
hold on
title('Averaged Energy vs Temperature')
xlabel("Temperature T (K)")
ylabel("Energy")
legend('4 spins','7 spins', '10 spins', '40 spins')
hold off
saveas(gcf, strcat(dir_name,'/',dat_str,num2str(N),'avgEvsT','.png'))

figure
plot(T, mean_nHS,"*-")
hold on
title('Averaged n_H_S vs Temperature')
xlabel("Temperature T (K)")
ylabel("n_H_S")
legend('4 spins','7 spins', '10 spins', '40 spins')
hold off
saveas(gcf, strcat(dir_name,'/',dat_str,num2str(N),'avgnHSvsT','.png'))

%%
function mean = aveS(N, spins)
sum_Si = sum(spins(2:N-1, 2:N-1), 'all');
mean = ((4*(N - 1)) + sum_Si)/(N*N);
end