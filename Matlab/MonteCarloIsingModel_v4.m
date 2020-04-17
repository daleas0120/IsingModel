% Monte Carlo Ising Model
% Ashley Dale
% Calls the following matlab files: initializeLattice.m,
% equilibrateSpins_H.m

%%
tic
clear;

N = 10; % square root of number of spins

k = 0.36:0.01:0.80;
k = [0.22 0.29 k 1 2 4];
%k = 0.44;

k_b = 8.617333262*10^-5;
J = 1; %coupling constant/exchange energy
mu = 1; %atomic magnetic moment

T = 0.44./k;%[0.001 0.01 0.1 10 100 1000];
H = 1; %external magnetic field
h = 0; %magnetic field

%Beta = 1./(k_b.*T);

evolution = 1e6;
frameRate = 1e6+1;



%Energy output variables
E = zeros(evolution, length(k));
E_img_name = strcat(dat_str, num2str(N),'spins_TotalEnergy');
Snn = zeros(1, length(k));

%Magnetism output variables
B = zeros(evolution, length(k));
B1_img_name = strcat(dat_str, num2str(N),'spins_TotalMagnetism_1.png');

p_name = {'a_', 'b_', 'c_', 'd_', 'e_', 'f_', 'g_', 'h_', 'i_', 'j_', 'k_',...
    'l_', 'm_', 'n_', 'o_', 'p_', 'q_', 'r_', 's_', 't_', 'u_', 'v_', 'w_', ...
    'x_', 'y_', 'z_', 'A_', 'B_', 'C_', 'D_'};

%results folder
t = datetime('now');
t.Format = "yyMMddHHmmss";
dat_str = string(t);
dir_name = strcat('..\..\',dat_str,'_',num2str(N),'spins');
mkdir(dir_name)
mkdir(dir_name,'frames')


for temp = 1:length(k)
    %create figure to view data
    figure;
    
    temp_name = num2str(k(temp));
    file_name = strcat(dir_name,'/',dat_str, num2str(N),'spins_k_', temp_name, '.txt');
    image_name = strcat(dir_name,'/',dat_str, num2str(N),'spins_k_', temp_name, '.png');
    m = 1;
    
    %initialize 2D lattice
    spins = initializeLattice(N);
    
    %copy spins for later comparison
    spins_last = spins;
    
    %let state reach equilibrium
    [spins, E(:, temp), B(:, temp)] = equilibrateSpins_H(...
        evolution, N, spins, k(temp), mu, h, H, J, frameRate, spins_last, dir_name);
    close all;
    
    %save spin matrix to text file
    writematrix(spins,file_name);
    
    toc
    
    %show new spins matrix
    figure;
    imagesc(spins)
    axis square;
    saveas(gcf, image_name)
    
    %add temp to legend for later plots
    plt_legend{temp} = num2str(k(temp));
    
end

%%
close all

figure
plot(T, abs(mean(B)./(N*N)),"*-")
hold on
title("Net Magnetism vs Temperature")
xlabel("Temperature T/T_c")
ylabel("Net Magnetism")
hold off
saveas(gcf, strcat(dir_name,'/','netMagvsT','_','.png'))

figure
plot(T, (0.44.*mean(E))./(N*N), "*-")
hold on
title("Net Energy vs Temperature")
xlabel("Temperature T/T_c")
ylabel("Net Energy")
hold off
saveas(gcf, strcat(dir_name,'/','netEvsT','_','.png'))
