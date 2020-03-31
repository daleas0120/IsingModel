% Monte Carlo Ising Model
tic

clear;

N = 200; % square root of number of spins

%k = 0.36:0.02:0.80;
%k = [0.22 0.29 k 1 2];
k = 0.44; 

%k_b = 1; %8.617333262*10^-5;
J = 1; %coupling constant/exchange energy
mu = 1; %atomic magnetic moment

T = 0.44./k;%[0.001 0.01 0.1 10 100 1000];
H = 1; %external magnetic field
h = 0; %magnetic field

%Beta = 1./(k_b.*T);

evolution = 5000;
frameRate = 1;

%results folder
dat_str = '200328_movie_';
dir_name = strcat(dat_str,num2str(N),'spins');
mkdir(dir_name)
mkdir(dir_name,'frames')

%Energy output variables
E = zeros(evolution, length(k));
E_img_name = strcat(dat_str, num2str(N),'spins_TotalEnergy');
Snn = zeros(1, length(k));

%Magnetism output variables
B = zeros(evolution, length(k));
B1_img_name = strcat(dat_str, num2str(N),'spins_TotalMagnetism_1.png');


for temp = 1:length(k)
    figure;
    
    temp_name = num2str(k(temp));
    file_name = strcat(dir_name,'/',dat_str, num2str(N),'spins_k_', temp_name, '.txt');
    image_name = strcat(dir_name,'/',dat_str, num2str(N),'spins_k_', temp_name, '.png');
    m = 1;
    
    %initialize 2D lattice
    spins = initializeLattice(N);
    
    spins_last = spins;
    
    %let state reach equilibrium
    [spins, E(:, temp), B(:, temp)] = equilibrateSpins_H(...
        evolution, N, spins, k(temp), mu, h, H, J, frameRate, spins_last, dir_name);
    close all;
    
    writematrix(spins,file_name);
    toc
    
    figure;
    imagesc(spins)
    axis square;
    saveas(gcf, image_name)
    
    plt_legend{temp} = num2str(k(temp));
    
end

%%
% figure
% plot(E)
% hold on
% %legend(plt_legend)
% title("Total Energy vs Time")
% xlabel("Evolutions (aka Time)")
% ylabel("Energy: -J\Sigma\sigma_i \sigma_j")
% hold off
% saveas(gcf, strcat(dir_name,'/','totE','_','.png'))
% writematrix(E, strcat(dir_name,'/','totE.txt'))
%
% figure
% plot(B1)
% hold on
% %legend(plt_legend)
% title("Spontaneous Magnetism Over Time")
% xlabel("Evolutions (aka Time)")
% ylabel("Magnetism: - \muH\Sigma\sigma_i")
% hold off
% saveas(gcf, strcat(dir_name,'/','totMag','_','.png'))
% writematrix(B1, strcat(dir_name,'/','totB.txt'))
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
