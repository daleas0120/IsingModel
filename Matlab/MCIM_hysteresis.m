% Monte Carlo Ising Model - Magnetism
tic
clear;

N = 25; % square root of number of spins

%results folder
dat_str = '200326e_';
dir_name = strcat('hysteresis/',dat_str,num2str(N),'spins');
mkdir(dir_name)
mkdir(dir_name,'frames')

k = 0.5;
%k = 0.36:0.02:0.80;
%k = [0.22 0.29 k 1 2];

%k_b = 1; %8.617333262*10^-5;
J = 1; %coupling constant/exchange energy
mu = 1; %atomic magnetic moment

T = 0.44./k;%[0.001 0.01 0.1 10 100 1000];
H = 1; %external magnetic field

h1 = 0:0.0025:0.15;
h2 = 0.15:-0.0025:-0.15;
h3 = -0.15:0.0025:0;
h = [h3 h1 h2]; %magnetic field

%Beta = 1./(k_b.*T);

evolution = 5000;
frameRate = 1000;

%Energy output variables
E = zeros(evolution, length(k));
E_img_name = strcat(dat_str, num2str(N),'spins_TotalEnergy');
Snn = zeros(1, length(k));

%Magnetism output variables
B = zeros(evolution, length(k));
B1_img_name = strcat(dat_str, num2str(N),'spins_TotalMagnetism_1.png');
%initialize 2D lattice
spins = initializeLattice(N);

for mag = 1:length(h)
    figure;
    
    mag_name = num2str(h(mag));
    file_name = strcat(dir_name,'/',dat_str, num2str(N),'spins_h_', mag_name, '.txt');
    image_name = strcat(dir_name,'/',dat_str, num2str(N),'spins_h_', mag_name, '.png');
    m = 1;
    
    %initialize 2D lattice
    %spins = initializeLattice(N);
    
    spins_last = spins;
    
    %let state reach equilibrium
    [spins, E(:, mag), B(:, mag)] = equilibrateSpins(...
        evolution, N, spins, k, mu, h(mag), H, J, frameRate, spins_last, dir_name);
    close all;
    
    writematrix(spins,file_name);
    toc
    
    figure;
    imagesc(spins)
    axis square;
    saveas(gcf, image_name)
    
    plt_legend{mag} = num2str(h(mag));
    
end
%%
B_val = (mean(B)./(N*N));
E_val = (0.44.*mean(E))./(N*N);

close all

figure
plot(h, B_val,"*-")
hold on
title("Net Magnetism vs Applied H Field")
xlabel("Magnetism H/H_c")
ylabel("Net Magnetism")
hold off
saveas(gcf, strcat(dir_name,'/','netMagvsT','_','.png'))

B_dat = [h' B_val'];
writematrix(B_dat, strcat(dir_name,'/','BvsH.txt'))

figure
plot(h, E_val, "*-")
hold on
title("Net Energy vs Applied H Field")
xlabel("Temperature H/H_c")
ylabel("Net Energy")
hold off
saveas(gcf, strcat(dir_name,'/','netEvsT','_','.png'))

E_dat = [h' E_val'];
writematrix(E_dat, strcat(dir_name,'/','EvsH.txt'))
