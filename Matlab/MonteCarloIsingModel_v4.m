% Monte Carlo Ising Model
% Ashley Dale
% Calls the following matlab files: initializeLattice.m,
% equilibrateSpins_H.m

%%
tic
clear;

N = 10; % square root of number of spins

k_up = 0.36:0.01:0.80;
k_down = 0.8:-0.01:0.36;

k = [0.22 0.29 k_up 1 2 4 4 2 1 k_down 0.29 0.22];
%k = 0.44;

k_b = 8.617333262*10^-5;
J = 1; %coupling constant/exchange energy
mu = 1; %atomic magnetic moment

T = 0.44./k;%[0.001 0.01 0.1 10 100 1000];
H = 1; %external magnetic field
h = 0; %magnetic field

%Beta = 1./(k_b.*T);

evolution = 1e4;
frameRate = 1e4+1;

%results folder
t = datetime('now');
t.Format = "yyMMddHHmmss";
dat_str = string(t);
dir_name = strcat('..\..\',dat_str,'_',num2str(N),'spins');
mkdir(dir_name)
mkdir(dir_name,'frames')


%Energy output variables
E = zeros(evolution, length(k));
E_img_name = strcat(dat_str, num2str(N),'spins_TotalEnergy');
Snn = zeros(1, length(k));

%Magnetism output variables
B = zeros(evolution, length(k));
B1_img_name = strcat(dat_str, num2str(N),'spins_TotalMagnetism_1.png');
%initialize 2D lattice
spins = initializeLattice(N);

for temp = 1:length(k)
    %create figure to view data
    figure;
    
    temp_name = num2str(k(temp));
    file_name = strcat(dir_name,'/',dat_str, num2str(N),'spins_k_', temp_name, '.txt');
    image_name = strcat(dir_name,'/',dat_str, num2str(N),'spins_k_', temp_name, '.png');
    m = 1;
    
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
    %figure;
    %imagesc(spins)
    %axis square;
    %saveas(gcf, image_name)
    
    %add temp to legend for later plots
    plt_legend{temp} = num2str(k(temp));
    
end

%%
%close all

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
