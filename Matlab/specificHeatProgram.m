%%Plot heat capacity from energy fluctuations C~ sE^2 K^2 where sE is 
%mean-square fluctuations in Energy

%% Get energy fluctuation data
%[file, path] = uigetfile("*.txt");
%file_path = strcat(path, file);

%E = readmatrix(file_path);
%k = [0.001 0.3 0.44 0.5 1 1000];
%k = [0.001 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.44 0.5 0.55 ...
%    0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 1];
k_c = 0.44;
k = 0.36:0.02:0.80;
k = [0.22 0.29 k 1 2]./k_c;

dat_str = '200325a_';
fig_name = strcat(dat_str, "Heat_Cap_vs_k.png");
%% find sE

E_sqrd = mean(E.*E);

E_avg = (mean(E)).^2;

sE = sqrt(E_sqrd - E_avg);

figure;

plot(k, sE, '-o')
hold on
title("Heat Capacity vs K/K_c")
xlabel("K")
ylabel("Heat Capacity C")
saveas(gcf, fig_name)