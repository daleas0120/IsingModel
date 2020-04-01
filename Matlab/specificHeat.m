%specificHeat.m
%Ashley Dale
%given a temperature vector k and an energy vector E(k) returns the
%specific heat, plots it, and saves the plot to the file

function specificHeat(k, E, dat_str)

fig_name = strcat(dat_str, "Heat_Cap_vs_k.png");

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
end