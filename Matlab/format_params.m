
bD_nom = num2str(big_delta_K);
J_nom = num2str(J_K);

J_ev = J_K*k_b; %coupling constant/exchange energy in eV
T_ev = T_K.*k_b;
bD_ev = big_delta_K*k_b;
k = J_ev./(k_b.*T_K); % dimensionless inverse temperature

big_delta = (k_b*big_delta_K)/J_ev;
T = (k_b.*T_K)./J_ev;
J = J_ev/J_ev;
G = (k_b*G)/J_ev;

T_inv = (J_ev.*T_K)./k_b;

%Energy output variables
E = zeros(1, length(T_K));
Snn = zeros(1, length(T_K));

%Magnetism output variables
B = zeros(1, length(T_K));

%Spin fraction output variables
nHS = zeros(length(T_K), dataPts);
nHS_evo = zeros(length(T_K), evo);