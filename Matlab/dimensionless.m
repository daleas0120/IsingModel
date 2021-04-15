k_b = 8.617333262*10^-5; %eV/K
mu = 1; %atomic magnetic moment

bD_nom1 = num2str(big_delta1);
J_nom1 = num2str(J1);

J_ev1 = J1*k_b; %coupling constant/exchange energy in eV
T_ev1 = T1.*k_b;
bD_ev1 = big_delta1*k_b;
G_ev1 = G1*k_b;

k1 = J_ev1./(k_b.*T1); % dimensionless inverse temperature

omega = (omega1)/abs(J1); %transverse field strength; declared in units of (K)
big_delta1 = (k_b*big_delta1)/abs(J_ev1);
T1 = (k_b.*T1)./abs(J_ev1);
J1 = J_ev1/abs(J_ev1);
G1 = G_ev1/abs(J_ev1);
T_inv1 = (abs(J_ev1).*T1)./k_b;


bD_nom2 = num2str(big_delta2);
J_nom2 = num2str(J2);

J_ev2 = J2*k_b; %coupling constant/exchange energy in eV
T_ev2 = T2.*k_b;
bD_ev2 = big_delta2*k_b;
G_ev2 = G2*k_b;

k2 = J_ev2./(k_b.*T2); % dimensionless inverse temperature

big_delta2 = (k_b*big_delta2)/abs(J_ev2);
T2 = (k_b.*T2)./abs(J_ev2);
J2 = J_ev2/abs(J_ev2);
G2 = G_ev2/abs(J_ev2);
T_inv2 = (abs(J_ev2).*T2)./k_b;
