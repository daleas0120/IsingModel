% PEEM Parameters

k_b = 8.617333262*10^-5;%eV/K 
mu = 1; %atomic magnetic moment
T_c = 4.515;

%% SIMULATION PARAMETERS

evo = 50e0; %number of MC steps to let the system burn in; this is discarded
dataPts = 50e0; %number of MC steps to evaluate the system
numTrials = 1; %number of times to repeat the experiment
frameRate = 1001; % provides a modulus to save snapshot of system
saveIntResults = false;% save intermediate results:

%% LATTICE PARAMETERS
boundCond = (0); %boundary condition
%L = [4, 7, 10, 40];
L = [1002];
D = 7;

%% MOLECULE PARAMETERS
bd = 1627.5;
omega = 1;

%% Way UP (LS to HS)
T_up = true;
J_inc =49;%
%T1 = 2200;
%T1 = [100:2:200];%K
T_inc = 165;
big_delta1 = bd;%K
%ln_g1 = 44.7/8.31; %ratio of degeneracy HS to LS
%ln_g1 = 67.5/8.31;
%ln_g1 = 48.2/8.31;
ln_g1 = 81.9/8.31;

G1 = 0;%K
H1 = 0; %external magnetic field

pLS1 = 0; %percentage of interior spins locked in LS
pHS1 = 0; %percentage of interior spins locked in HS
boundCond1 = (0); %boundary condition

%% WAY DOWN (HS to LS)
T_down = false;
%{
J2 = 20;%
%T2 = 0;
T2 = [200:-2:100];%K
big_delta2 = bd;%K
%ln_g2 = 47.4/8.31;
%ln_g2 = 49.8/8.31; %ratio of degeneracy HS to LS
%ln_g2 = 5.8002;
ln_g2 = ln_g1;
G2 = 0;%K
H2 = 0; %external magnetic field

pLS2 = 0; %percentage of interior spins locked in LS
pHS2 = 0; %percentage of interior spins locked in HS
boundCond2 = (0); %boundary condition
%}
