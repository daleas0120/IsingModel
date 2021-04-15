%% SIMULATION PARAMETERS

evo = 100e0; %number of MC steps to let the system burn in; this is discarded
dataPts = 100e0; %number of MC steps to evaluate the system
numTrials = 1; %number of times to repeat the experiment
frameRate = 101; % provides a modulus to save snapshot of system
saveIntResults = true;% save intermediate results:

%% LATTICE PARAMETERS
boundCond = (0); %boundary condition
%L = [4, 7, 10, 40];
L = [37];
D = 37;

%% MOLECULE PARAMETERS
bd = 1325;

%% Way UP (LS to HS) (Temp increasing)
J1 = 10;%
%T1 = 2200;
T1 = [0:25:500];%K
big_delta1 = bd;%K
%ln_g1 = 44.7/8.31; %ratio of degeneracy HS to LS
ln_g1 = 67.5/8.31;
G1 = 0;%K
H1 = 0; %external magnetic field

pLS1 = 0; %percentage of interior spins locked in LS
pHS1 = 0; %percentage of interior spins locked in HS

%% WAY DOWN (HS to LS) (Temp Decreasing)
J2 = 10;%
%T2 = 0;
T2 = [500:-25:0];%K
big_delta2 = bd;%K
%ln_g2 = 47.4/8.31;
ln_g2 = 75/8.31; %MOST RECENT VALUE; REVERT TO THIS
%ln_g2 = ln_g1;
G2 = 0;%K
H2 = 0; %external magnetic field

pLS2 = 0; %percentage of interior spins locked in LS
pHS2 = 0; %percentage of interior spins locked in HS

%% FIELD PARAMETERS
omega1 = J1/10;
omega2 = J2/10;

