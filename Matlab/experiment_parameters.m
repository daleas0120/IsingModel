%% Parameters for experiments:
%% 220228
%%{

bd = 3000;
k_b = 8.617333262*10^-5;%eV/K
weights = [1 0.7071 0 0];
L = [75]; % L <= (l mod 3 = 0)
D = 11; % D <= (d mod 3 = 0) + 2

substrateImg = 'D:\IsingModel\IsingModel\substrates\210427_1trialRunsflattenSpins.png';
substrateWeight = 1;

J_K = 47.5;%
T_K = 297;
big_delta_K = bd;%K
ln_g = 83.9/8.31;
G = 0;%K
H = 0; %external magnetic field

pLS = 0; %percentage of interior spins locked in LS
pHS = 0; %percentage of interior spins locked in HS
boundCond = (0); %boundary condition
omega = 1;

evo = 0; %number of MC steps to let the system burn in; this is discarded
dataPts = 1e3; %number of MC steps to evaluate the system
frameRate = 1; % provides a modulus to save snapshot of system
numTrials = 1; %number of times to repeat the experiment

% save intermediate results:
saveIntResults = false;

%image color
APSslideColor = [34/255 42/255 53/255];
%APSslideColor = [1 1 1];
%}
%% 170222
%{

bd = 3000;
k_b = 8.617333262*10^-5;%eV/K
weights = [1 0.7071 0 0];
L = [300]; % L <= (l mod 3 = 0)
D = 11; % D <= (d mod 3 = 0) + 2

substrateImg = 'D:\IsingModel\IsingModel\substrates\502spins_1.7437K_275img.png';
substrateWeight = 1;

J_K = 47.5;%
T_K = 297;
big_delta_K = bd;%K
ln_g = 83.9/8.31;
G = 0;%K
H = 0; %external magnetic field

pLS = 0; %percentage of interior spins locked in LS
pHS = 0; %percentage of interior spins locked in HS
boundCond = (0); %boundary condition
omega = 1;

evo = 0; %number of MC steps to let the system burn in; this is discarded
dataPts = 1e3; %number of MC steps to evaluate the system
frameRate = 1; % provides a modulus to save snapshot of system
numTrials = 1; %number of times to repeat the experiment

% save intermediate results:
saveIntResults = false;

%image color
APSslideColor = [34/255 42/255 53/255];
%APSslideColor = [1 1 1];
%}
%% February 15, 2022
% didn't finish because it was taking too long and the parameters were not
% correct
%{

bd = 3000.5;
k_b = 8.617333262*10^-5;%eV/K
weights = [1 0.7071 0.5 0];
L = [300]; % L <= (l mod 3 = 0)
D = 14; % D <= (d mod 3 = 0) + 2

substrateImg = 'D:\IsingModel\IsingModel\substrates\502spins_1.7437K_275img.png';
substrateWeight = 1;

J_K = 15;%
T_K = 297;
big_delta_K = bd;%K
ln_g = 83.9/8.31;
G = 0;%K
H = 0; %external magnetic field

pLS = 0; %percentage of interior spins locked in LS
pHS = 0; %percentage of interior spins locked in HS
boundCond = (0); %boundary condition
omega = 1;

evo = 0; %number of MC steps to let the system burn in; this is discarded
dataPts = 10e2; %number of MC steps to evaluate the system
frameRate = 1; % provides a modulus to save snapshot of system
numTrials = 1; %number of times to repeat the experiment

% save intermediate results:
saveIntResults = false;

%image color
APSslideColor = [34/255, 42/255, 53/255];
%APSslideColor = [1 1 1];
%}
%% February 14, 2022
%{

bd = 3000.5;
k_b = 8.617333262*10^-5;%eV/K
ln_g = 83.9/8.31;
J_K = 15;%
T_K = 297;
substrateImg = 'D:\IsingModel\IsingModel\substrates\502spins_1.7437K_275img.png';
substrateWeight = 1;
weights = [1 0.7071 0.5 0];
L = [30]; % L <= (l mod 3 = 0)
D = 14; % D <= (d mod 3 = 0) + 2


big_delta_K = bd;%K
G = 0;%K
H = 0; %external magnetic field
pLS = 0; %percentage of interior spins locked in LS
pHS = 0; %percentage of interior spins locked in HS
boundCond = (0); %boundary condition
omega = 1;

evo = 0; %number of MC steps to let the system burn in; this is discarded
dataPts = 10e2; %number of MC steps to evaluate the system
frameRate = 1; % provides a modulus to save snapshot of system
numTrials = 1; %number of times to repeat the experiment

% save intermediate results:
saveIntResults = false;

%image color
APSslideColor = [34/255, 42/255, 53/255];
%APSslideColor = [1 1 1];
%}
