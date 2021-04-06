% Monte Carlo Ising Model
% Ashley Dale
% Calls the following matlab files: initializeLattice.m,
% equilibrateSpins_H.m

%%
tic
clear;

PEEMparam()

%% CONVERSION TO DIMENSIONLESS UNITS (Formatting for program)

k_b = 8.617333262*10^-5;%eV/K
mu = 1; %atomic magnetic moment

if T_up == true
    bD_nom1 = num2str(big_delta1);
    J_nom1 = num2str(J_inc);
    
    J_ev1 = J_inc*k_b; %coupling constant/exchange energy in eV
    T_ev1 = T_inc.*k_b;
    bD_ev1 = big_delta1*k_b;
    G_ev1 = G1*k_b;
    
    k1 = J_ev1./(k_b.*T_inc); % dimensionless inverse temperature
    
    big_delta1 = (k_b*big_delta1)/abs(J_ev1);
    T1 = (k_b.*T_inc)./abs(J_ev1);
    J1 = J_ev1/abs(J_ev1);
    G1 = G_ev1/abs(J_ev1);
    T_inv1 = (abs(J_ev1).*T1)./k_b;
    
    %Magnetism output variables
    H1 = zeros(1, length(k1));
    
    %Spin fraction output variables
    n_HS1 = zeros(1, length(k1));
    nHS_by_layer1 = zeros(D, length(k1));
end

if T_down == true
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
    
    H2 = zeros(1, length(k2));
    n_HS2 = zeros(1, length(k2));
    nHS_by_layer2 = zeros(D, length(k2));
    
end

% naming system for the files and folders holding data from repeated trials
p_name = {'a_', 'b_', 'c_', 'd_', 'e_', 'f_', 'g_', 'h_', 'i_', 'j_', 'k_',...
    'l_', 'm_', 'n_', 'o_', 'p_', 'q_', 'r_', 's_', 't_', 'u_', 'v_', 'w_', ...
    'x_', 'y_', 'z_', 'A_', 'B_', 'C_', 'D_'};
%Energy output variables
E = zeros(1, length(k1));
Snn = zeros(1, length(k1));

%%
for p = 1:numTrials
    
    %%
    if numTrials>1 || ~saveIntResults
        % save all trials in a single directory at highest level
        t = datetime('now');
        t.Format = "yyMMdd";
        tryName = num2str(numTrials);
        dat_str0 = string(t);
        trial_dir = strcat('..\..\',dat_str0,'_',tryName,'_3DtrialRuns');
        mkdir(trial_dir)
    else
        % no group directory required
        trial_dir = '..\..';
    end
    
    %%
    for numSpins = 1:length(L)% square root of number of spins
        
        N = L(numSpins);
        
        %results folder for this particular data run
        if saveIntResults
            t = datetime('now');
            t.Format = "yyMMdd";
            dat_str = string(t);
            dir_name = strcat(trial_dir,'\',dat_str,p_name{p},'_',num2str(N),'spins3D');
            mkdir(dir_name)
            mkdir(dir_name,'frames')
        else
            dir_name = "";
        end
        
        %initialize 3D lattice
        
        %[spins, listLS] = initializeLattice3D_ones(...
        %    N, D, boundCond, pLS1, pHS1);
        
        [spins, listLS] = initializeLattice3D(...
            N, D, boundCond, pLS1, pHS1); %randomly initializes 3D lattice
        
        %[spins, listLS] = initializeLattice3D_pin(...
        %    N, D, boundCond, pLS1, pHS1, 2);
        
        origSpins = spins;
        
        % View initial lattice
        %%{
        figure
        spinVis(spins)
        axis equal
        pause(1)
        close
        %}
        figure;
        %%
        if T_up == true
            for temp = 1:length(k1)
                %copy spins for later comparison
                spins_last = spins;
                
                %let state reach equilibrium
                X = sprintf('Cooling %d x %d x %d spins to temp %f ....',...
                    N, N, D, T_inv1(temp));
                disp(X)
                
                [spins, ~, ~] = equilibrateSpins_3D(...
                    evo, spins, k1(temp), T1(temp), omega, H1, J1,...
                    big_delta1, ln_g1, listLS,...
                    frameRate, dir_name, saveIntResults);
                
                %take data
                fprintf("Taking Data\n")
                [spins, E(p, temp, numSpins), n_HS1(p, temp, numSpins)] = ...
                    equilibrateSpins_3D(...
                    dataPts, spins, k1(temp), T1(temp), omega, H1, J1, ...
                    big_delta1, ln_g1, listLS, ...
                    frameRate, dir_name, saveIntResults);
                
                nHS_by_layer1(:,temp) = nHS_layer(spins);
                H1(1, temp) = magnetism(spins(2:(L-1), 2:(L-1), 2:(D-1)));
                
                %close;
                %spinVis(spins)
                %pause(1)
                %axis equal
                toc
            end
        end
        %%
        if T_down == true
            %{
        for temp = 1:length(k2)
            %copy spins for later comparison
            spins_last = spins;
            
            %let state reach equilibrium
            X = sprintf('Cooling %d x %d x %d spins to temp %f ....',...
                N, N, D, T_inv2(temp));
            disp(X)
            
            [spins, ~, ~] = equilibrateSpins_3D(...
                evo, spins, k2(temp), T2(temp), mu, H2, J2, big_delta2, ln_g2, listLS, ...
                frameRate, dir_name, saveIntResults);
            
            %take data
            fprintf("Taking Data\n")
            [spins, E(p, temp, numSpins), n_HS2(p, temp, numSpins)] = ...
                equilibrateSpins_3D(...
                dataPts, spins, k2(temp), T2(temp), mu, H2, J2, ...
                big_delta2, ln_g2, listLS, ...
                frameRate, dir_name, saveIntResults);
            
            nHS_by_layer2(:, temp) = nHS_layer(spins);
            
            H2(1, temp) = magnetism(spins(2:(L-1), 2:(L-1), 2:(D-1)));
            
            %spinVis(spins)
            %axis equal
            %pause(1)
            
            toc
        end
            %}
        end
    end
end
%%

img = flattenSpins(spins, true, true, trial_dir);
binaryImg = binarizeSpins(img, 1, true, false, trial_dir);
nHS_img = n_HSfrac(img)
domains = bwboundaries(binaryImg);
detect_sizes = cellfun('size',domains,1);
[sz,index] = max(detect_sizes);
results = domains{index};
%results{i} = B;

%visual double check
figure;
imagesc(binaryImg(:,:,1));
axis square
colormap gray
colors=['b' 'g' 'r' 'c' 'm' 'y'
    ];
hold on
for k = 1:length(domains)
    boundary = domains{k};
    cidx = mod(k,length(colors))+1;
    plot(boundary(:,2),boundary(:,1),'g', 'LineWidth', 2)
end
hold off
