%ISING RENORMALIZATION STEP
clear
tic

snn_out_file = "snn_200311d_100spins_k_0pt44.txt";
snnn_out_file = "snnn_200311d_100spins_k_0pt44.txt";
img_out = "_renorm_200311d_100spins_k_0pt44";

cte = 1;

[file, path] = uigetfile("*.txt");
file_path = strcat(path, file);
spins(:,:,cte) = readmatrix(file_path);
file_list{cte, 1}= file_path;
cte = cte + 1;

while file ~= 0
    %import text file with spins
    [file, path] = uigetfile("*.txt");
    if file ~= 0
        file_path = strcat(path, file);
        %get file
        spins(:,:,cte) = readmatrix(file_path);
        file_list{cte, 1}= file_path;
        cte = cte + 1;
    end
end

[m, n, q] = size(spins);
maxRenorm = floor(log2(m));
S_nn = zeros(maxRenorm,q);
S_nnn = zeros(maxRenorm,q);

n = 0;

%snnn_fig = figure;
snn_fig = figure("InnerPosition", [10, 10, 1250, 400]);
for trial = 1:q
    name = strcat("set", num2str(trial),"/","set", num2str(trial), img_out,".png");
    n = 1;
    p = 1;
    %this gives us the zeroth order statistic
    S_nn(p, trial) = nearestN(spins(:,:,trial));
    S_nnn(p, trial) = nextNearestN(spins(:,:,trial));
    
    %plot(rows, columns, nth plot)
    
    subplot(1, maxRenorm, n)
    imagesc(spins(:,:, trial))
    axis square;
    
    % Renormalization
    S = spins(:,:, trial);
    
    for p = 2:maxRenorm %for each renormalization step
        
        n = n + 1

        spin_renorm = renormalize(S);
        
        S_nn(p, trial) = nearestN(spin_renorm);
        S_nnn(p, trial) = nextNearestN(spin_renorm);
        
        %plot(rows, columns, nth plot)
        subplot(1, maxRenorm, n)
        imagesc(spin_renorm)
        axis square
        
        S = spin_renorm;
    end
    saveas(gcf, name, 'png');
end

writematrix(S_nn, snn_out_file);
writematrix(S_nnn, snnn_out_file);

toc

function spin_renorm = renormalize(S)
len = length(S) - 2;
spin_renorm = zeros(floor(len/2 - 1));
idx = 1;
for i = 3:2:len
    jdx = 1;
    for j=3:2:len
        spin_renorm(idx, jdx) = ...
            sign(sum(sum(S(i:i+2,j:j+2)))+ 0.5);
        jdx = jdx + 1;
    end
    idx = idx + 1;
end
spin_renorm = padarray(spin_renorm,[2 2],0,'both');
end

function Snn = nearestN(spins)
Snn = 0;
for i = 3:length(spins) - 2
    for j = 3:length(spins) - 2
        Snn = Snn + spins(i,j)*(spins(i-1,j) + spins(i+1,j) +...
            spins(i,j-1) + spins(i,j+1));
    end
end
end

function Snnn = nextNearestN(spins)
Snnn = 0;
for i = 3:length(spins) - 2
    for j = 3:length(spins) - 2
        Snnn = Snnn + spins(i,j)*(spins(i-1,j-1) + spins(i+1,j-1) +...
            spins(i-1,j+1) + spins(i+1,j+1));
    end
end
end

