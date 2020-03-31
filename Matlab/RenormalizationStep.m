tic

%clear;

N = 200; % square root of number of spins
T = [0.1];%[0.001 0.01 0.1 10 100 1000];
%J = 1; %coupling constant
%k_b = 1; %8.617333262*10^-5;
K = 5; %number of renormalization steps
h = 0; %magnetic field
%Beta = 1./(k_b.*T);
k = [0.44 0.5 1 1000];
numTrials = 25;

y = zeros(K-1,length(T));
%S_nn = zeros(K,numTrials);
%S_nnn = zeros(K,numTrials);


for idx = 1:K-1
    
    a = mean(S_nn(idx,:))*mean(S_nn(idx+1,:)) - mean(S_nn(idx,:).*S_nn(idx,:));
    b = mean(S_nnn(idx,:))*mean(S_nn(idx+1,:)) - mean(S_nnn(idx,:).*S_nn(idx,:));
    c = mean(S_nn(idx,:))*mean(S_nnn(idx+1,:)) - mean(S_nn(idx,:).*S_nnn(idx,:));
    d = mean(S_nnn(idx,:))*mean(S_nnn(idx+1,:)) - mean(S_nnn(idx,:).*S_nnn(idx,:));
    
    alpha = mean(S_nn(idx,:))^2 - mean(S_nn(idx,:).^2);
    beta = mean(S_nn(idx,:))*mean(S_nnn(idx,:)) - mean(S_nn(idx,:).*S_nnn(idx,:));
    gamma = beta;
    delta = mean(S_nnn(idx,:))^2 - mean(S_nnn(idx,:).^2);
    
    %solve system of four equations
    big_A = [alpha 0 beta 0;...
        0 alpha 0 beta;...
        gamma 0 delta 0;...
        0 gamma 0 delta];
    
    little_b = [a b c d]';
    
    X = linsolve(big_A,little_b);
    
    A = X(1);
    B = X(2);
    C = X(3);
    D = X(4);
    
    matrix = [A B; C D];
    
    e_0 = eig(matrix);
    
    if e_0(1) > 0
        lambda = e_0(1);
    elseif e_0(2) > 0
        lambda = e_0(2);
    else
        fprintf('ERROR ')
        lambda = 1;
    end
    
    y(idx) = log(lambda)/log(4);
    
    fprintf('T = %f     lambda_1 = %f   lambda_2 =  %f \n', e_0(1), e_0(2))
end

for tmp = 1:length(k)
    leg(tmp) = {num2str(k(tmp))};
end
figure;
plot(y)
hold on
xlabel('Kth renormalization step')
ylabel('y value')
legend(leg)

toc
