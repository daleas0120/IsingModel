function S = sumNNNN(spins, i, j, weights, periodic)
%{

%}

[M, N] = size(spins);

sNN_wt = weights(1);
sNNN_wt = weights(2);
sNNNN_wt = weights(3);
sNNNNN_wt = weights(4);

%/ A B C D E
%V - n a o -
%W p b c d q
%X e f g h i
%Y r j k l s
%Z - t m u -

% (i-2) --> V
% (i-1) --> W
% (i-0) --> X
% (i+1) --> Y
% (i+2) --> Z

% (j-2) --> A
% (j-1) --> B
% (j-0) --> C
% (j+1) --> D
% (j+2) --> E


if periodic
    
    V = mod(M+i-2, M);
    if V == 0
        V = M;
    end
    W = mod(M+i-1, M);
    if W == 0
        W = M;
    end
    X = i;
    Y = mod(M+i+1, M);
    if Y == 0
        Y = M;
    end
    Z = mod(M+i+2, M);
    if Z == 0
        Z = M;
    end
    
    A = mod(N+j-2, N);
    if A == 0
        A = N;
    end
    B = mod(N+j-1, N);
    if B  == 0
        B = N;
    end
    
    C = j;
    D = mod(N+j+1, N);
    if D == 0
        D = N;
    end
    E = mod(N+j+2, N);
    if E == 0
        E = N;
    end
    
    sNN = spins(X, D)+...
        spins(Y, C)+...
        spins(W, C)+...
        spins(X, B);
    
    sNNN = spins(Y, D) +...
        spins(W, D)+...
        spins(Y, B)+...
        spins(W, B);
    
    sNNNN = spins(V, C)+...
        spins(Z, C)+...
        spins(X, A)+...
        spins(X, E);
    
    sNNNNN = spins(V, B)+...
        spins(V, D)+...
        spins(W, A)+...
        spins(W, E)+...
        spins(Y, A)+...
        spins(Y, E)+...
        spins(Z, B)+...
        spins(Z, D);
    
else
    
    sNN = spins(i, j+1)+...
        spins(i+1, j)+...
        spins(i-1, j)+...
        spins(i, j-1);
    
    sNNN = spins(i+1, j+1) +...
        spins(i-1, j+1)+...
        spins(i+1, j-1)+...
        spins(i-1, j-1);
    
    
    if i>2 && i<M-2 && j>2 && j<N-2
        sNNNN = spins(i-2, j)+...
            spins(i+2, j)+...
            spins(i, j-2)+...
            spins(i, j+2);
        
        sNNNNN = spins(i-2, j-1)+...
            spins(i-1, j-2)+...
            spins(i+1, j-2)+...
            spins(i+2, j-1)+...
            spins(i-2, j+1)+...
            spins(i-1, j+2)+...
            spins(i+1, j+2)+...
            spins(i+2, j+1);
            
    else
        sNNNN = 0;
        sNNNNN = 0;
    end
end

S = sNN_wt*sNN + sNNN_wt*sNNN + sNNNN_wt*sNNNN + sNNNNN_wt*sNNNNN;

end