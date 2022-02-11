function S = reduceLattice(spins, r)

% r is the reduction factor; e.g. 6 spins collapse to 1
[M, N] = size(spins);

newM = floor((M)/2);
newN = floor(2*(N)/3);

%S = zeros(newM+2, newN+2);

i = 2;
j = 2;

p = 1;
q = 1;
S = zeros();

while (i < M)
    while (j < N)
        
        if r == 6 %2x3 reduction
            S(q, p) = spins(i, j) +...
                spins(i+1, j)+...
                spins(i+1, j+1);
            
            j = j+1;
            
            S(q, p) = S(q, p) + spins(i, j) +...
                spins(i, j+1)+...
                spins(i+1, j+1);
            
            p = p+1;
            j = j+2;
        elseif r == 9 %3x3 reduction
            %             S(q, p) = spins(i, j) +...
            %                 spins(i - 1, j-1)+...
            %                 spins(i-1, j)+...
            %                 spins(i-1, j+1)+...
            %                 spins(i, j-1)+...
            %                 spins(i, j+1)+...
            %                 spins(i+1, j-1)+...
            %                 spins(i+1, j)+...
            %                 spins(i+1, j+1);
            
            S(q, p) = spins(i, j) + sumNNN(spins, i, j, [1, 1, 0], true);
            
            j = j+3;
            p = p+1;
        end
    end
    j = 2;
    i = i+3;
    q = q+1;
    p = 1;
    
    
end

S = S./r;

end


