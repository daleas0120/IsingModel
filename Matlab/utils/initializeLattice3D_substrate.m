function [spinD, listLS] = initializeLattice3D_substrate(...
    N,D,b,pLS,pHS,substrate)
%{
%initializeLattice3D.m
%Ashley Dale
%creates a randomly initialized NxN lattice

N: number of spins along an edge
D: number of stacked layers
b: [1, 0 -1]; lock edge spins into HS, open boundary, or LS
pLS: probability spin is locked in LS state
pHS: probability spin is locked in HS state
pin: locks the first (pin) layers in some orientation
%}

spinD = b.*ones(N); %bottom layer
listLS(1, :) = [0 0 0];

for kdx = 1:2
    spinD = cat(3, spinD, substrate);
    [cx, cy, cz] = ndgrid(2:N-1, 2:N-1, kdx+1);
    X = cx(:);
    Y = cy(:);
    Z = cz(:);
    
    pts = [X Y Z];
    
    listLS = [listLS; pts];
       
end



%creates additional pinned layers

for kdx = 3:(D - 2)
    
    spins = rand(N-2); %decide how many ones there are
    
    for idx = 1:(N-2)
        for jdx = 1:(N-2)
            
            if (spins(idx, jdx) > 0.5)
                
                spins(idx, jdx) = 1;
                
                if rand < 2*pHS %&& (kdx == 1)
                    listLS(pt, :) = [(idx+1) (jdx+1) (kdx+1)];
                    pt = pt + 1;
                end
                
            else
                spins(idx, jdx) = (-1);
                
                if rand < 2*pLS %&&(kdx == 1)
                    listLS(pt, :) = [(idx+1) (jdx+1) (kdx+1)];
                    pt = pt + 1;
                end
            end
        end
    end
    
    spins = padarray(spins,[1 1],b,'both');
    spinD = cat(3,spinD, spins);
end

spinD = cat(3, spinD, b.*ones(N));

end

