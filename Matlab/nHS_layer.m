function nHS = nHS_layer(spins)
%{
nHS_layer.m
Ashley Dale
takes as an input the entire lattice
returns a vector nHS consisting of the high spin fraction for each layer 

calls function n_HSfrac, which takes as input 2D matrix, ignores the
outermost border, and returns a nHS

%}

[~, ~, p] = size(spins);

%m is x direction
%n is y direction
%p is z direction

nHS = zeros(p, 1);

for i = 1:p
    nHS(i) = n_HSfrac(spins(:,:,i));
end

