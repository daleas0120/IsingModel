%{
spinVis.m
Ashley Dale
Visualizes a spin lattice
%}
function spinVis(spinD)
[m, n, p] = size(spinD);
for i = 1:m
    for j = 1:n
        for k = 1:p
            if spinD(i,j,k) == 1
                argC = 'r';
                argS = '^';
            elseif spinD(i,j,k) == (-1)
                argC = 'k';
                argS = 'v';
            else
                argC = 'b';
                argS = '*';
            end
            plot3(i,j,k,argS,'Color', argC, 'MarkerFaceColor',argC)
            hold on
        end
    end
end
grid on
axis equal
hold off

end