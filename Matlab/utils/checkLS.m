function checkLS(spins, listLS)
for idx = 1:max(size(listLS))
    if spins(listLS(idx)) ~= (-1)
        disp("ERROR\n")
    end
end

end