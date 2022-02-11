function saveSpinsImg(spins, fileName)

%remove boundary spins

S = spins(2:(end-1), 2:(end-1));

%convert to uint 255 scale

S = uint8(255.*S);

%write to file
imwrite(S, fileName)

end