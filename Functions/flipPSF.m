function outPSF = flipPSF(inPSF)
[Sy, Sx, Sz] = size(inPSF);

outPSF = zeros(Sy,Sx,Sz);
for i = 1:Sx
    for j = 1:Sy
        for k = 1:Sz
            outPSF(j,i,k) = inPSF(Sy-j+1,Sx-i+1,Sz-k+1);
        end
    end
end

dx = abs (mod( Sx , 2 ) - 1);
dy = abs (mod( Sy , 2 ) - 1);
dz = abs (mod( Sz , 2 ) - 1);

outPSF = circshift(outPSF, [dy, dx, dz]);

end
