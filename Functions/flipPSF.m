function outPSF = flipPSF(inPSF)
[Sy, Sx, Sz] = size(inPSF);

outPSF = inPSF(end:-1:1, end:-1:1, end:-1:1);

dx = abs (mod( Sx , 2 ) - 1);
dy = abs (mod( Sy , 2 ) - 1);
dz = abs (mod( Sz , 2 ) - 1);

outPSF = circshift(outPSF, [dy, dx, dz]);

end
