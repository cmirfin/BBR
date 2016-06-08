function outputImg = transformInterpolation2d(img,U,V)
img = double(img);
[m,n] = size(img);

[X,Y] = meshgrid(1:n,1:m);

Xp = X + U;
Yp = Y + V;

outputImg = griddata(Yp,Xp,img,Y,X);

end
