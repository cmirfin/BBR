% Christopher J. Mirfin 
% Sir Peter Mansfield Imaging Centre, University of Nottingham
% christopher.mirfin@dtc.ox.ac.uk
% 23/07/2016

function [outputImg,u2,v2] = weightedInterpolation(img,uBoundary,vBoundary,boundaryPoints,radius)
%WEIGHTEDINTERPOLATION is a custom interpolation scheme that interpolates
%a dense image transform from a sparse set of boundary transforms. The
%image background is discarded from the calculation.

%The interpolated value at position (x,y) is taken to be the maximum value
%transform to the local voxel neighbourhood (radius parameter regulates
%this). 

% background extraction
[m,n] = size(img);
background = reshape(kmeans(img(:),3),m,n);
backgroundClass = background(1,1); %K-means classifcation of background
background(background==backgroundClass) = 0; %background
background(background>0) = 1; %not background

background = background(:);
backgroundIdx = find(background == 0);

% use elastic transform
[u2,v2] = elasticSolver(uBoundary,vBoundary,boundaryPoints,m,n);
uF = u2(:);
vF = v2(:);
uF(backgroundIdx) = [];
vF(backgroundIdx) = [];
        

[Y,X] = meshgrid(1:n,1:m);
X = X(:);
Y = Y(:);
X0 = X;
Y0 = Y;

%remove background pixel positions
X0(backgroundIdx) = [];
Y0(backgroundIdx) = [];

%transform tissue pixels
Xp = X0 + uF;
Yp = Y0 + vF;
%these are the coordinates of the transformed intensity from each pixel in
%tissue

outputImg = zeros(m,n);

%create "boxes" loop over grid positions
for i = 1:m*n
    
    x = X(i);
    y = Y(i);
    
    ind = find(Xp <= x+radius & Xp >= x-radius & Yp <= y+radius & Yp >= y-radius);

    %Case 1: ind = [] (empty) - assign this point to background
    if isempty(ind) == 1 %TRUE
        outputImg(x,y) = 0;
    else
        xpoints = X0(ind);
        ypoints = Y0(ind);
        for j = 1:numel(xpoints)
            intensity = img(xpoints(j),ypoints(j));
        end
        %intensity = intensity/numel(xpoints);
        intensity = max(intensity);
        outputImg(x,y) = intensity;
    end
    
end

end

