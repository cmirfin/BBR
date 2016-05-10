function [boundary,normals] = boundaryNormal(img,Delta)

%function takes an image, finds boundaries and computes normal to these
%boundaries.
if nargin < 2
    Delta = 2;
else
    Delta = ceil(Delta)+1;
end

%preprocessing
% h = fspecial('gaussian',[5,5]);
% img = imfilter(img,h);

%img = medfilt2(img,[5,5]);

% thresh = 25;
% T = zeros(size(img));
% T(img>thresh) = 1;
% T(img<=thresh) = 0;

%find boundary points
[boundaryImg] = edge(img,'sobel');


%reduce 
newBoundaryImg = zeros(size(img));
newBoundaryImg(Delta:end-Delta,Delta:end-Delta) = boundaryImg(Delta:end-Delta,Delta:end-Delta);
[X,Y] = find(newBoundaryImg == 1);
boundary = [X(:),Y(:)];

[dx,dy] = gradient(img);

l=sqrt(dx.^2+dy.^2);
dx = dx./l; 
dy = dy./l;

nx = zeros(1,length(X));
ny = zeros(1,length(Y));

for i = 1:length(X)
    
    nx(i) = dx(X(i),Y(i));
    ny(i) = dy(X(i),Y(i));
end
normals = abs([nx;ny])';


%plot vectors
Nx = zeros(size(img));
Ny = zeros(size(img));

for i = 1:length(X)

Nx(X(i),Y(i)) = nx(i);
Ny(X(i),Y(i)) = ny(i);
end

%quiver(Nx,Ny);
end

