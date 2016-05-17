function [cost,phi] = boundaryCostNonRigid(phi,points,normals,fixedImage,Fx,Fy)

%parameters
M = 0.05;
DeltaIn = 2; %projection distance
DeltaOut = 2;

lambda = 0.1; %regularization
stepsize = 10;

[m,n] = size(fixedImage);
[X,Y] = meshgrid(1:n,1:m); %grid points of fixed image
maxIterations = 500;

for i = 1:maxIterations
    %transform boundary points
    projectedPoints = transformPointsNonRigid(points,phi,m,n,normals,DeltaIn,DeltaOut);
    N = size(projectedPoints,1); %number of transformed boundary points.

    x = [projectedPoints(:,1); projectedPoints(:,3)];
    y = [projectedPoints(:,2); projectedPoints(:,2)];

    %interpolate values
    fixedInterp = interp2(X,Y,fixedImage,y,x,'linear',1);
    valueIn = fixedInterp(1:N);
    valueOut = fixedInterp(N+1:end);

    FxInterp = interp2(X,Y,Fx,y,x,'linear',0); %x-gradient
    FyInterp = interp2(X,Y,Fy,y,x,'linear',0); %y-gradient

    Q = (M*(100*(valueIn-valueOut)./(0.5*(valueIn+valueOut))));

    J = 1 + tanh(Q);
    J = double(sum(J))/N;
    fprintf('Cost: %f \n',J);
    
    %gradients
    derivIn = [FxInterp(1:N),FyInterp(1:N)];
    derivOut = [FxInterp(N+1:end),FyInterp(N+1:end)];    
    
    grad = zeros(N,2);
    for j = 1:2
        grad(:,j) = (sech(Q).^2).*(valueOut.*derivIn(:,j) - valueIn.*derivOut(:,j))./((valueIn + valueOut).^2);
    end
    grad = double((400*M/N)*grad);

    h = -stepsize*grad;
    
    %update
    phi = phi + h; 
    cost(i) = J;
end

end

function rPrime = transformPointsNonRigid(r,phi,m,n,normals,DeltaIn,DeltaOut)
%output xIn, yIn, xOut, yOut
%NOTE: check if reversing of m,n is correct.
rIn = r - DeltaIn.*[normals(:,1),normals(:,2)];
rOut = r + DeltaOut.*[normals(:,1),normals(:,2)];

xuIn = rIn(:,1) + phi(:,1);
xuOut = rOut(:,1) + phi(:,1);
xuIn(xuIn<1) = 0;
xuIn(xuIn>m) = 0;
xuIn(isnan(xuIn))=0;

xuOut(xuOut<1) = 0;
xuOut(xuOut>m) = 0;
xuOut(isnan(xuOut))=0;

yvIn = rIn(:,2) + phi(:,2);
yvOut = rOut(:,2) + phi(:,2);
yvIn(yvIn<1) = 0;
yvIn(yvIn>n) = 0;
yvIn(isnan(yvIn))=0;

yvOut(yvOut<1) = 0;
yvOut(yvOut>n) = 0;
yvOut(isnan(yvOut))=0;

rPrime = [xuIn(:), yvIn(:), xuOut(:), yvOut(:)]; 

%remove points outside of image domain or NaN
[ind,~] = find(rPrime == 0);
rPrime(ind,:) = [];

end




