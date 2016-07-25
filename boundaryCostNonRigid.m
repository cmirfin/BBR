% Christopher J. Mirfin 
% Sir Peter Mansfield Imaging Centre, University of Nottingham
% christopher.mirfin@dtc.ox.ac.uk
% 23/07/2016

function [u1,v1] = boundaryCostNonRigid(points,normals,fixedImage,fixedPoints,DeltaIn,DeltaOut,alpha)
% BOUNDARYCOSTNONRIGID returns the displacement vectors based that
% transform the boundaries from one image to the fixedImage.
% This is non-rigid implementation of the BBR cost function proposed by
% Greve et al (2009).

%parameters
M = -0.5;
maxwarp = 500;

[m,n] = size(fixedImage);
[X,Y] = meshgrid(1:n,1:m); %grid points of fixed image


u1 = zeros(size(points,1),1);
v1 = u1;

% figure;
% plot(points(:,2),points(:,1),'.');
% [Fx, Fy] = gradient(fixedImage);

for i = 1:maxwarp

    %transform boundary points
    projectedPoints = transformPointsNonRigid(points,[u1,v1],m,n,normals,DeltaIn,DeltaOut);
    N = size(projectedPoints,1); %number of transformed boundary points.
    
        
    x = [projectedPoints(:,1); projectedPoints(:,3)];
    y = [projectedPoints(:,2); projectedPoints(:,2)];
    
    %interpolate values
    fixedInterp = interp2(X,Y,fixedImage,y,x,'linear',0);
    valueIn = fixedInterp(1:N);
    valueOut = fixedInterp(N+1:end);
    
    FxInterp = interp2(X,Y,Fx,y,x,'linear',0); %x-gradient
    FyInterp = interp2(X,Y,Fy,y,x,'linear',0); %y-gradient
    
    Q = (M*(100*(valueIn-valueOut)./(0.5*(valueIn+valueOut))));
    
    J = 1 + tanh(Q);
    J = double(sum(J))/N; % cost function
    
    % gradients
    derivIn = [FxInterp(1:N),FyInterp(1:N)];
    derivOut = [FxInterp(N+1:end),FyInterp(N+1:end)];
    
    grad = zeros(N,2);
    for j = 1:2
        grad(:,j) = (sech(Q).^2).*(valueOut.*derivIn(:,j) - valueIn.*derivOut(:,j))./((valueIn + valueOut).^2);
    end
    grad = double((400*M/N)*grad);
    Ix = grad(:,1); Iy = grad(:,2);
    [u1,v1]=solveFlow(Ix,Iy,u1,v1,alpha,fixedPoints);
    
    %regularization cost
%     [Du,Dv] = transformDerivatives(u1,v1,fixedPoints);
%     R = 0.5*alpha*sum(Du.^2 + Dv.^2);

    tpoints = round([u1,v1] + points);
    plot(tpoints(:,2),tpoints(:,1),'.');
    drawnow;
end


end

function rPrime = transformPointsNonRigid(r,phi,m,n,normals,DeltaIn,DeltaOut)
%output xIn, yIn, xOut, yOut
%NOTE: check if reversing of m,n is correct.
rIn = r - DeltaIn.*[normals(:,1),normals(:,2)];
rOut = r + DeltaOut.*[normals(:,1),normals(:,2)];

xuIn = rIn(:,1) + phi(:,1);
xuOut = rOut(:,1) + phi(:,1);
xuIn(xuIn<1) = rIn(xuIn<1,1);
xuIn(xuIn>m) = rIn(xuIn>m,1);
%xuIn(isnan(xuIn))=0;

xuOut(xuOut<1) = rOut(xuOut<1,1);
xuOut(xuOut>m) = rOut(xuOut>m,1);
%xuOut(isnan(xuOut))=0;

yvIn = rIn(:,2) + phi(:,2);
yvOut = rOut(:,2) + phi(:,2);
yvIn(yvIn<1) = rIn(yvIn<1,2);
yvIn(yvIn>n) = rIn(yvIn>n,2);
%yvIn(isnan(yvIn))=0;

yvOut(yvOut<1) = rOut(yvOut<1,2);
yvOut(yvOut>n) = rOut(yvOut>n,2);
%yvOut(isnan(yvOut))=0;

rPrime = [xuIn(:), yvIn(:), xuOut(:), yvOut(:)]; 

%remove points outside of image domain or NaN
% [ind,~] = find(rPrime == 0);
% rPrime(ind,:) = [];
% phi(ind,:) = [];
% normals(ind,:) = [];
end

function  [Du,Dv] = transformDerivatives(u,v,fixedPoints)
    %deriv=[-1,0,1]/2;
    deriv = [-1,1]/2;
    %D = imfilter([u,v],deriv','replicate'); %NOTE: transpose on deriv
    %Dv = imfilter(v,deriv','replicate'); 
    phi = [u,v];
    num = length(fixedPoints);
    startpoints = fixedPoints(1:num/2);
    endpoints = fixedPoints(num/2 + 1:end);
    for i = 1:num/2
         
        a = startpoints(i,1);
        b = endpoints(i,1);
        
        D(a:b,:) = (imfilter(phi(a:b,:),deriv',0));
          
    end
    Du = D(:,1);
    Dv = D(:,2);
end

function [u1,v1]=solveFlow(Ix,Iy,u0,v0,alpha,fixedPoints)
N = length(Ix);
num = length(fixedPoints);
s = fixedPoints(1:num/2);
e = fixedPoints(num/2 + 1:end);
%laplace operator
x1=ones(N,1).*alpha; x1(1)=0;
x2=ones(N,1).*alpha; x2(end)=0;
S=spdiags([x1(:),x2(:)],[1,-1],N,N);

C=sum(S,1);
%don't regularise between boundaries (e.g. between end/start points)
% x1(fixedPoints(2:end-1) + 1) = 0;
% x2(fixedPoints(2:end-1) - 1) = 0;

x1(s(2:end-1) + 1) = 0;
x2(s(2:end-1) - 1) = 0;
x1(e(2:end-1) + 1) = 0;
x2(e(2:end-1) - 1) = 0;

C(fixedPoints) = 1;

L=spdiags([x1(:),x2(:),-C(:)],[1,-1,0],N,N); 
%replicate L for u,v
A=[L,sparse(N,N);sparse(N,N),L];

b=double([Ix(:);Iy(:)]); %only incremental

%solve Gauss-Newton update-step with Diffusion Regularisation
uv1=A\b;

u1=reshape(uv1(1:N),N,1);
v1=reshape(uv1(N+1:end),N,1);

u1=u1+u0;
v1=v1+v0;
end



