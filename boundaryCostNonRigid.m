function [u1,v1] = boundaryCostNonRigid(points,normals,fixedImage,Fx,Fy,endpoints)



%parameters
M = 0.05;
DeltaIn = 2; %projection distance
DeltaOut = 2;

alpha = 0.3; %regularization

[m,n] = size(fixedImage);
[X,Y] = meshgrid(1:n,1:m); %grid points of fixed image
maxwarp = 200;

u1 = zeros(size(points,1),1);
v1 = u1;

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
    Ix = grad(:,1); Iy = grad(:,2);
    [u1,v1]=solveFlowMldivide(Ix,Iy,u1,v1,alpha,endpoints);
    
    plot(i,J,'.');
    hold on;
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
% [ind,~] = find(rPrime == 0);
% rPrime(ind,:) = [];
% phi(ind,:) = [];
% normals(ind,:) = [];
end

function  [D_phi, D2_phi] = transformDerivatives(phi)
    deriv=[-1,0,1]/2;
    D_phi = imfilter(phi,deriv','replicate'); %NOTE: transpose on deriv
    
    deriv2=[1,-2,1]/4;
    D2_phi = imfilter(phi,deriv2','replicate'); %NOTE: transpose on deriv2
    
    
end

function [u1,v1]=solveFlowMldivide(Ix,Iy,u0,v0,alpha,endpoints)
N = length(Ix);

%laplace operator
x1=ones(N,1).*alpha; x1(1)=0;
x2=ones(N,1).*alpha; x2(end)=0;
S=spdiags([x1(:),x2(:)],[1,-1],N,N);

C=sum(S,1);
x1(endpoints + 1) = 0;
x2(endpoints) = 0;
C(endpoints) = 0;
C(1) = 0;
C(end) = 0;

L=spdiags([x1(:),x2(:),-C(:)],[1,-1,0],N,N); 
%replicate L for u,v

A=[L,sparse(N,N);sparse(N,N),L];

%set-up linear system of equations

%b=A*double([u0(:);v0(:)])-double([Ix(:);Iy(:)]);
b=double([Ix(:);Iy(:)]); %only incremental

tic;
%solve Gauss-Newton update-step with Diffusion Regularisation
uv1=-A\b;

%[uv1,flagpcg]=pcg(A,b,[],20);
%[uv1,flagbig]=bicgstab(A,b,1E-2,20,[],[],double([u0(:);v0(:)]));
%[uv1,flag,res,nit]=sor(A',b,1.9,25,1E-4,double([u0(:);v0(:)]));
tOp=toc;
u1=reshape(uv1(1:N),N,1);
v1=reshape(uv1(N+1:end),N,1);


u1=u1+u0;
v1=v1+v0;
end

