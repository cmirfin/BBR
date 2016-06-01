function [u1,v1] = boundaryCostNonRigid(points,normals,fixedImage,Fx,Fy,endpoints)

%parameters
M = -0.5;
DeltaIn = 2; %projection distance
DeltaOut = 2;

alpha = 0.6; %regularization

[m,n] = size(fixedImage);
[X,Y] = meshgrid(1:n,1:m); %grid points of fixed image
maxwarp = 300;

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
    
    Q = M*(100*(valueIn-valueOut)./(0.5*(valueIn+valueOut)));
    
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
    [u1,v1]=solveFlow(Ix,Iy,u1,v1,alpha,endpoints);
    
    %regularization cost
    %[Du,Dv] = transformDerivatives(u1,v1,endpoints);
    %R = 0.5*alpha*sum(Du.^2 + Dv.^2);
       
    plot(i,J,'.b');
    hold on;
    %plot(i,R,'.r')
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

function  [Du,Dv] = transformDerivatives(u,v,endpoints)
    deriv=[-1,0,1]/2;
    %D = imfilter([u,v],deriv','replicate'); %NOTE: transpose on deriv
    %Dv = imfilter(v,deriv','replicate'); 
    phi = [u,v];
    D(1:endpoints(1),:) = imfilter(phi(1:endpoints(1),:),deriv',0);
    for i = 1:numel(endpoints)-1
         
        a = endpoints(i);
        b = endpoints(i+1);
        
        D(a:b,:) = imfilter(phi(a:b,:),deriv',0);
                
    end
    D(endpoints(end):length(u),:) = imfilter(phi(endpoints(end):length(u),:),deriv',0);
    
    D(D > 0.1) = [];
    
    Du = D(:,1);
    Dv = D(:,2);
end

function [u1,v1]=solveFlow(Ix,Iy,u0,v0,alpha,endpoints)
N = length(Ix);

%laplace operator
x1=ones(N,1).*alpha; x1(1)=0;
x2=ones(N,1).*alpha; x2(end)=0;
S=spdiags([x1(:),x2(:)],[1,-1],N,N);

C=sum(S,1);
x1(endpoints + 1) = 0;
x2(endpoints(endpoints > 1) - 1) = 0;
C(endpoints) = 1;
C(1) = 1;
C(end) = 1;

L=spdiags([x1(:),x2(:),-C(:)],[1,-1,0],N,N); 
%replicate L for u,v
A=[L,sparse(N,N);sparse(N,N),L];

%set-up linear system of equations

% b=A*double([u0(:);v0(:)])-double([Ix(:);Iy(:)]);
b=double([Ix(:);Iy(:)]); %only incremental

tic;
%solve Gauss-Newton update-step with Diffusion Regularisation
uv1=A\b;

%[uv1,flagpcg]=pcg(A,b,[],20);
%[uv1,flagbig]=bicgstab(A,b,1E-2,20,[],[],double([u0(:);v0(:)]));
%[uv1,flag,res,nit]=sor(A',b,1.9,25,1E-4,double([u0(:);v0(:)]));
tOp=toc;
u1=reshape(uv1(1:N),N,1);
v1=reshape(uv1(N+1:end),N,1);

u1=u1+u0;
v1=v1+v0;
end



