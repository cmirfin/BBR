function [uFull,vFull] = elasticSolver(u,v, pixelPositions,u0,v0)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
mu = 1e-9;

[m,n] = size(u0);

x1=mu*ones(m,n); x1(1,:)=0;
x2=mu*ones(m,n); x2(end,:)=0;
y=mu*ones(m,n);
S=spdiags([x1(:),x2(:),y(:),y(:)],[1,-1,m,-m],m*n,m*n);
C=sum(S,1);
% L=spdiags([x1(:),x2(:),y(:),y(:),-C(:)],[1,-1,m,-m,0],m*n,m*n); 
% %replicate L for u,v
% A=[L,sparse(m*n,m*n);sparse(m*n,m*n),L];

%set-up linear system of equations

bu = u0(:);
bv = v0(:);
p = pixelPositions(:,1) + (pixelPositions(:,2)-1)*m; %convert to positions in image vector

bu(p) = u(:);
bv(p) = v(:);

b = [bu(:);bv(:)];


D = zeros(m*n,1);
D(p) = 1;

C = D(:)+C(:);
L=spdiags([x1(:),x2(:),y(:),y(:),-C(:)],[1,-1,m,-m,0],m*n,m*n); 
%replicate L for u,v
A=[L,sparse(m*n,m*n);sparse(m*n,m*n),L];

%b=-double([Ixt(:);Iyt(:)]); %only incremental

tic;
%solve Gauss-Newton update-step with Diffusion Regularisation
uv1=-A\b;
%[uv1,flagpcg]=pcg(A,b,[],20);
%[uv1,flagbig]=bicgstab(A,b,1E-2,20,[],[],double([u0(:);v0(:)]));
%[uv1,flag,res,nit]=sor(A',b,1.9,25,1E-4,double([u0(:);v0(:)]));
tOp=toc;
uFull=reshape(uv1(1:m*n),m,n);
vFull=reshape(uv1(m*n+1:end),m,n);

%u1=u1+u0;
%v1=v1+v0;

end

