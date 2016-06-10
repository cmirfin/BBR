function [uFull,vFull] = elasticSolver(u,v,pixelPositions,m,n)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
mu = 1;
%[m,n] = size(u0);
u0 = zeros(m,n);
v0 = u0;

x1=mu*ones(m,n); x1(1,:)=0;
x2=mu*ones(m,n); x2(end,:)=0;
y=mu*ones(m,n);
S=spdiags([x1(:),x2(:),y(:),y(:)],[1,-1,m,-m],m*n,m*n);
C=sum(S,1);

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


%solve Gauss-Newton update-step with Diffusion Regularisation
uv1=-A\b;

uFull=reshape(uv1(1:m*n),m,n);
vFull=reshape(uv1(m*n+1:end),m,n);

end

