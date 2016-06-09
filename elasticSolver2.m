function [uFull,vFull] = elasticSolver2(u,v,pixelPositions,backgroundIdx,img)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
mu = 1e-9;

[m,n] = size(img);
u0 = zeros(m,n);
v0 = u0;

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



% 
% background = background(:);
% backgroundIdx = fbackgroundIdx(background == 0);
bu(p) = u(:);
bv(p) = v(:);
bu(backgroundIdx) = [];
bv(backgroundIdx) = [];

b = [bu(:);bv(:)];


D = zeros(m*n,1);
D(p) = 1;
D(backgroundIdx) = [];
C(backgroundIdx) = [];
x1(backgroundIdx) = [];
x2(backgroundIdx) = [];
y(backgroundIdx) = [];

C = D(:)+C(:);
N = m*n-length(backgroundIdx);
L=spdiags([x1(:),x2(:),y(:),y(:),-C(:)],[1,-1,m,-m,0],N,N); 
%replicate L for u,v
A=[L,sparse(N,N);sparse(N,N),L];

uv1=-A\b;


% uFull=reshape(uv1(1:m*n),m,n);
% vFull=reshape(uv1(m*n+1:end),m,n);

uFull=uv1(1:N);
vFull=uv1(N+1:end);

%interpolation
% img = double(img(:));
% img(backgroundIdx) = []; %remove background pixels
% % 
% [X,Y] = meshgrid(1:n,1:m);
% % 
% X0 = X(:);
% Y0 = Y(:);
% % 
% X0(backgroundIdx) = [];
% Y0(backgroundIdx) = [];
% % 
% Xp = X0 + uFull;
% Yp = Y0 + vFull;
% 
% %convert to 2d coordinates
% 
% 
% %outputImg = griddata(Yp,Xp,img,Y,X);
% F = scatteredInterpolant(Yp,Xp,img,'linear','linear');
% outputImg = F(Y,X);

end

