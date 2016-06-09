function [uFull,vFull] = vectorFieldExtrapolation(uConstrained,vConstrained, pixelPositions,m,n)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

x1=ones(m,n); x1(1,:)=0;
x2=ones(m,n); x2(end,:)=0;
y=ones(m,n);
S=spdiags([x1(:),x2(:),y(:),y(:)],[1,-1,m,-m],m*n,m*n);
C=sum(S,1);
L=spdiags([x1(:),x2(:),y(:),y(:),-C(:)],[1,-1,m,-m,0],m*n,m*n); 

pc = pixelPositions(:,1) + (pixelPositions(:,2)-1)*m; %convert to positions in image vector
pv = 1:m*n;
mask = pv(:);
mask(pc) = 0;
mask(mask > 0) = 1;
pv(pc) = [];

Lc = L(:,pc);
Lv = L(:,pv);

% r = qr(Lv,0);
% inverseLv = r\(r'\Lv');
uVariable = -Lv\(Lc*uConstrained);
vVariable = -Lv\(Lc*vConstrained);

%put in original vectors and reshape
uFull = zeros(m*n,1);
vFull = uFull;

uFull(mask == 0) = uConstrained;
uFull(mask == 1) = uVariable;

vFull(mask == 0) = vConstrained;
vFull(mask == 1) = vVariable;

uFull = reshape(uFull,m,n);
vFull = reshape(vFull,m,n);
end

