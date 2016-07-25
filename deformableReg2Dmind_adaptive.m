% Code by MP Heinrich

function [u1,v1,deformed]=deformableReg2Dmind_adaptive(im1,im2,alpha)

if nargin<4
    flag2=false;
end

im1=single(im1); im2=single(im2);

levels=[4,2,1]; % resolution levels


u1=zeros(2,2,'single'); v1=u1;u2=u1; v2=u1;
% initialise flow fields with 0

timeOpt=0;
warps=[20,20,10]; % number of warps per level (increase if needed)


count=1;
figure;
for j=1:length(levels)
    
    hs=1;%fspecial('gaussian',[15,15],levels(j)/2);
    
    %pre-compute mind
%     if(flag2==1)
%         im2f=imresize(imfilter(im2,hs,'replicate'),size(im1)./levels(j),'antialiasing',false);
%         im1f=imresize(imfilter(im1,hs,'replicate'),size(im1)./levels(j),'antialiasing',false);
%         mind1=MIND_descriptor2D(im1f,0,1);
%         mind2o=MIND_descriptor2D(im2f,0,1);
%     end
    
    maxwarp=warps(j);
    
    % resize images for current level
    [u1,v1]=resizeFlow2d(u1,v1,ceil(size(im1)./levels(j)));
    % upsample flow to current level
    
    for i=1:maxwarp
        count=count+1;
        [ux,vx]=resizeFlow2d(u1,v1,size(im1));
        warped1=imWarp(ux,vx,im2);
        
        warped1f=imresize(imfilter(warped1,hs,'replicate'),size(im1)./levels(j),'antialiasing',false);

        if(flag2==1)
            mind2=imWarp(u1,v1,mind2o);
            
        else
            
            im1f=imresize(imfilter(im1,hs,'replicate'),size(im1)./levels(j),'antialiasing',false);
            % transform images to intermediate space
            
            % extract modality independent neighbourhood descriptor
            mind1=MIND_descriptor2D(im1f,1,1);
            mind2=MIND_descriptor2D(warped1f,1,1);
        end
        % calculate derivates of MIND
        [Sx,Sy,S]=MINDgrad2s(mind1,mind2);%warped2f,warped1f);%
        
        % solve Euler-Lagrange equations using mldivide or iterative
        % solvers (like conjugate gradient, SOR)
        % includes (compositive) diffusion regularisation in cost term
        
        [u1,v1,tOp1]=solveFlowMldivide(im1f,Sx,Sy,S,u1,v1,alpha);
        % median filer to remove noisy estimates (w/o smoothing over edges)
        u1=medfilt2(u1,[3,3]);
        u2=medfilt2(u2,[3,3]);
        
        %end
        timeOpt=timeOpt+tOp1;
        
        % visualise result
        subplot(1,2,1); showMap(ux,vx);
        subplot(1,2,2); imagesc(im1(:,:,1)-warped1(:,:,1)); axis image; axis off; colormap gray(256); drawnow;
    end
    
end


if(levels(end)>1)
    [u1,v1]=resizeFlow2d(u1,v1,size(im1));
end
deformed=imWarp(u1,v1,im2);
% generate final output volumes

timeOpt
end

function [Sx,Sy,S,Sxi,Syi]=MINDgrad2s(mind1,mind2)
%symmetric calculation of derivatives of MIND
[m,n,p]=size(mind1);


%deriv=[-1,8,0,-8,1]/12; %5-point derivative filter (as recommended by W. Yu, Francetti 2008).
deriv=[-1,0,1]/2;

Sx=imfilter(mind2,deriv,'replicate');
Sy=imfilter(mind2,deriv','replicate'); %NOTE transpose on deriv

Sxi=imfilter(mind1,deriv,'replicate');
Syi=imfilter(mind1,deriv','replicate');

S=mind2-mind1;

Sx(:,[1,n],:)=0;
Sy([1,m],:,:)=0;

Sxi(:,[1,n],:)=0;
Syi([1,m],:,:)=0;


end

function [u1,v1,tOp]=solveFlowMldivide(I,Ix,Iy,It,u0,v0,alpha)
[m,n,o]=size(Ix);

sigma2=2*8.^2;

%laplace operator that penalises deviations for
%edges with low intensity difference (adaptive)
% define Laplacian matrix
w_d=alpha*exp(-(I(2:end,:)-I(1:end-1,:)).^2/sigma2);
w_d(end+1,:)=0;
w_r=alpha*exp(-(I(:,2:end)-I(:,1:end-1)).^2/sigma2);
w_r(:,end+1)=0;
w_l=alpha*exp(-(I(:,2:end)-I(:,1:end-1)).^2/sigma2);
w_l=[zeros(m,1),w_l];
w_u=alpha*exp(-(I(2:end,:)-I(1:end-1,:)).^2/sigma2);
w_u=[zeros(1,n);w_u];
d=w_l+w_r+w_u+w_d+1e-6;
L=spdiags([-d(:),w_d(:),w_u(:),w_r(:),w_l(:)],[0,-1,1,-m,m],m*n,m*n);

%replicate L for u,v
M=[L,sparse(m*n,m*n);sparse(m*n,m*n),L];

%structure tensor
hs=1;%fspecial('gaussian',[7,7],0.5);
Ix=single(Ix); Iy=single(Iy); It=single(It); u0=single(u0); v0=single(v0);
%local smoothing (similar to Bruhn et al. 2005 IJCV)
Ixx=mean(imfilter(Ix.^2,hs,'replicate'),3);
Ixy=mean(imfilter(Ix.*Iy,hs,'replicate'),3);
Ixt=mean(imfilter(Ix.*It,hs,'replicate'),3);
Iyy=mean(imfilter(Iy.^2,hs,'replicate'),3);
Iyt=mean(imfilter(Iy.*It,hs,'replicate'),3);

dIxx=spdiags(Ixx(:),0,m*n,m*n);
dIxy=spdiags(Ixy(:),0,m*n,m*n);
dIyy=spdiags(Iyy(:),0,m*n,m*n);

%set-up linear system of equations
A=double([dIxx,dIxy;dIxy,dIyy])-M;
b=M*double([u0(:);v0(:)])-double([Ixt(:);Iyt(:)]);
%b=-double([Ixt(:);Iyt(:)]); %only incremental

tic;
%solve Gauss-Newton update-step with Diffusion Regularisation
uv1=A\b;
%[uv1,flagpcg]=pcg(A,b,[],20);
%[uv1,flagbig]=bicgstab(A,b,1E-2,20,[],[],double([u0(:);v0(:)]));
%[uv1,flag,res,nit]=sor(A',b,1.9,25,1E-4,double([u0(:);v0(:)]));
tOp=toc;
u1=reshape(uv1(1:m*n),m,n);
v1=reshape(uv1(m*n+1:end),m,n);

u1=u1+u0;
v1=v1+v0;
end

function [u2,v2]=resizeFlow2d(u1,v1,size1)

size2=size(u1);
scales=(size1./size2);

[x,y]=meshgrid(0:size1(2)-1,0:size1(1)-1);
x=single(x);
y=single(y);

x=max(min(x./scales(2)+1,size2(2)),1);
y=max(min(y./scales(1)+1,size2(1)),1);


u2=interp2(single(u1).*scales(2),x,y,'*linear');
v2=interp2(single(v1).*scales(1),x,y,'*linear');

end

function [ B ] = imWarp( flowHor, flowVer, Bin, method )
if nargin<4
    method='linear';
end
%This function warps B towards A
[m,n,o]=size(Bin);

[x y] = meshgrid(1:size(Bin,2),1:size(Bin,1));
xu=x+flowHor;
xu(xu<1)=1;
xu(xu>n)=n;

xu(isnan(xu))=0;
yv=y+flowVer;
yv(yv<1)=1;
yv(yv>m)=m;
yv(isnan(yv))=0;

for i=1:size(Bin,3)
    B(:,:,i) = interp2(Bin(:,:,i), xu,yv, method);
end
B(isnan(B)) = Bin(isnan(B));
end

function [u1s,v1s,u2s,v2s]=symmetricMapping(u1,v1,u2,v2,it)
if nargin<5
    it=10;
end

[u2i,v2i]=fastInverse(u2./2,v2./2,it);
[u1i,v1i]=fastInverse(u1./2,v1./2,it);


[u1s,v1s]=combineDeformation(u1./2,v1./2,u2i,v2i,'compositive');
[u2s,v2s]=combineDeformation(u2./2,v2./2,u1i,v1i,'compositive');

end

function [u_combined,v_combined]=combineDeformation(u1st,v1st,u2nd,v2nd,method)

if nargin<5
    method='compositive';
end
u2nd(isnan(u2nd))=0;
u1st(isnan(u1st))=0;
v2nd(isnan(v2nd))=0;
v1st(isnan(v1st))=0;

if strcmp(method,'compositive')
    u_combined=imWarp(u2nd,v2nd,u1st)+u2nd;
    v_combined=imWarp(u2nd,v2nd,v1st)+v2nd;
else
    u_combined=u1st+u2nd;
    v_combined=v1st+v2nd;
end

end

function [u2,v2]=fastInverse(u1,v1,it)
if nargin<3
    it=10;
end
[m,n]=size(u1);
[x,y]=meshgrid(1:n,1:m);
u2=zeros(size(u1));
v2=zeros(size(v1));

for i=1:it
    u2n=-interp2(u1,min(max(x+u2,1),n),min(max(y+v2,1),m),'linear');
    v2n=-interp2(v1,min(max(x+u2,1),n),min(max(y+v2,1),m),'linear');
    u2=u2n;
    v2=v2n;
end

end

function rgbim=showMap(x,y)

hsvim=ones([size(x),3]);
[theta,rho]=cart2pol(-x,y);
[hist1,ind1]=hist(rho(:),100);
max90=ind1(find(cumsum(hist1)>(0.9.*sum(hist1).*ones(1,length(hist1))),1));
rho2=min(rho,max90)./max90;

[x2,y2]=pol2cart(theta,rho2);
[theta,rho2]=cart2pol(x2,y2);

theta2=(-theta+pi)./(2*pi);

hsvim(:,:,1)=theta2;
hsvim(:,:,2)=rho2;
%hsvim(:,:,3)=rho2./2+0.5;
rgbim=hsv2rgb(hsvim);

imagesc(rgbim); axis image; axis off;

end