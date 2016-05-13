function [J, grad] = boundaryCost(theta,movingImg,fixedImg)

J = 0;
grad = zeros(size(theta));
hessian = zeros(size(theta,1),size(theta,1));

movingImg = imguidedfilter(movingImg,'NeighborhoodSize',[4,4]);
fixedImg = imguidedfilter(fixedImg,'NeighborhoodSize',[4,4]);

%parameters
M = 0.01;
DeltaW = 2; %white matter projection distance
DeltaG = 2; %gray matter projection distance

[m,n] = size(fixedImg);

tImg = transformImage(movingImg,theta); %transform movingImg by theta.

[Fx,Fy] = gradient(fixedImg);


[boundary,normals] = boundaryNormal(tImg,max(DeltaW,DeltaG));%get normals to boundary
N = length(boundary); %find number of boundary points

%get x,y coordinaties of boundary pixels
rPointsIn = boundary - round(DeltaG.*[normals(:,1),normals(:,2)]);
rPointsOut = boundary + round(DeltaW.*[normals(:,1),normals(:,2)]);

%plotBoundaries(rPointsIn,rPointsOut,boundary,fixedImg);

for i = 1:N
        
    xIn = rPointsIn(i,1);
    yIn = rPointsIn(i,2);
    xOut = rPointsOut(i,1);
    yOut = rPointsOut(i,2);
    
    %evaluate points either side of boundary
    valueIn = fixedImg(xIn,yIn);
    valueOut = fixedImg(xOut,yOut);
    
%     thetaJacobianIn = [xIn 0 yIn 0 1 0; 0 xIn 0 yIn 0 1];
%     thetaJacobianOut = [xOut 0 yOut 0 1 0; 0 xOut 0 yOut 0 1];

    if size(theta,1) < 2
        thetaJacobianIn = [1 ;1];
        thetaJacobianOut = [1;1];
    else
        thetaJacobianIn = [1 0; 0 1];
        thetaJacobianOut = [1 0; 0 1];
    end
    
    derivIn = [Fx(xIn,yIn), Fy(xIn,yIn)]*thetaJacobianIn;
    derivOut = [Fx(xOut,yOut), Fy(xOut,yOut)]*thetaJacobianOut;
    
    %Update cost function and gradient
    Q = abs(M*(100*(valueIn-valueOut)/(0.5*(valueIn+valueOut))));
    J = J + 1 - tanh(Q);
    grad = grad + (M*(sech(Q).^2)*2*(abs(valueOut.*derivIn) - abs(valueIn.*derivOut))./(valueIn+valueOut).^2);
    %grad = grad + (M*2*(abs(valueOut.*derivIn) - abs(valueIn.*derivOut))./(valueIn+valueOut).^2);
 
end

%J = double(sum((tImg(:) - fixedImage(:)).^2)); %SSD
J = double(J./N);
grad = (200/N)*(double(grad));

%fprintf('Cost: %f\n', J);
%disp(grad);

end

function tImg= transformImage(img, t)
p = zeros(1,6);
if length(t) < 2
    p(6) = t(1);
elseif length(t) < 3
    p(5:6) = t(1:2);
else
    p = t;
end
% p(6) = t;
A = [1+p(1) p(3) p(5); p(2) 1+p(4) p(6); 0 0 1]';

% A = [1 t(3) t(1);t(4) 1 t(2); 0 0 1]';

tform = affine2d(A);

tImg = imwarp_same(img, tform,'bilinear');

end

function plotBoundaries(pointsIn, pointsOut,boundary,img)

array = zeros(size(img));

len = length(pointsIn);

for i=1:len
    x = boundary(i,1);
    y = boundary(i,2);
    array(x,y) = 1;
    x = pointsOut(i,1);
    y = pointsOut(i,2);
    array(x,y) = 0.5;
    x = pointsIn(i,1);
    y = pointsIn(i,2);
    array(x,y) = 0.8;
  
end

    imagesc(array);
    drawnow;

end

