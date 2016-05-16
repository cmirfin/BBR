%script to test cost function
clear cost; clear grad;
[shiftX, shiftY] = meshgrid(-5:1:5);

for i = 1:size(shiftX,1)
    for j =1:size(shiftY,2)
    
        theta = [shiftX(i,j),shiftY(i,j)];
        [cost(i,j), gradtmp] = boundaryCost2(theta,rIn,rOut,fixedImage,Fx,Fy);
        gx(i,j) = gradtmp(1);
        gy(i,j) = gradtmp(2);
    end
end

figure;
surf(shiftX,shiftY,cost);
%contour(shiftX,shiftY,cost);
% figure;
% scatter(shift,cost);
% 
% title('Cost as function of theta')
% ylabel('Cost')
% xlabel('x shift')
% 
% figure;
% scatter(shift,grad);
% title('Gradient as function of theta')
% ylabel('Gradient')
% xlabel('Theta(1)')

