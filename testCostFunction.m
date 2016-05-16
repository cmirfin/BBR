%script to test cost function
clear cost; clear grad;
shiftX = [-10:0.1:10];
for i = 1:length(shiftX)
    theta = shiftX(i);
    [cost(i),grad(i)] = boundaryCost2(theta,rIn,rOut,fixedImage,Fx,Fy);

end

figure;
scatter(shiftX,cost);

title('Cost as function of theta')
ylabel('Cost')
xlabel('x shift')

figure;
scatter(shiftX,grad);
title('Gradient as function of theta')
ylabel('Gradient')
xlabel('Theta(1)')


