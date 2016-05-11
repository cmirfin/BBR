%script to test cost function

[shiftX, shiftY] = meshgrid(-4:0.5:10);

for i = 1:size(shiftX,1)
    for j =1:size(shiftY,2)
    
        theta = [shiftX(i,j),shiftY(i,j)];
        [cost(i,j)] = boundaryCost(theta,im1,im2);
    end
end

figure;
%surf(shiftX,shiftY,cost);
contour(shiftX,shiftY,cost);
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

