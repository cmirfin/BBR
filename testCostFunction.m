%script to test cost function

shift = [1:0.5:10];
shift = [-shift, 0, shift];

for i = 1:length(shift)
    
    
    [cost(i),grad(i)] = boundaryCost(shift(i),im1,im1);
    
end

figure;
scatter(shift,cost);

title('Cost as function of theta')
ylabel('Cost')
xlabel('x shift')

figure;
scatter(shift,grad);
title('Gradient as function of theta')
ylabel('Gradient')
xlabel('Theta(1)')

