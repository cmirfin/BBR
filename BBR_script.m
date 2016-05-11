
addpath('~/Documents/ONBI-Project1/HistoRegModified/RegCode/NIfTI_20140122/');
addpath('~/Documents/ONBI-Project1/HistoRegModified/RegCode/OpticalFlowMIND/');

%% load images
path1 =  '~/Documents/ONBI-Project1/HistoRegModified/ExampleData/mge3d.nii.gz';
path2 = '~/Documents/ONBI-Project1/HistoRegModified/ExampleData/PLI/Transmittance_CC01.tif';

[im1, im2] = load_images(path1,path2);

%% pre-processing

im1 = medfilt2(im1,[5,5]);

%% create artificial image

A = [1 0 4; 0 1 6; 0 0 1]';
tform = affine2d(A);
im2 = imwarp_same(im1,tform);

%% Optimization of cost function using local solver

%  Set options for fminunc
options = optimset('GradObj', 'on', 'MaxIter', 500,'Display','iter');
% initial_theta = [0,0,0,0,0,4.5];
initial_theta = [7, 7];

%  Run fminunc to obtain the optimal theta
%  This function will return theta and the cost 
[theta, cost,exitflag,output] = fminunc(@(t)(boundaryCost(t, im1, im2)), initial_theta, options);

%options = optimset('OutputFcn', @outfun);
%[theta, cost] = fminsearch(@(t)(boundaryCost(t, im1, im2)), initial_theta);

% Print theta to screen
fprintf('Cost at theta found by fminunc: %f\n', cost);
fprintf('theta: \n');
fprintf(' %f \n', theta);   

%% Global solver approach
tic;
opts = optimoptions(@fmincon,'Algorithm','interior-point');
problem = createOptimProblem('fmincon','objective',@(x) (boundaryCost(x,im1,im2)),'x0',[3,3],'lb',[0,0],'ub',[8,8],'options',opts);
gs = GlobalSearch;
[theta,cost] = run(gs,problem)
toc

%% peform transformation
% A = [1 + theta(1) theta(3) theta(5); theta(2) 1 + theta(4) theta(6); 0 0 1]'; 
A = [1 0 theta(1); 0 1 theta(2); 0 0 1]'; 
tform = affine2d(A);
registeredIm = imwarp_same(im1, tform);
figure;
imagesc(registeredIm);