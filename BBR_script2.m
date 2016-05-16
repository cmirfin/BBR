%boundary script

addpath('~/Documents/ONBI-Project1/HistoRegModified/RegCode/NIfTI_20140122/');
addpath('~/Documents/ONBI-Project1/HistoRegModified/RegCode/OpticalFlowMIND/');

%% load images
path1 =  '~/Documents/ONBI-Project1/HistoRegModified/ExampleData/mge3d.nii.gz';
path2 = '~/Documents/ONBI-Project1/HistoRegModified/ExampleData/PLI/Transmittance_CC01.tif';

[im1, im2] = load_images(path1,path2);

%% pre-processing

fixedImage = medfilt2(im1,[5,5]);
%% create artificial image

A = [1 0 0; 0 1 0; 0 0 1]';
tform = affine2d(A);
%movingImage = imwarp(fixedImage,tform,'OutputView');
movingImage = imwarp_same(fixedImage,tform);
%% find boundary and projection points in moving image (im1)
DeltaIn = 2; %projection distance
DeltaOut = 2;

[boundaryPoints, normals] = boundaryNormal(movingImage,max(DeltaIn,DeltaOut)); %NOTE - normals are given as absolute values.
rIn = boundaryPoints - round(DeltaIn.*[normals(:,1),normals(:,2)]);
rOut = boundaryPoints + round(DeltaOut.*[normals(:,1),normals(:,2)]);

[Fx,Fy] = gradient(fixedImage);
%% optimization

%  Set options for fminunc
options = optimset('GradObj', 'on', 'MaxIter', 500,'Display','iter');
% initial_theta = [0,0,0,0,0,4.5];
initial_theta = [3,2];

%  Run fminunc to obtain the optimal theta
%  This function will return theta and the cost 
[theta, cost,exitflag,output] = fminunc(@(t)(boundaryCost2(t,rIn,rOut,fixedImage,Fx,Fy)), initial_theta, options);

%options = optimset('OutputFcn', @outfun);
%[theta, cost] = fminsearch(@(t)(boundaryCost(t, im1, im2)), initial_theta);

% Print theta to screen
fprintf('Cost at theta found by fminunc: %f\n', cost);
fprintf('theta: \n');
fprintf(' %f \n', theta);   
