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

A = [1 0 8; 0 1 0; 0 0 1]';
tform = affine2d(A);

movingImage = imwarp_same(fixedImage,tform);
%% find boundary and projection points in moving image (im1)
DeltaIn = 2; %projection distance
DeltaOut = 2;

[boundaryPoints, normals] = boundaryNormal(movingImage,max(DeltaIn,DeltaOut)); %NOTE - normals are given as absolute values.

[Fx,Fy] = gradient(fixedImage);
%% optimization

initialTransform = 6*ones(length(boundaryPoints),2);

[cost,phi] = boundaryCostNonRigid(initialTransform,boundaryPoints,normals,fixedImage,Fx,Fy);
scatter(1:length(cost),cost);

%% visualization

tpoints = round(phi + boundaryPoints);
boundaryImg = zeros(size(fixedImage));
for i = 1:length(tpoints)
    boundaryImg(tpoints(i,1),tpoints(i,2)) = 1;
    boundaryImg(boundaryPoints(i,1),boundaryPoints(i,2)) = 0.5;
end
imagesc(boundaryImg)