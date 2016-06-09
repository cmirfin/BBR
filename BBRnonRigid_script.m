%boundary script

addpath('~/Documents/ONBI-Project1/HistoRegModified/RegCode/NIfTI_20140122/');
addpath('~/Documents/ONBI-Project1/HistoRegModified/RegCode/OpticalFlowMIND/');

%% load images
%path1 =  '~/Documents/ONBI-Project1/HistoRegModified/ExampleData/mge3d.nii.gz';
%path2 = '~/Documents/ONBI-Project1/HistoRegModified/ExampleData/PLI/Transmittance_CC01.tif';

%[im1, im2] = load_images(path1,path2);

load('fixedImage.mat'); %loads unaltered PLI fixed image
load('inputMIND.mat'); %loads moving image from output of MIND analysis 
%% pre-processing
fixedImage(fixedImage == 0) = 255; %get rid off patched image edges.
fixedImage = medfilt2(fixedImage,[5,5]);
fixedImage(:,end-5:end) = 255;
fixedImage(:,1:5) = 255;

% [m,n] = size(movingImage);
% backgroundMovingImage = reshape(kmeans(movingImage(:),3),m,n);
% backgroundClass = backgroundMovingImage(1,1); %K-means classifcation of background
% backgroundMovingImage(backgroundMovingImage==backgroundClass) = 0; %background
% backgroundMovingImage(backgroundMovingImage>0) = 1; %not background

%movingImage = medfilt2(movingImage,[5,5]);
%% find boundary and projection points in moving image (im1)
DeltaIn = 2; %projection distance
DeltaOut = 2;

[boundaryPoints, normals] = boundaryNormal(movingImage,max(DeltaIn,DeltaOut)); %NOTE: normals are unit vectors with direction
%[boundaryPoints, normals] = tissueBoundary(backgroundMovingImage,backgroundMovingImage,max(DeltaIn,DeltaOut));
[boundaryPoints,normals,endpoints,startpoints] = sortBoundaries(boundaryPoints,normals); %NOTE: add third input to graphically display boundaries
 %remove start and final point

fixedPoints = [startpoints;endpoints];
%% solve

[u,v] = boundaryCostNonRigid(boundaryPoints,normals,fixedImage,fixedPoints,0.5);

%% boundary visualization

tpoints = round([u,v] + boundaryPoints);
%boundaryImg = zeros(size(movingImage));
boundaryImg = fixedImage;
for i = 1:length(boundaryPoints)
    boundaryImg(boundaryPoints(i,1),boundaryPoints(i,2)) = 100;
    boundaryImg(tpoints(i,1),tpoints(i,2)) = 255;
    
end
figure;
imagesc(boundaryImg)

%% interpolating vector field

%intialize flow fields
[m,n] = size(movingImage);
[uF,vF] = elasticSolver(u,v,boundaryPoints,m,n);
figure;
quiver(uF,vF);
outputImage = transformInterpolation2d(movingImage,uF,vF);
figure;
imagesc(outputImage);
%outputImage = medfilt2(outputImage,[5,5]);