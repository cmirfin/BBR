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
%fixedImage = 255 - fixedImage;
%fixedImage(:,end-5:end) = 1e-2;
fixedImage(:,end-5:end) = 255;
movingImage = medfilt2(movingImage,[5,5]);
%% find boundary and projection points in moving image (im1)
DeltaIn = 2; %projection distance
DeltaOut = 2;

[boundaryPoints, normals] = boundaryNormal(movingImage,max(DeltaIn,DeltaOut)); %NOTE: normals are unit vectors with direction
[order,endpoints] = sortManyBoundaryPoints(boundaryPoints);
%[order,endpoints,chainLengths] = sortBoundaries(boundaryPoints);
boundaryPoints = boundaryPoints(order,:);
normals = normals(order,:);

[Fx,Fy] = gradient(fixedImage);
%% solve

[u,v] = boundaryCostNonRigid(boundaryPoints,normals,fixedImage,Fx,Fy,endpoints);


%% boundary visualization

tpoints = round([u,v] + boundaryPoints);
boundaryImg = zeros(size(movingImage));
for i = 1:length(boundaryPoints)
    boundaryImg(boundaryPoints(i,1),boundaryPoints(i,2)) = 100;
    boundaryImg(tpoints(i,1),tpoints(i,2)) = 200;
    
end
figure;
imagesc(boundaryImg)

%% interpolating vector field

%intialize flow fields
u0 = zeros(size(movingImage));
v0 = u0;

[uF,vF] = elasticSolver(u,v,boundaryPoints,u0,v0);
figure;
quiver(uF,vF);
outputImage = transformInterpolation2d(movingImage,uF,vF);
figure;
imagesc(outputImage);
%outputImage = medfilt(outputImage,[5,5]);