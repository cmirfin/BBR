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
% fixedImage = 255 - fixedImage;
% fixedImage(:,end-5:end) = 0;
fixedImage(:,end-5:end) = 255;
movingImage = medfilt2(movingImage,[5,5]);
%% find boundary and projection points in moving image (im1)
DeltaIn = 2; %projection distance
DeltaOut = 2;

[boundaryPoints, normals] = boundaryNormal(movingImage,max(DeltaIn,DeltaOut)); %NOTE - normals are given as absolute values.
%[order,endpoints] = sortBoundaryPoints(boundaryPoints);
[order,endpoints] = sortManyBoundaryPoints(boundaryPoints);
boundaryPoints = boundaryPoints(order,:);
normals = normals(order,:);

[Fx,Fy] = gradient(fixedImage);
%% solve

[u,v] = boundaryCostNonRigid(boundaryPoints,normals,fixedImage,Fx,Fy,endpoints);


%% visualization

tpoints = round([u,v] + boundaryPoints);
boundaryImg = zeros(size(fixedImage));
for i = 1:length(boundaryPoints)
    boundaryImg(boundaryPoints(i,1),boundaryPoints(i,2)) = 50;
    boundaryImg(tpoints(i,1),tpoints(i,2)) = 100;
    
end
figure;
imagesc(boundaryImg)