%% addpaths

addpath('~/Documents/ONBI-Project1/HistoRegModified/RegCode/NIfTI_20140122/');
addpath('~/Documents/ONBI-Project1/HistoRegModified/RegCode/OpticalFlowMIND/');

filename1 = '~/Documents/ONBI-Project1/HistoRegModified/ExampleData/mge3d.nii.gz';
filename2 = '~/Documents/ONBI-Project1/HistoRegModified/ExampleData/PLI/Transmittance_CC02.tif';

%% read in MRI image (and average MGE)

vol1 = load_nii(filename1);
vol1=mean(vol1.img,4);
im1_orig=reshape(vol1(96+1-71,:,:),128,256);
inputImage = imresize(im1_orig(8:117,37:223)/530,[220,355]);
%inputImage = padarray(inputImage,[20,20],'both');

clear vol1; clear im1_orig;

%% read in PLI image and preprocess
im2_orig = imread(filename2);
fixedImage = padarray(imresize(single(flipud(im2_orig)),[210,355]),[5,0],255);
fixedImage(fixedImage <= 0) = 255; %get rid off patched image edges.
fixedImage(fixedImage > 255) = 255;
fixedImage = straightLineRemover(fixedImage,25,255); %vertical direction
fixedImage = straightLineRemover(fixedImage',25,255)'; %NOTE transposes - horizontal direction
fixedImage = medfilt2(fixedImage,[3,3]);
% fixedImage(:,end-5:end) = 255;
% fixedImage(:,1:5) = 255;
%fixedImage = padarray(fixedImage,[20,20],'both');
%fixedImage(fixedImage <=0) = 255;

clear im2_orig;
%% perform registration using locally adaptive MIND

%[u1,v1,movingImage]=deformableReg2Dmind_adaptive(fixedImage,inputImage,.3);
[u1,v1,movingImage]=deformableReg2Dmind_asym(fixedImage,inputImage,.3);
rgb_disc=gray2blue(movingImage)+gray2orange(fixedImage);
rgb_before=gray2blue(inputImage)+gray2orange(fixedImage);

magnitude=grey2parula(min(sqrt(u1.^2+v1.^2),20)./20);
%% write output images
% imwrite(magnitude,'mri_hist_magn_disc.png');
% imwrite(rgb_before,'mri_hist_overlay_before.png');
% imwrite(rgb_disc,'mri_hist_overlay_disc.png');
% imwrite(im1/255,'mri_hist_mri.png');
% imwrite(im2/255,'mri_hist_hist.png');
% imwrite(movingImage/255,'warped_mri.png');

%% find boundary and projection points in moving image (im1)
DeltaIn = 3; %projection distance
DeltaOut = 3;

[boundaryPoints, normals] = boundaryNormal(movingImage,max(DeltaIn,DeltaOut)); %NOTE: normals are unit vectors with direction
%[boundaryPoints, normals] = tissueBoundary(backgroundMovingImage,backgroundMovingImage,max(DeltaIn,DeltaOut));
[boundaryPoints,normals,endpoints,startpoints] = sortBoundaries(boundaryPoints,normals); %NOTE: add third input to graphically display boundaries
 %remove start and final point

fixedPoints = [startpoints;endpoints];

%% solve

[u,v] = boundaryCostNonRigid(boundaryPoints,normals,fixedImage,fixedPoints,DeltaIn,DeltaOut,0.3);

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

%% Extrapolation of vector field and interpolation of image

outputImage = weightedInterpolation(movingImage,u,v,boundaryPoints,'elasticSolver',1); %NOTE: specify extrapolaton scheme

figure;
imagesc(outputImage);
%outputImage = medfilt2(outputImage,[5,5]);