% Christopher J. Mirfin 
% Sir Peter Mansfield Imaging Centre, University of Nottingham
% christopher.mirfin@dtc.ox.ac.uk
% 23/07/2016

% For MIND code (stage 1)- Code by MP Heinrich
% MP Heinrich, M Jenkinson, M Bhushan, T Matin, F Gleeson, M Brady, JA Schnabel.
% MIND: Modality Independent Neighbourhood Descriptor for Multi-modal Deformable Registration 
% Medical Image Analysis. vol. 16(7) 2012, pp. 1423?1435

%% 
clear;
close all;
clc;

%% Load Images
% NOTE - requires nifti extraction tools

filename1 = ''; % MRI 
filename2 = ''; % PLI 

% read in MRI image (and average MGE)

vol1 = load_nii(filename1);
vol1=mean(vol1.img,4);
im1_orig=reshape(vol1(96+1-71,:,:),128,256);
inputImage = imresize(im1_orig(8:117,37:223)/530,[220,355]);
%inputImage = padarray(inputImage,[20,20],'both');

clear vol1; clear im1_orig;
% read in PLI image and preprocess
im2_orig = imread(filename2);
fixedImage = padarray(imresize(single(flipud(im2_orig)),[210,355]),[5,0],255);
fixedImage(fixedImage <= 0) = 255; %get rid off patched image edges.
fixedImage(fixedImage > 255) = 255;
fixedImage = straightLineRemover(fixedImage,25,255); % vertical direction
fixedImage = straightLineRemover(fixedImage',25,255)'; %NOTE transposes - horizontal direction
fixedImage = medfilt2(fixedImage,[3,3]);

clear im2_orig;

%% Stage 1 - MIND (standard asymmetric or with locally adaptive weighting)
[u1,v1,movingImage]=deformableReg2Dmind_asym(fixedImage,inputImage,.3);
%[u1,v1,movingImage]=deformableReg2Dmind_adaptive(fixedImage,inputImage,.3);

%% Stage 2 - BBR

DeltaIn = 3; %projection distance
DeltaOut = 3;

[boundaryPoints, normals] = boundaryNormal(movingImage,max(DeltaIn,DeltaOut)); %NOTE: normals are unit vectors with direction
[boundaryPoints,normals,endpoints,startpoints] = sortBoundaries(boundaryPoints,normals); %NOTE: add third input to graphically display boundaries

fixedPoints = [startpoints;endpoints];

% solve
[u,v] = boundaryCostNonRigid(boundaryPoints,normals,fixedImage,fixedPoints,DeltaIn,DeltaOut,0.3);

% Extrapolation of vector field and interpolation of image
[outputImage,u2,v2] = weightedInterpolation(movingImage,u,v,boundaryPoints,1);


%% Visualization 

figure;
subplot(2,2,1); imagesc(inputImage); title('input MRI'); colormap(gray);
subplot(2,2,2); imagesc(fixedImage); title('fixed PLI'); colormap(gray);
subplot(2,2,3); imagesc(movingImage); title('Stage 1 - MIND'); colormap(gray);
subplot(2,2,4); imagesc(outputImage); title('Stage 2 - MIND + BBR'); colormap(gray);


