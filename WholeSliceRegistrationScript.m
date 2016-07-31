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

%% addpaths
% NOTE - requires nifti extraction tools

filename1 = 'FILEPATH/x_CC/MRI/3DGRE_to_PLI.nii.gz'; % MRI
filename2 = 'FILEPATH/x_CC/PLI/Inclination.tif'; % PLI

%% read in MRI image (and average MGE)

vol1 = load_nii(filename1);
vol1=mean(vol1.img,4);
inputImage = squeeze(vol1(:,259,:));
inputImage = inputImage(192:316,201:237);
inputImage = imresize(inputImage,[89,36]);
clear vol1;

%% read in PLI image and preprocess
fixedImage = single((imread(filename2)));
fixedImage(fixedImage <= 0) = 255; %get rid off patched image edges.
fixedImage(fixedImage > 255) = 255;
fixedImage = straightLineRemover(fixedImage,25,255); %vertical direction
fixedImage = straightLineRemover(fixedImage',25,255)'; %NOTE transposes - horizontal direction
fixedImage = imresize(fixedImage,0.1); % downsample
fixedImage(fixedImage > 250) = 240;
fixedImage = rot90(fixedImage,-1); % rotates image
fixedImage = ArtificialEdgeCrop(fixedImage,230,20); % crops edges that appear to be artificial cuts
fixedImage = imresize(fixedImage,0.1); % downsample

%% perform registration using locally adaptive MIND

%[u1,v1,movingImage]=deformableReg2Dmind_adaptive(fixedImage,inputImage,.3);
[u1,v1,movingImage]=deformableReg2Dmind_asym(fixedImage,inputImage,.3);


%% find boundary and projection points in moving image (im1)
DeltaIn = 3; %projection distance
DeltaOut = 3;

[boundaryPoints, normals] = boundaryNormal(movingImage,max(DeltaIn,DeltaOut)); %NOTE: normals are unit vectors with direction
[boundaryPoints,normals,endpoints,startpoints] = sortBoundaries(boundaryPoints,normals); %NOTE: add third input to graphically display boundaries
 %remove start and final point

fixedPoints = [startpoints;endpoints];

%% solve

[u,v] = boundaryCostNonRigid(boundaryPoints,normals,fixedImage,fixedPoints,DeltaIn,DeltaOut,0.3);

%% Extrapolation of vector field and interpolation of image

outputImage = weightedInterpolation(movingImage,u,v,boundaryPoints,1); %NOTE: specify extrapolaton scheme

%% Visualization 

figure;
subplot(2,2,1); imagesc(inputImage); title('input MRI'); colormap(gray);
subplot(2,2,2); imagesc(fixedImage); title('fixed PLI'); colormap(gray);
subplot(2,2,3); imagesc(movingImage); title('Stage 1 - MIND'); colormap(gray);
subplot(2,2,4); imagesc(outputImage); title('Stage 2 - MIND + BBR'); colormap(gray);
