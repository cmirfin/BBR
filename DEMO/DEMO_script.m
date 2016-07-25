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
addpath('../') 
%% Load Images (NOTE - these images have already been preprocessed from raw images)


filename1 = 'inputImage.tif';
filename2 = 'fixedImage.tif';

inputImage = single(imread(filename1));
fixedImage = single(imread(filename2));


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


