% MOVABLEROI displays an image with a predefined rectangle ROI of size (m,n). 
% The ROI can be dragged around the screen with pixel locations displayed.
% Double-clicking crops the image to the given ROI.

% Christopher J. Mirfin 
% Sir Peter Mansfield Imaging Centre, University of Nottingham
% christopher.mirfin@dtc.ox.ac.uk
% 12/08/2016

function croppedImage = movableROI(Image,m,n)

S = [1 1 n-1 m-1]; % initialise box in top-left hand corner
figure;
imagesc(Image);
title('Please drag box over ROI - double click to finish');
h = imrect(gca, S);

addNewPositionCallback(h,@(p) title(mat2str(p,3)));
fcn = makeConstrainToRectFcn('imrect',get(gca,'XLim'),get(gca,'YLim'));
setPositionConstraintFcn(h,fcn)

position = wait(h);
croppedImage = imcrop(Image,position);
imagesc(croppedImage); title('Cropped MRI image');   % the output image of your ROI

end
