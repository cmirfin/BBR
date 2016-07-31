% Christopher J. Mirfin 
% Sir Peter Mansfield Imaging Centre, University of Nottingham
% christopher.mirfin@dtc.ox.ac.uk
% 23/07/2016

function croppedImg = ArtificialEdgeCrop(img,threshold, edgeTolerance)
%ARTIFICIALEDGECROP Crops edges of input image by first binarizing image
%with user supplied threshold (default 230). Straight line edges are then
%found, and if within edgeTolerance of image, they will be cropped. 

if nargin < 2
    threshold = 230;
    edgeTolerance = 100;
elseif nargin < 3
    edgeTolerance = 100;
end

[m,n] = size(img);
%binarize image with threshold
BW = img;
BW(BW<threshold) = 0;
BW(BW>0) = 1;

BW = medfilt2(BW,[9,9]); %median filter to smooth salt/pepper noise

BW = edge(BW,'sobel'); %find all edges
points = [];

[H,T,R] = hough(BW);
P  = houghpeaks(H,5,'threshold',ceil(0.3*max(H(:))));
lines = houghlines(BW,T,R,P,'FillGap',50,'MinLength',5);


figure, imagesc(img), hold on

for k = 1:length(lines)
   xy = [lines(k).point1; lines(k).point2];
   plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','black');

   % Plot beginnings and ends of lines
   plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
   plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');
   
   points = [points;xy];

end
fprintf('Please inspect quality of artificial edge detection \n');
reply = input('Do you wish to continue to crop this image Y/N [Y]: ', 's');

if reply == 'Y';

    val1 = max(points(find(points(:,2)<edgeTolerance),2));
    val2 = max(points(find(points(:,1)<edgeTolerance),1));
    val3 = min(points(find(points(:,2)> m- edgeTolerance),2));
    val4 = min(points(find(points(:,1)> n- edgeTolerance),1));
    
    if isempty(val1) == 1
        val1 = 1;
    end
    if isempty(val2) == 1
        val2 = 1;
    end
    if isempty(val3) == 1
        val3 = m;
    end
    if isempty(val4) == 1
        val4 = n;
    end
    fprintf('Image cropped \n');
    croppedImg = img(val1:val3,val2:val4);
    figure;
    imagesc(croppedImg);
else
    fprintf('Image not cropped \n');
    croppedImg = img;
    figure;
    imagesc(croppedImg);
end
end

