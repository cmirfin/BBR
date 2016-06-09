function outputImg = straightLineRemover(img,minLength,fillValue)
%STRAIGHTLINEREMOVER removes vertical lines greater than minLength from
%image input. Fills with fillValue.

[m,n] = size(img);
img = img(:);
outputImg = img;
img = round(img);

dif = [1;diff(img)];
%img(dif == 0 & img<255) = NaN;
BW = zeros(size(img));
BW(dif == 0 & img<max(img)) = 1;
dif2 = [0; diff(BW)];
indplus = find(dif2==1);
indneg = find(dif2== -1);

lns = indneg - indplus;
idx = find(lns > minLength);
for i=1:numel(idx)
   outputImg(indplus(idx(i)):indneg(idx(i))) = fillValue; 
end
    
outputImg = reshape(outputImg,m,n);