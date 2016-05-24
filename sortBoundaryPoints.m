%sort boundary points
function [result,position] = sortBoundaryPoints(points)

%data = [2,2 ; 2,3 ; 1,2 ; 1,3 ; 2,1 ; 1,1 ; 3,2 ; 3,3 ; 3 ,1];
dist = pdist2(points,points); %calculates Euclidian distance between points

% N = size(points,1);
% result = NaN(1,N);
% result(1) = 1; % first point is first row in data matrix
% 
% for ii=2:N
%     dist(:,result(ii-1)) = Inf;
%     [~, closest_idx] = min(dist(result(ii-1),:));
% 
%     result(ii) = closest_idx;
% 
%    
% end
% 
% sortedPoints = points(result,:);

N = size(points,1);
result = NaN(1,N);
result(1) = 1; % first point is first row in data matrix

position = [];

for ii=2:N
    dist(:,result(ii-1)) = Inf; %ensures that we do not reconnect to previous point
    [val, closest_idx] = min(dist(result(ii-1),:)); %find next nearest neighbour
    
    result(ii) = closest_idx;
    if val > 2
       %flag as new line segment if greater than 5 pixels apart
        position = [position, ii];
        
    end
           
end

%sortedPoints = points(result,:); %sort points in line order
