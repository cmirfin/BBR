%sort boundary points
function [newOrder,newPosition] = sortManyBoundaryPoints(points)

%data = [2,2 ; 2,3 ; 1,2 ; 1,3 ; 2,1 ; 1,1 ; 3,2 ; 3,3 ; 3 ,1];
dist = pdist2(points,points); %calculates Euclidian distance between points

N = size(points,1);
order = NaN(1,N);
order(1) = 1; % first point is first row in data matrix

position = [];

for ii=2:N
    dist(:,order(ii-1)) = Inf; %ensures that we do not reconnect to previous point
    [val, closest_idx] = min(dist(order(ii-1),:)); %find next nearest neighbour
    
    order(ii) = closest_idx;
    if val > sqrt(2);
       %flag as new line segment if greater than 2 pixels apart
        position = [position, ii];
        
    end
           
end

%sortedPoints = points(order,:); %sort points in line order

newOrder = [];
newPosition = [];

boundaryLengths = diff(position);
for k = 1:numel(boundaryLengths)
    if boundaryLengths(k) > 10
        %only include boundary segments that have above certain length
        %newPoints = [newPoints; sortedPoints(position(k):position(k+1)-1,:)];
        
        newOrder = [newOrder, order(position(k):position(k+1)-1)];
        newPosition = [newPosition; squeeze(length(newOrder))];
    end
end

%remove end condition:
newPosition(end) = [];