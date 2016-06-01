%sort boundary points
function [order,endpoint,chainLengths] = sortBoundaries(points)

tolerance = sqrt(2);
%find start point of chain
d = sqrt(points(:,1).^2 + points(:,2).^2); %Euclidian distance from origin
[~, ind] = sort(d,'ascend');
points = points(ind,:); %sort points based on distance from (0,0)

dist = pdist2(points,points); %calculates Euclidian distance between points

N = size(points,1);
%order = NaN(1,N);
%order(1) = 1; % first point is first row in data matrix

endpoint = [];
chainLengths = [];
figure;
plot(points(:,1),points(:,2),'.');
hold on;
count = 0;
idx = 1;
for i = 1:N-1

    dist(:,idx) = Inf;
    previous_idx = idx;
    [val, idx] = min(dist(idx,:));
    if val > tolerance
        endpoint = [endpoint, previous_idx];
        d(order) = 0;
        idx = find(d > 0,1,'first');
        chainLengths = [chainLengths, count];
        count = 0;
    end
    order(i) = idx;
    count = count+1;
    plot(points(idx,1),points(idx,2),'.r');
    drawnow;
end
endpoint = endpoint(:);
plot(points(endpoint,1),points(endpoint,2),'.g','MarkerSize',10);
order(end) = [];
order = [1;order(:)];
chainLengths = chainLengths(:);


%delete boundaries less than minLength

% sortedPoints = points(order,:); %sort points in line order

% newOrder = [];
% newPosition = [];
% 
% % boundaryLengths = diff(endpoint);
% for k = 1:numel(chainLengths)-1
%     if chainLengths(k) > minLength
%         %only include boundary segments that have above certain length
%         %newPoints = [newPoints; sortedPoints(position(k):position(k+1)-1,:)];
%         
%         newOrder = [newOrder; order(endpoint(k)+1:endpoint(k)+chainLengths(k))];
%         newPosition = [newPosition; chainLengths(k)];
%     end
% end

%remove end condition:
% newPosition(end) = [];

