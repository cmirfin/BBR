%sort boundary points
function [sortedPoints,sortedNormals,endpoint,chainLengths,startpoint] = sortBoundaries(points,normals)

tolerance = sqrt(2);
%find start point of chain
d = sqrt(points(:,1).^2 + points(:,2).^2); %Euclidian distance from origin
[d, ind] = sort(d,'ascend');
points = points(ind,:); %sort points based on distance from (0,0)
normals = normals(ind,:);


dist = pdist2(points,points); %calculates Euclidian distance between points

N = size(points,1);
%order = NaN(1,N);
%order(1) = 1; % first point is first row in data matrix

endpoint = [];
chainLengths = [];
startpoint = 1;
figure;
plot(points(:,1),points(:,2),'.');
hold on;
count = 0;
idx = 1;
for i = 1:N

    dist(:,idx) = Inf;
    previous_idx = idx;
    [val, idx] = min(dist(idx,:));
    if val > tolerance
%        endpoint = [endpoint, previous_idx];
        d(order) = 0;
        idx = find(d > 0,1,'first'); %NOTE: this line may be unnecessary
%        startpoint = [startpoint,idx];
%        chainLengths = [chainLengths, count];
        endpoint = [endpoint, count];
       fprintf('Count: %d \n',count);
    end
    order(i) = idx;
    count = count+1;
   plot(points(idx,1),points(idx,2),'.r');
    drawnow;
end
endpoint = [endpoint(:);count];
chainLengths = chainLengths(:);
%endpoint = cumsum(chainLengths);
startpoint = [1;endpoint(1:end-1)+1];
%startpoint = startpoint(:);

%plot(points(startpoint,1),points(startpoint,2),'.g','MarkerSize',20);
%order(end) = [];
order = order(:);


sortedPoints = points(order,:);
sortedNormals = normals(order,:);

plot(sortedPoints(endpoint,1),sortedPoints(endpoint,2),'.b','MarkerSize',20);
plot(sortedPoints(startpoint,1),sortedPoints(startpoint,2),'.g','MarkerSize',20);