%sort boundary points
function [sortedPoints,sortedNormals,endpoint,startpoint] = sortBoundaries(points,normals,flag)

if nargin < 3
    %no plotting of boundaries
    fprintf('Display boundaries switched off \n');
    flag = false;
else
    fprintf('Display boundaries switched on \n');
    flag = true;
end

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
order = [];
endpoint = [];

if flag == true
    
    figure;
    plot(points(:,1),points(:,2),'.');
    hold on;
end
flagempty = 0;
count = 0;
idx = 1;
for i = 1:N

    dist(:,idx) = Inf;
    [val, idx] = min(dist(idx,:));
    if val > tolerance
        if length(order) == 0
            flagempty = 1;
        end
        d(order) = 0;
        idx = find(d > 0,1,'first'); %NOTE: this line may be unnecessary
        endpoint = [endpoint, count];
       fprintf('Count: %d \n',count);
    end
    order(i) = idx;
    count = count+1;
    if flag == true
        plot(points(idx,1),points(idx,2),'.r');
        drawnow;
    end
end
endpoint = [endpoint(:);count];
if flagempty == 1
    endpoint(1) = [];
end
%chainLengths = chainLengths(:);
%endpoint = cumsum(chainLengths);
startpoint = [1;endpoint(1:end-1)+1];
%startpoint = startpoint(:);

%plot(points(startpoint,1),points(startpoint,2),'.g','MarkerSize',20);
%order(end) = [];
order = order(:);


sortedPoints = points(order,:);
sortedNormals = normals(order,:);

if flag == true
    plot(sortedPoints(endpoint,1),sortedPoints(endpoint,2),'.b','MarkerSize',20);
    plot(sortedPoints(startpoint,1),sortedPoints(startpoint,2),'.g','MarkerSize',20);
end