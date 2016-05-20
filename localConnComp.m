function array = localConnComp(boundaryPoints, backgroundSegmentation)

numPoints = size(boundaryPoints,1);
R = 9;
[m,n] = size(backgroundSegmentation);
count = 0;
centre = ((R-1)/2) + 1;
for i = 1:numPoints
    
    x0 = boundaryPoints(i,1); %neighbourhood centre point
    y0 = boundaryPoints(i,2);
    xlb = x0 - (R-1)/2;
    xub = x0 + (R-1)/2;
    ylb = y0 - (R-1)/2;
    yub = y0 + (R-1)/2;
    
    xlb(xlb < 1) = 1;
    xub(xub > m) = m;
    ylb(ylb < 1) = 1;
    yub(yub > n) = n;
    
    neighbourhood = backgroundSegmentation(xlb:xub, ylb:yub);
    %should insert check here that we are not going outside of array
    
    %check for other boundary points in neighbourhood
    [idx,~] = find(boundaryPoints(:,1) >= xlb & boundaryPoints(:,1) < xub & boundaryPoints(:,2) >= ylb & boundaryPoints(:,2) < yub);
    
    L = length(idx); %NOTE: minumum will be L = 1
    
    %convert to local indices.
    
    
    if L == 1
        %if there are no other boundary points in local neighbourhood
        return
    else
        %there are > 1 boundary point in neighbourhood.
        %run connected component analysis
        
        labels = bwlabel(neighbourhood); %connected component analysis
        specialLabel = labels(centre,centre); %label of our centre boundary point
        
        for j = 1:L
            
            localIndX = boundaryPoints(idx(j),1) - x0 + centre;
            localIndY = boundaryPoints(idx(j),2) - y0 + centre;
            
            if labels(localIndX,localIndY) == specialLabel
                %boundary point in locally connected. Add regularisation
                %fprintf('Connected boundary point! \n');
                count = count + 1;
            end
        end
        
        
    end
    
end

end

    
    