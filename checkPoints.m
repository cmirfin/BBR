function [xCheck,yCheck] = checkPoints(X,Y,m,n)

    dig = 4;
    %%function to check boundaries and if points overlap
    N = length(X)/2;
    
    X1 = X(1:N); 
    X2 = X(N+1:end);
    Y1 = Y(1:N);
    Y2 = Y(N+1:end);
    
    zIn = [round(X1,dig),round(Y1,dig)];
    [z1new, idxrow] = unique(zIn,'rows','stable');
    
    zOut = [round(X2,dig),round(Y2,dig)];
    z2new = zOut(idxrow,:);

    [z2new, idxrow] = unique(z2new,'rows','stable');
    z1new = z1new(idxrow,:);
%     
%     l1 = ismember(z1new,z2new);
%     
%     z1 = z1new(l1 == 0,:);
%     z2 = z2new(l1 == 0,:);
%     
%     zfinal = [z1;z2];
%     
%     xCheck = zfinal(:,1);
%     yCheck = zfinal(:,2);
    z = [z1new;z2new];
    N2 = size(z,1); %reduced length
    [~, idxrow] = unique(z,'rows','stable');
    duplicateIdx = setdiff(1:N2,idxrow);
    
    idx = duplicateIdx - (N2)/2;
    
    z1new(idx,:) = [];
    z2new(idx,:) = [];
    
    xCheck = [z1new(:,1);z2new(:,1)];
    yCheck = [z1new(:,2);z2new(:,2)];
    
    %remove points outside of fixed image domain
    

    