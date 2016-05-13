function [J, grad] = boundaryCost2(theta,rIn,rOut,fixedImage,Fx,Fy)

% J = 0;
% grad = zeros(size(theta));

%parameters
M = 0.5;

[m,n] = size(fixedImage);

%transform rIn, rOut
rPrimeIn = transformPoints(rIn,theta,m,n);
rPrimeOut = transformPoints(rOut,theta,m,n);

%check points
PointsX = [rPrimeIn(:,1); rPrimeOut(:,1)];
PointsY = [rPrimeIn(:,2); rPrimeOut(:,2)];
[xPrime,yPrime] = checkPoints(PointsX,PointsY);

N = length(xPrime)/2; %number of boundary points 
[X,Y] = meshgrid(1:n,1:m);

%interpolate values
fixedInterp = interp2(X,Y,fixedImage,xPrime,yPrime,'spline');
valueIn = fixedInterp(1:N);
valueOut = fixedInterp(N+1:end);

FxInterp = interp2(X,Y,Fx,xPrime,yPrime,'spline');
FxInterpIn = FxInterp(1:N);
FxInterpOut = FxInterp(N+1:end);

FyInterp = interp2(X,Y,Fy,xPrime,yPrime,'spline');
FyInterpIn = FyInterp(1:N);
FyInterpOut = FyInterp(N+1:end);

Q = abs(M*(100*(valueIn-valueOut)./(0.5*(valueIn+valueOut))));

J = 1 - tanh(Q);
J = double(sum(J)./N);
fprintf('Cost: %f \n',J);

thetaJacobianIn = [1 ;1];
thetaJacobianOut = [1;1];
derivIn = [FxInterpIn, FyInterpIn]*thetaJacobianIn;
derivOut = [FxInterpOut, FyInterpOut]*thetaJacobianOut;

grad = (M*(sech(Q).^2).*2.*((valueOut.*derivIn) - (valueIn.*derivOut))./(valueIn+valueOut).^2);

grad = (200/N)*(double(sum(grad,1)));

end

function rPrime = transformPoints(r,t,m,n)

%Build affine matrix
%A = [1+t(1), t(3), t(5); t(2), 1+t(4), t(6); 0,0,1];
A = [1, 0, t(1); 0, 1, 0; 0,0,1];
num = size(r,1);
xPrime = zeros(num,1);
yPrime = zeros(num,1);

%Apply to each point
for i = 1:num
   vec = A*[r(i,1); r(i,2); 1];
%    
%    %do not accept points outside of fixed image domain
%    if vec(1) > 1 && vec(1) < m && vec(2) > 1 && vec(2) < n
%        xPrime(i) = vec(1);
%        yPrime(i) = vec(2);
%    end
    rPrime(i,1) = vec(1);
    rPrime(i,2) = vec(2);
end

end

function [xCheck,yCheck] = checkPoints(X,Y,m,n)

    %%function to check boundaries and if points overlap
    N = length(X)/2;
    
    X1 = X(1:N); 
    X2 = X(N+1:end);
    Y1 = Y(1:N);
    Y2 = Y(N+1:end);
    
    zIn = [X1,Y1];
    [z1new, idxrow] = unique(zIn,'rows','stable');
    
    zOut = [X2,Y2];
    z2new = zOut(idxrow,:);

    [z2new, idxrow] = unique(z2new,'rows','stable');
    z1new = z1new(idxrow,:);
    
    l1 = ismember(z1new,z2new);
    
    z1 = z1new(l1 == 0,:);
    z2 = z2new(l1 == 0,:);
    
    zfinal = [z1;z2];
    
    xCheck = zfinal(:,1);
    yCheck = zfinal(:,2);
end


