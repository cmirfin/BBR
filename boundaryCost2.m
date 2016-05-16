function [J, grad] = boundaryCost2(theta,rIn,rOut,fixedImage,Fx,Fy)

%parameters
M = 0.005;

[m,n] = size(fixedImage);

%transform rIn, rOut
rPrimeIn = transformPoints(rIn,theta,m,n);
rPrimeOut = transformPoints(rOut,theta,m,n);

xPrime = [rPrimeIn(:,1); rPrimeOut(:,1)];
yPrime = [rPrimeIn(:,2); rPrimeOut(:,2)];
[xPrime,yPrime] = checkPoints(xPrime,yPrime); %rmove overlapping points

N = length(xPrime)/2; %number of boundary points 
[X,Y] = meshgrid(1:n,1:m); %grid points of fixed image

%interpolate values
fixedInterp = interp2(X,Y,fixedImage,yPrime,xPrime,'linear',1);
valueIn = fixedInterp(1:N);
valueOut = fixedInterp(N+1:end);

FxInterp = interp2(X,Y,Fx,yPrime,xPrime,'linear',0); %x-gradient
FyInterp = interp2(X,Y,Fy,yPrime,xPrime,'linear',0); %y-gradient

Q = (M*(100*(valueIn-valueOut)./(0.5*(valueIn+valueOut))));

J = 1 + tanh(Q);
J = double(sum(J))/N;
fprintf('Cost: %f \n',J);

%transform derivatives with Jacobian dT/dtheta
derivIn = transformDerivatives(theta,xPrime(1:N),yPrime(1:N),FxInterp(1:N),FyInterp(1:N));
derivOut = transformDerivatives(theta,xPrime(N+1:end),yPrime(N+1,end),FxInterp(N+1:end),FyInterp(N+1:end));

grad = zeros(1,length(theta));
for i = 1:length(theta)
    grad(i) = sum((sech(Q).^2).*(valueOut.*derivIn(:,i) - valueIn.*derivOut(:,i))./((valueIn + valueOut).^2));
end

grad = double((400*M/N)*grad);

end

function rPrime = transformPoints(r,t,m,n)

numParameters = length(t);
switch numParameters
    case 1
        A = [1, 0, t; 0, 1, 0; 0,0,1]';
    case 2
        A = [1, 0, t(1); 0, 1, t(2); 0,0,1]';
    otherwise
        A = [1+t(1), t(3), t(5); t(2), 1+t(4), t(6); 0,0,1]';
end

tform = affine2d(A);
rPrime = transformPointsForward(tform,r);

end

function Jacobian = transformDerivatives(theta,x,y,dx,dy)
len = length(x);
dx = dx(:);
dy = dy(:);
numParameters = length(theta);
switch numParameters
    case 1
        %Jacobian = [1;1];
        Jacobian = dx+dy; %NOTE: this should be dx or dy.
    case 2
        %Jacobian = [1 0; 0 1];
        Jacobian = [dx,dy];
    otherwise
        Jacobian = zeros(len,6);
        Jacobian(:,1) = x.*dx;
        Jacobian(:,2) = x.*dy;
        Jacobian(:,3) = y.*dx;
        Jacobian(:,4) = y.*dy;
        Jacobian(:,5) = dx;
        Jacobian(:,6) = dy;
end

end



