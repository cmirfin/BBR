function [CS, norms] = splineBoundary(points)

% 
% n = length (x);
% s = zeros (n,1);
% for i = 2:n
%      s(i) = s(i-1) + sqrt((x(i)-x(i-1))^2+(y(i)-y(i-1))^2);
% end
% sout = linspace(0, s(n), m)';
% xout = spline(s, x, sout);
% yout = spline(s, y, sout);
% figure;
% plot (x, y, 'ro', xout, yout, 'b');
% end

CS = cscvn(points'); %NOTE transpose on points
der = fnder(CS);
dydx = abs(ppval(der,CS.breaks)'); %NOTE transpose

dy = dydx(:,1);
dx = dydx(:,2);

l=sqrt(dx.^2+dy.^2);
dx = dx./l; 
dy = dy./l;

norms = [dx(:),dy(:)];