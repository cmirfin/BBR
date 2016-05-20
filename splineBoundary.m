function splineBoundary(x,y,m)

% t = linspace(0,1,numel(x) ); % define the parameter t
% fitX = fit( t, x, 'cubicinterp'); % fit x as a function of parameter t
% fitY = fit( t, y, 'cubicinterp'); % fit y as a function of parameter t
% plot( fitX, fitY, 'r-' ); % plot the parametric curve


n = length (x);
s = zeros (n,1);
for i = 2:n
     s(i) = s(i-1) + sqrt((x(i)-x(i-1))^2+(y(i)-y(i-1))^2);
end
sout = linspace(0, s(n), m)';
xout = spline(s, x, sout);
yout = spline(s, y, sout);
figure;
plot (x, y, 'ro', xout, yout, 'b');
end
