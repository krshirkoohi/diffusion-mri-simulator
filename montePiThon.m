n = 1000;
r = 1; % radius
x = rand(n,1)*2*r-r;
y = rand(n,1)*2*r-r;
mask = x.^2 + y.^2 < r^2;
figure;
ax = axes();
hold(ax);
scatter(x(mask),y(mask),50, 'r.')
scatter(x(~mask),y(~mask),50, 'b.')
axis equal
ax.XLim = [-1 1];
ax.YLim = [-1 1];