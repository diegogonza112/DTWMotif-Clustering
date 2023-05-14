t = 0:.1:2*pi;
v = 0:.1:2.2*pi;
x0 = sin(t) + rand(size(t)) * .1;
% x1 = sin(.9*t) + rand(size(t)) * .1;
x1 = sin(.9*v) + rand(size(v)) * .1;
w=2;

figure
% plot(t, x0, t, x1);
plot(t, x0, v, x1);
hold on

if (exist('w','var'))
    [cost, path] = warping_path_w(x0,x1,w);
    title(sprintf('Cost: %f; Warp limit: %d',cost(end),w))
else
    [cost, path] = warping_path_w(x0,x1);
    title(sprintf('Cost: %f; no warp limit',cost(end)))
end
% plot(v, x1(path));
plot(v(1:length(path)),x1(path))
legend({'x_0', 'x_1', 'x_1 warped to x_0'});
clear('w');