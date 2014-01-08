% draws a gaussian into a figure that must already be open. this allows to draw several gaussian into a figure
% NOTE: this function takes standard deviation not FWHM!
function drawGaussian(p, color)
% % adding 0.5 is required to adjust to the coordinate system of imagesc
% mx = p(1) + 0.5;
% my = p(2) + 0.5;
mx = p(1);
my = p(2);
sx = p(3);
sy = p(4);
r = p(5);
U = [[sx*sx r*sx*sy]; [r*sx*sy sy*sy]];
[V, D] = eig(U);

plot(mx, my, [color '.']);

t = 0:0.1:2*pi;
t = [t 0];
circ = [cos(t); sin(t)];

circ = V * (sqrt(D) * circ);
plot(circ(1,:) + mx, circ(2,:) + my, [color '-']);
end
