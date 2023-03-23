% % Matched Filter
function [gx,gy]=matchedfilter(IM,sigma)
epsilon=1e-2;
halfsize = ceil(sigma*sqrt(-2*log(sqrt(2*pi)*sigma*epsilon)));
size= 2*halfsize+1;
t=3;
p = normcdf([-3*sigma 3*sigma]);
s = (diff(p)/(2*t*sigma));
%generate a 2-D Gaussian kernel along x direction
for i=1:size
for j=1:size
 u=[i-halfsize-1 j-halfsize-1];
hx(i,j)=gauss(u(1),sigma)- s;
end
end
hx=hx/sqrt(sum(sum(abs(hx).*abs(hx))));
%generate a 2-D Gaussian kernel along y direction
hy=hx';
%2-D filtering
gx=imfilter(IM,hx,'replicate','conv');
gy=imfilter(IM,hy,'replicate','conv');

function y = gauss(x,sigma)
%Gaussian
y = exp(-x^2/(2*sigma^2)) / (sigma*sqrt(2*pi));

function y = dgauss(x,sigma)
%first order derivative of Gaussian
y = -x * gauss(x,sigma) / sigma^2;