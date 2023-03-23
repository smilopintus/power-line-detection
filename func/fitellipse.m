%=========================================================================================================================================================
% Reference: http://homepages.inf.ed.ac.uk/rbf/CVonline/LOCAL_COPIES/FITZGIBBON/ELLIPSE/
%  Copyright (c) 1999, 2005, Andrew Fitzgibbon, Maurizio Pilu, Bob Fisher
%
% Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING

function [uCentre, vCentre, Ru, Rv, thetarad] = fitellipse(X,Y)
% FITELLIPSE Least-squares fit of ellipse to 2D points.
% A = FITELLIPSE(X,Y) returns the parameters of the best-fit
% ellipse to 2D points (X,Y).
% The returned vector A contains the center, radii, and orientation
% of the ellipse, stored as (Cx, Cy, Rx, Ry, theta_radians)

% normalize data
mx = mean(X);
my = mean(Y);
sx = (max(X)-min(X))/2;
sy = (max(Y)-min(Y))/2;
x = (X-mx)/sx;
y = (Y-my)/sy;
% Force to column vectors
x = x(:);
y = y(:);
% Build design matrix
D = [ x.*x x.*y y.*y x y ones(size(x)) ];
% Build scatter matrix
S = D'*D;
% Build 6x6 constraint matrix
C(6,6) = 0; C(1,3) = -2; C(2,2) = 1; C(3,1) = -2;
% Solve eigensystem
[gevec, geval] = eig(S,C);
% Find the negative eigenvalue
I = find(real(diag(geval)) < 1e-8 & ~isinf(diag(geval)));
% Extract eigenvector corresponding to negative eigenvalue
A = real(gevec(:,I));
% unnormalize
par = [
A(1)*sy*sy, ...
A(2)*sx*sy, ...
A(3)*sx*sx, ...
-2*A(1)*sy*sy*mx - A(2)*sx*sy*my + A(4)*sx*sy*sy, ...
-A(2)*sx*sy*mx - 2*A(3)*sx*sx*my + A(5)*sx*sx*sy, ...
A(1)*sy*sy*mx*mx + A(2)*sx*sy*mx*my + A(3)*sx*sx*my*my ...
- A(4)*sx*sy*sy*mx - A(5)*sx*sx*sy*my ...
+ A(6)*sx*sx*sy*sy ...
]';
% Convert to geometric radii, and centers
thetarad = 0.5*atan2(par(2),par(1) - par(3));
cost = cos(thetarad);
sint = sin(thetarad);
sin_squared = sint.*sint;
cos_squared = cost.*cost;
cos_sin = sint .* cost;
Ao = par(6);
Au = par(4) .* cost + par(5) .* sint;
Av = - par(4) .* sint + par(5) .* cost;
Auu = par(1) .* cos_squared + par(3) .* sin_squared + par(2) .* cos_sin;
Avv = par(1) .* sin_squared + par(3) .* cos_squared - par(2) .* cos_sin;
% ROTATED = [Ao Au Av Auu Avv]
tuCentre = - Au./(2.*Auu);
tvCentre = - Av./(2.*Avv);
wCentre = Ao - Auu.*tuCentre.*tuCentre - Avv.*tvCentre.*tvCentre;
uCentre = tuCentre .* cost - tvCentre .* sint;
vCentre = tuCentre .* sint + tvCentre .* cost;
Ru = -wCentre./Auu;
Rv = -wCentre./Avv;
Ru = sqrt(abs(Ru)).*sign(Ru);
Rv = sqrt(abs(Rv)).*sign(Rv);

end

