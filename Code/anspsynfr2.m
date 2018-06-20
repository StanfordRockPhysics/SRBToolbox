function [ysm,y] = anspsynfr2(beta,n,theta,a)
%function [ysm,y] = anspsynfr2(beta,n,theta,a)
% ANSPSYNFR2 simulates an anisotropic 2-d random fractal field,
% by spectral synthesis, and by rotation and scaling.
% beta is the spectral exponent, the power spectrum falls off as frequency^BETA.
% ysm is output matrix of size 2*n+1 by 2*n+1
% theta is the angle from x-axis of the major axis of the correlation function
% a is the anisotropic factor > 1, can take only integer values
%
% See also SPSYNGS2 and SPSYNFRAC2

%Written by M. Sengupta, T. Mukerji, 1996

theta=90-theta;
y=spsynfrac2(beta,2*n);
y=imrotate(y,theta); 
[r,c]=size(y);
y=resample(y(:),a,1);
y=reshape(y,r*a,c);
y=imrotate(y,360-theta);
[r1,c1]=size(y);
ysm=y(floor(r1/2)-n:floor(r1/2)+n+1,floor(c1/2)-n:floor(c1/2)+n+1);


