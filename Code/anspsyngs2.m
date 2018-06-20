function [ysm,y] = anspsyngs2(beta,n,theta,a)
%function [ysm,y] = anspsyngs2(beta,n,theta,a)
% ANSPSYNGS2 simulates an anisotropic 2-d random gaussian field,
% by spectral synthesis, and by rotation and scaling.
% beta is the correlation length
% ysm is output matrix of size 2*n+1 by 2*n+1
% theta is the angle from x-axis of the major axis of the correlation function
% a is the anisotropic factor > 1, can take only integer values
%
% See also SPSYNGS2 and SPSYNFRAC2

%Written by M. Sengupta and T. Mukerji, 1996

theta=90-theta;
y=spsyngs2(beta,2*n);
y=imrotate(y,theta);
[r,c]=size(y);
y=resample(y(:),a,1);
y=reshape(y,r*a,c);
y=imrotate(y,360-theta);
[r1,c1]=size(y);
ysm=y(floor(r1/2)-n:floor(r1/2)+n+1,floor(c1/2)-n:floor(c1/2)+n+1);


