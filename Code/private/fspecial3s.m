function h = fspecial3(type,P1,P2)
%function h = fspecial3(type,P1,P2)
%FSPECIAL3 Create predefined filters.
%   H = FSPECIAL3(TYPE) creates a three-dimensional filter H of the
%   specified type. (FSPECIAL3 returns H as a computational
%   molecule, which is the appropriate form to use with FILTER.)
%   TYPE is a string having one of these values:
%
%        'gaussian'  for a Gaussian lowpass filter
%        'average'   for an averaging filter
%
%   Depending on TYPE, FSPECIAL3 can take additional parameters
%   which you can supply.  These parameters all have default
%   values. 
%
%   H = FSPECIAL3('gaussian',N,SIGMA) returns a rotationally symmetric
%   Gaussian lowpass filter with standard deviation SIGMA (in pixels).
%   N specifies the number of rows in H. The correlation coefficient of
%   each axis is zero. 
%   Either SIGMA and N can be a vector with length 3, each term 
%   corresponding to the standard deviation/filter length in 
%   1st/2nd/3rd axis. 
%   If SIGMA is a 3*3 matrix, SIGMA will be considered as the covariance
%   matrix of the output filter. (Remember to normalize the covariance
%   matrix to the grid in pixels.) Covariance matrix may be computed with
%   cov.m
%   If you do not specify the parameters, FSPECIAL3 uses the default 
%   values of N = 3 and 0.5 for SIGMA.
%
%   H = FSPECIAL3('average',N) returns an averaging filter. N specifies
%   the number of rows in H. N can be a vector with length 3, each term 
%   corresponding to the filter length in each axis.
%   If you do not specify N, FSPECIAL3 uses the default value of N = 3.
%
%   See also CONV, FILTER, FSPECIAL, FSPECIAL1.

%   Clay M. Thompson 11-17-92
%   Copyright 1993-1998 The MathWorks, Inc.  All Rights Reserved.
%   $Revision: 5.9 $  $Date: 1997/11/24 15:34:46 $

%   Modified by Isao Takahashi 11/23/99

if nargin==0, error('Not enough input arguments.'); end
type = [type,'  '];
code = lower(type(1:2));
if nargin<2, P1 = 3; end
if nargin<3, P2 = .5; end

if length(P1)==3,
  siz = P1;
elseif length(P1)==1,
  siz = P1*ones(1,3);
else
  error('The second parameter must be a scalar or 3-term vector.');
end
if length(P2)==3,
  std = P2;
elseif length(P2)==1 ,
  std = P2*ones(1,3);
else
  error('The third parameter must be a scalar or 3-term vector.');
end

if all(code=='ga'), % Gaussian filter
  [x,y,z] = meshgrid(-(siz(1)-1)/2:(siz(1)-1)/2,-(siz(2)-1)/2:(siz(2)-1)/2,-(siz(3)-1)/2:(siz(3)-1)/2);
if min(size(P2))==1
  h = exp(-(x.*x/(2*std(1)*std(1)) + y.*y/(2*std(2)*std(2)) + z.*z/(2*std(3)*std(3))));
  h = h/sum(sum(sum(h)));

elseif size(P2)==[3 3];
C=P2;
for i=1:siz(1)
for j=1:siz(2)
for k=1:siz(3)
X=[x(i,j,k) y(i,j,k) z(i,j,k)];
  h(i,j,k) = exp(-.5*(X*inv(C)*X'));
end
end
end
  h = h/sum(sum(sum(h)));
end

elseif all(code=='av'), % Smoothing filter
  h = ones(siz)/prod(siz);

else
  error('Unknown filter type.');

end
