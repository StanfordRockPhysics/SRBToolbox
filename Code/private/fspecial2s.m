function h = fspecial2s(type,P1,P2,P3)
%FSPECIAL Create predefined filters.
%   H = FSPECIAL(TYPE) creates a two-dimensional filter H of the
%   specified type. (FSPECIAL returns H as a computational
%   molecule, which is the appropriate form to use with FILTER2.)
%   TYPE is a string having one of these values:
%
%        'gaussian'  for a Gaussian lowpass filter
%        'sobel'     for a Sobel horizontal edge-emphasizing
%                       filter 
%        'prewitt'   for a Prewitt horizontal edge-emphasizing
%                       filter 
%        'laplacian' for a filter approximating the
%                       two-dimensional Laplacian operator
%        'log'       for a Laplacian of Gaussian filter
%        'average'   for an averaging filter
%        'unsharp'   for an unsharp contrast enhancement filter
%
%   Depending on TYPE, FSPECIAL can take additional parameters
%   which you can supply.  These parameters all have default
%   values. 
%
%   H = FSPECIAL('gaussian',N,SIGMA) returns a rotationally
%   symmetric Gaussian lowpass filter with standard deviation
%   SIGMA (in pixels). N and SIGMA is a 1-by-2 vector specifying 
%   the number of rows and columns in H. (N can also be a scalar, 
%   in which case H is N-by-N.) 
%   SIGMA can be a 2-by-2 matrix, in which case SIGMA is interpreted
%   as the covariance matrix of the filter. (Remember to normalize it
%   to the grid in pixels.) Covariance matrix may be computed with 
%   cov.m
%   If you do not specify the parameters,FSPECIAL uses the default 
%   values of [3 3] for N and 0.5 for SIGMA.
%
%   H = FSPECIAL('sobel') returns this 3-by-3 horizontal edge
%   finding and y-derivative approximation filter:
%
%       [1 2 1;0 0 0;-1 -2 -1].
%
%   To find vertical edges, or for x-derivates, use -h'.
%
%   H = FSPECIAL('prewitt') returns this 3-by-3 horizontal edge
%   finding and y-derivative approximation filter:
%
%       [1 1 1;0 0 0;-1 -1 -1].
%
%   To find vertical edges, or for x-derivates, use -h'.
%
%   H = FSPECIAL('laplacian',ALPHA) returns a 3-by-3 filter
%   approximating the shape of the two-dimensional Laplacian
%   operator. The parameter ALPHA controls the shape of the
%   Laplacian and must be in the range 0.0 to 1.0. FSPECIAL uses
%   the default value of 0.2 if you do not specify ALPHA.
%
%   H = FSPECIAL('log',N,SIGMA) returns a rotationally symmetric
%   Laplacian of Gaussian filter with standard deviation SIGMA
%   (in pixels). N is a 1-by-2 vector specifying the number of
%   rows and columns in H. (N can also be a scalar, in which case
%   H is N-by-N.) If you do not specify the parameters, FSPECIAL
%   uses the default values of [5 5] for N and 0.5 for SIGMA.
%
%   H = FSPECIAL('average',N) returns an averaging filter. N is a
%   1-by-2 vector specifying the number of rows and columns in
%   H. (N can also be a scalar, in which case H is N-by-N.) If
%   you do not specify N, FSPECIAL uses the default value of 
%   [3 3].
%
%   H = FSPECIAL('unsharp',ALPHA) returns a 3-by-3 unsharp
%   contrast enhancement filter. FSPECIAL creates the unsharp
%   filter from the negative of the Laplacian filter with
%   parameter ALPHA. ALPHA controls the shape of the Laplacian
%   and must be in the range 0.0 to 1.0. FSPECIAL uses the
%   default value of 0.2 if you do not specify ALPHA.
%
%   Example
%   -------
%       I = imread('saturn.tif');
%       h = fspecial('unsharp',0.5);
%       I2 = filter2(h,I)/255;
%       imshow(I), figure, imshow(I2)
%
%   See also CONV2, EDGE, FILTER2, FSAMP2, FWIND1, FWIND2.

%   Clay M. Thompson 11-17-92
%   Copyright 1993-1998 The MathWorks, Inc.  All Rights Reserved.
%   $Revision: 5.9 $  $Date: 1997/11/24 15:34:46 $

if nargin==0, error('Not enough input arguments.'); end
type = [type,'  '];
code = lower(type(1:2));
if nargin>1,
  if ~(all(size(P1)==[1 1]) | all(size(P1)==[1 2])),
     error('The second parameter must be a scalar or a 1-by-2 size vector.');
  end
  if length(P1)==1, siz = [P1 P1]; else siz = P1; end
end

if all(code=='ga'), % Gaussian filter
  if nargin<2, siz = [3 3]; end
  if nargin<3, std = .5; else std = P2; end
  if length(std)==1, std=[std std]; end
  if nargin<4, rho = 0; else rho = P3; end
  [x,y] = meshgrid(-(siz(2)-1)/2:(siz(2)-1)/2,-(siz(1)-1)/2:(siz(1)-1)/2);
if min(size(P2))==1
  x_norm=x/std(1); y_norm=y/std(2);
  h = exp(-(x_norm.*x_norm + y_norm.*y_norm - 2*rho*x_norm.*y_norm)/(2*(1 - rho^2)));
elseif size(P2)==[2 2];
C=P2;
for i=1:siz(1)
for j=1:siz(2)
X=[x(i,j) y(i,j)];
  h(i,j) = exp(-.5*(X*inv(C)*X'));
end
end
end
  h = h/sum(sum(h));

elseif all(code=='so'), % Sobel filter
  h = [1 2 1;0 0 0;-1 -2 -1];

elseif all(code=='pr'), % Prewitt filter
  h = [1 1 1;0 0 0;-1 -1 -1];

elseif all(code=='la'), % Laplacian filter
  if nargin<2, alpha = 1/5; else alpha = P1; end
  alpha = max(0,min(alpha,1));
  h1 = alpha/(alpha+1); h2 = (1-alpha)/(alpha+1);
  h = [h1 h2 h1;h2 -4/(alpha+1) h2;h1 h2 h1];

elseif all(code=='lo'), % Laplacian of Gaussian
  if nargin<2, siz = [5 5]; end
  if nargin<3, std = .5; else std = P2; end
  [x,y] = meshgrid(-(siz(2)-1)/2:(siz(2)-1)/2,-(siz(1)-1)/2:(siz(1)-1)/2);
  std2 = std*std;
  h1 = exp(-(x.*x + y.*y)/(2*std2));
  h = h1.*(x.*x + y.*y - 2*std2)/(2*pi*(std^6));
  h = h - sum(h(:))/prod(size(h)); % make the filter sum to zero
  
elseif all(code=='av'), % Smoothing filter
  if nargin<2, siz = [3 3]; end
  h = ones(siz)/prod(siz);

elseif all(code=='un'), % Unsharp filter
  if nargin<2, alpha = 1/5; else alpha = P1; end
  h = [0 0 0;0 1 0;0 0 0] - fspecials('laplacian',alpha);

else
  error('Unknown filter type.');

end
