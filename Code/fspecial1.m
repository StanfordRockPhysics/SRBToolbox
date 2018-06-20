function h = fspecial1(type,P1,P2)
%FSPECIAL1 Create one dimenstional predefined filters.
%   H = FSPECIAL1(TYPE) creates a one-dimensional filter H of the
%   specified type. (FSPECIAL1 returns H as a computational
%   molecule, which is the appropriate form to use with FILTER.)
%   TYPE is a string having one of these values:
%
%        'gaussian'  for a Gaussian lowpass filter
%        'average'   for an averaging filter
%
%   Depending on TYPE, FSPECIAL1 can take additional parameters
%   which you can supply.  These parameters all have default
%   values. 
%
%   H = FSPECIAL1('gaussian',N,SIGMA) returns a symmetric 
%   Gaussian lowpass filter with standard deviation 
%   SIGMA (in pixels). N specifies the number of rows in H.
%   If you do not specify the parameters, FSPECIAL1 uses the 
%   default values of N = 3 and 0.5 for SIGMA.
%
%   H = FSPECIAL1('average',N) returns an averaging filter. N specifies
%   the number of rows in H. If you do not specify N, FSPECIAL1 uses 
%   the default value of N = 3.
%
%   See also CONV, FILTER, FSPECIAL, FSPECIAL3.

%   Clay M. Thompson 11-17-92
%   Copyright 1993-1998 The MathWorks, Inc.  All Rights Reserved.
%   $Revision: 5.9 $  $Date: 1997/11/24 15:34:46 $

%   Modified by Isao Takahashi 

if nargin==0, error('Not enough input arguments.'); end
type = [type,'  '];
code = lower(type(1:2));
if nargin>1,
  if ~all(size(P1)==[1 1]),
     error('The second parameter must be a scalar.');
  end
  siz = P1;
end

if all(code=='ga'), % Gaussian filter
  if nargin<2, siz = 3; end
  if nargin<3, std = .5; else std = P2; end
  x = (-(siz-1)/2:(siz-1)/2);
  h = exp(-(x.*x)/(2*std*std));
  h = h/sum(h);

elseif all(code=='av'), % Smoothing filter
  if nargin<2, siz = 3; end
  h = ones(1,siz)/siz;
else
  error('Unknown filter type.');

end
