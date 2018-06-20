function [no,xo1] = hist1d(y,x1,weight)
%HIST1D  1 Dimensional Histogram.
%   N = HIST1d(Y) bins the elements of Y into 15 equally spaced containers
%   and returns the number of elements in each container.  
%   
%   Y must be a one-column matrix.
%
%   N = HIST1D(Y,X1), where X1 is a scalar, uses X1 bins.
%
%   N = HIST1D(Y,X1), where X1 is a vector, returns the distribution of Y
%	among bins with centers specified by X1.
%
%   N = HIST1D(Y,X1,W), where W is a vector with the same length as Y,
%	returns N as the sum of the weights specified by W, in each bin.
%
%   [N,X1] = HIST2D(...) also returns the position of the bin centers in X1.
%
%   HIST1D(...) without output arguments produces a histgram bar plot of
%   the results.
%
%   See also HIST, HIST2D.

%   J.N. Little 2-06-86
%   Revised 10-29-87, 12-29-88 LS
%   Revised 8-13-91 by cmt, 2-3-92 by ls.
%   Copyright (c) 1984-98 by The MathWorks, Inc.
%   $Revision: 5.13 $  $Date: 1997/12/02 19:27:08 $

%   Modified by Isao.Takahashi 1998/6/30
%   Modified by T. Mukerji 1998/7/12

if nargin == 0
   error('Requires one or two input arguments.')
end
if (nargin == 1), x1 = 10; weight = 1; end
if (nargin == 2), weight = 1; end;
[m,n]=size(y);
if min(m,n)~=1
   error('Requires vector input')
end 
if isstr(y) | isstr(x1) 
   error('Input arguments must be numeric.')
end
   miny1 = min(y(:,1)); maxy1 = max(y(:,1));


if length(x1) == 1
   nbin1=x1;
   if miny1 == maxy1,
        miny1 = miny1 - floor(x1/2) - 0.5;
        maxy1 = maxy1 + ceil(x1/2) - 0.5;
   end
   binwidth1 = (maxy1 - miny1) ./ x1;
   xx1 = miny1 + binwidth1*(0:x1);
   xx1(length(xx1)) = maxy1;
   x1 = [xx1(1:length(xx1)-1) + binwidth1/2]';

   y1 = y(:,1);
   y1 = ceil((y1 - miny1)/binwidth1);
   y1 = y1 + (y1==0) - (y1>nbin1);

else
   xx1 = x1(:)';
   binwidth1 = [diff(xx1) 0];
   xx1 = [xx1(1)-binwidth1(1)/2 xx1+binwidth1/2];
   xx1(1) = miny1; xx1(length(xx1)) = maxy1;

   nbin1 = length(xx1)-1;
   y1 = y(:,1);
   for k=1:length(y1);
	y1(k) = sum(y1(k)>=xx1);
   end
   y1 = y1 + (y1==0) - (y1>nbin1);

end

S  = sparse(y1,1:length(y1),weight,nbin1,length(y1));
nn = full(sum(S'));

if nargout == 0
   bar(x1,nn,'hist');
else
   no = nn; xo1 = x1;
end;
