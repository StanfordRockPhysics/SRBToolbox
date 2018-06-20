function [no,xo1,xo2] = hist2d(y,x1,x2,weight)
%HIST2D  2 Dimensional Histogram.
%   N = HIST2d(Y) bins the elements of Y into 15 equally spaced containers
%   and returns the number of elements in each container.  
%   
%   Y must be a two-column matrix.
%
%   N = HIST2D(Y,X1,X2), where X1 and X2 are scalars, uses X1*X2 bins.
%
%   N = HIST2D(Y,X1,X2), where X1 and X2 are vectors, returns the 
%       distribution of Y among bins with centers specified by X1 and X2.
%
%   N = HIST2D(Y,X1) does the same operation as N = HIST(Y,X1,X1)
%
%   N = HIST2D(Y,X1,X2,W), where W is a vector with the same length as Y,
%       returns N as the sum of the weights, specified by W, in each bin.
%	Every element of W must be between 0 and 1. Sum(N) is equal to Sum(N).
%
%   [N,X1,X2] = HIST2D(...) also returns the position of the bin centers 
%               in X1 and X2.
%
%   HIST2D(...) without output arguments produces a grayscale plot of
%   the results.
%
%   See also HIST, HIST1D.

%   J.N. Little 2-06-86
%   Revised 10-29-87, 12-29-88 LS
%   Revised 8-13-91 by cmt, 2-3-92 by ls.
%   Copyright (c) 1984-98 by The MathWorks, Inc.
%   $Revision: 5.13 $  $Date: 1997/12/02 19:27:08 $

%   Modified by Isao.Takahashi 1998/6/30
%   Modified by T. Mukerji 1998/7/12
%   Modified by Isao.Takahashi 1999/4/19

if nargin == 0
   error('Requires one or more input arguments.')
end
if (nargin == 1) x1 = 15; x2 = 15; weight = 1; end
if (nargin == 2) x2 = x1; weight = 1; end;
if (nargin == 3) weight = 1; end;
if (max(weight)>1)|(min(weight)<0), 
   error('Weight must be between 0 and 1')
end
[m,n]=size(y);
if n==1
   if nargout == 0
        hist(y,x1);
   else
        [no, xo1] = hist(y,x1);
        xo1 = xo1'; xo2 = xo1; 
   end
elseif n~=2
   error('Requires matrix input with two columns')
else 
if isstr(y) | isstr(x1) | isstr(x2) 
   error('Input arguments must be numeric.')
end
   miny1 = min(y(:,1)); maxy1 = max(y(:,1));
   miny2 = min(y(:,2)); maxy2 = max(y(:,2));

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
if length(x2) == 1
   nbin2=x2;
   if miny2 == maxy2,
        miny2 = miny2 - floor(x2/2) - 0.5; 
        maxy2 = maxy2 + ceil(x2/2) - 0.5;
   end
   binwidth2 = (maxy2 - miny2) ./ x2;
   xx2 = miny2 + binwidth2*(0:x2);
   xx2(length(xx2)) = maxy2;
   x2 = [xx2(1:length(xx2)-1) + binwidth2/2]';

   y2 = y(:,2);
   y2 = ceil((y2 - miny2)/binwidth2);
   y2 = y2 + (y2==0) - (y2>nbin2);

else
   xx2 = x2(:)';
   binwidth2 = [diff(xx2) 0];
   xx2 = [xx2(1)-binwidth2(1)/2 xx2+binwidth2/2];
   xx2(1) = miny2; xx2(length(xx2)) = maxy2;

   nbin2 = length(xx2)-1;
   y2 = y(:,2);
   for k=1:length(y2);
	y2(k) = sum(y2(k)>=xx2);
   end
   y2 = y2 + (y2==0) - (y2>nbin2);

end

ytemp = (y2-1)*nbin1+y1;
S  = sparse(ytemp,1:length(ytemp),weight,nbin1*nbin2,length(ytemp));
%size(S), nbin1, nbin2
nn = reshape(full(sum(S')),nbin1,nbin2);


%The following is the algorithm used in HIST (1-d histogram)
%nn = zeros(nbin1,nbin2);
%for i=2:nbin1
%   for j=2:nbin2
%        nn(i,j) = sum((y(:,1) <= xx1(i)).*(y(:,2) <= xx2(j)) );
%   end
%end
%nn = nn(2:nbin1,:) - nn(1:nbin1-1,:);
%nn = nn(:,2:nbin2) - nn(:,1:nbin2-1);
%end;

if nargout == 0
   imagesc(x1,x2,nn'); axis xy; colormap(1-gray);
else
   no = nn; xo1 = x1; xo2 = x2;
end;
end;
