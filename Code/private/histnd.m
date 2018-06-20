function [no,xo] = histnd(y,x,varargin)

if ~iscell(x), x={x}; end;

if nargout~=0
switch length(x)
case 1
[no,xo1]=hist1d(y,x{1},varargin{:});
xo={xo1};
case 2
[no,xo1,xo2]=hist2ds(y,x{1},x{2},varargin{:});
xo={xo1,xo2};
case 3
[no,xo1,xo2,xo3]=hist3d(y,x{1},x{2},x{3},varargin{:});
xo={xo1,xo2,xo3};
end

else
switch length(x)
case 1
hist1d(y,x{1},varargin{:});
case 2
hist2ds(y,x{1},x{2},varargin{:});
case 3
hist3d(y,x{1},x{2},x{3},varargin{:});
end
end
