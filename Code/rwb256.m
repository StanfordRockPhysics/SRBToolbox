function map=rwb
%function map=rwb
%red-white-blue colormap of size 256 by 3.

map=[1 0 0;1 1 1;0 0 1];             
map=interp1([0;127;255],map,[0:255]);

