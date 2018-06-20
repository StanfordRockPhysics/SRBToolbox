function map=rwb
%function map=rwb
%red-white-blue colormap of size 64 by 3.

map=[1 0 0;1 1 1;0 0 1];             
map=interp1([0;31;63],map,[0:63]);

