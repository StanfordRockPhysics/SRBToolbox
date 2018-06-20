function map=rwk
%function map=rwk
%red-white-black colormap of size 64 by 3.

map=[1 0 0;1 1 1;0 0 0];             
map=interp1([0;31;63],map,[0:63]);

