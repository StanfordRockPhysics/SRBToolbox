function map=seisrev(m)
% colormap seisrev
%
% Specifies a colormap for seismic displays. Blue for positive; red for
% negative -- the reverse as that specified in seis.m
% This version defines equally spaced lt blue, dark blue, light gray, red, and yellow;
% then interpolates.
% This version uses 256 colors
% see also seis.m

% Written by Gary Mavko

    x=[1 16 48 84 128 176 192 208 248 272];
    r0=[.6 .6 .1 .1 1 .8  1  1 1 .9];
    g0=[.6 .6 .2 .2 1 .3 .1 .2 1 .9]; 
    b0=[ 1  1  1  1 1 .3 .1 .2 0  0];

    r=interp1(x,r0,[1:256]);
    g=interp1(x,g0,[1:256]);
    b=interp1(x,b0,[1:256]);
    map = [r' g' b'];
    map = flipud(map);