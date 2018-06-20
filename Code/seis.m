function map=seis(m)
% colormap seis
%
% Specifies a colormap for seismic displays.  Blue for positive; red for negative
% This version defines equally spaced light blue, dark blue, light gray, red, and yellow;
% then interpolates.
% This version uses 256 colors
% see also seisrev.m

% written by Gary Mavko

    x=[1 16 48 84 128 176 192 208 248 272];
    r0=[.6 .6 .1 .1 1 .8  1  1 1 .9];
    g0=[.6 .6 .2 .2 1 .3 .1 .2 1 .9]; 
    b0=[ 1  1  1  1 1 .3 .1 .2 0  0];

    r=interp1(x,r0,[1:256]);
    g=interp1(x,g0,[1:256]);
    b=interp1(x,b0,[1:256]);
    map = [r' g' b'];