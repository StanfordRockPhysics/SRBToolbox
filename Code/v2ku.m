function [k,u,nu]=v2ku(vp,vs,rho)
% [K,U,NU]=V2KU(VP,VS,RHO)
% get bulk and shear moduli (K,U) and Poisson's ratio (NU)
% from velocities (VP,VS) and density (RHO)
%
% input:    VP,VS,RHO
% output:   K,U,NU
% note:	    for isotropic media only;

% Written by Frank Liu

vp2=vp.*vp;
vs2=vs.*vs;

u=rho.*vs2;
k=rho.*vp2-4*u/3;

nu=(vp2-2*vs2)./(vp2-vs2)/2;

