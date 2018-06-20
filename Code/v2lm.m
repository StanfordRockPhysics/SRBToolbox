function [lambda,mu,nu]=v2lm(vp,vs,rho)
%[LAMBDA,MU,NU]=V2LM(VP,VS,RHO)
% convert (VP,VS,RHO) to LAMBDA and MU 
% input:   VP, VS, RHO (P and S velocities, and density) 
% output:  LAMBDA, MU  (elastic constants)
%          NU Poisson's ratio

% Written by Frank Liu

lambda = rho.*(vp.^2 - 2*vs.^2); mu= rho.*vs.^2;
nu=(vp.^2-2*vs.^2)./(vp.^2-vs.^2)/2;

