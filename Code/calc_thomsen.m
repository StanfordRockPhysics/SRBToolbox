function [e,g,d]=calc_thomsen(c)
% [e,g,d]=calc_thomsen(c)
% Calculates Thomsen's parameters, epsilon, gamma and delta from 6x6
% stiffness matrix for a VTI medium
% see also c2anis, c2vti

% Written by: Kaushik Bandyopadhyay, 2008

e=(c(1,1)-c(3,3))/(2*c(3,3));
g=(c(6,6)-c(4,4))/(2*c(4,4));
d =((c(1,3)+c(4,4))^2 - (c(3,3) - c(4,4))^2)/(2*c(3,3)*(c(3,3)-c(4,4)));