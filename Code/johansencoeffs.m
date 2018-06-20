function [w000, w200, w400] = johansencoeffs(a);
% [w000, w200, w400] = johansencoeffs(a);
%
% Computes coefficients to legendre polynomial expansion for Johansen's
% orientation distribution function for shales:
%
%      odf = a^2./(x.^2 + a^2*(1-x.^2))^(3/2) ;   where x = cos(theta)
%
% Coefficients are obtained from interpolation from previously computed
% w200 vs. a, and w400 vs. a  found by integrating
%
%      integral( p2*odf)   from -1:1
%      integral( p4*odf)   from -1:1
%
% in function legendrecoeffs.m
%
% reference:  Johansen et al, 2004, Geophysical Prospecting, 52, 133-149

% written by G. Mavko, July 2005

a0=[1:.5:20];
w2 = [   0.0000    0.0068    0.0116    0.0153    0.0181    0.0204   0.0222    0.0238    0.0251    0.0262 ...
          0.0271    0.0279    0.0287    0.0293    0.0299    0.0304  0.0309    0.0313    0.0317    0.0320 ...
          0.0324    0.0327    0.0329    0.0332    0.0334    0.0337  0.0339    0.0341    0.0343    0.0344 ...
          0.0346    0.0347    0.0349    0.0350    0.0352    0.0353  0.0354    0.0355    0.0356];
w4 = [-0.0000    0.0018    0.0051    0.0086    0.0119    0.0150     0.0177    0.0200    0.0222    0.0240 ...
        0.0257    0.0272    0.0286    0.0298    0.0310    0.0320    0.0329    0.0338    0.0346    0.0353 ...
        0.0360    0.0366    0.0372    0.0378    0.0383    0.0388    0.0392    0.0397    0.0401    0.0404 ...
        0.0408    0.0411    0.0415    0.0418    0.0421    0.0423    0.0426    0.0429    0.0431];

w000 = 1/(sqrt(32)*pi^2);
w200 = interp1(a0, w2, a);
w400 = interp1(a0, w4, a);
