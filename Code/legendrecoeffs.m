function [w000, w200, w400, a] = legendrecoeffs
% [w000, w200, w400, a] = legendrecoeffs
%
% Computes coefficients to legendre polynomial expansion for Johansen's
% orientation distribution function for shales:
%
%      odf = a^2./(x.^2 + a^2*(1-x.^2))^(3/2) ;   where x = cos(theta)
%
% where we let the compaction parameter a range from 1 to 5
% Computes coefficients w200 and w400 by integrating
%
%      integral( p2*odf)   from -1:1
%      integral( p4*odf)   from -1:1
%
% reference:  Johansen et al, 2004, Geophysical Prospecting, 52, 133-149

% written by G. Mavko, July 2005

p2 = '((3*x.^2-1)/2)';
p4 = '((35*x.^4-30*x.^2+3)/8)';

a=[1:.5:20];

for k=1:length(a),
     a2 = num2str(a(k).^2);
     odf = ['((' a2 './(x.^2 + ' a2 '.*(1-x.^2)).^(3/2))/(8*pi*pi))'];
     w200(k) = sqrt(5/2)*quad(inline([p2 '.*' odf]),-1,1);
     w400(k) = sqrt(9/2)*quad(inline([p4 '.*' odf]),-1,1);
end;

w000 = 1/(sqrt(32)*pi^2);
phi=1-a.*(1-.4);

% figure; plot(a,w200,'r'); xlabel('a'); ylabel('W200 and W400');
% hold on; plot(a,w400); legend('W200', 'W400')
% figure; plot(w200,w400); xlabel('w200'); ylabel('w400');
% figure; plot(phi, a); xlabel('Phi'); ylabel('a'); xlim([0 .8])
% figure; plot(phi, w200,'r', phi, w400, 'b'); xlabel('Phi'); ylabel('W'); xlim([0 .8])
% 
% 
% 
