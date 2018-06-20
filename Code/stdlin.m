function [v,qinv,m,w]=stdlin(v0,vinf,ro,wcr,w)
% [V,QINV,M,W]=STDLIN(V0,VINF,RO,WCR,W)
% Calculate and plot velocity dispersion and attenuation (1/Q) for 
% the standard linear viscoelastic solid 
%
% Inputs: 
%   V0:   Low (zero) frequency velocity
%   VINF: High (infinite) frequency velocity
%   RO:   Density
%   WCR:  Critical frequency in radians (WCR = 2*pi*fc)
%   W:    Frequencies (in radians = 2*pi*f) at which to calculate 
%         Velocity and attenuation
%         Vector should be preferably logspaced over wide frequency range
% Outputs:
%   V:    Velocity (vector) at frequencies W        
%   QINV: Attenuation (1/Q) at frequencies W
%   M:    Complex modulus as frequencies W
%   W:    Frequency vector

%Written by T. Mukerji

m0=ro*v0^2; minf=ro*vinf^2;
%w=logspace(d1,d2,100);
m=minf*(m0+i*sqrt(m0*minf)*(w/wcr))./(minf+i*sqrt(m0*minf)*(w/wcr));
v=sqrt(real(m)/ro); qinv=imag(m)./real(m);
if nargout==0
subplot(1,2,1), semilogx(w,v,'r'); xlabel('w=2*pi*f '); ylabel('V');
subplot(1,2,2), semilogx(w,qinv,'c'); xlabel('w=2*pi*f '); ylabel('1/Q');
end
