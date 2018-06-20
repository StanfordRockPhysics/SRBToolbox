function [Rxz,Ryz]=Rorsym(a1,b1,e11,d11,e12,d12,g1,rho1,a2,b2,e21,d21,e22,d22,g2,rho2,the)
% function [Rxz,Ryz]=Rorsym(a1,b1,e11,d11,e12,d12,g1,rho1,a2,b2,e21,d21,e22,d22,g2,rho2,the)
% Rorsym- calculates the reflectivity in the symmetry plane for interfaces between 2 orthorhombic media  
%
% a1, b1, a2, b2: P- and S-wave vertical velocities of upper medium (1) and lower medium (2) 
% e11, d11, e12, d12: equivalent epsilon and delta Thomsen's anisotropic parameters in the two symmetry 
% planes of the orthorhombic medium for the upper medium (first index indicates the upper medium (1), 
%second index indicates the plane of symmetry (1 - plane perpendicular to x, 2 - plane perpendicular to y);
%
% g1 the vertical shear wave splitting parameter for the upper medium (1);
%
% e21, d21, e22, d22: equivalent epsilon and delta Thomsen's anisotropic parameters in the two symmetry 
% planes of the orthorhombic medium for the lower and upper media (first index indicates the lower (2) medium, 
%second index indicates the plane of symmetry (1 - plane perpendicular to x, 2 - plane perpendicular to y);
%
% g2 the vertical shear wave splitting parameter for the lower medium (2);
% the, phi: the incidence angle and azimuth of seismic waves
%
% Rxz  the PP reflectivity as a function of angle of incidence in xz plane (13)
% Ryz the PP reflectivity as a function of angle of incidence in yz plane (23)
%
% reference:
% Ruger, A., 1998, Variation of P-wave reflectivity coefficients with
% offset and azimuth in anisotropic media
% Geophysics, Vol 63, No 3, p935

%  Diana Sava 04/2001


Z1=rho1.*a1;
Z2=rho2.*a2;
G1=rho1.*b1.*b1;
G2=rho2.*b2.*b2;
G=(G1+G2)./2;
dG=G2-G1;

a=(a1+a2)./2;
da=a2-a1;
b=(b1+b2)./2;
db=b2-b1;

the=the.*pi./180;

sthe2=sin(the).^2;
tthe2=tan(the).^2;

f=(2.*b./a).^2;
Rxz=(Z2-Z1)./(Z2+Z1)+0.5.*(da./a-f.*(dG./G-2*(g2-g1))+d22-d12).*sthe2+.5*(da./a+e22-e12).*sthe2.*tthe2;  
Ryz=(Z2-Z1)./(Z2+Z1)+0.5.*(da./a-f.*dG./G+d21-d11).*sthe2+.5*(da./a+e21-e11).*sthe2.*tthe2;  

