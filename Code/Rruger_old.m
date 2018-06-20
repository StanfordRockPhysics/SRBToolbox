function R=Rruger(a01,b01,e1,g1,d1,rho1,a02,b02,e2,g2,d2,rho2,the,phi)
% function R=Rruger(a01,b01,e1,g1,d1,rho1,a02,b02,e2,g2,d2,rho2,the,phi)
% Rruger - calculate the reflectivity in weakly anisotropic media using Ruger's
%  approximation
% a01, b01: P- and S-wave velocities of upper media
% e1, g1, d1: epsilon gamma delta, thomsen wealy anisotropic parameters
% a02, b02, e2, g2, d2: velocities and aniso parameters of lower medium
% the, phi: the incidence angle and azimuth of seismic waves
% phi=0 is the azimuth perpendicular to the fracture plane
%

% reference:
% Ruger, A., 1997, P-wave reflection coefficients for transversely
%    isotropic models with verical and horizontal axis of symmetry,
%    Geophysics, Vol 62, No 3, p713
%

% written by Li Teng 07/20/98
%


a1=a01.*sqrt(1+2.*e1);
a2=a02.*sqrt(1+2.*e2);
b1=b01.*sqrt(1+2.*g1);
b2=b02.*sqrt(1+2.*g2);
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

ev1=-e1./(1+2.*e1);
ev2=-e2./(1+2.*e2);
f1=1-(b01./a01).^2;
f2=1-(b02./a02).^2;
dv1=(d1-2.*e1.*(1+e1./f1))./((1+2.*e1).*(1+2.*e1./f1));
dv2=(d2-2.*e2.*(1+e2./f2))./((1+2.*e2).*(1+2.*e2./f2));

the=the.*pi./180;
phi=phi.*pi./180;

[the, phi]=meshgrid(the,phi);

cphi2=cos(phi).^2;
sphi2=1-cphi2;
cphi4=cphi2.^2;
sthe2=sin(the).^2;
tthe2=tan(the).^2;

f=(2.*b./a).^2;
R=(Z2-Z1)./(Z2+Z1)+0.5.*(da./a-f.*dG./G+((dv2-dv1)+2.*f.*(g2-g1)).*cphi2).*sthe2+0.5.*(da./a+(ev2-ev1).*cphi4+(dv2-dv1).*sthe2.*cphi2).*sthe2.*tthe2;
