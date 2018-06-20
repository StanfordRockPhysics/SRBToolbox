function R=Rruger(a1,b1,e1,g1,d1,rho1,a2,b2,e2,g2,d2,rho2,the,phi)
% function R=Rruger(a1,b1,e1,g1,d1,rho1,a2,b2,e2,g2,d2,rho2,the,phi)
% Rruger - calculate the reflectivity in weakly anisotropic HTI media using Ruger's
%  approximation
%
% a1, b1: P- and S-wave velocities of upper medium perpendicular to
% symmetry axis (vertical velocities for HTI)
% e1, g1, d1: epsilon gamma delta, Thomsen weakly anisotropic parameters of the uper medium
% a2, b2: P- and S-wave velocities perpendicular to symmetry axis in the lower medium
% e3,g2,d2: thomsen's aniso parameters of lower medium
% the, phi: the incidence angle and azimuth of seismic waves
% phi=0 is the azimuth perpendicular to the fracture plane
%
% R reflectivity as a function of incidence angle and azimuth
%
% reference:
% Ruger, A., 1998, Variation of P-wave reflectivity coefficients with
% offset and azimuth in anisotropic media
%    Geophysics, Vol 63, No 3, p935
%

% written by Li Teng 07/20/98 / modified by Diana Sava 09/2000
%



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
a01=a1./sqrt(1+2*e1);
b01=b1./sqrt(1+2*g1);
a02=a2./sqrt(1+2*e2);
b02=b2./sqrt(1+2*g2);
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
R=(Z2-Z1)./(Z2+Z1)+0.5.*(da./a-f.*dG./G+((dv2-dv1)+2.*f.*(g2-g1)).*cphi2).*sthe2+0.5.*(da./a+(ev2-ev1).*cphi4+(dv2-dv1).*sphi2.*cphi2).*sthe2.*tthe2;

