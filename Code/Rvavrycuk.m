function Rpp=Rvavrycuk(C1,C2,ro1,ro2,the,phi);
%function Rpp=Rvavrycuk(C1,C2,ro1,ro2,the,phi);
%calculates the PP reflectivity for the interface between 2 weakly,
%but arbitrarily anisotropic media;
%
%input parameters: 
%C1, C2  stiffnesses of the 2 media
%ro1, ro2  densities of the 2 media
%the: angle of incidence
%phi: azimuth (for HTI azimuth=0 is perpendicular to the fracture plane)
%
%output parameter: Rpp -reflectivity as a function of angle and azimuth
%References: Vavrycuk and Psencik, 1998 Geophysics 63, 2129-2141;

%written by Diana Sava, 04/2001

A1=C1./ro1;
A2=C2./ro2;
A11_1=A1(1,1);
A12_1=A1(1,2);
A13_1=A1(1,3);
A22_1=A1(2,2);
A23_1=A1(2,3);
A33_1=A1(3,3);
A44_1=A1(4,4);
A55_1=A1(5,5);
A66_1=A1(6,6);
A16_1=A1(1,6);
A26_1=A1(2,6);
A36_1=A1(3,6);
A45_1=A1(4,5);
A11_2=A2(1,1);
A12_2=A2(1,2);
A13_2=A2(1,3);
A22_2=A2(2,2);
A23_2=A2(2,3);
A33_2=A2(3,3);
A44_2=A2(4,4);
A55_2=A2(5,5);
A66_2=A2(6,6);
A16_2=A2(1,6);
A26_2=A2(2,6);
A36_2=A2(3,6);
A45_2=A2(4,5);
a1=sqrt(A33_1);
a2=sqrt(A33_2);
b1=sqrt(A55_1);
b2=sqrt(A55_2);
Z1=ro1.*a1;
Z2=ro2.*a2;
G1=ro1.*A55_1;
G2=ro2.*A55_2;
da=a2-a1;
aa=(a2+a1)./2;
dZ=Z2-Z1;
aZ=(Z2+Z1)./2;
dG=G2-G1;
aG=(G1+G2)./2;
ab=(b2+b1)./2;
db=b2-b1;
dro=ro2-ro1;
aro=(ro2+ro1)./2;
the=the.*pi./180;
phi=phi.*pi./180;
[the,phi]=meshgrid(the,phi);
cphi2=cos(phi).^2;
sphi2=1-cphi2;
cphi4=cphi2.^2;
sthe2=sin(the).^2;
tthe2=tan(the).^2;
sphi4=sin(phi).^4;
cphi=cos(phi);
sphi=sin(phi);
sthe=sin(the);
cphi3=cos(phi).^3;
sphi3=sin(phi).^3;
dA33=A33_2-A33_1;
dA362A45_1=A36_1+2*A45_1;
dA362A45_2=A36_2+2*A45_2;
dA362A45=dA362A45_2-dA362A45_1;
dA44=A44_2-A44_1;
dA55=A55_2-A55_1;
dA45=A45_2-A45_1;
dA11A33_1=A11_1-A33_1;
dA11A33_2=A11_2-A33_2;
dA11A33=dA11A33_2-dA11A33_1;
dA22A33_1=A22_1-A33_1;
dA22A33_2=A22_2-A33_2;
dA22A33=dA22A33_2-dA22A33_1;
dA11A33=dA11A33_2-dA11A33_1;
dA16=A16_2-A16_1;
dA26=A26_2-A26_1;
d1=(A13_1+2*A55_1-A33_1);
d2=(A13_2+2*A55_2-A33_2);
d=d2-d1;
f1=(A23_1+2*A44_1-A33_1);
f2=(A23_2+2*A44_2-A33_2);
f=f2-f1;
g1=(A44_1-A55_1)/2./A33_1;
g2=(A44_2-A55_2)/2./A33_2;
g=g2-g1;
k1=(A36_1+2*A45_1)./A33_1;
k2=(A36_2+2*A45_2)./A33_2;
k=k2-k1;
l1=A45_1./A33_1;
l2=A45_2./A33_2;
l=l2-l1;
m1=(A11_1-A33_1)./2./A33_1;
m2=(A11_2-A33_2)./2./A33_2;
m=m2-m1;
n1=(A22_1-A33_1)./2./A33_1;
n2=(A22_2-A33_2)./2./A33_2;
n=n2-n1;
o1=(A12_1+2*A66_1-A33_1);
o2=(A12_2+2*A66_2-A33_2);
o=o2-o1;
p1=A16_1./A33_1;
p2=A16_2./A33_2;
p=p2-p1;
q1=A26_1./A33_1;
q2=A26_2./A33_2;
q=q2-q1;

Rpp1=(aro.*dA33+2*aa.^2.*dro)./4./aro./aa.^2;
Rpp2=.5*(d./aa.^2.*cphi2+f./aa.^2.*sphi2+2*dA362A45/aa.^2.*cphi.*sphi-4*dA55./aa.^2.*cphi2-8*dA45./aa.^2.*cphi.*sphi-4*dA44./aa.^2.*sphi2-4*ab.^2.*dro./aro./aa.^2+dA33./2./aa.^2).*sthe2;
Rpp3=.5*(dA33./2./aa.^2+dA11A33./2./aa.^2.*cphi4+dA22A33./2./aa^2.*sphi4+o./aa.^2.*cphi2.*sphi2+2*(dA16.*cphi2+dA26.*sphi2).*sphi.*cphi./aa.^2).*sthe2.*tthe2;
Rpp=Rpp1+Rpp2+Rpp3;
