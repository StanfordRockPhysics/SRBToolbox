function [cc,anis]=v2cti(vv,r)
% [CC,ANIS]=V2CTI(VV,R)
% Get elastic stiffnesses (Cij) from velocities (vp,vs) at 0,45,90 degrees from
% the symmetry axis for a TI medium 
%
% input:    VV=[vp0 vp45 vp90 vs0 vsh90]; R-density.
% output:   CC=[c11 c33 c44 c66 c13] 
%           ANIS=[epsilon gamma delta deltasv]; Thomsen's parameters
% note:     vp,vs and r can all be vectors (so vv, cc, anis are matrices).

% Written by Frank Liu

vp0=vv(:,1);
vp45=vv(:,2);
vp90=vv(:,3);
vs0=vv(:,4);
vsh90=vv(:,5);

r2=r.*r;
v2=vp45.*vp45;
v4=v2.*v2;

c11=r.*vp90.*vp90;
c12=c11-2*r.*vsh90.*vsh90;
c33=r.*vp0.*vp0;
c44=r.*vs0.*vs0;
c1=c11+c33+2*c44;
c2=(c11+c44).*(c33+c44);
c13=-c44+sqrt(4*r2.*v4-2*r.*v2.*c1+c2);
c66=(c11-c12)/2;

cc=[c11 c33 c44 c66 c13];

e=(c11-c33)./c33/2;
g=(c66-c44)./c44/2;
c3=c13+c44;
c4=c33-c44;
c5=c11-c44;
d=(c3.*c3-c4.*c4)./c33./c4/2;
dsv=(c5.*c4-c3.*c3)./c44./c4/2;

anis=[e g d dsv];

