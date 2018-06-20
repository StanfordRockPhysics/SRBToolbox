function [Cw sigma11 sigma22 sigma33 vp1 vp3 vs3 vs1 rhob]=walton_aniso(k, mu, por, rhomin,  n, E)
% [Cw sigma11 sigma22 sigma33 vp1 vp3 vs3 vs1 rhob]=walton_aniso(k, mu, por, rhomin,  n, E)
% e.g., E=[0 0 0; 0 0 0; 0 0 .01]; [Cwalton sigma11 sigma22 sigma33 vp1 vp3 vs3 vs1 rhob]=walton_aniso(36*10^9,45*10^9,.38, 2650, 8, E)
% kim, gmin - Bulk and shear modulus of grain
% phi - porosity, 
% n - Co-ordination number
% E - strain tensor. 
% Reference: Walton, K, 1987, The Effective Elastic Moduli of a Random
% Packing of Spheres, J. Mech. Phys. Solids
%
% Written by: Kaushik Bandyopadhyay, 2008

lam= k - 2*mu/3;
B = (1/(4*pi))*(1./mu + 1./(mu+lam));
C = (1/(4*pi))*(1./mu - 1./(mu+lam));
F=(3*(1-por)*n)/(4*pi*pi*B*(2*B+C));
alpha=((1-por)*n*sqrt(E(3,3)))./(32*pi*pi*B);
beta=((1-por)*n*sqrt(E(3,3)))./(32*pi*pi*(2*B+C));
exact_c44=alpha+7*beta;
c44_test2=2*alpha+5*beta;
X_I1sq = dblquad(@(theta,phi)int_X_I1sq(theta,phi,E),0,pi,0,2*pi);
X_I1pow4 = dblquad(@(theta,phi)int_X_I1pow4(theta,phi,E),0,pi,0,2*pi);
X_I2sq = dblquad(@(theta,phi)int_X_I2sq(theta,phi,E),0,pi,0,2*pi);
X_I2pow4 = dblquad(@(theta,phi)int_X_I2pow4(theta,phi,E),0,pi,0,2*pi);
X_I3sq = dblquad(@(theta,phi)int_X_I3sq(theta,phi,E),0,pi,0,2*pi);
X_I3pow4 = dblquad(@(theta,phi)int_X_I3pow4(theta,phi,E),0,pi,0,2*pi);
X_I1sq_I2sq = dblquad(@(theta,phi)int_X_I1sq_I2sq(theta,phi,E),0,pi,0,2*pi);
X_I1sq_I3sq = dblquad(@(theta,phi)int_X_I1sq_I3sq(theta,phi,E),0,pi,0,2*pi);
X_I2sq_I3sq = dblquad(@(theta,phi)int_X_I2sq_I3sq(theta,phi,E),0,pi,0,2*pi);


c1111=F*(4*B*X_I1sq + 2*C*X_I1pow4);
c2222=F*(4*B*X_I2sq + 2*C*X_I2pow4);
c3333=F*(4*B*X_I3sq + 2*C*X_I3pow4);
c1122=F*(2*C*X_I1sq_I2sq);

c1133=F*(2*C*X_I1sq_I3sq);
c2233=F*(2*C*X_I2sq_I3sq);
c2323=F*(B*(X_I2sq+X_I3sq) + 2*C*X_I2sq_I3sq);
[X_I1sq 1/8];
[X_I3sq 1/4];
[X_I1sq_I3sq 1/24];
% test=(B*(X_I1sq+X_I3sq) + 2*C*X_I1sq_I3sq)
test2=((1-por)*n./(32*pi*pi*B))+7*((1-por)*n./(32*pi*pi*(2*B+C)));
test3=(3*(1-por)*n./(4*pi*pi*B*(2*B+C)))*((B/8)+(B/4)+(2*C/24));
test4=((1-por)*n/(32*pi*pi))*(9*B+C)./(B*(2*B+C));

c1313_2=F*B*(X_I1sq)+ F*B*(X_I3sq) + F*2*C*X_I1sq_I3sq;
c1313=F*(B*(X_I1sq)+B*(X_I3sq) + 2*C*X_I1sq_I3sq);
c1212=F*(B*(X_I1sq+X_I2sq) + 2*C*X_I1sq_I2sq);

Cw=zeros(6,6);
Cw(1,1)=c1111; Cw(2,2)=c2222; Cw(3,3)=c3333;
Cw(1,2)=c1122; Cw(1,3)=c1133; Cw(2,1)=Cw(1,2);Cw(3,1)=Cw(1,3);
Cw(2,3)=c2233; Cw(3,2)=Cw(2,3); 
Cw(4,4)=c2323; Cw(5,5)=c1313; Cw(6,6)=c1212;

F_sig=((1-por)*n)/(pi*pi*B*(2*B+C));
Y1_int = dblquad(@(theta,phi)int_Y1(theta,phi,E),0,pi,0,2*pi);
Y2_int = dblquad(@(theta,phi)int_Y2(theta,phi,E),0,pi,0,2*pi);

Y3_int = dblquad(@(theta,phi)int_Y3(theta,phi,E),0,pi,0,2*pi);
Y4_int = dblquad(@(theta,phi)int_Y4(theta,phi,E),0,pi,0,2*pi);

Y5_int = dblquad(@(theta,phi)int_Y5(theta,phi,E),0,pi,0,2*pi);
Y6_int = dblquad(@(theta,phi)int_Y6(theta,phi,E),0,pi,0,2*pi);

sigma11=F_sig*(B*Y1_int - C*Y2_int);
sigma22=F_sig*(B*Y3_int - C*Y4_int);
sigma33=F_sig*(B*Y5_int - C*Y6_int);

rhob = rhomin*(1-por);
vp1 = sqrt(Cw(1,1)/rhob);
vp3 = sqrt(Cw(3,3)/rhob);
vs3 = sqrt(Cw(4,4)/rhob);
vs1 = sqrt(Cw(6,6)/rhob);

