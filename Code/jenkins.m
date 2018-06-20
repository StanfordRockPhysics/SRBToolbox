function [VP_Jen,VS_Jen,K_eff,G_eff] = jenkins(dia,Kmin,Gmin,Rhomin,Pressure,phi,friction,k)
% [VP_Jen,VS_Jen,K_eff,G_eff] = jenkins(dia,Kmin,Gmin,Rhomin,Pressure,phi,friction,k)
% Computes the Vp, Vs, bulk and shear moduli  of isotropic, random aggregate of identical and 
% frictionless spheres using the Jenkins et al. (2005)model. 
%Example
% [VP_Jen,VS_Jen,K_eff,G_eff]= jenkins (0.1995*10^-3,36.6*10^9,45*10^9,2.65*10^3,40*10^6,0.4,0,6);
%
% INPUTS
%
%dia        Grain diameter (m)
%Kmin       Mineral bulk-modulus(Pa)
%Gmin       Mineral shear-modulus (Pa)
%Rhomin     Mineral  density (kg/m3)
%Pressure   Effective pressure (Pa)
%phi        Porosity (in fraction,e.g., 0.4)
%friction   '0' for frictionless spheres (original Jenkins paper), 
%           '1' for frictional spheres (an adhoc multiplication term from Walton's model)                  
%k     	     Coordination number (the model is not valid for k <  5.1)
%
% OUTPUTS
%K_eff      Effective Bulk modulus (MPa)
%G_eff      Effective Shear modulus (MPa)
%VP_Jen     P-wave Velocity (m/sec)
%VS_Jen     S-wave Velocity (m/sec)

%Written by Tanima Dutta July'2005, last modified  June'2009

%Reference: Jenkins, J., Johnson, D., Ragione, L. La., and  Makse, H, 2005,
% Fluctuations and the effective moduli of an isotropic, random aggregate 
%of identical, frictionless spheres: Journal of the Mechanics and Physics of Solids, 53, 197 – 225

%If k is not specified, use Murphy's relation of coordination number and
%porosity
if nargin==7,
ctemp=[14.007 12.336 10.843 9.5078 8.3147 7.2517 6.3108 5.4878 4.7826 4.1988 3.7440];
por=[.2:.05:.7];
k=interp1(por,ctemp,phi);
end;
    
%Define coefficients based on co-ordination no
%(equation C.1)
g1 = -(0.52.* (k-2).*(k-4) + 0.10.*k.*(k-2) - 0.13.*k.*(k-4) - 0.01.*k.*k); g1= g1 / (16*pi);

g2 = (0.44.* (k-2).*(k-4) - 0.24.*k.*(k-2) - 0.11.*k.*(k-4) - 0.14.*k.*k); g2= g2 / (16*pi);

g3 = -(0.44.* (k-2).*(k-4) - 0.42.*k.*(k-2) - 0.11.*k.*(k-4) + 0.04.*k.*k); g3= g3 / (16*pi);
%(equation C.3)

rho1= (1.96 .* (k-2).* (k-4) + 3.30 .* k.* (k-2) + 0.49 .*k.* (k-4) + 0.32.*k.*k);rho1= rho1 / (16*pi);

rho2= - (2.16 .* (k-2).* (k-4) + 2.30 .* k.* (k-2) + 0.54 .*k.* (k-4) - 0.06.*k.*k);rho2= rho2 / (16*pi);

alpha1= (19*k -22) / 48; alpha2= (22 - 3*k) / 16;
w_1=(166-11*k)/128; w_2=-(k+14)/128;
alpha2_tilda= (18-9*k)/48; 
w1_tilda= (38-11*k)/128; 
%(adding C.1 and C.2)
a1=w1_tilda +g1; a2= w_2+ g2; a3= w_2+ g3;

%(adding C.3 and C.4)
b1=rho1 +alpha1; b2= rho2+ alpha2_tilda; 
%(equation 26 and C.2 nd C.1)
kappa_1 = (a1 ) -(alpha1.* w1_tilda + w1_tilda.* alpha2_tilda +2.*w_2.*alpha2_tilda);
kappa_2 = (a2 ) -(w_2.*alpha1);
kappa_3 = (a3 ) -(w_2.* alpha2_tilda +w_2.*alpha1);
%(equation 27 and C.3 nd C.4)
n_1=b1 - (alpha1.*alpha1);
n_2=b2 - (2.*alpha1.*alpha2_tilda + alpha2_tilda.* alpha2_tilda);


e_1= n_1.*w_1 + n_2.*w_1 + 2.*n_2.*w_2;
e_2= n_1.*w_2;
e_3= n_2.*w_2 + n_1.*w_2;

%calculate PR,v_solid,delta and K_n 

PR= (3*Kmin - 2*Gmin)/ (2*(3*Kmin + Gmin));
v_solid= 1-phi;


a=(3*pi/2);
b=(1- PR)./ (v_solid.*k);
c=Pressure/Gmin;
d = (a.*b.*c).^(2/3);
delta=dia .*d;

d2=sqrt(delta);
K_n= (Gmin .* ((dia).^0.5).*d2)./ (1-PR);

%calculate Effective moduli (multiplication term is added from Walton's
%model)

G_eff= ((k .* v_solid) ./ (5*pi*dia)) .* K_n .* ( 1 - 2.*(((k/3).^(-1)) .* (w_1+2*w_2) - ...
        ((k/3).^(-2)) .* (kappa_1+2*kappa_2) + ((k/3).^(-3)) .* (e_1+2*e_2))) ...
        .* ((2-PR + 3*friction.*(1-PR))./(2-PR));
    
G_eff=G_eff/ 10^6;   

lambda_eff= ((k .* v_solid) ./ (5*pi*dia)) .* K_n .* (1-  2* ((k/3).^(-1)) .* (w_1+7*w_2) + ...
        2*((k/3).^(-2)) .* (kappa_1+2*kappa_2+5*kappa_3) - 2*((k/3).^(-3)) .* (e_1+2*e_2+5*e_3))...
        .* ((2-PR -2*friction*(1-PR))./(2-PR));  
lambda_eff=lambda_eff/ 10^6;  
  
K_eff= lambda_eff + (2/3)*G_eff;
    
%calculate velocities
RHOB=(1-phi).*Rhomin;
[VP_Jen,VS_Jen]=ku2v(K_eff*10^6,G_eff*10^6,RHOB);