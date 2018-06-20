function [Ksat,Gsat]=SolGassmn(Kd,Kg,Kif,Gd,Gg,Gif,phi)
%[KSAT,GSAT]=SolGassmn(Kd,Kg,Kif,Gd,Gg,Uif,phi)
%
% Generalization of Brown and Korringa’s and Gassmann equations
%    for a solid infill of the pore space by Ciz and Shapiro
% INPUTS:
% Kd: bulk modulus of original dry rock
% Kg: bulk modulus of grain material
% Kif: bulk modulus of infill material
% Gd: shear modulus of original dry rock
% Gg: shear modulus of grain material
% Gif: shear modulus of infill material
% phi: porosity
% OUTPUTS:
% Ksat: saturated bulk modulus
% Gsat: saturated shear modulus

%Written by K. Wolf

kd=1./Kd;
kg=1./Kg;
kif=1./Kif;
gd=1./Gd;
gg=1./Gg;
gif=1./Gif;

Ksat = 1./( kd - ((kd-kg).^2)./(phi*(kif-kg)+(kd-kg)) );
Gsat = 1./( gd - ((gd-gg).^2)./(phi*(gif-gg)+(gd-gg)) );
