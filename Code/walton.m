function [k,mu]=walton(k0,mu0,p,phi)
%[K,MU]=WALTON(K0,MU0,P,PHI)
%calculates bulk and shear moduli (K, MU) of sphere packs using Walton's model.
%inputs: mineral moduli K0, MU0; pressure P, porosity PHI.

%Written by T. Mukerji

  lambda0=k0-(2/3)*mu0; c=9;
  a=(1/(4*pi))*(1./mu0 - 1./(mu0+lambda0));
  b=(1/(4*pi))*(1./mu0 + 1./(mu0+lambda0));
  k=(1/6)*(3*c^2*(1-phi).^2.*p./(pi^4*b.^2)).^(1/3); 
  mu=(3/5)*k.*(5*b+a)./(2*b+a);
