function [vp,vs]=waltonv(sr,kmin,gmin,rhomin,pres,phi,C)
%function [VP,VS]=WALTON(KMIN,GMIN,RHOMIN,C,PHI,PRES,SR)
%calculates Vp and Vs of dry sphere packs at hydrostatic pressure condition
%using Walton's model.
%
%inputs:
%	SR	:	'r' for Walton's rough model, 's' for smooth model.
%	KMIN	:	mineral bulk modulus.
%	GMIN	:	mineral shear modulus.
%	RHOMIN	:	mineral density.
%	PHI	:	porosity.
%	PRES	:	pressure in Pa.
%	C	:	coordination number.
%outputs:
%	VP	:	Vp in m/s.
%	VS	:	Vs in m/s.
%
%	The program use the porosity-coordination number relation in page150 of 
%	the Rock Physics Handbook, if C is omitted.
%
%	From chapter 1 (Equations 3.32, 3.33, 3.46, 3.46, 3.48) in 
%	Seisic and acoustic velocities in reservoir rocks, vol.2
%	by Wang and Nur. 
%	Original paper is Walton, K., 1987, 
%	"The effective elastic moduli of a random packing of spheres"
%	J. Mech. Phys. Solids, vol.35, p213-226.

%Written by Isao Takahashi 4/12/00

if nargin<=6
ctemp=[14.007 12.336 10.843 9.5078 8.3147 7.2517 6.3108 5.4878 4.7826 4.1988 3.7440]';
por=[.2:.05:.7]';
C=interp1(por,ctemp,phi);
end;

lmin=kmin-2/3*gmin;

n_f=C.^2./(1-phi);
b=(1./gmin+1./(lmin+gmin))./(4*pi);
c=(1./gmin-1./(lmin+gmin))./(4*pi);
a=(3*n_f.*pres./(pi^4*b.^2)).^(1/3);

if sr=='r'
vp=sqrt(a.*(10*b+3*c)./(10*rhomin.*(2*b+c)));
vs=sqrt(a.*(5*b+c)./(10*rhomin.*(2*b+c)));
elseif sr=='s'
vp=sqrt(a.*3./(10*rhomin));
vs=sqrt(a.*1./(10*rhomin));
end
