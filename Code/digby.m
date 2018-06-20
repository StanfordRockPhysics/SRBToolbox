function [Keff, mueff]=digby(Ks, Mus, Pr, phi, C, a, R)
% [Keff, mueff]=digby(Ks, Mus, Pr, phi, C, a, R)
% Computes effective, elastic properties for a random pack of identical 
% spheres with an intial contact radius using Dibgy (1981).
% The output of this program will match with the smooth model of 
% Walton (1987), When the initial radius is zero.
% Example:
% [Keff, mueff]=digby(35*10^9, 45*10^9, 10*10^6, .36, 6, 5e-6, 5e-4); 
% 
% Inputs:
% Ks, Mus: grain moduli (in Pa)
% phi: porosity (in fraction)
% C: coordination number
% PR: pressure (in Pa)
% a:initial contact radius (in meter)
% R: grain radius (in meter)
%         
% Outputs:
% Keff: Effective bulk modulus (Pa)
% Mueff: Effective shear modulus (Pa)

% written by Kaushik Bandyopadhyay and Tanima Dutta: March, 2009

% Internally calls function: srb_cubicfcn()

nu=(3*Ks-2*Mus)./(2*(3*Ks+Mus));

for l=1:length(phi)
    for m=1:length(Pr)
        % a*x^3 + b*x^2 + c*x + d = 0
        [x,nroot]=srb_cubicfcn(1, 0, (3/2)*(a/R)^2, -(3*pi.*(1-nu).*Pr(m))./(2*C.*(1-phi(l)).*Mus));

        d=x(1); %Real root
        b_by_R=sqrt(d.^2 + (a/R).^2);

        Sn_by_R=(4*Mus.*b_by_R)./(1-nu);
        St_by_R=(8*Mus.*a/R)./(2-nu);

        Keff(l,m)=C.*(1-phi(l)).*Sn_by_R./(12*pi);
        mueff(l,m)=C.*(1-phi(l)).*(Sn_by_R+1.5*St_by_R)./(20*pi);
    end
end
