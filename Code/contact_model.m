function [keff mueff]=contact_model (Ks, Mus,C,phi, pr, b, R,sel)
%[keff mueff]=contact_model (Ks, Mus,C,phi, pr, b, R,sel)
%Calculates effective elastic moduli predicted by different contact 
%models using the same functional form
%Inputs are by input dialog box, if called without input arguments. 
%
%Input
%
%Ks,Mus     bulk and shear moduli of solid
%C          coordination number
%phi,pr     porosity and effective pressure
%b          initial adhered radius
%R          grain radius
%sel        '1' for Hertz-Mindlin and Walton-rough, 
%           '2' for Walton-smooth, '3' for Digby
%Output
%
%keff,mueff effective bulk and shear moduli
%Example:
%
% [keff mueff]=contact_model(36, 45, 6, .36, 10, 0, 5e-4,1);
% [keff mueff]=contact_model(36, 45, 6, .36, 10, 5e-6, 5e-4,3);

% Reference: Norris and Johnson, 1997, Nonlinear Elasticity of Granular
% Media, Journal of Applied Mechanics

% Written by Tanima Dutta and Kaushik Bandyopadhyay, 2009

% internally calls function radiodlg


if nargin==0
    [sel] = radiodlg({'Select contact model:'}, {'Hetrz-Mindlin/Walton-Rough','Walton-Smooth','Digby'});
    prompt1={'Ks(GPa)','Gs (GPa)','Coord.#','Phi','Pressure (MPa)', 'Initial contact radius (For Digby)', 'Grain Radius (For Digby)'};
    defans1={'36','45','6','.36','10', '0', '5e-4'};

    getpar=inputdlg(prompt1,'Hetrz-Mindlin/Walton-Rough',1,defans1);
    for k=1:length(getpar), param(k)=str2num(getpar{k}); end;
    Ks=param(1)*10^9; Mus=param(2)*10^9; C=param(3); phi=param(4); 
    pr=param(5)*10^6; b=param(6); R=param(7); 
    if b>R
        errordlg('Contact radius can not be greater than grain radius','Error');   
    end
else
    Ks=Ks*10^9; Mus=Mus*10^9; pr=pr*10^6; 
end
format short

nu=(3*Ks-2*Mus)/(2*(3*Ks+Mus));

e=0.000001:.0001:.04; %hydrostatic strain

Cn=(4*Mus)/(1-nu);
Ct=(8*Mus)./(2-nu);
V0=(4/3)*pi*R^3;

if sel==1
    a_n=sqrt(R.*R.*e);
    a_t=a_n;
    A_n=(2/3).*(R.*e).*a_n;
    elseif sel==2
        a_n=sqrt(R.*R.*e);
        a_t=0;
        A_n=(2/3).*(R.*e).*a_n;
    elseif sel==3
        a_n=sqrt(sqrt((R^2).*((e.*R).^2)+(b^4)/4)+(b^2)/2);
        a_t=b;
        A_n=(2/3)*(e.*R).*(a_n+(b^2)./(2*a_n));
end;


lam=(1-phi)*C.*(Cn.*a_n-Ct.*a_t)./(20*pi*R);
mu_nor=(1-phi)*C.*(Cn.*a_n+(3/2)*Ct.*a_t)./(20*pi*R);
K_nor=(1-phi)*C.*((5/3)*Cn.*a_n)./(20*pi*R);
p_nor=(1-phi).*C*R.*Cn.*A_n./(3*V0); %Pressure in Pa

keff = interp1(p_nor,K_nor,pr);
mueff = interp1(p_nor,mu_nor,pr);
