function [c11,c33,c44,c12,c13,v1_frac, v2_frac, Ceff]=sca_aniso(C1,C2,a,v1_frac, v2_frac)
% [c11,c33,c44,c12,c13,v1_frac, v2_frac, Ceff]=sca_aniso(C1,C2,a,v1_frac, v2_frac)
% C1, C2 - inclusion stiffnesses: 6x6 Voigt matrix, a - aspect ratio,
% Example: 
% [c11,c33,c44,c12,c13,v1_frac, v2_frac, Ceff]=sca_aniso(C1,C2,.1,.5, .5);
% Interchanging the matrix and inclusion will give you stiffest and softest moduli.

%Kaushik Bandyopadhyay Feb 07, 2008
%Kaushik Bandyopadhyay May 25, 2009: Corrected Ceff to output as Voigt
%matrix instead of Kelvin's

I=eye(6,6);     %Sixth rank identity tensor
%Pre and post multiply pp to a TI Voigt matrix to take it to Kelvin's
pp=[1 0 0 0 0 0; ...
    0 1 0 0 0 0; ...
    0 0 1 0 0 0; ...
    0 0 0 sqrt(2) 0 0; ...
    0 0 0 0 sqrt(2) 0; ...
    0 0 0 0 0 sqrt(2)];

%Initial guess for the effective stiffness
C1_k=pp*C1*pp;        %C1 in Kelvin
C2_k=pp*C2*pp;        %C2 in Kelvin
Ceff=(C1_k+C2_k)./2;

err_eff=norm(Ceff-C2_k);
i=1;
c11(i)=Ceff(1,1); c33(i)=Ceff(3,3); c44(i)=0.5*Ceff(4,4); c12(i)=Ceff(1,2); c13(i)=Ceff(1,3);

while (err_eff>.001)
    
    Ceff_old=Ceff;
    [G, P] = calc_PandG2(Ceff(1,1),Ceff(3,3),.5*Ceff(4,4),Ceff(1,2),Ceff(1,3),1/a);     
    %Note I send 1/a as aspect ratio. Eshelby has a different definition of aspect ratio.
       
    Q1=inv(I+P*(C1_k-Ceff));  %Q for first component
    Q2=inv(I+P*(C2_k-Ceff));  %Q for second component
    
    Ceff=(v1_frac*C1_k*Q1+v2_frac*C2_k*Q2)*inv(v1_frac*Q1 + v2_frac*Q2);
    c11(i)=Ceff(1,1); c33(i)=Ceff(3,3); c44(i)=0.5*Ceff(4,4); c12(i)=Ceff(1,2); c13(i)=Ceff(1,3);
    
    err_eff=norm(Ceff-Ceff_old);
end
Ceff(4,4)=0.5*Ceff(4,4);
Ceff(5,5)=0.5*Ceff(5,5);
Ceff(6,6)=0.5*Ceff(6,6);
% figure, plot(v2_frac,dc11,'r'), hold on,plot(v2_frac,dc33,'b'), plot(phi,c11,'r:'), plot(phi,c33,'b:')
% figure, plot(v2_frac,(dc11-dc12)/2,'r'), hold on, plot(v2_frac,dc44,'b'), plot(phi,(c11-c12)/2,'r:'), plot(phi,c44,'b:')
% figure, plot(v2_frac,dc13,'r'), hold on, plot(phi,c13,'b')    
