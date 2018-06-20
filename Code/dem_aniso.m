function [c11,c33,c44,c12,c13,v1_frac, v2_frac]=dem_aniso(C1,C2,a,max_C2_frac)
% [c11,c33,c44,c12,c13,v1_frac, v2_frac]=dem_aniso(C1,C2,a,max_C2_frac);
% C1-background stiffness, C2 - inclusion stiffness: 6x6 Voigt matrix, 
% a - aspect ratio,
% max_C2_frac - inclusion fraction goes from zero to (max_C2_frac-1)
% percentage
% Example:
% [dc11,dc33,dc44,dc12,dc13,v1_frac,v2_frac]=dem_aniso(C1,C2,.1,101)
% Interchanging the matrix and inclusion will give you stiff and soft
% estimations.

% Kaushik Bandyopadhyay Feb 07, 2008

max_C2_frac=round(max_C2_frac);
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

Ceff=C1_k;
i=1;
c11(i)=Ceff(1,1); c33(i)=Ceff(3,3); c44(i)=0.5*Ceff(4,4); c12(i)=Ceff(1,2); c13(i)=Ceff(1,3);
dv=.01;
v2_frac(i)=0;
v1_frac(i)=1-v2_frac(i);

for i=2:max_C2_frac
    [G, P] = calc_PandG2(Ceff(1,1),Ceff(3,3),.5*Ceff(4,4),Ceff(1,2),Ceff(1,3),1/a);     
    %Eshelby has a different definition of aspect ratio.
    Q1=inv(I+P*(C1_k-Ceff));  %Q for first component
    Q2=inv(I+P*(C2_k-Ceff));  %Q for second component
    dC=(dv/(1-v2_frac(i-1)))*(C2_k-Ceff)*Q2;
    v2_frac(i)=v2_frac(i-1)+dv;
    v1_frac(i)=v1_frac(i-1)-dv;
    Ceff=Ceff+dC;
    c11(i)=Ceff(1,1); c33(i)=Ceff(3,3); c44(i)=0.5*Ceff(4,4); c12(i)=Ceff(1,2); c13(i)=Ceff(1,3);
end
% figure, plot(v2_frac,dc11,'r'), hold on,plot(v2_frac,dc33,'b'), plot(phi,c11,'r:'), plot(phi,c33,'b:')
% figure, plot(v2_frac,(dc11-dc12)/2,'r'), hold on, plot(v2_frac,dc44,'b'), plot(phi,(c11-c12)/2,'r:'), plot(phi,c44,'b:')
% figure, plot(v2_frac,dc13,'r'), hold on, plot(phi,c13,'b')    
