function C = toemodel(toe, c_ref, pr_ref, sigma)
% Calculates the stiffness matrix at a particular stress using the third
% order elastic (TOE) coefficients.
%Inputs
%toe = isotropic third order elastic constants;
%c_ref = stiffness at the reference stress (VTI); send a vector [c11 c33 c44 c66 c13]
%at the reference stress ; In GPa
%pr_ref = pressure at the reference stress; In MPa
%pr= pressures at which you have measurements : [5 10 15 20]; in MPa
%sigma=[s1 s2 s3 0 0 0]; %Stress at which youwant to predict: MPa
%Outputs
%C = stiffness matrix

% Reference: Prioul et al., 2004, Nonlinear rock physics model for estimation of 3D subsurface stress
% in anisotropic formations: Theory and laboratory verification, Geophysics

% Written by: Kaushik Bandyopadhyay - 2009

c011=c_ref(1); c033=c_ref(2); c044=c_ref(3); c066=c_ref(4); c013=c_ref(5); c012=c011-2*c066; 
c_ref=makec([c011 (c011-2*c066) c013(1) c033 c044]);

T_ref=[-pr_ref; -pr_ref; -pr_ref;  0; 0; 0;].*.001;
pr=sigma.*.001;
E_ref=inv(c_ref)*T_ref;

c111=toe(1); c112=toe(2); c123=toe(3); 

    T=-pr;
    E=inv(c_ref)*T;
    dE11=E(1)-E_ref(1);
    dE22=E(2)-E_ref(2);
    dE33=E(3)-E_ref(3);

    c144=(c112-c123)/2;
    c155=(c111-c112)/4;

c11=c011 + c111.*dE11 + c112.*(dE22+dE33);
c22=c011 + c111.*dE22 + c112.*(dE11+dE33);
c33=c033 + c111.*dE33 + c112.*(dE11+dE22);
c44=c044 + c144.*dE11 + c155.*(dE22+dE33);
c55=c044 + c144.*dE22 + c155.*(dE11+dE33);
c66=c066 + c144.*dE33 + c155.*(dE11+dE22);
c12=c012 + c112.*(dE11+dE22) + c123.*dE33;
c13=c013 + c112.*(dE11+dE33) + c123.*dE22;
c23=c013 + c112.*(dE22+dE33) + c123.*dE11;
% C= [c11; c33; c44; c66; c13];
C=[c11 c12 c13 0 0 0; c12 c22 c23 0 0 0; c13 c23 c33 0 0 0; 0 0 0 c44 0 0; 0 0 0 0 c55 0; 0 0 0 0 0 c66];
% Check the anisotropy parameters
% [epsi_x, epsi_y, gam_x, gam_y, del_x, del_y, del_3] = thomp_ortho(C);
% [epsi_x, epsi_y, gam_x, gam_y, del_x, del_y, del_3]