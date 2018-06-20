function [eff_c eff_s, eff_c_frm_s]=sch_muir_bkus(c1,c2,v1,v2)
% [eff_c eff_s, eff_c_frm_s]=sch_muir_bkus(c1,c2,v1,v2)
% Effective stiffness and compliance matrix for layered medium using
% Schoenberg & Muir theory.
%c1 and c2 the stiffness matrix for individual layers
%v1 and v2 are thir thickness (in fraction)
%Outputs are the effective stiffness and compliance matrices
%eff_c: Effective stiffness 
%eff_s: Effective compliance 
%eff_c_frm_s: Effective stiffness inverted from the effective comliance
%matrix
%A calculus for finely layered anisotropic media: Schoenberg & Muir, 1989, Geophysics
%Written by: Kaushik Bandyopadhyay - 2008



Cnn1=[c1(3,3) c1(3,4) c1(3,5); c1(3,4) c1(4,4) c1(4,5); c1(3,5) c1(4,5) c1(5,5)]
Cnn2=[c2(3,3) c2(3,4) c2(3,5); c2(3,4) c2(4,4) c2(4,5); c2(3,5) c2(4,5) c2(5,5)];

Ctn1=[c1(1,3) c1(1,4) c1(1,5); c1(2,3) c1(2,4) c1(2,5); c1(3,6) c1(4,6) c1(5,6)];
Ctn2=[c2(1,3) c2(1,4) c2(1,5); c2(2,3) c2(2,4) c2(2,5); c2(3,6) c2(4,6) c2(5,6)];

Ctt1=[c1(1,1) c1(1,2) c1(1,6); c1(1,2) c1(2,2) c1(2,6); c1(1,6) c1(2,6) c1(6,6)];
Ctt2=[c2(1,1) c2(1,2) c2(1,6); c2(1,2) c2(2,2) c2(2,6); c2(1,6) c2(2,6) c2(6,6)];

Cnn=inv(v1*inv(Cnn1) + v2*inv(Cnn2))
Ctn=(v1*Ctn1*inv(Cnn1) + v2*Ctn2*inv(Cnn2))*Cnn;
Ctt=(v1*Ctt1+v2*Ctt2) - (v1*Ctn1*inv(Cnn1)*Ctn1' + v2*Ctn2*inv(Cnn2)*Ctn2')...
    + (v1*Ctn1*inv(Cnn1) + v2*Ctn2*inv(Cnn2))*Cnn*(v1*inv(Cnn1)*Ctn1'+v2*inv(Cnn2)*Ctn2');

eff_c=[Ctt(1,1) Ctt(1,2) Ctn(1,1) Ctn(1,2) Ctn(1,3) Ctt(1,3);
        Ctt(2,1) Ctt(2,2) Ctn(2,1) Ctn(2,2) Ctn(2,3) Ctt(2,3);
        Ctn(1,1) Ctn(2,1) Cnn(1,1) Cnn(2,1) Cnn(1,3) Ctn(3,1);
        Ctn(1,2) Ctn(2,2) Cnn(2,1) Cnn(2,2) Cnn(2,3) Ctn(3,2);
        Ctn(1,3) Ctn(2,3) Cnn(3,1) Cnn(3,2) Cnn(3,3) Ctn(3,3);
        Ctt(1,3) Ctt(2,3) Ctn(3,1) Ctn(3,2) Ctn(3,3) Ctt(3,3)];
    
%Compute the effective compliance matrix
%Nichols, Muir and Schoenberg, SEG expanded abstract, 1989
s1=inv(c1);
s2=inv(c2);

Snn1=[s1(3,3) s1(3,4) s1(3,5); s1(3,4) s1(4,4) s1(4,5); s1(3,5) s1(4,5) s1(5,5)];
Snn2=[s2(3,3) s2(3,4) s2(3,5); s2(3,4) s2(4,4) s2(4,5); s2(3,5) s2(4,5) s2(5,5)];

Stn1=[s1(1,3) s1(1,4) s1(1,5); s1(2,3) s1(2,4) s1(2,5); s1(3,6) s1(4,6) s1(5,6)];
Stn2=[s2(1,3) s2(1,4) s2(1,5); s2(2,3) s2(2,4) s2(2,5); s2(3,6) s2(4,6) s2(5,6)];

Snt1=Stn1';
Snt2=Stn2';

Stt1=[s1(1,1) s1(1,2) s1(1,6); s1(1,2) s1(2,2) s1(2,6); s1(1,6) s1(2,6) s1(6,6)];
Stt2=[s2(1,1) s2(1,2) s2(1,6); s2(1,2) s2(2,2) s2(2,6); s2(1,6) s2(2,6) s2(6,6)];



Snn=(v1*Snn1+v2*Snn2) - (v1*Snt1*inv(Stt1)*Stn1 + v2*Snt2*inv(Stt2)*Stn2) ...
    + (v1*Snt1*inv(Stt1) + v2*Snt2*inv(Stt2))*inv(v1*inv(Stt1)+v2*inv(Stt2))*(v1*inv(Stt1)*Stn1 + v2*inv(Stt2)*Stn2);

Stn=(inv(v1*inv(Stt1) + v2*inv(Stt2)))*(v1*inv(Stt2)*Stn2 + v2*inv(Stt2)*Stn2);

Stt=inv(v1*inv(Stt1) + v2*inv(Stt2));

eff_s = [Stt(1,1) Stt(1,2) Stn(1,1) Stn(1,2) Stn(1,3) Stt(1,3);
        Stt(2,1) Stt(2,2) Stn(2,1) Stn(2,2) Stn(2,3) Stt(2,3);
        Stn(1,1) Stn(2,1) Snn(1,1) Snn(2,1) Snn(1,3) Stn(3,1);
        Stn(1,2) Stn(2,2) Snn(2,1) Snn(2,2) Snn(2,3) Stn(3,2);
        Stn(1,3) Stn(2,3) Snn(3,1) Snn(3,2) Snn(3,3) Stn(3,3);
        Stt(1,3) Stt(2,3) Stn(3,1) Stn(3,2) Stn(3,3) Stt(3,3)];
    
eff_c_frm_s = inv(eff_s);        