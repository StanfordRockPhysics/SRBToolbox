function [C_bkus]=aniso_bkus(c1,c2,v1,v2)
%[C_bkus]=aniso_bkus(c1,c2,v1,v2)
% backus average when individual layers are VTI
%Input: 
%c1, c2: Stiffness matrix of the two layers in Voigt notation
%v1, v2: Volume fractions
%Output: 
%C_bkus=Effective stiffness of the layered medium
%See also bkusc, bkuslog
%Kaushik Bandyopadhyay - 2007
c1_11=c1(1,1);
c1_33=c1(3,3);
c1_12=c1(1,2);
c1_13=c1(1,3);
c1_44=c1(4,4);
c1_66=c1(6,6);


c2_11=c2(1,1);
c2_33=c2(3,3);
c2_12=c2(1,2);
c2_13=c2(1,3);
c2_44=c2(4,4);
c2_66=c2(6,6);

c11=(((v1*(c1_13/c1_33)+v2*(c2_13/c2_33))^2)/((v1/c1_33)+(v2/c2_33))) - (v1*(c1_13^2/c1_33)+v2*(c2_13^2/c2_33)) + (v1*c1_11+v2*c2_11);
c12=c11- (v1*c1_11+v2*c2_11) + (v1*c1_12+v2*c2_12);
c13=(v1*(c1_13/c1_33)+v2*(c2_13/c2_33))/((v1/c1_33)+(v2/c2_33));
c33=1/((v1/c1_33)+(v2/c2_33));
c44=1/((v1/c1_44)+(v2/c2_44));
c66=((v1*c1_66)+(v2*c2_66));

C_bkus=zeros(6,6);
C_bkus(1,1)=c11;
C_bkus(2,2)=c11;
C_bkus(3,3)=c33;
C_bkus(4,4)=c44;
C_bkus(5,5)=c44;
C_bkus(6,6)=c66;
C_bkus(1,2)=c12;
C_bkus(1,3)=c13;
C_bkus(2,1)=c12;
C_bkus(3,1)=c13;
C_bkus(2,3)=c13;
C_bkus(3,2)=c13;