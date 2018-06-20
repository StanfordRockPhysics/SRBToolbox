function C=makec(cc)
% C=makec(cc)
% Creates a 6x6 stiffness matrix from five independent moduli of a VTI medium
% input:  cc=[c11,c12,c13,c33,c44];
% output: C=6x6 stiffness matrix
% Kaushik Bandyopahdyay - 2007
C=zeros(6,6);
C(1,1)=cc(1);
C(2,2)=C(1,1);
C(3,3)=cc(4);
C(4,4)=cc(5);
C(5,5)=C(4,4);
C(6,6)=0.5*(cc(1)-cc(2));
C(1,2)=cc(2);
C(2,1)=C(1,2);
C(1,3)=cc(3);
C(2,3)=C(1,3);
C(3,2)=C(2,3);
C(3,1)=C(1,3);

