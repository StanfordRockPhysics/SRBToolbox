function [F_sc,POR]=berrysc_res(Rmin,Rw,a11,a12,a21,a22)
%[F_sc,POR]=berrysc_res(Rmin,Rw,a11,a12,a21,a22)
%Effective formation factor using Berryman (1995)'s Self-Consistent
% method for spheroidal pores and spheroidal grains.
% Input
%	Rmin,Rw: resistivities of the two constituent
%		     phases (min:mineral phase, w: pore fluid)
%   a11, a12: aspect ratios for spheroids of phase 1 (mineral)
%   a21, a22: aspect ratios for spheroids of phase 2 (pore fluid)
% Output
%	F_sc:  	Effective formation factor (Rt/Rw)
%	POR:	Porosity (fraction of phase 2)
% Example: Assuming quartz (aspect ratios=1), and a brine with resistivity 0.1 (aspect ratios =1), we have:
% Rmin=10^15; Rw=0.1;[F_sc,POR]=berrysc_res(Rmin,Rw,1,1,1,1);
% figure;semilogy(POR,F_sc);ylim([1 10^4]); 
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note of caution: Numerical unstability occurs as we approach zero
% porosity or the percolation threshold.
% Written by C. Gomez 2008
%
%******* Berryman SC oblate and prolate spheroids *****************
cond1=1/Rmin;cond2=1/Rw;
cond_br=[];
j=0;POR=(0:0.01:1);
epsilon1=1e-4;epsilon2=1e-4;
cond_new=0.;

tol_c=1e-6;

%prolate spheroids:a1>a2=a3
if a11>a12;
ee=sqrt(1-(a12/a11)^2);
L1=((1-ee^2)/(2*ee^3))*(log((1+ee)/(1-ee))-2*ee);
end
if a11<a12
% %oblate spheroids:a1<a2=a3
ee=sqrt(1-(a11/a12)^2);
L1=(1/ee^2)*(1-(sqrt(1-ee^2)*asin(ee))/ee);
end
if a11==a12
    L1=1/3;%spherical inclusions
end

%prolate spheroids:a1>a2=a3
if a21>a22;
ee=sqrt(1-(a21/a22)^2);
L2=((1-ee^2)/(2*ee^3))*(log((1+ee)/(1-ee))-2*ee);
end
if a21<a22
% %oblate spheroids:a1<a2=a3
ee=sqrt(1-(a21/a22)^2);
L2=(1/ee^2)*(1-(sqrt(1-ee^2)*asin(ee))/ee);
end
if a21==a22
    L2=1/3;%spherical inclusions
end

La=L1;Lb=(1-L1)/2;Lc=Lb;
La2=L2;Lb2=(1-L2)/2;Lc2=Lb2;k=1;
xi=0.0:0.01:1;
for x1=0:0.01:1
if x1==0, x1=x1+epsilon1;end;
if x1==1, x1=x1-epsilon2;end;
 x2=1-x1;
 j=j+1;
condsc=1/3400;
del_c=abs(condsc-cond_new);
condo=1;
niter=0;  %number of iterations
while ((del_c > abs(tol_c)) & (niter<30000))
    aa=log(cond1/condo);bb=log(cond2/condo);cc=log(condsc/condo);
    t1=(La2*exp(bb)+(1-La2)*exp(cc));t2=(Lb2*exp(bb)+(1-Lb2)*exp(cc));t3=(Lc2*exp(bb)+(1-Lc2)*exp(cc));
    t4=(La*exp(aa)+(1-La)*exp(cc));t5=(Lb*exp(aa)+(1-Lb)*exp(cc));t6=(Lc*exp(aa)+(1-Lc)*exp(cc));
    tt=(1/t1+1/t2+1/t3)/(1/t4+1/t5+1/t6);
    ecc=((x1*exp(aa))/(x2*tt)+exp(bb))/(1+x1/(x2*tt));
    ecc_sc=condsc/condo;
    del_c=abs(ecc_sc-ecc);
 	condsc=ecc*condo;
 	niter=niter+1;
end		
	cond_br=[cond_br,condsc]; 
end
lxi=length(xi);
for j=0:lxi-1;
jj=j+1;
F_sc(jj)=1./(Rw.*cond_br(lxi-j));
end
