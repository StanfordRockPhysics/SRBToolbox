function [vphfsat,vshfsat,kurf,muurf,dataout]=squirt(vdata,vphp,vshp,rodry,rofl,kfl,k0,por)
%[VPHF,VSHF,KURF,MUURF,DATAOUT]=SQUIRT(VDATA,VPHP,VSHP,RODRY,ROFL,KFL,K0,POR)
%
% Implements Mavko's isotropic squirt model to calculate high-frequency
% saturated Vp and Vs from dry velocity-pressure data. Biot dispersion
% is not taken into account.
% Plots the measured (dry and saturated) and predicted (saturated low and
% high frequency) velocities versus effective pressure.
%
% Inputs: 
%  VDATA: matrix of 5 columns containing dry and saturated measurements.
%         column 1= effective pressure,  column 2= Vp dry, column 3= Vs dry
%         column 4= Vp sat. measured, column 5= Vs sat. measured.
%         The saturated measurement columns may be set to Nan if measurements
%         are not available.
%  VPHP, VSHP:  Vp, Vs at high effective stress.       
%  RODRY: Dry rock density.
%  ROFL, KFL: Fluid density and bulk modulus.
%  K0, POR: Mineral bulk modulus and rock porosity.
% 
% Outputs:
%  VPHF, VSHF: High-frequency saturated Vp, Vs.
%  KURF, MUURF:Unrelaxed, wet-frame bulk and shear moduli.
%  DATAOUT: matrix of 9 columns containing measured and predicted velocities.
%           col.1= eff. pressure, col.2= Vpdry, col.3= Vsdry
%           col.4= Vp_sat measured, col.5= Vs_sat. measured
%           col.6= Vp_sat from Gassmann, col.7= Vs_sat from Gassmann
%           col.8= Vp_sat high frequency from squirt
%           col.9= Vs_sat high frequency from squirt
%           

%Written by T. Mukerji

p=vdata(:,1); vpdry=vdata(:,2); vsdry=vdata(:,3); vpsatm=vdata(:,4);
vssatm=vdata(:,5); rosat=rodry+por*rofl;
mudry=rodry.*vsdry.^2; kdry=rodry.*(vpdry.^2-(4/3)*vsdry.^2);
khp=rodry*(vphp^2-(4/3)*vshp^2); muhp=rodry*vshp^2;
kurf=khp; khfsat=gassmnk(kurf,0.,kfl,k0,por);
muhfsat=15*kurf*mudry.*kdry./(15*kurf*kdry-4*mudry.*(kurf-kdry));
muurf=muhfsat;
vphfsat=sqrt((khfsat+(4/3)*muhfsat)./rosat);
vshfsat=sqrt(muhfsat./rosat);
[vpgass,vsgass]=gassmnv(vpdry,vsdry,rodry,0.,0.,rofl,kfl,k0,por);
mugass=rosat.*vsgass.^2; kgass=rosat.*(vpgass.^2-(4/3)*vsgass.^2);
dataout=[p,vpdry,vsdry,vpsatm,vssatm,vpgass,vsgass,vphfsat,vshfsat];

subplot(2,1,1)
h=plot(p,vpdry,'o',p,vpsatm,'*',p,vpgass,'--',p,vphfsat,'-');
set(h,'markersize',6), set(h,'linewidth',2)
%plot(p,kdry,'o',p,kgass,'x',p,khfsat*ones(size(p)),'*');
xlabel('Effective Pressure (MPa)'), ylabel('Vp (m/s)');
subplot(2,1,2)
h=plot(p,vsdry,'o',p,vssatm,'*',p,vsgass,'--',p,vshfsat,'-');
set(h,'markersize',6), set(h,'linewidth',2)
xlabel('Effective Pressure (MPa)'), ylabel('Vs (m/s)');
legend(gca,'dry','sat.','Gassmann','squirt');
