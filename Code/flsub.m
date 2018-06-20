function varargout=flsub(velp1,vels1,poro,dens)
% [velp2,vels2,ss1,ss2,ttu]=flsub(velp1,vels1,poro,dens)
%
% Interactive fluid substitution allowing mouse-based selection of region of
% interest within velocity image VELP1 on which to do fluid substitution.
% Fluid and mineral properties, saturations, and other parameters
% are obtained from dialog boxes. Calculates and plots the initial and
% final velocity images, and normal incidence seismic sections.
% Required inputs: VELP1 (matrix), initial P velocity section;
%                  VELS1 (matrix), initial S velocity, give 0 if unavailable
%                  PORO  (matrix or scalar) porosity
% Outputs (all optional): VELP2, VELS2, P and S velocity after fluid 
% substitution; SS1, SS2, before and after seismic sections with time
% axis given by TTU vector. Optional input DENS (scalar or matrix same
% size as VELP1) is the density. Default=2200 kg/m3. This is required
% for seismogram calculation. Within the region of fluid substitution
% density is calculated from porosity, mineral modulus, fluid saturation
% and fluid densities.
%
% Fluid substitution is done using Gassmann's equations, with the uniform 
% effective fluid (Reuss average) saturation model.
% Fluid properties are from Batzle-Wang empirical relations.
% Requires reservoir fluid properties such as API, GOR, Salinity, etc. 
% Normal incidence seismograms are calculated from reflectivity and
% horizontal Fresnel zone calculations (see EZSEIS).
%
% See also FLSUBK, FLPROPUI, GASSMNV, GASSMNK, EZSEIS, EZSEIS2

%Written by T. Mukerji, 1997

if nargin<4, dens=2200; end; 
if prod(size(dens))==1,dens=dens.*ones(size(velp1));end;

imagesc(velp1),title('velocity'), colormap pink, colorbar, hold on, drawnow;
h=helpdlg({'Select polygonal region of interest using the mouse. Click twice or hit return when done'},'Select region'); drawnow; pause(2);

[bw,x,y]=roipoly; plot(x,y,'-w','linewidth',3); drawnow; hold off; 
if ishandle(h),delete(h); end;

indx=find(bw); vp1=velp1(indx); 
if prod(size(vels1))~=1, vs1=vels1(indx); else vs1=vels1; end;
if prod(size(poro))~=1, phi=poro(indx); else phi=poro; end;

defans={min(vp1), max(vp1), min(vs1), max(vs1), min(phi), max(phi)};
for k=1:length(defans),defans{k}=num2str(defans{k});end;
prompt={'min Vp','max Vp','min Vs','max Vs','min porosity','max porosity'};
limits=inputdlg(prompt,'Select Vp, Vs, and porosity limits',1,defans);

vpmin=str2num(limits{1}); vpmax=str2num(limits{2});
vsmin=str2num(limits{3}); vsmax=str2num(limits{4});
phimin=str2num(limits{5}); phimax=str2num(limits{6});
dxvp=find((vp1>=vpmin)&(vp1<=vpmax)); indx1=dxvp;
dxvs=find((vs1>=vsmin)&(vs1<=vsmax));
dxphi=find((phi>=phimin)&(phi<=phimax));
if size(vels1,1)*size(vels1,2)~=1, indx1=intersect(indx1,dxvs); end;
if size(poro,1)*size(poro,2)~=1, indx1=intersect(indx1,dxphi); end;
if length(indx1)<1, warndlg('no points satisfy limits'), drawnow, end;
vp1=vp1(indx1);
if size(vels1,1)*size(vels1,2)~=1, vs1=vs1(indx1); end;
if size(poro,1)*size(poro,2)~=1, phi=phi(indx1); end;

defans={2650,36,80,25,20,0.8,30000,300,0.9,0,0.1,0};
for k=1:length(defans),defans{k}=num2str(defans{k});end;
ssdefans={25, 0, 10, 10, 100};
for k=1:length(ssdefans),ssdefans{k}=num2str(ssdefans{k});end;
prmt={'Mineral density (kg/m3)', ...
      'Mineral bulk modulus (GPa)', ...
      'Temperature (Celsius)', ...
      'Pressure (MPa)',   ...
      'Oil API',        ...
      'Gas gravity',    ...
      'Brine salinity (ppm)', ...
      'GOR (L/L)', ...
      'S_oil initial', ...
      'S_gas initial', ...
      'S_oil final', ...
      'S_gas final'};
ssprmt={'Seismic frequency (Hz)', ...
        'Noise to Signal energy ratio', ...
        'Horizontal grid spacing dx (m)', ...
        'Vertical grid spacing dz (m)', ...
        'Depth to top of image (m)'};

button=1;
while 1
switch button
 case 0
break
 otherwise
flinpt=flsubdlg(prmt,'Mineral & Fluids',1,defans);
defans=flinpt; button=prod(size(flinpt));

if button ~= 0
ro0=str2num(flinpt{1}); k0=str2num(flinpt{2}).*1e9;
tmpr=str2num(flinpt{3}); press=str2num(flinpt{4});
oilapi=str2num(flinpt{5}); gg=str2num(flinpt{6});
sal=str2num(flinpt{7}); gor=str2num(flinpt{8});
so1=str2num(flinpt{9}); sg1=str2num(flinpt{10}); sw1=max([0,1-so1-sg1]);
so1=so1./(so1+sg1+sw1); sg1=sg1./(so1+sg1+sw1);
so2=str2num(flinpt{11}); sg2=str2num(flinpt{12}); sw2=max([0,1-so2-sg2]);
so2=so2./(so2+sg2+sw2); sg2=sg2./(so2+sg2+sw2);

[kfl1,rofl1]=batzle(2,sal,oilapi,gg,gor,0,0.1,press,tmpr,so1,sg1);
[kfl2,rofl2]=batzle(2,sal,oilapi,gg,gor,0,0.1,press,tmpr,so2,sg2);
kfl1=kfl1.*1e9; kfl2=kfl2.*1e9; rofl1=rofl1.*1e3; rofl2=rofl2.*1e3;
ro1=(1-phi).*ro0+rofl1.*phi; ro2=(1-phi).*ro0+rofl2.*phi;
[vp2,vs2]=gassmnv(vp1,vs1,ro1,rofl1,kfl1,rofl2,kfl2,k0,phi);

if length(vp1>0), 
%hist(vp1); hold on; hist(vp2); drawnow;
npt=length(vp1); mnvp1=mean(vp1); mnvp2=mean(vp2);
dvp=vp2-vp1; maxdvp=max(dvp); mindvp=min(dvp); 
rmsdvp=100*sqrt(mean(dvp.^2))/mean(vp1);
statprmt={'No. of points selected',...
          'K_fluid initial (GPa)',...
          'K_fluid final (GPa)',...
          'Mean Velocity initial (m/s)',...
          'Mean Velocity final (m/s)', ...
          'Min. change (final - initial) (m/s)', ...
          'Max. change (final - initial) (m/s)', ...
	  'RMS change (percent)'};
statdefans={npt,kfl1./1e9,kfl2./1e9,mnvp1,mnvp2,mindvp,maxdvp,rmsdvp};
outputdlg(statprmt,'Fluid sub stats',1,statdefans); drawnow;
end;

velp2=velp1; velp2(indx(indx1))=vp2; 
vels2=vels1; vels2(indx(indx1))=vs2;
dens1=dens; dens1(indx(indx1))=ro1; dens2=dens; dens2(indx(indx1))=ro2;

ssinpt=inputdlg(ssprmt, 'Seismogram',1,ssdefans);
ssdefans=ssinpt;
freq=str2num(ssinpt{1}); snr=str2num(ssinpt{2});
dx=str2num(ssinpt{3}); dz=str2num(ssinpt{4}); top=str2num(ssinpt{5});

h1=figure('units','normalized','position',[0.02 0.35 0.48 0.56]);
zax=[top:dz:dz*(size(velp1,1)-1)]; xax=[0:dx:dx*(size(velp1,2)-1)];
caxv=[ min([velp1(:);velp2(:)]),max([velp1(:);velp2(:)]) ];
subplot(311), imagesc(xax,zax,velp1); caxis(caxv);  colorbar;
ylabel('depth (m)'); title('Velocity before');
subplot(312), imagesc(xax,zax,velp2); caxis(caxv);  colorbar;
ylabel('depth (m)'); title('Velocity after');
subplot(313), imagesc(xax,zax,velp2-velp1); colorbar;
xlabel('distance'); ylabel('depth (m)'); title('Velocity after-before');
colormap pink; zoom(h1,'on'); drawnow; 

h2=figure('units','normalized','position',[0.51 0.35 0.48 0.56]);
[ss1,ss2,ttu]=ezseis2(velp1,dens1,velp2,dens2,dx,dz,top,freq,snr); clf;
subplot(311), seisrwb(ss1,ttu(1),ttu(2)-ttu(1),0,dx,0,ss1);
xlabel(''); ylabel('time'); title('Seismograms before');
subplot(312), seisrwb(ss2,ttu(1),ttu(2)-ttu(1),0,dx,0,ss1);
xlabel(''); ylabel('time'); title('Seismograms after');
subplot(313), seisrwb(ss2-ss1,ttu(1),ttu(2)-ttu(1),0,dx,0,ss1);
xlabel('distance'); ylabel('time'); title('Seismograms after-before');
zoom(h2,'on'); drawnow; pause(2);

end %if button~=0
end %switch
end %while

outputcell={velp2,vels2,ss1,ss2,ttu};
for k=1:nargout, varargout(k)=outputcell(k); end;

