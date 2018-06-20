function [velp2,vels2]=fldsubimg(velp1,vels1,ro1,rofl1,kfl1,rofl2,kfl2,k0,poro)
%[velp2,vels2]=fldsubimg(velp1,vels1,ro1,rofl1,kfl1,rofl2,kfl2,k0,poro)

imagesc(velp1),colorbar, hold on, drawnow;
h=helpdlg({'Select polygonal region of interest using the mouse. Click twice or hit return when done'},'Select region'); drawnow; pause(2);

[bw,x,y]=roipoly; plot(x,y,'-w','linewidth',3); drawnow; hold off; 
if ishandle(h),delete(h); end;

dx=find(bw); vp1=velp1(dx); 
if size(vels1,1)*size(vels1,2)~=1, vs1=vels1(dx); else vs1=vels1; end;
if size(poro,1)*size(poro,2)~=1, phi=poro(dx); else phi=poro; end;

defans={min(vp1), max(vp1), min(vs1), max(vs1), min(phi), max(phi)};
prompt={'min Vp','max Vp','min Vs','max Vs','min porosity','max porosity'};
limits=inputdlg(prompt,'Select Vp, Vs, and porosity limits',1,defans);

vpmin=str2num(limits{1}); vpmax=str2num(limits{2});
vsmin=str2num(limits{3}); vsmax=str2num(limits{4});
phimin=str2num(limits{5}); phimax=str2num(limits{6});
dxvp=find((vp1>=vpmin)&(vp1<=vpmax)); dx1=dxvp;
dxvs=find((vs1>=vsmin)&(vs1<=vsmax));
dxphi=find((phi>=phimin)&(phi<=phimax));
if size(vels1,1)*size(vels1,2)~=1, dx1=intersect(dx1,dxvs); end;
if size(poro,1)*size(poro,2)~=1, dx1=intersect(dx1,dxphi); end;
if length(dx1)<1, warndlg('no points satisfy limits'), drawnow, end;
vp1=vp1(dx1);
if size(vels1,1)*size(vels1,2)~=1, vs1=vs1(dx1); end;
if size(poro,1)*size(poro,2)~=1, phi=phi(dx1); end;
if size(ro1,1)*size(ro1,2)~=1, ro1=ro1(dx1); end;
if size(rofl1,1)*size(rofl1,2)~=1, rofl1=rofl1(dx1); end;
if size(kfl1,1)*size(kfl1,2)~=1, kfl1=kfl1(dx1); end;
if size(rofl2,1)*size(rofl2,2)~=1, rofl2=rofl2(dx1); end;
if size(kfl2,1)*size(kfl2,2)~=1, kfl2=kfl2(dx1); end;
if size(k0,1)*size(k0,2)~=1, k0=k0(dx1); end;

[vp2,vs2]=gassmnv(vp1,vs1,ro1,rofl1,kfl1,rofl2,kfl2,k0,phi);
clf; if length(vp1>0), hist(vp1); hold on; hist(vp2); drawnow; pause(3); end;
velp2=velp1; velp2(dx(dx1))=vp2; vels2=vels1; vels2(dx(dx1))=vs2;
cax=[ min([velp1(:);velp2(:)]),max([velp1(:);velp2(:)]) ];

subplot(2,1,1), imagesc(velp1), caxis(cax), title('before'), colorbar;
subplot(2,1,2), imagesc(velp2), caxis(cax), title('after'), colorbar;
